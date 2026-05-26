//! FASTA / monomer-TSV I/O.
//!
//! All public functions return `anyhow::Result<T>` and attach `with_context`
//! at every fallible step. The CLI dispatchers in `cmd/*` typically `?` or
//! `.expect()` the result, so when something is corrupt the user sees the
//! full chain — file path, line number, parse cause — rather than a bare
//! `unwrap()` panic without provenance.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::fs::File;

use anyhow::{anyhow, Context, Result};

use crate::monomer::Monomer;

/// Read a FASTA file into `(name, sequence_bytes)` pairs. Sequences are
/// uppercased on read so downstream code can assume ASCII-upper ACGTN.
pub fn read_fasta(path: &str) -> Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(path).with_context(|| format!("opening FASTA {}", path))?;
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);
    let mut arrays: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| {
            format!("reading FASTA {} at line {}", path, lineno + 1)
        })?;
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() && !seq.is_empty() {
                arrays.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
            }
            name = rest.trim().to_string();
            seq.clear();
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !name.is_empty() && !seq.is_empty() {
        arrays.push((name, seq));
    }
    Ok(arrays)
}

/// Same as [`read_fasta`] but returns the sequence as a `String` — convenient
/// for callers that work with `&str` slices.
pub fn read_fasta_strings(path: &str) -> Result<Vec<(String, String)>> {
    read_fasta(path)?
        .into_iter()
        .map(|(n, s)| {
            String::from_utf8(s)
                .map(|s| (n.clone(), s))
                .with_context(|| format!("non-UTF8 bytes in FASTA {} for record {}", path, n))
        })
        .collect()
}

/// Read ArraySplitter .monomers.tsv and filter alpha satellite base monomers.
/// Returns (monomers grouped by array_id, header_indices)
pub fn read_monomers_tsv(
    path: &str,
    min_length: u32,
    max_length: u32,
) -> Result<(HashMap<String, Vec<Monomer>>, Stats)> {
    let file = File::open(path).with_context(|| format!("opening monomers TSV {}", path))?;
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);

    let mut arrays: HashMap<String, Vec<Monomer>> = HashMap::new();
    let mut stats = Stats::default();

    let mut lines = reader.lines();

    // Parse header
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow!("monomers TSV {} is empty (no header line)", path))?
        .with_context(|| format!("reading header of monomers TSV {}", path))?;
    let headers: Vec<&str> = header_line.split('\t').collect();
    let col = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| *h == name)
            .ok_or_else(|| anyhow!("column '{}' not found in header of {}", name, path))
    };

    let col_array_id = col("array_id")?;
    let col_type = col("type")?;
    let col_idx = col("idx")?;
    let col_length = col("length")?;
    let col_period = col("period")?;
    let col_source = col("source")?;
    let col_ed_prev = col("ed_prev")?;
    let col_ed_next = col("ed_next")?;
    let col_sequence = col("sequence")?;

    for (lineno, line_result) in lines.enumerate() {
        let line = line_result.with_context(|| {
            format!("reading monomers TSV {} at line {}", path, lineno + 2)
        })?;
        stats.total_rows += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= col_sequence {
            continue;
        }

        // Filter: base_monomer only
        if fields[col_type] != "base_monomer" {
            stats.filtered_type += 1;
            continue;
        }

        // Filter: period 171 or 172
        let period: u16 = match fields[col_period].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        if period != 171 && period != 172 {
            stats.filtered_period += 1;
            continue;
        }

        // Filter: length
        let length: u32 = match fields[col_length].parse() {
            Ok(l) => l,
            Err(_) => continue,
        };
        if length < min_length || length > max_length {
            stats.filtered_length += 1;
            continue;
        }

        let seq = fields[col_sequence];
        if seq.is_empty() || seq == "-" {
            continue;
        }

        let monomer = Monomer {
            array_id: fields[col_array_id].to_string(),
            idx: fields[col_idx].parse().unwrap_or(0),
            length,
            period,
            source: fields[col_source].to_string(),
            ed_prev: fields[col_ed_prev].to_string(),
            ed_next: fields[col_ed_next].to_string(),
            sequence: seq.to_string(),
            strand_flipped: false,
            rotation_offset: 0,
            anchor_confidence: 0.0,
        };

        arrays.entry(monomer.array_id.clone()).or_default().push(monomer);
    }

    // Sort monomers within each array by idx
    for monomers in arrays.values_mut() {
        monomers.sort_by_key(|m| m.idx);
    }

    stats.n_monomers = arrays.values().map(|v| v.len() as u64).sum();
    stats.n_arrays = arrays.len() as u64;

    Ok((arrays, stats))
}

#[derive(Default, Debug, serde::Serialize)]
pub struct Stats {
    pub total_rows: u64,
    pub filtered_type: u64,
    pub filtered_period: u64,
    pub filtered_length: u64,
    pub n_monomers: u64,
    pub n_arrays: u64,
    pub arrays_kept_fwd: u64,
    pub arrays_flipped: u64,
    pub seed_length: usize,
}

/// Write canonical monomers to TSV
pub fn write_monomers_tsv(path: &str, arrays: &HashMap<String, Vec<Monomer>>) -> Result<()> {
    use std::io::Write;
    let file = File::create(path).with_context(|| format!("creating monomers TSV {}", path))?;
    let mut writer = std::io::BufWriter::with_capacity(64 * 1024 * 1024, file);

    writeln!(writer, "array_id\tidx\tlength\tperiod\tsource\ted_prev\ted_next\tstrand_flipped\trotation_offset\tanchor_confidence\tsequence")
        .with_context(|| format!("writing header to {}", path))?;

    // Collect and sort all monomers
    let mut all: Vec<&Monomer> = arrays.values().flat_map(|v| v.iter()).collect();
    all.sort_by(|a, b| a.array_id.cmp(&b.array_id).then(a.idx.cmp(&b.idx)));

    for m in all {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
            m.array_id, m.idx, m.length, m.period, m.source,
            m.ed_prev, m.ed_next, m.strand_flipped, m.rotation_offset,
            m.anchor_confidence, m.sequence
        )
        .with_context(|| format!("writing row to {}", path))?;
    }
    Ok(())
}
