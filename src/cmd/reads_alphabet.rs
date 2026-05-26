use std::io::{BufRead, BufReader, Write};
use crate::monomer::hpc;

/// Pass 2: Convert satellite reads to alphabet notation.
/// Input: FASTA of satellite reads (from reads_extract).
/// For each read:
///   1. Homopolymer compress
///   2. Scan for chain anchor sites (both strands)
///   3. Assign letters between consecutive anchors
///   4. Output: read_id, length, n_monomers, letter_string, site_positions
///
/// Usage: reads_alphabet <satellite_reads.fasta> <chains.json>
pub fn run_from_args(args: Vec<String>) {
    if args.len() < 3 {
        eprintln!("Usage: reads_alphabet <satellite_reads.fasta> <chains.json>");
        std::process::exit(1);
    }
    let input = &args[1];
    let chains_path = &args[2];

    let sites = load_chain_sites(chains_path);
    let site_names: Vec<String> = (0..sites.len()).map(|i| {
        (b'a' + (i as u8)) as char
    }).map(|c| c.to_string()).collect();

    eprintln!("Loaded {} sites", sites.len());

    // Prepare HPC versions
    let sites_hpc: Vec<Vec<u8>> = sites.iter().map(|s| hpc(s)).collect();

    // Read FASTA, process each read
    let file = std::fs::File::open(input).expect("Cannot open input");
    let reader = BufReader::with_capacity(16 * 1024 * 1024, file);

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());

    writeln!(out, "read_id\tread_len\tn_sites\tn_monomers\tletter_string\tsite_positions").unwrap();

    let mut name = String::new();
    let mut seq = Vec::new();
    let mut total_reads = 0usize;
    let mut total_monomers = 0usize;

    for line in reader.lines() {
        let line = line.unwrap();
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() && !seq.is_empty() {
                process_read(&name, &seq, &sites, &sites_hpc, &site_names, &mut out,
                    &mut total_reads, &mut total_monomers);
            }
            name = rest.split_whitespace().next().unwrap_or(rest).to_string();
            seq.clear();
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !name.is_empty() && !seq.is_empty() {
        process_read(&name, &seq, &sites, &sites_hpc, &site_names, &mut out,
            &mut total_reads, &mut total_monomers);
    }

    out.flush().unwrap();
    eprintln!("\nDone: {} reads, {} monomers", total_reads, total_monomers);
}

fn process_read(
    name: &str,
    seq: &[u8],
    sites: &[Vec<u8>],
    _sites_hpc: &[Vec<u8>],
    site_names: &[String],
    out: &mut impl Write,
    total_reads: &mut usize,
    total_monomers: &mut usize,
) {
    let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();

    // Find all site occurrences using exact match (not HPC — too many FP)
    let mut hits: Vec<(usize, usize, usize)> = Vec::new(); // (pos, pos, site_idx)

    for (si, site) in sites.iter().enumerate() {
        if site.len() > seq_upper.len() { continue; }
        for i in 0..=seq_upper.len() - site.len() {
            if seq_upper[i..i + site.len()] == site[..] {
                hits.push((i, i, si));
            }
        }
    }

    // Sort by position
    hits.sort_by_key(|h| h.0);

    // Dedup: keep best hit within 5bp window
    let mut deduped: Vec<(usize, usize, usize)> = Vec::new();
    for h in &hits {
        if deduped.last().is_none_or(|last: &(usize, usize, usize)| h.0.abs_diff(last.0) > 3) {
            deduped.push(*h);
        }
    }

    if deduped.is_empty() {
        return;
    }

    // Build site_order string and letter string
    let _site_order: String = deduped.iter()
        .map(|(_, _, si)| site_names[*si].as_str())
        .collect::<Vec<_>>()
        .join("");

    // Build letter string: cut at first site of each period
    // Find the chain period in HPC space from consecutive identical sites
    let mut letters = Vec::new();
    let first_site = deduped[0].2;
    let mut mono_start_positions = Vec::new();

    for (i, h) in deduped.iter().enumerate() {
        if h.2 == first_site {
            mono_start_positions.push(i);
        }
    }

    // For each monomer (between consecutive first-site occurrences)
    for w in mono_start_positions.windows(2) {
        let start_idx = w[0];
        let end_idx = w[1];
        let sites_in_mono: String = deduped[start_idx..end_idx].iter()
            .map(|(_, _, si)| site_names[*si].as_str())
            .collect::<Vec<_>>()
            .join("");
        letters.push(sites_in_mono);
    }

    let n_monomers = letters.len();
    let letter_string = letters.join("|");

    // Site positions (original coords)
    let site_pos_str: String = deduped.iter()
        .map(|(_, orig, si)| format!("{}:{}", site_names[*si], orig))
        .collect::<Vec<_>>()
        .join(",");

    *total_reads += 1;
    *total_monomers += n_monomers;

    writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}",
        name, seq.len(), deduped.len(), n_monomers, letter_string, site_pos_str).unwrap();

    if (*total_reads).is_multiple_of(10000) {
        eprintln!("\r  {} reads, {} monomers", total_reads, total_monomers);
    }
}

fn load_chain_sites(path: &str) -> Vec<Vec<u8>> {
    let content = std::fs::read_to_string(path).expect("Cannot read chains file");
    let mut sites = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        for key in &["\"sequence\"", "\"seq\""] {
            if let Some(start) = line.find(key) {
                if let Some(colon) = line[start..].find(':') {
                    let after_colon = &line[start + colon + 1..];
                    if let Some(q1) = after_colon.find('"') {
                        if let Some(q2) = after_colon[q1 + 1..].find('"') {
                            let seq = &after_colon[q1 + 1..q1 + 1 + q2];
                            let seq_bytes = seq.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<_>>();
                            if seq_bytes.len() >= 6 && !sites.contains(&seq_bytes) {
                                sites.push(seq_bytes);
                            }
                        }
                    }
                }
            }
        }
    }
    sites
}
