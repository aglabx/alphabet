use std::io::{BufRead, BufReader, Write, Read as IoRead};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;
use crate::monomer::{hpc, revcomp};

/// Pass 1: Extract reads containing satellite DNA by chain grammar check.
/// Not just motif presence — checks ORDER and DISTANCES between motifs.
/// A read passes only if motifs appear in chain order at expected spacing.
///
/// Usage: reads_extract <input.fastq.gz> <chains.json> [min_chain_hits=3]
//
// The decompressor child is intentionally not `.wait()`ed: it's a single
// streaming pipe whose lifetime equals this whole subcommand's stdin reader.
// When `run_from_args` returns, the child gets SIGPIPE and exits; the OS
// reaps it before this process exits. Refactoring to track the Child handle
// would mean threading it through the `Box<dyn IoRead>` reader abstraction.
#[allow(clippy::zombie_processes)]
pub fn run_from_args(args: Vec<String>) {
    if args.len() < 3 {
        eprintln!("Usage: reads_extract <input.fastq.gz> <chains.json> [min_chain_hits=3]");
        std::process::exit(1);
    }
    let input = &args[1];
    let chains_path = &args[2];
    let min_chain_hits: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);
    let n_threads: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(8);

    // Load chain: motifs with positions
    let chain = load_chain(chains_path);
    let period = chain.period;
    eprintln!("Chain: period={}bp, {} motifs", period, chain.motifs.len());
    for m in &chain.motifs {
        eprintln!("  {} pos={} seq={} ({}bp) hpc={}",
            m.name, m.position, String::from_utf8_lossy(&m.seq), m.seq.len(),
            String::from_utf8_lossy(&m.seq_hpc));
    }
    eprintln!("Min chain hits: {}", min_chain_hits);

    // Expected pairwise distances between consecutive motifs
    let mut expected_dists: Vec<(usize, usize, usize)> = Vec::new(); // (motif_i, motif_j, expected_dist)
    for i in 0..chain.motifs.len() {
        let j = (i + 1) % chain.motifs.len();
        let dist = if chain.motifs[j].position > chain.motifs[i].position {
            chain.motifs[j].position - chain.motifs[i].position
        } else {
            period - chain.motifs[i].position + chain.motifs[j].position
        };
        expected_dists.push((i, j, dist));
        eprintln!("  {} -> {}: expected {}bp", chain.motifs[i].name, chain.motifs[j].name, dist);
    }

    let is_gzipped = input.ends_with(".gz");
    let is_fastq = input.contains(".fastq") || input.contains(".fq");

    let reader: Box<dyn IoRead + Send> = if is_gzipped {
        let decompress = if which_exists("pigz") { "pigz" } else { "gzip" };
        eprintln!("Streaming: {} -dc {}", decompress, input);
        let child = Command::new(decompress)
            .args(["-dc", input])
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("Failed to start decompressor");
        Box::new(child.stdout.unwrap())
    } else {
        Box::new(std::fs::File::open(input).expect("Cannot open file"))
    };

    let buf_reader = BufReader::with_capacity(32 * 1024 * 1024, reader);
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::with_capacity(16 * 1024 * 1024, stdout.lock());

    let mut lines = buf_reader.lines();
    let mut total_reads = 0usize;
    let mut passed_reads = 0usize;
    let mut total_bp = 0usize;
    let mut passed_bp = 0usize;

    let chain = Arc::new(chain);
    let expected_dists = Arc::new(expected_dists);
    let batch_size = 10000;

    loop {
        // Read batch
        let mut batch: Vec<(String, String)> = Vec::with_capacity(batch_size);
        for _ in 0..batch_size {
            let header = match lines.next() {
                Some(Ok(l)) => l,
                _ => break,
            };
            let seq = match lines.next() {
                Some(Ok(l)) => l,
                _ => break,
            };
            if is_fastq {
                let _ = lines.next();
                let _ = lines.next();
            }
            total_bp += seq.len();
            total_reads += 1;
            batch.push((header, seq));
        }

        if batch.is_empty() { break; }

        // Process batch in parallel
        let chunk_size = batch.len().div_ceil(n_threads);
        let batch_arc = Arc::new(batch);
        let mut handles = Vec::new();
        let results: Arc<Mutex<Vec<(String, String, usize)>>> = Arc::new(Mutex::new(Vec::new()));

        for t in 0..n_threads {
            let batch_arc = Arc::clone(&batch_arc);
            let chain = Arc::clone(&chain);
            let expected_dists = Arc::clone(&expected_dists);
            let results = Arc::clone(&results);
            let start = t * chunk_size;
            let end = ((t + 1) * chunk_size).min(batch_arc.len());

            handles.push(thread::spawn(move || {
                let mut local_results = Vec::new();
                for idx in start..end {
                    let (ref header, ref seq) = batch_arc[idx];
                    let seq_upper: Vec<u8> = seq.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
                    let score_fwd = check_chain_grammar(&seq_upper, &chain, &expected_dists, period);
                    let best = if score_fwd >= min_chain_hits {
                        score_fwd
                    } else {
                        let rc = revcomp(&seq_upper);
                        let score_rc = check_chain_grammar(&rc, &chain, &expected_dists, period);
                        score_fwd.max(score_rc)
                    };
                    if best >= min_chain_hits {
                        local_results.push((header.clone(), seq.clone(), best));
                    }
                }
                results.lock().unwrap().extend(local_results);
            }));
        }
        for h in handles { h.join().unwrap(); }

        // Write results
        let res = results.lock().unwrap();
        for (header, seq, score) in res.iter() {
            passed_reads += 1;
            passed_bp += seq.len();
            let name = if header.starts_with('@') || header.starts_with('>') {
                &header[1..]
            } else {
                header
            };
            writeln!(out, ">{} chain_score={}", name.split_whitespace().next().unwrap_or(name), score).unwrap();
            writeln!(out, "{}", seq).unwrap();
        }

        if total_reads % 100000 < batch_size {
            eprintln!("  {} reads, {} passed ({:.1}%), {:.1}Gb",
                total_reads, passed_reads,
                if total_reads > 0 { passed_reads as f64 / total_reads as f64 * 100.0 } else { 0.0 },
                total_bp as f64 / 1e9);
        }
    }

    out.flush().unwrap();

    eprintln!("\nDone:");
    eprintln!("  Total: {} reads, {:.2}Gb", total_reads, total_bp as f64 / 1e9);
    eprintln!("  Passed: {} reads ({:.2}%), {:.2}Gb",
        passed_reads,
        if total_reads > 0 { passed_reads as f64 / total_reads as f64 * 100.0 } else { 0.0 },
        passed_bp as f64 / 1e9);
}

/// Check if sequence contains chain motifs in correct order with correct distances.
/// Returns number of chain-consistent motif pairs found.
/// Optimized: early exit once min_required reached, only check nearby hits.
fn check_chain_grammar(
    seq: &[u8],
    chain: &Chain,
    expected_dists: &[(usize, usize, usize)],
    period: usize,
) -> usize {
    // Find all motif hits using ORIGINAL sequences (not HPC — too many FP with short HPC motifs)
    // ONT accuracy ~95%, 11bp exact match hits ~55% of true sites — enough for chain grammar
    let mut hits: Vec<(usize, usize)> = Vec::new(); // (pos, motif_idx)

    for (mi, motif) in chain.motifs.iter().enumerate() {
        if motif.seq.len() > seq.len() { continue; }
        for i in 0..=seq.len() - motif.seq.len() {
            if seq[i..i + motif.seq.len()] == motif.seq[..] {
                hits.push((i, mi));
            }
        }
    }

    if hits.len() < 2 { return 0; }
    hits.sort_by_key(|h| h.0);

    // Only check consecutive hits (not all pairs) — O(N) not O(N²)
    let tolerance = (period as f64 * 0.3) as usize;
    let mut consistent_pairs = 0usize;
    let max_look_ahead = 5; // only check next 5 hits, not all

    for w in 0..hits.len() - 1 {
        let (pos_a, mi_a) = hits[w];
        let end = (w + max_look_ahead).min(hits.len());
        for v in w + 1..end {
            let (pos_b, mi_b) = hits[v];
            if pos_b - pos_a > period * 2 { break; }

            for &(exp_i, exp_j, exp_dist) in expected_dists {
                if mi_a == exp_i && mi_b == exp_j {
                    let actual_dist = pos_b - pos_a;
                    if actual_dist.abs_diff(exp_dist) <= tolerance {
                        consistent_pairs += 1;
                        if consistent_pairs >= 3 { return consistent_pairs; } // early exit: it's satellite
                    }
                }
            }
        }
    }

    consistent_pairs
}

fn which_exists(cmd: &str) -> bool {
    Command::new("which").arg(cmd).output().map(|o| o.status.success()).unwrap_or(false)
}

struct Motif {
    name: String,
    seq: Vec<u8>,
    seq_hpc: Vec<u8>,
    position: usize,
}

struct Chain {
    period: usize,
    motifs: Vec<Motif>,
}

fn load_chain(path: &str) -> Chain {
    let content = std::fs::read_to_string(path).expect("Cannot read chains file");

    // Extract period
    let mut period = 171usize; // default
    for line in content.lines() {
        if line.contains("\"target_period\"") {
            if let Some(colon) = line.find(':') {
                let val = line[colon + 1..].trim().trim_matches(|c: char| !c.is_numeric());
                if let Ok(p) = val.parse::<usize>() {
                    if p > 0 { period = p; }
                }
            }
        }
    }

    // Extract motifs
    let mut motifs = Vec::new();
    let mut i = 0;
    let lines: Vec<&str> = content.lines().collect();

    while i < lines.len() {
        let line = lines[i].trim();

        // Look for motif objects with "sequence" and "position_in_monomer"
        if line.contains("\"sequence\"") && line.contains(":") {
            let mut name = String::new();
            let mut seq = String::new();
            let mut pos = 0usize;

            // Scan nearby lines for name, sequence, position
            let start = i.saturating_sub(5);
            let end = (i + 5).min(lines.len());
            for j in start..end {
                let l = lines[j].trim();
                if l.contains("\"name\"") {
                    name = extract_string_value(l);
                }
                if l.contains("\"sequence\"") || l.contains("\"seq\"") {
                    seq = extract_string_value(l);
                }
                if l.contains("\"position_in_monomer\"") || l.contains("\"position\"") {
                    pos = extract_number_value(l);
                }
            }

            if !seq.is_empty() && seq.len() >= 6 {
                let seq_bytes: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
                let seq_hpc = hpc(&seq_bytes);
                if !motifs.iter().any(|m: &Motif| m.seq == seq_bytes) {
                    motifs.push(Motif {
                        name: if name.is_empty() { format!("M{}", motifs.len()) } else { name },
                        seq: seq_bytes,
                        seq_hpc,
                        position: pos,
                    });
                }
            }
        }
        i += 1;
    }

    // Sort by position
    motifs.sort_by_key(|m| m.position);

    Chain { period, motifs }
}

fn extract_string_value(line: &str) -> String {
    if let Some(colon) = line.find(':') {
        let after = &line[colon + 1..];
        if let Some(q1) = after.find('"') {
            if let Some(q2) = after[q1 + 1..].find('"') {
                return after[q1 + 1..q1 + 1 + q2].to_string();
            }
        }
    }
    String::new()
}

fn extract_number_value(line: &str) -> usize {
    if let Some(colon) = line.find(':') {
        let after = line[colon + 1..].trim().trim_matches(|c: char| !c.is_numeric() && c != '.');
        if let Ok(n) = after.parse::<f64>() {
            return n as usize;
        }
    }
    0
}
