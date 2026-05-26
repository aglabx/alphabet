use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read as IoRead};
use std::sync::{Arc, Mutex};
use std::thread;
use std::process::{Command, Stdio};
use crate::monomer::revcomp;

/// Streaming CENP-B/box scanner for FASTQ/FASTA (plain or gzipped).
/// Reads one record at a time via streaming, never loads all into memory.
/// For .gz files uses pigz/gzip via pipe — no Rust gz dependency needed.
///
/// Usage: scan_reads <input.fastq.gz|.fasta|.fastq> [threads] [--pattern NAME:PATTERN ...]
//
// See reads_extract.rs for the rationale on not `.wait()`ing the decompressor.
#[allow(clippy::zombie_processes)]
pub fn run_from_args(args: Vec<String>) {
    if args.len() < 2 {
        eprintln!("Usage: scan_reads <input.fastq.gz|.fasta|.fastq> [threads] [--pattern NAME:PATTERN ...]");
        std::process::exit(1);
    }
    let input = &args[1];
    let n_threads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);

    // Built-in patterns: chr-specific CENP-B boxes
    let mut patterns: Vec<(String, Vec<u8>, u32)> = vec![
        // (name, pattern, min_score) — N=any, letters=must match
        // Full CENP-B with score threshold
        ("CENPB_7of9".into(), b"NTTCGNNNNANNCGGGN".to_vec(), 7),
        ("CENPB_6of9".into(), b"NTTCGNNNNANNCGGGN".to_vec(), 6),
        // Chr5-specific (exact 17bp)
        ("chr5".into(), b"CTTCGTTGGAAACGGGT".to_vec(), 17),
        // Chr11-specific (exact 17bp)
        ("chr11".into(), b"TTTCGTTGGAAACGGGT".to_vec(), 17),
    ];

    // Parse --pattern NAME:SEQ[:MIN_SCORE]
    let mut i = 3;
    while i < args.len() {
        if args[i] == "--pattern" && i + 1 < args.len() {
            let spec = &args[i + 1];
            let parts: Vec<&str> = spec.splitn(3, ':').collect();
            if parts.len() >= 2 {
                let name = parts[0].to_string();
                let pat = parts[1].to_uppercase().into_bytes();
                let min_score = parts.get(2).and_then(|s| s.parse().ok())
                    .unwrap_or(pat.iter().filter(|&&b| b != b'N').count() as u32);
                patterns.push((name, pat, min_score));
            }
            i += 2;
        } else {
            i += 1;
        }
    }

    let is_gzipped = input.ends_with(".gz");
    let is_fastq = input.contains(".fastq") || input.contains(".fq");

    eprintln!("Input: {} ({}{})", input,
        if is_fastq { "FASTQ" } else { "FASTA" },
        if is_gzipped { ", gzipped" } else { "" });
    eprintln!("Threads: {}", n_threads);
    for (name, pat, ms) in &patterns {
        let n_fixed = pat.iter().filter(|&&b| b != b'N').count();
        eprintln!("  {}: {} ({}/{}bp, min_score={})", name, String::from_utf8_lossy(pat), n_fixed, pat.len(), ms);
    }

    // Open input stream (pigz for .gz, direct for plain)
    let reader: Box<dyn IoRead + Send> = if is_gzipped {
        let decompress = if which_exists("pigz") { "pigz" } else { "gzip" };
        eprintln!("Decompressing with: {} -dc", decompress);
        let child = Command::new(decompress)
            .args(["-dc", input])
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .expect("Failed to start decompressor");
        Box::new(child.stdout.unwrap())
    } else {
        let file = std::fs::File::open(input).expect("Cannot open file");
        Box::new(file)
    };

    let buf_reader = BufReader::with_capacity(32 * 1024 * 1024, reader);
    let patterns = Arc::new(patterns);
    let counters: Arc<Mutex<HashMap<String, usize>>> = Arc::new(Mutex::new(HashMap::new()));

    let mut lines = buf_reader.lines();
    let batch_size = 50000;
    let mut read_count = 0usize;
    let mut bp_count = 0usize;

    loop {
        // Read batch of sequences
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(batch_size);

        for _ in 0..batch_size {
            let seq_str = if is_fastq {
                let header = match lines.next() {
                    Some(Ok(l)) if l.starts_with('@') => l,
                    Some(Ok(_)) => continue,
                    _ => break,
                };
                let _ = header;
                let seq = match lines.next() {
                    Some(Ok(l)) => l,
                    _ => break,
                };
                let _ = lines.next(); // +
                let _ = lines.next(); // qual
                seq
            } else {
                let header = match lines.next() {
                    Some(Ok(l)) if l.starts_with('>') => l,
                    _ => break,
                };
                let _ = header;
                match lines.next() {
                    Some(Ok(l)) => l,
                    _ => break,
                }
            };

            bp_count += seq_str.len();
            read_count += 1;
            batch.push(seq_str.into_bytes());
        }

        if batch.is_empty() {
            break;
        }

        // Process batch in parallel
        let chunk_size = batch.len().div_ceil(n_threads);
        let batch_arc = Arc::new(batch);
        let mut handles = Vec::new();

        for t in 0..n_threads {
            let batch_arc = Arc::clone(&batch_arc);
            let patterns = Arc::clone(&patterns);
            let counters = Arc::clone(&counters);
            let start = t * chunk_size;
            let end = ((t + 1) * chunk_size).min(batch_arc.len());

            handles.push(thread::spawn(move || {
                let mut local: HashMap<String, usize> = HashMap::new();

                for idx in start..end {
                    let seq = &batch_arc[idx];
                    let rc = revcomp(seq);

                    for (name, pat, min_score) in patterns.iter() {
                        let plen = pat.len();
                        if seq.len() < plen { continue; }
                        let n_fixed = pat.iter().filter(|&&b| b != b'N').count() as u32;

                        for strand_seq in [&seq[..], &rc[..]] {
                            for i in 0..=strand_seq.len() - plen {
                                let mut score = 0u32;
                                for j in 0..plen {
                                    if pat[j] != b'N' && strand_seq[i + j].to_ascii_uppercase() == pat[j] {
                                        score += 1;
                                    }
                                }
                                // For exact patterns (min_score == n_fixed), require all match
                                // For CENP-B style, require >= min_score
                                if *min_score == n_fixed {
                                    if score == n_fixed {
                                        *local.entry(name.clone()).or_insert(0) += 1;
                                    }
                                } else if score >= *min_score {
                                    *local.entry(name.clone()).or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }

                let mut g = counters.lock().unwrap();
                for (k, v) in local {
                    *g.entry(k).or_insert(0) += v;
                }
            }));
        }

        for h in handles { h.join().unwrap(); }

        // Progress every 100K reads
        if read_count % 100000 < batch_size {
            let c = counters.lock().unwrap();
            eprint!("\r  {:.1}M reads, {:.1}Gb", read_count as f64 / 1e6, bp_count as f64 / 1e9);
            for (name, _, _) in patterns.iter() {
                eprint!("  {}:{}", name, c.get(name).copied().unwrap_or(0));
            }
        }
    }

    eprintln!("\n\nDone: {} reads, {:.2}Gb", read_count, bp_count as f64 / 1e9);

    // Output
    let c = counters.lock().unwrap();
    let total_cenpb = c.get("CENPB_7of9").copied().unwrap_or(1).max(1);

    println!("# scan_reads results");
    println!("# input: {}", input);
    println!("# reads: {}", read_count);
    println!("# total_bp: {}", bp_count);
    println!("# total_Gb: {:.2}", bp_count as f64 / 1e9);
    println!("pattern\thits\thits_per_Mb\tpct_of_CENPB");

    for (name, _, _) in patterns.iter() {
        let cnt = c.get(name).copied().unwrap_or(0);
        let per_mb = cnt as f64 / (bp_count as f64 / 1e6);
        let pct = cnt as f64 / total_cenpb as f64 * 100.0;
        println!("{}\t{}\t{:.4}\t{:.2}%", name, cnt, per_mb, pct);
    }
}

fn which_exists(cmd: &str) -> bool {
    Command::new("which").arg(cmd).output().map(|o| o.status.success()).unwrap_or(false)
}
