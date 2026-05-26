use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};
use std::thread;
use crate::monomer::hpc;

/// Classify ONT/HiFi reads using existing CHM13 alphabet in HPC space.
///
/// 1. Loads CHM13 chains.json (motifs) + v4_alphabet.json (letter definitions)
/// 2. HPC-compresses motif sequences
/// 3. For each read: HPC → find chain sites → site_order → lookup letter
/// 4. Outputs per-read and per-monomer letter assignments
///
/// Usage: classify_reads <reads.fasta> <chains.json> <alphabet.json> [threads]
pub fn run_from_args(args: Vec<String>) {
    if args.len() < 4 {
        eprintln!("Usage: classify_reads <reads.fasta> <chains.json> <alphabet.json> [threads]");
        std::process::exit(1);
    }
    let input = &args[1];
    let chains_path = &args[2];
    let alphabet_path = &args[3];
    let n_threads: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(8);

    // Load motifs: use the 5 chain sites from motif_cut v4 (not the 12 v1 motifs)
    // These are the sites that define the alphabet in v4_alphabet.json
    let motifs = [Motif { name: "a".into(), seq: b"ACATCACAAAG".to_vec(), position: 0 },
        Motif { name: "b".into(), seq: b"AGAATGCTTCT".to_vec(), position: 19 },
        Motif { name: "c".into(), seq: b"GAAGATATTTC".to_vec(), position: 32 },
        Motif { name: "d".into(), seq: b"TCCACTTGCAG".to_vec(), position: 54 },
        Motif { name: "e".into(), seq: b"AAAGAGTGTTT".to_vec(), position: 111 }];
    let _ = load_motifs(chains_path); // still parse for validation
    let n_motifs = motifs.len();
    eprintln!("Loaded {} motifs", n_motifs);

    // HPC compress motifs
    let motifs_hpc: Vec<Vec<u8>> = motifs.iter().map(|m| hpc(&m.seq)).collect();
    let motifs_orig: Vec<Vec<u8>> = motifs.iter().map(|m| m.seq.clone()).collect();
    for (i, m) in motifs.iter().enumerate() {
        eprintln!("  {}: {} ({}bp) → HPC: {} ({}bp)",
            m.name, String::from_utf8_lossy(&m.seq), m.seq.len(),
            String::from_utf8_lossy(&motifs_hpc[i]), motifs_hpc[i].len());
    }

    // Load alphabet: key -> letter name
    let alphabet = load_alphabet(alphabet_path, n_motifs);
    eprintln!("Loaded {} letter definitions", alphabet.len());

    // Expected HPC period (from motif positions)
    let orig_period = 171usize;
    let hpc_ratio = 0.727; // from CHM13 data
    let hpc_period = (orig_period as f64 * hpc_ratio) as usize;
    eprintln!("HPC period estimate: {}bp", hpc_period);

    // Expected distances between consecutive motifs in HPC space
    let mut expected_dists: Vec<usize> = Vec::new();
    for i in 0..motifs.len() {
        let j = (i + 1) % motifs.len();
        let dist = if motifs[j].position > motifs[i].position {
            motifs[j].position - motifs[i].position
        } else {
            orig_period - motifs[i].position + motifs[j].position
        };
        let hpc_dist = (dist as f64 * hpc_ratio) as usize;
        expected_dists.push(hpc_dist);
    }

    // Read FASTA, process in batches
    let file = std::fs::File::open(input).unwrap();
    let reader = BufReader::with_capacity(16 * 1024 * 1024, file);

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, "read_id\tread_len\thpc_len\tn_monomers\tletters\tletter_counts").unwrap();

    let motifs_hpc = Arc::new(motifs_hpc);
    let motifs_orig = Arc::new(motifs_orig);
    let alphabet = Arc::new(alphabet);
    let expected_dists = Arc::new(expected_dists);

    let mut lines = reader.lines();
    let batch_size = 5000;
    let mut total_reads = 0usize;
    let mut total_monomers = 0usize;
    let letter_counts_global: Arc<Mutex<HashMap<String, usize>>> = Arc::new(Mutex::new(HashMap::new()));

    loop {
        let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(batch_size);

        // Read batch
        loop {
            let header = match lines.next() {
                Some(Ok(l)) if l.starts_with('>') => l,
                _ => break,
            };
            let seq = match lines.next() {
                Some(Ok(l)) => l.into_bytes(),
                _ => break,
            };
            let name = header[1..].split_whitespace().next().unwrap_or(&header[1..]).to_string();
            batch.push((name, seq));
            if batch.len() >= batch_size { break; }
        }

        if batch.is_empty() { break; }
        total_reads += batch.len();

        // Process batch in parallel
        let chunk_size = batch.len().div_ceil(n_threads);
        let batch = Arc::new(batch);
        let results: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let mut handles = Vec::new();

        for t in 0..n_threads {
            let batch = Arc::clone(&batch);
            let motifs_hpc = Arc::clone(&motifs_hpc);
            let motifs_orig = Arc::clone(&motifs_orig);
            let alphabet = Arc::clone(&alphabet);
            let _expected_dists = Arc::clone(&expected_dists);
            let results = Arc::clone(&results);
            let letter_counts_global = Arc::clone(&letter_counts_global);
            let start = t * chunk_size;
            let end = ((t + 1) * chunk_size).min(batch.len());

            handles.push(thread::spawn(move || {
                let mut local_lines = Vec::new();
                let mut local_counts: HashMap<String, usize> = HashMap::new();

                for idx in start..end {
                    let (ref name, ref seq) = batch[idx];
                    let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
                    let seq_hpc = hpc(&seq_upper);

                    // Find motifs in HPC read (try both original and HPC motifs)
                    let mut hits: Vec<(usize, usize)> = Vec::new(); // (pos_in_hpc, motif_idx)

                    // Primary: search HPC motifs in HPC sequence
                    for (mi, motif_hpc) in motifs_hpc.iter().enumerate() {
                        if motif_hpc.len() > seq_hpc.len() { continue; }
                        for i in 0..=seq_hpc.len() - motif_hpc.len() {
                            if seq_hpc[i..i + motif_hpc.len()] == motif_hpc[..] {
                                hits.push((i, mi));
                            }
                        }
                    }

                    // Also search original motifs in original sequence for better coverage
                    let mut orig_hits: Vec<(usize, usize)> = Vec::new();
                    for (mi, motif_orig) in motifs_orig.iter().enumerate() {
                        if motif_orig.len() > seq_upper.len() { continue; }
                        for i in 0..=seq_upper.len() - motif_orig.len() {
                            if seq_upper[i..i + motif_orig.len()] == motif_orig[..] {
                                // Map to HPC position
                                let hpc_pos = orig_to_hpc_pos(&seq_upper, i);
                                orig_hits.push((hpc_pos, mi));
                            }
                        }
                    }

                    // Merge hits (prefer orig if available, add HPC-only)
                    for h in &orig_hits {
                        if !hits.iter().any(|existing| existing.1 == h.1 && existing.0.abs_diff(h.0) < 3) {
                            hits.push(*h);
                        }
                    }

                    hits.sort_by_key(|h| h.0);

                    // Dedup within 3bp HPC window
                    let mut deduped: Vec<(usize, usize)> = Vec::new();
                    for h in &hits {
                        if deduped.last().is_none_or(|last| h.0.abs_diff(last.0) > 2) {
                            deduped.push(*h);
                        }
                    }

                    if deduped.len() < 2 { continue; }

                    // Cut monomers: find first-motif (motif 0) occurrences
                    let first_motif_idx = 0;
                    let cut_positions: Vec<usize> = deduped.iter()
                        .filter(|h| h.1 == first_motif_idx)
                        .map(|h| h.0)
                        .collect();

                    if cut_positions.len() < 2 { continue; }

                    // For each monomer, determine which motifs are present
                    let mut letters = Vec::new();
                    let n_motifs = motifs_hpc.len();

                    for w in 0..cut_positions.len() - 1 {
                        let mono_start = cut_positions[w];
                        let mono_end = cut_positions[w + 1];
                        let mono_len = mono_end - mono_start;

                        // Skip if too short or too long
                        if mono_len < hpc_period / 2 || mono_len > hpc_period * 3 { continue; }

                        // Which motifs are in this monomer?
                        let mut present = vec![false; n_motifs];
                        present[first_motif_idx] = true;
                        let mut dists = Vec::new();

                        let mono_hits: Vec<&(usize, usize)> = deduped.iter()
                            .filter(|h| h.0 >= mono_start && h.0 < mono_end)
                            .collect();

                        for h in &mono_hits {
                            present[h.1] = true;
                        }

                        // Compute distances between consecutive present sites
                        let mut prev_pos = mono_start;
                        for h in &mono_hits {
                            if h.0 > prev_pos {
                                dists.push(h.0 - prev_pos);
                            }
                            prev_pos = h.0;
                        }

                        // Build key: presence pattern + binned distances (v4 format: bin=20bp)
                        let presence: String = present.iter().map(|&p| if p { '1' } else { '0' }).collect();
                        // Distances between consecutive PRESENT sites, binned to 20bp
                        let mut key_dists = Vec::new();
                        let mut prev_pos_opt: Option<usize> = None;
                        for h in &mono_hits {
                            if let Some(pp) = prev_pos_opt {
                                let d = h.0 - pp;
                                // Scale back from HPC to original space for binning
                                let d_orig = (d as f64 / hpc_ratio) as usize;
                                let binned = (d_orig / 20) * 20;
                                key_dists.push(format!("{}", binned));
                            }
                            prev_pos_opt = Some(h.0);
                        }
                        let binned_str = key_dists.join(",");
                        let key = format!("{}:{}", presence, binned_str);

                        // Lookup letter: exact match first, then closest by presence pattern
                        let letter = alphabet.get(&key)
                            .map(|s| s.as_str())
                            .unwrap_or_else(|| {
                                // Find closest: match by presence pattern, then closest distances
                                find_closest_letter(&presence, &dists, &alphabet)
                            });

                        letters.push(letter.to_string());
                        *local_counts.entry(letter.to_string()).or_insert(0) += 1;
                    }

                    if letters.is_empty() { continue; }

                    let letter_str = letters.join("|");
                    local_lines.push(format!("{}\t{}\t{}\t{}\t{}",
                        name, seq.len(), seq_hpc.len(), letters.len(), letter_str));
                }

                results.lock().unwrap().extend(local_lines);
                let mut g = letter_counts_global.lock().unwrap();
                for (k, v) in local_counts { *g.entry(k).or_insert(0) += v; }
            }));
        }

        for h in handles { h.join().unwrap(); }

        // Write results
        let res = results.lock().unwrap();
        for line in res.iter() {
            writeln!(out, "{}", line).unwrap();
            total_monomers += line.matches('|').count() + 1;
        }

        if total_reads % 10000 < batch_size {
            eprintln!("  {} reads, {} monomers", total_reads, total_monomers);
        }
    }

    out.flush().unwrap();

    // Summary
    let counts = letter_counts_global.lock().unwrap();
    let mut sorted: Vec<_> = counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));

    eprintln!("\nDone: {} reads, {} monomers", total_reads, total_monomers);
    eprintln!("\nLetter distribution:");
    let total_letters: usize = sorted.iter().map(|(_, c)| **c).sum();
    for (letter, count) in sorted.iter().take(20) {
        let pct = **count as f64 / total_letters as f64 * 100.0;
        eprintln!("  {}: {} ({:.1}%)", letter, count, pct);
    }
    eprintln!("  ... {} total distinct letters", sorted.len());
}

fn find_closest_letter<'a>(presence: &str, _dists: &[usize], alphabet: &'a HashMap<String, String>) -> &'a str {
    // Find letter with same presence pattern, or closest one
    // First: exact presence match (any distances)
    let mut best_letter = "?";
    let mut best_score = 0usize;

    for (key, letter) in alphabet {
        let parts: Vec<&str> = key.splitn(2, ':').collect();
        if parts.is_empty() { continue; }
        let key_presence = parts[0];

        // Score = number of matching presence bits
        let score: usize = presence.chars().zip(key_presence.chars())
            .filter(|(a, b)| a == b)
            .count();

        // Prefer exact presence match, then most similar
        if key_presence == presence {
            return letter; // exact presence match — take first one
        }
        if score > best_score {
            best_score = score;
            best_letter = letter;
        }
    }
    best_letter
}

fn orig_to_hpc_pos(orig: &[u8], orig_pos: usize) -> usize {
    let mut hpc_pos = 0usize;
    let mut prev = 0u8;
    for (i, &b) in orig.iter().enumerate() {
        if b != prev {
            if i >= orig_pos { return hpc_pos; }
            hpc_pos += 1;
            prev = b;
        }
    }
    hpc_pos
}

struct Motif {
    name: String,
    seq: Vec<u8>,
    position: usize,
}

fn load_motifs(path: &str) -> Vec<Motif> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut motifs = Vec::new();
    let lines: Vec<&str> = content.lines().collect();

    for i in 0..lines.len() {
        let line = lines[i].trim();
        if !line.contains("\"sequence\"") { continue; }

        let mut name = String::new();
        let mut seq = String::new();
        let mut pos = 0usize;

        let start = i.saturating_sub(5);
        let end = (i + 5).min(lines.len());
        for j in start..end {
            let l = lines[j].trim();
            if l.contains("\"name\"") { name = extract_str(l); }
            if l.contains("\"sequence\"") { seq = extract_str(l); }
            if l.contains("\"position_in_monomer\"") { pos = extract_num(l); }
        }

        if seq.len() >= 6 {
            let seq_bytes = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
            if !motifs.iter().any(|m: &Motif| m.seq == seq_bytes) {
                motifs.push(Motif {
                    name: if name.is_empty() { format!("M{}", motifs.len()) } else { name },
                    seq: seq_bytes,
                    position: pos,
                });
            }
        }
    }
    motifs.sort_by_key(|m| m.position);
    motifs
}

fn load_alphabet(path: &str, _n_motifs: usize) -> HashMap<String, String> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut alphabet = HashMap::new();

    // Parse JSON manually: look for "key" and "name" pairs
    let lines: Vec<&str> = content.lines().collect();
    let mut current_key = String::new();

    for line in &lines {
        let line = line.trim();
        if line.contains("\"key\"") {
            current_key = extract_str(line);
        }
        if line.contains("\"name\"") && !current_key.is_empty() {
            let current_name = extract_str(line);
            if !current_name.is_empty() && !current_key.is_empty() {
                // Expand key to full n_motifs presence pattern if needed
                // Key format: "11111:20,20,40,20" (5 sites) → need to expand to n_motifs sites
                // For now, store as-is and match by presence prefix
                alphabet.insert(current_key.clone(), current_name);
            }
            current_key.clear();
        }
    }
    alphabet
}

fn extract_str(line: &str) -> String {
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

fn extract_num(line: &str) -> usize {
    if let Some(colon) = line.find(':') {
        let after = line[colon + 1..].trim().trim_matches(|c: char| !c.is_numeric());
        if let Ok(n) = after.parse::<f64>() { return n as usize; }
    }
    0
}
