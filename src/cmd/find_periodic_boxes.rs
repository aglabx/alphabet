use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::thread;
use crate::io::read_fasta_strings;

pub fn run_from_args(args: Vec<String>) {
    if args.len() < 2 {
        eprintln!("Usage: find_periodic_boxes <input.fasta> [threads]");
        std::process::exit(1);
    }
    let input = &args[1];
    let n_threads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(32);

    let boxlen: usize = 16;
    let min_spacing: usize = 330;
    let max_spacing: usize = 360;
    let min_shared: usize = 8;
    let min_fixed: usize = 7;
    let max_fixed: usize = 12;
    let min_pairs: usize = 50;

    let arrays = Arc::new(read_fasta_strings(input).unwrap_or_else(|e| panic!("{:?}", e)));
    let total_bp: usize = arrays.iter().map(|(_, s)| s.len()).sum();
    let n_mon = total_bp as f64 / 717.0;
    eprintln!("{} arrays, {:.2}Mb, ~{:.0} monomers", arrays.len(), total_bp as f64 / 1e6, n_mon);

    // Phase 1: collect degenerate patterns from 330-360bp pairs
    eprintln!("Phase 1: collecting degenerate patterns from {}-{}bp pairs ({} threads)...", min_spacing, max_spacing, n_threads);

    let pattern_counts = Arc::new(Mutex::new(HashMap::<String, usize>::new()));
    let n_arrays = arrays.len();

    let mut handles = Vec::new();
    for t in 0..n_threads {
        let arrays = Arc::clone(&arrays);
        let pattern_counts = Arc::clone(&pattern_counts);
        handles.push(thread::spawn(move || {
            let mut local: HashMap<String, usize> = HashMap::new();
            let mut idx = t;
            while idx < n_arrays {
                let seq = arrays[idx].1.as_bytes();
                for i in 0..seq.len().saturating_sub(max_spacing + boxlen) {
                    let k1 = &seq[i..i + boxlen];
                    if k1.contains(&b'N') { continue; }

                    for offset in min_spacing..=max_spacing {
                        let j = i + offset;
                        if j + boxlen > seq.len() { break; }
                        let k2 = &seq[j..j + boxlen];
                        if k2.contains(&b'N') { continue; }

                        let mut shared = 0usize;
                        for p in 0..boxlen {
                            if k1[p] == k2[p] { shared += 1; }
                        }
                        if shared < min_shared { continue; }

                        let mut pat = Vec::with_capacity(boxlen);
                        for p in 0..boxlen {
                            if k1[p] == k2[p] {
                                pat.push(k1[p]);
                            } else {
                                pat.push(b'N');
                            }
                        }
                        let n_fixed = pat.iter().filter(|&&b| b != b'N').count();
                        if n_fixed >= min_fixed && n_fixed <= max_fixed {
                            let key = unsafe { String::from_utf8_unchecked(pat) };
                            *local.entry(key).or_insert(0) += 1;
                        }
                    }
                }
                idx += n_threads;
            }
            // Merge into global
            let mut global = pattern_counts.lock().unwrap();
            for (k, v) in local {
                *global.entry(k).or_insert(0) += v;
            }
        }));
    }
    for h in handles { h.join().unwrap(); }

    let pattern_counts = Arc::try_unwrap(pattern_counts).unwrap().into_inner().unwrap();
    eprintln!("{} distinct degenerate patterns found", pattern_counts.len());

    let mut sorted_patterns: Vec<_> = pattern_counts.into_iter()
        .filter(|(_, c)| *c >= min_pairs)
        .collect();
    sorted_patterns.sort_by_key(|b| std::cmp::Reverse(b.1));
    sorted_patterns.truncate(200);

    eprintln!("{} patterns with >= {} pairs", sorted_patterns.len(), min_pairs);

    // Phase 2: score periodicity
    eprintln!("Phase 2: scoring periodicity of top patterns...");

    #[derive(Debug)]
    struct PatResult {
        pattern: String,
        n_fixed: usize,
        pair_count: usize,
        total_hits: usize,
        on_target: usize,
        total_spacings: usize,
        tigd4_overlap: usize,
        spacing_dist: HashMap<usize, usize>,
    }

    let results_lock = Arc::new(Mutex::new(Vec::<PatResult>::new()));
    let sorted_patterns = Arc::new(sorted_patterns);

    let mut handles = Vec::new();
    let n_pats = sorted_patterns.len();
    for t in 0..n_threads.min(n_pats) {
        let arrays = Arc::clone(&arrays);
        let sorted_patterns = Arc::clone(&sorted_patterns);
        let results_lock = Arc::clone(&results_lock);
        handles.push(thread::spawn(move || {
            let mut pi = t;
            while pi < n_pats {
                let (ref pat, pair_count) = sorted_patterns[pi];
                let pat_bytes = pat.as_bytes();
                let mut all_positions: Vec<Vec<usize>> = Vec::new();
                let mut total_hits = 0usize;

                for (_, seq) in arrays.iter() {
                    let seq = seq.as_bytes();
                    let mut positions = Vec::new();
                    for i in 0..seq.len().saturating_sub(boxlen - 1) {
                        let mut ok = true;
                        for p in 0..boxlen {
                            if pat_bytes[p] != b'N' && seq[i + p] != pat_bytes[p] {
                                ok = false;
                                break;
                            }
                        }
                        if ok {
                            positions.push(i);
                            total_hits += 1;
                        }
                    }
                    if !positions.is_empty() {
                        all_positions.push(positions);
                    }
                }

                let mut spacing_dist: HashMap<usize, usize> = HashMap::new();
                let mut on_target = 0usize;
                let mut total_spacings = 0usize;

                for positions in &all_positions {
                    for i in 0..positions.len().saturating_sub(1) {
                        let d = positions[i + 1] - positions[i];
                        let bin = d / 10 * 10;
                        *spacing_dist.entry(bin).or_insert(0) += 1;
                        total_spacings += 1;
                        if d >= min_spacing && d <= max_spacing {
                            on_target += 1;
                        }
                    }
                }

                let n_fixed = pat_bytes.iter().filter(|&&b| b != b'N').count();
                let tigd4 = b"TTNNGGNNANGGNNGG";
                let mut overlap = 0;
                for p in 0..boxlen.min(16) {
                    if pat_bytes[p] != b'N' && tigd4[p] != b'N' && pat_bytes[p] == tigd4[p] {
                        overlap += 1;
                    }
                }

                results_lock.lock().unwrap().push(PatResult {
                    pattern: pat.clone(),
                    n_fixed,
                    pair_count,
                    total_hits,
                    on_target,
                    total_spacings,
                    tigd4_overlap: overlap,
                    spacing_dist,
                });

                pi += n_threads.min(n_pats);
            }
        }));
    }
    for h in handles { h.join().unwrap(); }

    let mut results = Arc::try_unwrap(results_lock).unwrap().into_inner().unwrap();
    results.sort_by_key(|b| std::cmp::Reverse(b.on_target));

    println!("# Periodic box candidates");
    println!("# Target spacing: {}-{}bp", min_spacing, max_spacing);
    println!("# Reference: TTNNGGNNANGGNNGG (TIGD4, 9/16 fixed)");
    println!("# {} arrays, {:.2}Mb, ~{:.0} monomers", arrays.len(), total_bp as f64 / 1e6, n_mon);
    println!("# rank\tpattern\tfixed\tpairs\thits\thits/mon\ton_target\ttotal_sp\tfrac%\ttag[overlap]\ttop_spacings");

    for (i, r) in results.iter().enumerate().take(50) {
        let frac = if r.total_spacings > 0 { r.on_target as f64 / r.total_spacings as f64 * 100.0 } else { 0.0 };
        let per_mon = r.total_hits as f64 / n_mon;
        let tag = if r.tigd4_overlap >= 4 { "TIGD4" } else { "NEW" };

        let mut sorted_sp: Vec<_> = r.spacing_dist.iter().collect();
        sorted_sp.sort_by(|a, b| b.1.cmp(a.1));
        let top_sp: String = sorted_sp.iter().take(5)
            .map(|(d, c)| format!("{}:{}", d, c))
            .collect::<Vec<_>>().join(" ");

        println!("{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{:.1}\t{}[{}]\t{}",
            i + 1, r.pattern, r.n_fixed, r.pair_count,
            r.total_hits, per_mon, r.on_target, r.total_spacings,
            frac, tag, r.tigd4_overlap, top_sp);
    }

    let tigd4_count = results.iter().filter(|r| r.tigd4_overlap >= 4).count();
    let new_count = results.iter().filter(|r| r.tigd4_overlap < 4).count();
    eprintln!("\n=== SUMMARY ===");
    eprintln!("TIGD4-related: {}", tigd4_count);
    eprintln!("NEW candidates: {}", new_count);

    if new_count > 0 {
        eprintln!("\nTop NEW candidates:");
        for r in results.iter().filter(|r| r.tigd4_overlap < 4).take(10) {
            let frac = if r.total_spacings > 0 { r.on_target as f64 / r.total_spacings as f64 * 100.0 } else { 0.0 };
            eprintln!("  {} ({}/16 fixed) hits={} periodic={:.0}%", r.pattern, r.n_fixed, r.total_hits, frac);
        }
    }
}

