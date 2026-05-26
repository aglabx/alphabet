use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use crate::io::read_fasta;

#[derive(Parser)]
#[command(name = "find_phase", about = "Find optimal monomer phase by entropy minimization + CENP-B box stability")]
struct Args {
    /// Input: .10kb.fasta
    input: String,

    /// Output: phase_score.tsv
    #[arg(short, long, default_value = "phase_score.tsv")]
    output: String,

    /// Target period
    #[arg(long, default_value_t = 171)]
    period: usize,

    /// Period tolerance for array selection
    #[arg(long, default_value_t = 5)]
    period_tol: usize,

    /// Min array length (bp)
    #[arg(long, default_value_t = 50000)]
    min_array_len: usize,

    /// Max arrays to process
    #[arg(long, default_value_t = 20)]
    max_arrays: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

const CENP_B_POSITIONS: [(usize, u8); 9] = [
    (1, b'T'), (2, b'T'), (3, b'C'), (4, b'G'),
    (9, b'A'), (12, b'C'), (13, b'G'), (14, b'G'), (15, b'G'),
];

#[derive(Debug, Serialize)]
struct PhaseScore {
    array_id: String,
    shift: usize,
    n_monomers: usize,
    mean_entropy: f64,
    low_entropy_positions: usize,  // positions with H < 0.8
    entropy_variance: f64,
    cenpb_count: usize,           // monomers with CENP-B score >= 7
    cenpb_mean_position: f64,
    cenpb_position_sd: f64,
    total_score: f64,
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let p = args.period;

    // Read arrays
    eprintln!("Reading {}...", args.input);
    let all_arrays = read_fasta(&args.input).unwrap_or_else(|e| panic!("{:?}", e));
    eprintln!("  {} total arrays", all_arrays.len());

    // Filter: period-matched, large enough
    let mut arrays: Vec<&(String, Vec<u8>)> = all_arrays.iter().filter(|(name, seq)| {
        let fields: Vec<&str> = name.split('_').collect();
        let period: usize = fields.last().and_then(|s| s.parse().ok()).unwrap_or(0);
        period.abs_diff(p) <= args.period_tol && seq.len() >= args.min_array_len
    }).collect();
    arrays.sort_by_key(|b| std::cmp::Reverse(b.1.len()));
    arrays.truncate(args.max_arrays);
    eprintln!("  {} arrays selected (period ~{}, len >= {})", arrays.len(), p, args.min_array_len);

    // For each array, scan all 171 phases
    let all_results: Vec<Vec<PhaseScore>> = arrays.par_iter().map(|(name, seq)| {
        let mut results = Vec::with_capacity(p);

        for shift in 0..p {
            // Slice into monomers at this phase
            let n_monomers = (seq.len().saturating_sub(shift)) / p;
            if n_monomers < 10 { continue; }

            let monomers: Vec<&[u8]> = (0..n_monomers)
                .map(|i| &seq[shift + i * p..shift + (i + 1) * p])
                .collect();

            // Metric A: column entropy
            let mut entropies = vec![0.0f64; p];
            for pos in 0..p {
                let mut counts = [0u32; 4];
                let mut valid = 0u32;
                for mono in &monomers {
                    let idx = match mono[pos] {
                        b'A' | b'a' => 0, b'C' | b'c' => 1,
                        b'G' | b'g' => 2, b'T' | b't' => 3, _ => continue,
                    };
                    counts[idx] += 1;
                    valid += 1;
                }
                if valid == 0 { entropies[pos] = 2.0; continue; }
                let mut h = 0.0f64;
                for &c in &counts {
                    if c > 0 {
                        let freq = c as f64 / valid as f64;
                        h -= freq * freq.log2();
                    }
                }
                entropies[pos] = h;
            }

            let mean_entropy = entropies.iter().sum::<f64>() / p as f64;

            // Metric B: contrast
            let low_entropy_count = entropies.iter().filter(|&&h| h < 0.8).count();
            let entropy_mean = mean_entropy;
            let entropy_var = entropies.iter()
                .map(|&h| (h - entropy_mean).powi(2))
                .sum::<f64>() / p as f64;

            // Metric C: CENP-B box position compactness
            let mut cenpb_positions: Vec<usize> = Vec::new();
            for mono in &monomers {
                let mut best_score = 0u32;
                let mut best_pos = 0usize;
                for i in 0..p.saturating_sub(16) {
                    let score: u32 = CENP_B_POSITIONS.iter()
                        .filter(|&&(off, base)| mono[i + off].to_ascii_uppercase() == base)
                        .count() as u32;
                    if score > best_score {
                        best_score = score;
                        best_pos = i;
                    }
                }
                if best_score >= 7 {
                    cenpb_positions.push(best_pos);
                }
            }

            let cenpb_count = cenpb_positions.len();
            let cenpb_mean = if cenpb_count > 0 {
                cenpb_positions.iter().sum::<usize>() as f64 / cenpb_count as f64
            } else { 0.0 };
            let cenpb_sd = if cenpb_count > 1 {
                let var = cenpb_positions.iter()
                    .map(|&p| (p as f64 - cenpb_mean).powi(2))
                    .sum::<f64>() / cenpb_count as f64;
                var.sqrt()
            } else { 171.0 }; // penalty if no box found

            // Total score: lower = better
            // Combine: low entropy (good), high contrast (good), low cenpb_sd (good)
            let score = mean_entropy
                - (low_entropy_count as f64 / p as f64) * 0.5  // bonus for conserved positions
                - entropy_var * 0.3                              // bonus for contrast
                + (cenpb_sd / p as f64) * 0.5;                  // penalty for box spread

            results.push(PhaseScore {
                array_id: name.clone(),
                shift,
                n_monomers,
                mean_entropy,
                low_entropy_positions: low_entropy_count,
                entropy_variance: entropy_var,
                cenpb_count,
                cenpb_mean_position: cenpb_mean,
                cenpb_position_sd: cenpb_sd,
                total_score: score,
            });
        }

        results
    }).collect();

    // Write all results
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.output).unwrap());
        writeln!(f, "array_id\tshift\tn_monomers\tmean_entropy\tlow_entropy_pos\tentropy_variance\tcenpb_count\tcenpb_mean_pos\tcenpb_pos_sd\ttotal_score").unwrap();
        for results in &all_results {
            for r in results {
                writeln!(f, "{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{:.1}\t{:.1}\t{:.4}",
                    r.array_id, r.shift, r.n_monomers, r.mean_entropy,
                    r.low_entropy_positions, r.entropy_variance,
                    r.cenpb_count, r.cenpb_mean_position, r.cenpb_position_sd,
                    r.total_score).unwrap();
            }
        }
    }

    // Per-array best phase
    eprintln!("\n=== BEST PHASE PER ARRAY ===");
    eprintln!("{:<55} {:>5} {:>6} {:>6} {:>6} {:>6} {:>8} {:>6}",
        "array", "shift", "H_mean", "lowH", "H_var", "cenB", "cenB_sd", "score");

    let mut all_best_shifts: Vec<usize> = Vec::new();
    for results in &all_results {
        if results.is_empty() { continue; }
        let best = results.iter().min_by(|a, b| a.total_score.partial_cmp(&b.total_score).unwrap()).unwrap();
        let short_name = &best.array_id[best.array_id.len().saturating_sub(55)..];
        eprintln!("{:<55} {:>5} {:>6.3} {:>6} {:>6.3} {:>6} {:>8.1} {:>6.3}",
            short_name, best.shift, best.mean_entropy, best.low_entropy_positions,
            best.entropy_variance, best.cenpb_count, best.cenpb_position_sd, best.total_score);
        all_best_shifts.push(best.shift);
    }

    // Family canonical phase = mode of best shifts
    if !all_best_shifts.is_empty() {
        let mut shift_counts: HashMap<usize, usize> = HashMap::new();
        for &s in &all_best_shifts {
            *shift_counts.entry(s).or_insert(0) += 1;
        }
        let mut sorted: Vec<(usize, usize)> = shift_counts.into_iter().collect();
        sorted.sort_by_key(|b| std::cmp::Reverse(b.1));

        eprintln!("\n=== FAMILY CANONICAL PHASE ===");
        eprintln!("Top shifts:");
        for (shift, count) in sorted.iter().take(10) {
            eprintln!("  shift={}: {} arrays", shift, count);
        }
        let canonical_shift = sorted[0].0;
        eprintln!("\nCanonical phase: shift={}", canonical_shift);

        // CENP-B box position at canonical phase
        eprintln!("\n=== CENP-B BOX AT CANONICAL PHASE ===");
        for results in &all_results {
            if let Some(r) = results.iter().find(|r| r.shift == canonical_shift) {
                if r.cenpb_count > 0 {
                    eprintln!("  {} : {} B+ monomers, box at position {:.0} ± {:.1}",
                        &r.array_id[r.array_id.len().saturating_sub(40)..],
                        r.cenpb_count, r.cenpb_mean_position, r.cenpb_position_sd);
                }
            }
        }
    }

    eprintln!("\nWritten to {}", args.output);
}

