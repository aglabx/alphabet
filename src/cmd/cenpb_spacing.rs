use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::thread;
use crate::monomer::revcomp;
use crate::io::read_fasta_strings;

pub fn run_from_args(args: Vec<String>) {
    if args.len() < 2 {
        eprintln!("Usage: cenpb_spacing <input.fasta> [min_score] [threads]");
        std::process::exit(1);
    }
    let input = &args[1];
    let min_score: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(7);
    let n_threads: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(32);

    // CENP-B box: nTTCGnnnnAnnCGGGn (17bp)
    // Fixed positions (0-indexed): 1=T, 2=T, 3=C, 4=G, 9=A, 12=C, 13=G, 14=G, 15=G
    let fixed: [(usize, u8); 9] = [
        (1, b'T'), (2, b'T'), (3, b'C'), (4, b'G'),
        (9, b'A'),
        (12, b'C'), (13, b'G'), (14, b'G'), (15, b'G'),
    ];

    let arrays = Arc::new(read_fasta_strings(input).unwrap_or_else(|e| panic!("{:?}", e)));
    let total_bp: usize = arrays.iter().map(|(_, s)| s.len()).sum();
    eprintln!("{} arrays, {:.1}Mb", arrays.len(), total_bp as f64 / 1e6);

    let spacings = Arc::new(Mutex::new(HashMap::<i64, usize>::new()));
    let total_hits = Arc::new(Mutex::new(0usize));
    let hits_fwd = Arc::new(Mutex::new(0usize));
    let hits_rc = Arc::new(Mutex::new(0usize));
    let n_arrays = arrays.len();

    let mut handles = Vec::new();
    for t in 0..n_threads {
        let arrays = Arc::clone(&arrays);
        let spacings = Arc::clone(&spacings);
        let total_hits = Arc::clone(&total_hits);
        let hits_fwd = Arc::clone(&hits_fwd);
        let hits_rc = Arc::clone(&hits_rc);
        handles.push(thread::spawn(move || {
            let mut local_spacings = HashMap::<i64, usize>::new();
            let mut local_hits = 0usize;
            let mut local_fwd = 0usize;
            let mut local_rc = 0usize;

            let mut idx = t;
            while idx < n_arrays {
                let seq = arrays[idx].1.as_bytes();
                let rc = revcomp(seq);
                let mut hits: Vec<(usize, char)> = Vec::new();

                for (strand_seq, strand) in [(seq, '+'), (&rc[..], '-')] {
                    if strand_seq.len() < 17 { continue; }
                    for i in 0..strand_seq.len() - 16 {
                        let mut score = 0u32;
                        for &(pos, base) in &fixed {
                            if strand_seq[i + pos].to_ascii_uppercase() == base {
                                score += 1;
                            }
                        }
                        if score >= min_score {
                            let real_pos = if strand == '+' { i } else { seq.len() - i - 17 };
                            hits.push((real_pos, strand));
                            local_hits += 1;
                            if strand == '+' { local_fwd += 1; } else { local_rc += 1; }
                        }
                    }
                }

                // Dedup close hits
                hits.sort_by_key(|h| h.0);
                let mut deduped: Vec<(usize, char)> = Vec::new();
                for h in &hits {
                    if deduped.last().is_none_or(|last| h.0.abs_diff(last.0) > 5) {
                        deduped.push(*h);
                    }
                }

                // Spacings
                for i in 0..deduped.len().saturating_sub(1) {
                    let d = deduped[i + 1].0 as i64 - deduped[i].0 as i64;
                    *local_spacings.entry(d).or_insert(0) += 1;
                }

                idx += n_threads;
            }

            // Merge
            {
                let mut g = spacings.lock().unwrap();
                for (k, v) in local_spacings { *g.entry(k).or_insert(0) += v; }
            }
            *total_hits.lock().unwrap() += local_hits;
            *hits_fwd.lock().unwrap() += local_fwd;
            *hits_rc.lock().unwrap() += local_rc;
        }));
    }
    for h in handles { h.join().unwrap(); }

    let spacings = Arc::try_unwrap(spacings).unwrap().into_inner().unwrap();
    let total_hits = *total_hits.lock().unwrap();
    let hits_fwd = *hits_fwd.lock().unwrap();
    let hits_rc = *hits_rc.lock().unwrap();

    let total_sp: usize = spacings.values().sum();
    let _period_guess = if total_bp > 0 { total_bp as f64 / (total_hits as f64 / 2.0) } else { 0.0 };

    println!("# CENP-B box spacing (score >= {})", min_score);
    println!("# {} arrays, {:.1}Mb", arrays.len(), total_bp as f64 / 1e6);
    println!("# Hits: {} ({:.2}/monomer at any period)", total_hits, total_hits as f64 / (total_bp as f64 / 120.0));
    println!("# Strand: fwd={} rc={}", hits_fwd, hits_rc);
    println!("# spacing\tcount\tpercent");

    let mut sorted_sp: Vec<_> = spacings.into_iter().filter(|&(d, _)| d > 0).collect();
    sorted_sp.sort_by_key(|&(d, _)| d);

    for (d, c) in &sorted_sp {
        let pct = *c as f64 * 100.0 / total_sp as f64;
        if pct >= 0.2 {
            let bar: String = "#".repeat((pct * 2.0).min(60.0) as usize);
            println!("{}\t{}\t{:.1}%\t{}", d, c, pct, bar);
        }
    }
}

