use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use crate::monomer::revcomp;
use crate::io::read_fasta;

#[derive(Parser)]
#[command(name = "motif_graph", about = "Scan alpha arrays for anchor motifs, build transition graph")]
struct Args {
    /// Input: .10kb.fasta file
    input: String,

    /// Output: graph edges TSV
    #[arg(short, long, default_value = "motif_graph_edges.tsv")]
    output: String,

    /// Output: motif hits TSV (per-array motif positions)
    #[arg(long, default_value = "motif_hits.tsv")]
    hits: String,

    /// Output: graph summary JSON
    #[arg(long, default_value = "motif_graph.json")]
    report: String,

    /// Max hamming distance for fuzzy motif match
    #[arg(long, default_value_t = 2)]
    max_mismatch: u32,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

// The 5 anchor motifs found ab initio in Pass 2
const ANCHORS: [(&str, &str); 5] = [
    ("M1", "ACATCACAAAG"),  // pos 21-32, conservation 87.2%
    ("M2", "AGAATGCTTCT"),  // pos 41-52, conservation 92.6%
    ("M3", "GAAGATATTTC"),  // pos 67-78, conservation 87.4%
    ("M4", "TCCACTTGCAG"),  // pos 112-123, conservation 86.4%
    ("M5", "AAAGAGTGTTT"),  // pos 132-143, conservation 91.4%
];

#[derive(Debug, Clone)]
struct MotifHit {
    motif_idx: usize,  // 0-4
    position: usize,   // position in array
    mismatches: u32,
    strand: char,       // '+' or '-'
}

#[derive(Debug, Clone, Hash, Eq, PartialEq, Serialize)]
struct Edge {
    from_motif: String,
    to_motif: String,
    linker_length: i32,
}

#[derive(Debug, Serialize)]
struct EdgeStats {
    from_motif: String,
    to_motif: String,
    linker_length: i32,
    count: u64,
    fraction: f64,
}

#[derive(Debug, Serialize)]
struct GraphReport {
    n_arrays: usize,
    n_arrays_with_hits: usize,
    total_motif_hits: usize,
    hits_per_motif: Vec<(String, usize)>,
    total_edges: usize,
    unique_edges: usize,
    top_edges: Vec<EdgeStats>,
    linker_length_summary: LinkerSummary,
}

#[derive(Debug, Serialize)]
struct LinkerSummary {
    /// For each consecutive motif pair in expected order (M1->M2, M2->M3, etc)
    canonical_pairs: Vec<CanonicalPair>,
}

#[derive(Debug, Serialize)]
struct CanonicalPair {
    from: String,
    to: String,
    expected_distance: i32,
    median_observed: i32,
    count: u64,
    length_distribution: Vec<(i32, u64)>,
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    // Build motif sequences and their reverse complements
    let motifs: Vec<(String, Vec<u8>, Vec<u8>)> = ANCHORS.iter().map(|(name, seq)| {
        let fwd = seq.as_bytes().to_vec();
        let rc = revcomp(&fwd);
        (name.to_string(), fwd, rc)
    }).collect();

    // --- Step 1: Read FASTA, filter alpha arrays ---
    eprintln!("Reading {}...", args.input);
    let arrays = read_alpha_arrays(&args.input);
    eprintln!("  {} alpha satellite arrays (period 171/172)", arrays.len());

    // --- Step 2: Scan each array for motif hits ---
    eprintln!("Scanning for {} motifs (max {} mismatches)...", motifs.len(), args.max_mismatch);

    let array_hits: Vec<(String, usize, Vec<MotifHit>)> = arrays.par_iter().map(|(name, seq)| {
        let hits = scan_motifs(seq, &motifs, args.max_mismatch);
        (name.clone(), seq.len(), hits)
    }).collect();

    let total_hits: usize = array_hits.iter().map(|(_, _, h)| h.len()).sum();
    let arrays_with_hits = array_hits.iter().filter(|(_, _, h)| !h.is_empty()).count();
    eprintln!("  Total motif hits: {}", total_hits);
    eprintln!("  Arrays with hits: {}/{}", arrays_with_hits, arrays.len());

    // Count per motif
    let mut per_motif = [0usize; 5];
    for (_, _, hits) in &array_hits {
        for h in hits {
            per_motif[h.motif_idx] += 1;
        }
    }
    for (i, (name, _)) in ANCHORS.iter().enumerate() {
        eprintln!("    {}: {} hits", name, per_motif[i]);
    }

    // --- Step 3: Build edges from consecutive hits ---
    eprintln!("Building transition graph...");

    let _all_edges: Vec<Edge> = Vec::new();
    let mut edge_counts: HashMap<Edge, u64> = HashMap::new();

    for (_array_name, _array_len, hits) in &array_hits {
        if hits.len() < 2 { continue; }

        // Hits are already sorted by position from scan_motifs
        for i in 0..hits.len() - 1 {
            let from = &hits[i];
            let to = &hits[i + 1];

            let motif_len = ANCHORS[from.motif_idx].1.len();
            let linker_len = (to.position as i32) - (from.position as i32) - (motif_len as i32);

            let edge = Edge {
                from_motif: ANCHORS[from.motif_idx].0.to_string(),
                to_motif: ANCHORS[to.motif_idx].0.to_string(),
                linker_length: linker_len,
            };

            *edge_counts.entry(edge).or_insert(0) += 1;
        }
    }

    let total_edge_count: u64 = edge_counts.values().sum();
    eprintln!("  Unique edge types: {}", edge_counts.len());
    eprintln!("  Total edge instances: {}", total_edge_count);

    // --- Step 4: Sort and report ---
    let mut edge_stats: Vec<EdgeStats> = edge_counts.iter().map(|(e, &c)| {
        EdgeStats {
            from_motif: e.from_motif.clone(),
            to_motif: e.to_motif.clone(),
            linker_length: e.linker_length,
            count: c,
            fraction: c as f64 / total_edge_count as f64,
        }
    }).collect();
    edge_stats.sort_by(|a, b| b.count.cmp(&a.count));

    eprintln!("\n=== TOP 30 EDGES ===");
    eprintln!("{:<5} {:<5} {:>8} {:>8} {:>6}", "from", "to", "linker", "count", "%");
    for e in edge_stats.iter().take(30) {
        eprintln!("{:<5} {:<5} {:>8} {:>8} {:>5.1}%",
            e.from_motif, e.to_motif, e.linker_length, e.count, e.fraction * 100.0);
    }

    // Canonical pairs analysis (expected order: M1->M2->M3->M4->M5->M1)
    let canonical_pairs_names = [
        ("M1", "M2"), ("M2", "M3"), ("M3", "M4"), ("M4", "M5"), ("M5", "M1"),
    ];

    let mut canonical_pairs: Vec<CanonicalPair> = Vec::new();
    for &(from, to) in &canonical_pairs_names {
        let mut lengths: Vec<(i32, u64)> = edge_stats.iter()
            .filter(|e| e.from_motif == from && e.to_motif == to)
            .map(|e| (e.linker_length, e.count))
            .collect();
        lengths.sort_by_key(|&(l, _)| l);

        let total: u64 = lengths.iter().map(|(_, c)| c).sum();
        let median = if !lengths.is_empty() {
            let mut cum = 0u64;
            let half = total / 2;
            lengths.iter().find(|&&(_l, c)| { cum += c; cum >= half }).map(|&(l, _)| l).unwrap_or(0)
        } else { 0 };

        // Expected distance based on anchor positions in canonical monomer
        let from_pos = ANCHORS.iter().position(|(n, _)| *n == from).unwrap();
        let to_pos = ANCHORS.iter().position(|(n, _)| *n == to).unwrap();
        let from_start = [21, 41, 67, 112, 132][from_pos];
        let to_start = [21, 41, 67, 112, 132][to_pos];
        let expected = if to_start > from_start {
            to_start - from_start - 11 // motif length
        } else {
            171 - from_start + to_start - 11
        };

        canonical_pairs.push(CanonicalPair {
            from: from.to_string(),
            to: to.to_string(),
            expected_distance: expected,
            median_observed: median,
            count: total,
            length_distribution: lengths.iter().take(20).cloned().collect(),
        });
    }

    eprintln!("\n=== CANONICAL PAIRS (expected monomer order) ===");
    for cp in &canonical_pairs {
        eprintln!("  {}->{}:  expected={:>3}, observed_median={:>3}, count={}",
            cp.from, cp.to, cp.expected_distance, cp.median_observed, cp.count);
    }

    // --- Step 5: Write outputs ---
    eprintln!("\nWriting outputs...");

    // Motif hits TSV
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.hits).unwrap());
        writeln!(f, "array_id\tarray_len\tmotif\tposition\tmismatches\tstrand").unwrap();
        for (name, len, hits) in &array_hits {
            for h in hits {
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}",
                    name, len, ANCHORS[h.motif_idx].0, h.position, h.mismatches, h.strand).unwrap();
            }
        }
    }

    // Edges TSV
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.output).unwrap());
        writeln!(f, "from_motif\tto_motif\tlinker_length\tcount\tfraction").unwrap();
        for e in &edge_stats {
            writeln!(f, "{}\t{}\t{}\t{}\t{:.6}",
                e.from_motif, e.to_motif, e.linker_length, e.count, e.fraction).unwrap();
        }
    }

    // Report JSON
    let report = GraphReport {
        n_arrays: arrays.len(),
        n_arrays_with_hits: arrays_with_hits,
        total_motif_hits: total_hits,
        hits_per_motif: ANCHORS.iter().enumerate().map(|(i, (n, _))| (n.to_string(), per_motif[i])).collect(),
        total_edges: total_edge_count as usize,
        unique_edges: edge_counts.len(),
        top_edges: edge_stats.into_iter().take(100).collect(),
        linker_length_summary: LinkerSummary { canonical_pairs },
    };
    std::fs::write(&args.report, serde_json::to_string_pretty(&report).unwrap()).unwrap();

    eprintln!("Done.");
}

// ============ MOTIF SCANNING ============

fn scan_motifs(seq: &[u8], motifs: &[(String, Vec<u8>, Vec<u8>)], max_mismatch: u32) -> Vec<MotifHit> {
    let mut hits: Vec<MotifHit> = Vec::new();

    for (midx, (_, fwd, rc)) in motifs.iter().enumerate() {
        let mlen = fwd.len();
        if seq.len() < mlen { continue; }

        for i in 0..=(seq.len() - mlen) {
            let window = &seq[i..i + mlen];

            // Check forward
            let fwd_mm = hamming_bytes(window, fwd);
            if fwd_mm <= max_mismatch {
                hits.push(MotifHit {
                    motif_idx: midx,
                    position: i,
                    mismatches: fwd_mm,
                    strand: '+',
                });
                continue; // don't double-count
            }

            // Check reverse complement
            let rc_mm = hamming_bytes(window, rc);
            if rc_mm <= max_mismatch {
                hits.push(MotifHit {
                    motif_idx: midx,
                    position: i,
                    mismatches: rc_mm,
                    strand: '-',
                });
            }
        }
    }

    // Sort by position
    hits.sort_by_key(|h| h.position);

    // Deduplicate: if multiple motifs hit overlapping positions, keep best
    dedup_hits(&mut hits, motifs[0].1.len());

    hits
}

fn dedup_hits(hits: &mut Vec<MotifHit>, motif_len: usize) {
    if hits.len() < 2 { return; }

    let mut keep = vec![true; hits.len()];
    for i in 0..hits.len() {
        if !keep[i] { continue; }
        for j in (i + 1)..hits.len() {
            if !keep[j] { continue; }
            // If overlapping (same motif, nearby position)
            if hits[j].position < hits[i].position + motif_len {
                // Keep the one with fewer mismatches
                if hits[j].mismatches < hits[i].mismatches {
                    keep[i] = false;
                    break;
                } else {
                    keep[j] = false;
                }
            } else {
                break; // sorted, no more overlaps
            }
        }
    }

    let mut write = 0;
    for read in 0..hits.len() {
        if keep[read] {
            hits.swap(write, read);
            write += 1;
        }
    }
    hits.truncate(write);
}

fn hamming_bytes(a: &[u8], b: &[u8]) -> u32 {
    let mut d = 0u32;
    for i in 0..a.len().min(b.len()) {
        if !a[i].eq_ignore_ascii_case(&b[i]) { d += 1; }
    }
    d
}

// ============ FASTA READING ============

fn read_alpha_arrays(path: &str) -> Vec<(String, Vec<u8>)> {
    read_fasta(path)
        .unwrap_or_else(|e| panic!("{:?}", e))
        .into_iter()
        .filter(|(name, _)| {
            let fields: Vec<&str> = name.split('_').collect();
            fields
                .last()
                .and_then(|s| s.parse::<u32>().ok())
                .map(|p| p == 171 || p == 172)
                .unwrap_or(false)
        })
        .collect()
}
