use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::Serialize;
use crate::monomer::kmer_hash;
use crate::io::read_fasta;

#[derive(Parser)]
#[command(name = "discover_chains", about = "Chain-first discovery: find periodic anchor chains in satellite arrays")]
struct Args {
    /// Input: .10kb.fasta file
    input: String,

    /// Output: chains JSON
    #[arg(short, long, default_value = "chains.json")]
    output: String,

    /// Target period
    #[arg(long, default_value_t = 0)]
    period: usize,

    /// Period tolerance
    #[arg(long, default_value_t = 30)]
    period_tol: usize,

    /// Starting k-mer size for anchor search
    #[arg(long, default_value_t = 8)]
    k_start: usize,

    /// Max EM iterations (discover → enrich → recalculate)
    #[arg(long, default_value_t = 3)]
    max_iter: usize,

    /// Max suprachromosomal families to discover
    #[arg(long, default_value_t = 5)]
    max_sf: usize,

    /// Min family support (fraction of arrays)
    #[arg(long, default_value_t = 0.5)]
    min_family_support: f64,

    /// Min positional concentration (fraction of hits in ±5bp of mode)
    #[arg(long, default_value_t = 0.5)]
    min_positional_concentration: f64,

    /// Max anchors to keep after filtering
    #[arg(long, default_value_t = 200)]
    max_anchors: usize,

    /// Max spacer variance for link (bp)
    #[arg(long, default_value_t = 5)]
    max_spacer_sd: usize,

    /// Min link support (fraction of arrays with both anchors at correct spacing)
    #[arg(long, default_value_t = 0.3)]
    min_link_support: f64,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Clone, Serialize)]
struct Site {
    id: usize,
    sequence: String,
    length: usize,
    // INVARIANT: canonical position within the monomer, always in [0, period).
    // All mutations (e.g. leftward growth in grow_sites) must go through
    // modular arithmetic (rem_euclid) to preserve this.
    position_mod_p: usize,
    family_support: f64,
    copy_support_median: f64,
    positional_concentration: f64,
    n_arrays: usize,
}

#[derive(Debug, Clone, Serialize)]
struct Link {
    from_site: usize,
    to_site: usize,
    spacer_mean: f64,
    spacer_sd: f64,
    support: f64,
    n_arrays: usize,
}

#[derive(Debug, Clone, Serialize)]
struct Chain {
    sites: Vec<usize>,     // site IDs in order
    links: Vec<usize>,     // link IDs in order
    period_coverage: f64,  // fraction of period covered by sites+spacers
    family_support: f64,   // fraction of arrays where this chain is found
    n_arrays: usize,
}

#[derive(Debug, Serialize)]
struct ChainReport {
    target_period: usize,
    n_arrays: usize,
    n_candidate_anchors: usize,
    n_sites: usize,
    n_links: usize,
    n_chains: usize,
    n_enriched_arrays: usize,
    n_enriched_families: usize,
    sites: Vec<Site>,
    links: Vec<Link>,
    chains: Vec<Chain>,
    /// Flat list of all motifs (backward-compat: motif_cut used to read this).
    /// New consumers should use `families` below, which carries per-SF
    /// motif bundles AND the array membership each SF has claimed.
    motifs: Vec<MotifCompat>,
    families: Vec<SfFamily>,
    provisional_families: Vec<ProvisionalFamily>,
}

#[derive(Debug, Clone, Serialize)]
struct SfFamily {
    /// Stable id: "P<period>SF<sf_index>"
    id: String,
    period: usize,
    sf_index: usize,
    motifs: Vec<MotifCompat>,
    /// Original array names (no _rc) claimed by this SF.
    arrays: Vec<String>,
    n_arrays: usize,
}

#[derive(Debug, Clone, Serialize)]
struct MotifCompat {
    name: String,
    sequence: String,
    position_in_monomer: usize,
    conservation: f64,
}

#[derive(Debug, Serialize)]
struct ProvisionalFamily {
    id: i32,
    canonical_order: String,
    n_arrays: usize,
    total_bp: usize,
    core_sites: Vec<bool>,
    core_links: Vec<bool>,
    n_core_sites: usize,
    n_core_links: usize,
    periods: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize)]
struct ArraySignature {
    name: String,
    length: usize,
    period: usize,
    chr: String,
    best_strand: char,
    n_sites_present: usize,
    n_chain_links_confirmed: usize,
    chain_compatible: bool,
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    // --- Read arrays ---
    eprintln!("Reading {}...", args.input);
    let all_arrays = read_fasta(&args.input).unwrap_or_else(|e| panic!("{:?}", e));
    eprintln!("  {} total arrays", all_arrays.len());

    // --- Detect periods ---
    // Count all periods, sort by frequency
    let mut period_counts: HashMap<usize, usize> = HashMap::new();
    for (name, _) in &all_arrays {
        let fields: Vec<&str> = name.split('_').collect();
        if let Some(p) = fields.last().and_then(|s| s.parse::<usize>().ok()) {
            *period_counts.entry(p).or_insert(0) += 1;
        }
    }
    let mut sorted_periods: Vec<(usize, usize)> = period_counts.into_iter().collect();
    // Tie-break by period ascending so the iteration order is deterministic across runs
    // (HashMap iteration is randomized; without the second key, equal-count periods would
    // shuffle and the per-period claim order would change run-to-run, propagating to
    // different chains.json output for low-seed-count families).
    sorted_periods.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    // If specific period requested, use only that. Otherwise, process all with ≥5 arrays.
    let periods_to_process: Vec<usize> = if args.period > 0 {
        vec![args.period]
    } else {
        sorted_periods.iter()
            .map(|(period, _)| *period)
            .collect()
    };

    eprintln!("  Periods to process: {:?}", periods_to_process);
    eprintln!("  All periods: {:?}", &sorted_periods[..sorted_periods.len().min(15)]);

    // Track which arrays are claimed globally
    let mut claimed_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut all_sf_motifs: Vec<(String, Vec<MotifCompat>)> = Vec::new(); // (label, motifs)
    let mut sf_families: Vec<SfFamily> = Vec::new();
    let mut global_sf_id = 0usize;

    let k = args.k_start;

    // ========== OUTER PERIOD LOOP ==========
    for &target_period in &periods_to_process {
    let p = target_period;

    // Seed arrays: exact period match, not yet claimed
    let period_arrays: Vec<(String, Vec<u8>)> = all_arrays.iter().filter(|(name, _)| {
        if claimed_names.contains(name) { return false; }
        let fields: Vec<&str> = name.split('_').collect();
        let period: usize = fields.last().and_then(|s| s.parse().ok()).unwrap_or(0);
        period == target_period
    }).cloned().collect();

    if period_arrays.is_empty() {
        eprintln!("\n  Period {}: no unclaimed arrays, skipping.", target_period);
        continue;
    }

    eprintln!("\n{}", "~".repeat(70));
    eprintln!("~~~ PERIOD {} : {} seed arrays ~~~", target_period, period_arrays.len());
    eprintln!("{}", "~".repeat(70));

    // ========== OUTER SF LOOP ==========
    for sf_id in 0..args.max_sf {

    // Seed arrays for this SF: period-matched arrays NOT yet claimed
    let sf_seed_arrays: Vec<(String, Vec<u8>)> = if sf_id == 0 {
        period_arrays.clone()
    } else {
        period_arrays.iter()
            .filter(|(name, _)| !claimed_names.contains(name))
            .cloned()
            .collect()
    };

    if sf_seed_arrays.is_empty() {
        eprintln!("\n  SF{}: no unclaimed seed arrays, stopping.", sf_id);
        break;
    }

    eprintln!("\n{}", "#".repeat(70));
    eprintln!("### P{}:SF{} (global SF{}) — {} seed arrays ###", target_period, sf_id, global_sf_id, sf_seed_arrays.len());
    eprintln!("{}", "#".repeat(70));

    // ========== ITERATIVE EM LOOP (per SF) ==========
    let mut working_arrays: Vec<(String, Vec<u8>)> = sf_seed_arrays;
    let mut prev_n_arrays = 0usize;
    let mut final_sites: Vec<Site> = Vec::new();
    let mut final_chains: Vec<Chain> = Vec::new();
    let mut final_n_compatible = 0usize;
    let mut final_enriched: Vec<ArraySignature> = Vec::new();

    for iteration in 0..args.max_iter {
        let n_arrays = working_arrays.len();
        eprintln!("\n{}", "=".repeat(60));
        eprintln!("=== ITERATION {} : {} arrays ===", iteration, n_arrays);
        eprintln!("{}", "=".repeat(60));

        if n_arrays == prev_n_arrays && iteration > 0 {
            eprintln!("  Converged — no new arrays added.");
            break;
        }
        prev_n_arrays = n_arrays;

    // === STEP 1: Scan candidate anchors ===
    eprintln!("\n=== Step 1: Scanning {}-mers across {} arrays ===", k, n_arrays);

    // For each k-mer: which arrays contain it, at what absolute positions
    let per_array_kmers: Vec<HashMap<u64, Vec<usize>>> = working_arrays.par_iter().map(|(_, seq)| {
        let mut kmer_positions: HashMap<u64, Vec<usize>> = HashMap::new();
        for i in 0..seq.len().saturating_sub(k - 1) {
            if let Some(hash) = kmer_hash(&seq[i..i + k], k) {
                kmer_positions.entry(hash).or_default().push(i);
            }
        }
        kmer_positions
    }).collect();

    // For each k-mer: compute inter-occurrence distances → check periodicity
    let mut global_kmers: HashMap<u64, AnchorCandidate> = HashMap::new();
    for arr_kmers in &per_array_kmers {
        for (&hash, positions) in arr_kmers {
            let entry = global_kmers.entry(hash).or_insert_with(|| AnchorCandidate {
                hash,
                n_arrays: 0,
                mod_p_positions: Vec::new(),
                total_hits: 0,
                positional_concentration: 0.0,
                mode_position: 0,
                inter_occurrence_distances: Vec::new(),
            });
            entry.n_arrays += 1;
            entry.total_hits += positions.len();

            // Compute inter-occurrence distances
            if positions.len() >= 2 {
                for i in 0..positions.len() - 1 {
                    let d = positions[i + 1] - positions[i];
                    entry.inter_occurrence_distances.push(d);
                }
            }
        }
    }

    eprintln!("  {} unique {}-mers found", global_kmers.len(), k);

    // Filter by family support
    let min_arrays = (n_arrays as f64 * args.min_family_support) as usize;
    let mut candidates: Vec<AnchorCandidate> = global_kmers.into_values()
        .filter(|c| c.n_arrays >= min_arrays)
        .collect();
    eprintln!("  {} with family support >= {:.0}% ({} arrays)",
        candidates.len(), args.min_family_support * 100.0, min_arrays);

    // Score by periodicity: what fraction of inter-occurrence distances are close to period P?
    let period_tol = p / 10; // ±10%
    for c in &mut candidates {
        if c.inter_occurrence_distances.is_empty() {
            c.positional_concentration = 0.0;
            continue;
        }
        let periodic_count = c.inter_occurrence_distances.iter().filter(|&&d| {
            // Check if d is close to any multiple of P (1x, 2x, 3x...)
            if d == 0 { return false; }
            let nearest_multiple = ((d + p / 2) / p) * p;
            d.abs_diff(nearest_multiple) <= period_tol
        }).count();
        c.positional_concentration = periodic_count as f64 / c.inter_occurrence_distances.len() as f64;

        // mode_position: use first occurrence in each array mod P, find mode
        c.mod_p_positions.clear();
        for arr_kmers in &per_array_kmers {
            if let Some(positions) = arr_kmers.get(&c.hash) {
                if let Some(&first) = positions.first() {
                    c.mod_p_positions.push(first % p);
                }
            }
        }
        if !c.mod_p_positions.is_empty() {
            // Find mode with ±5bp window
            let window = 5;
            let mut best_count = 0;
            let mut best_pos = 0;
            for center in 0..p {
                let count = c.mod_p_positions.iter().filter(|&&pos| {
                    let d = pos.abs_diff(center);
                    let d = d.min(p - d);
                    d <= window
                }).count();
                if count > best_count { best_count = count; best_pos = center; }
            }
            c.mode_position = best_pos;
        }
    }

    // Filter by periodicity
    candidates.retain(|c| c.positional_concentration >= args.min_positional_concentration);
    eprintln!("  {} with periodicity >= {:.0}%",
        candidates.len(), args.min_positional_concentration * 100.0);

    // Filter out low-complexity k-mers (max 2 of same base)
    candidates.retain(|c| {
        let mut base_counts = [0u32; 4];
        let mut h = c.hash;
        for _ in 0..k {
            base_counts[(h & 3) as usize] += 1;
            h >>= 2;
        }
        let max_base = *base_counts.iter().max().unwrap();
        max_base as usize <= k * 2 / 3 // no base > 2/3 of k-mer
    });
    eprintln!("  {} after low-complexity filter", candidates.len());

    // Sort by score (support * concentration) and take top
    candidates.sort_by(|a, b| b.score().partial_cmp(&a.score()).unwrap().then_with(|| a.hash.cmp(&b.hash)));
    candidates.truncate(args.max_anchors);
    eprintln!("  Kept top {} candidates", candidates.len());

    // Decode sequences for top candidates (need to scan arrays again)
    eprintln!("  Decoding sequences...");
    let candidate_hashes: Vec<u64> = candidates.iter().map(|c| c.hash).collect();
    let hash_to_seq = decode_kmer_sequences(&working_arrays, k, &candidate_hashes);

    // Build Site objects
    let sites: Vec<Site> = candidates.iter().enumerate().map(|(i, c)| {
        let seq = hash_to_seq.get(&c.hash).cloned().unwrap_or_else(|| format!("?k{}", c.hash));
        Site {
            id: i,
            sequence: seq,
            length: k,
            position_mod_p: c.mode_position,
            family_support: c.n_arrays as f64 / n_arrays as f64,
            copy_support_median: c.total_hits as f64 / c.n_arrays as f64,
            positional_concentration: c.positional_concentration,
            n_arrays: c.n_arrays,
        }
    }).collect();

    // Cluster overlapping sites (same mod-P position ± k/2, similar sequence)
    eprintln!("  Clustering overlapping sites...");
    let sites = cluster_sites(sites, k, p);
    eprintln!("  {} sites after clustering", sites.len());

    // Save pre-growth sites for link discovery (exact hash matching)
    let pre_growth_sites = sites.clone();

    // === STEP 1b: Grow site boundaries ===
    eprintln!("\n=== Step 1b: Growing site boundaries ===");
    let sites = grow_sites(sites, &working_arrays, n_arrays, 0.7, p);

    // === STEP 1c: Dedup and enforce spacer gaps ===
    eprintln!("\n=== Step 1c: Dedup + enforce spacers ===");
    let sites = enforce_spacers(sites, p);

    // === STEP 2: Discover links (using exact 8-mer positions from Step 1) ===
    eprintln!("\n=== Step 2: Discovering links between {} sites ===", pre_growth_sites.len());

    let site_hashes: Vec<u64> = candidates.iter().take(pre_growth_sites.len()).map(|c| c.hash).collect();

    let links = discover_links(
        &working_arrays, &pre_growth_sites, &site_hashes, k, p, n_arrays,
        args.max_spacer_sd, args.min_link_support,
    );
    eprintln!("  {} links discovered", links.len());

    if !links.is_empty() {
        eprintln!("\n  Top 20 links:");
        for l in links.iter().take(20) {
            eprintln!("    S{} -> S{}: spacer={:.0}±{:.1}, support={:.1}%, arrays={}",
                l.from_site, l.to_site, l.spacer_mean, l.spacer_sd,
                l.support * 100.0, l.n_arrays);
        }
    }

    // === STEP 3: Assemble chains ===
    eprintln!("\n=== Step 3: Assembling chains ===");

    let chains = assemble_chains(&sites, &links, p);

    if chains.is_empty() {
        eprintln!("  No chains found. Falling back to top sites as independent anchors.");
    } else {
        for (i, chain) in chains.iter().enumerate() {
            let site_names: Vec<String> = chain.sites.iter()
                .map(|&si| format!("S{}({})", si, &sites[si].sequence[..sites[si].sequence.len().min(8)]))
                .collect();
            eprintln!("  Chain {}: {} sites, coverage={:.0}%, support={:.0}%",
                i, chain.sites.len(), chain.period_coverage * 100.0, chain.family_support * 100.0);
            eprintln!("    {}", site_names.join(" → "));
        }
    }

    // === STEP 4: Enrichment — chain grammar compatibility ===
    eprintln!("\n=== Step 4: Enrichment — chain grammar check on ALL {} arrays ===", all_arrays.len());

    let site_seqs: Vec<Vec<u8>> = sites.iter().map(|s| s.sequence.as_bytes().to_vec()).collect();
    let site_rcs: Vec<Vec<u8>> = site_seqs.iter().map(|s| s.iter().rev().map(|&b| match b {
        b'A'|b'a'=>b'T', b'T'|b't'=>b'A', b'C'|b'c'=>b'G', b'G'|b'g'=>b'C', o=>o
    }).collect()).collect();
    let n_sites = sites.len();

    // Build expected chain order and distances from the chain
    let chain_order: Vec<usize> = if !chains.is_empty() {
        chains[0].sites.clone()
    } else {
        (0..n_sites).collect()
    };
    // Expected distances between consecutive chain sites (from links)
    let expected_distances: Vec<Option<f64>> = (0..chain_order.len()).map(|i| {
        let from = chain_order[i];
        let to = chain_order[(i + 1) % chain_order.len()];
        links.iter().find(|l| l.from_site == from && l.to_site == to)
            .map(|l| l.spacer_mean)
    }).collect();
    let distance_tolerance = 0.3; // ±30% of expected distance

    let enriched: Vec<ArraySignature> = all_arrays.par_iter().map(|(name, seq)| {
        // Find all site positions on both strands
        let mut fwd_positions: Vec<Vec<usize>> = vec![Vec::new(); n_sites];
        let mut rc_positions: Vec<Vec<usize>> = vec![Vec::new(); n_sites];

        for (si, (fwd, rc)) in site_seqs.iter().zip(site_rcs.iter()).enumerate() {
            let slen = fwd.len();
            if seq.len() < slen { continue; }
            // Exact match only — fuzzy is wrong, grammar check handles specificity
            for i in 0..=(seq.len() - slen) {
                let w = &seq[i..i + slen];
                if w == fwd.as_slice()
                    && fwd_positions[si].last().is_none_or(|&last| i - last >= slen) {
                        fwd_positions[si].push(i);
                    }
                if w == rc.as_slice()
                    && rc_positions[si].last().is_none_or(|&last| i - last >= slen) {
                        rc_positions[si].push(i);
                    }
            }
        }

        let fwd_total: usize = fwd_positions.iter().map(|v| v.len()).sum();
        let rc_total: usize = rc_positions.iter().map(|v| v.len()).sum();
        let (positions, strand) = if fwd_total >= rc_total {
            (fwd_positions, '+')
        } else {
            (rc_positions, '-')
        };

        // Site is "present" only if it has sufficient occupancy (≥50% of expected monomers)
        let expected_monomers = seq.len().checked_div(p).unwrap_or(1).max(1);
        let min_hits = (expected_monomers as f64 * 0.3).max(3.0) as usize; // ≥30% occupancy, min 3 hits
        let site_present: Vec<bool> = positions.iter().map(|v| v.len() >= min_hits).collect();
        let n_sites_present = site_present.iter().filter(|&&b| b).count();
        // Zero out positions for non-present sites (to avoid false grammar matches)
        let positions: Vec<Vec<usize>> = positions.into_iter().zip(site_present.iter())
            .map(|(p, &present)| if present { p } else { Vec::new() })
            .collect();

        // Chain grammar check:
        // For each pair of found sites (A, B) where A < B in chain order:
        //   expected distance = sum of link distances from A to B through chain
        //   check if observed distance ≈ expected (± tolerance)
        // Gaps (missing sites) are OK — sum through them.

        // Precompute cumulative distance along chain from site 0
        // cum_dist[i] = distance from chain_order[0] to chain_order[i]
        let mut cum_dist: Vec<f64> = vec![0.0; chain_order.len()];
        for ci in 1..chain_order.len() {
            let prev_si = chain_order[ci - 1];
            let prev_len = site_seqs[prev_si].len() as f64;
            let link_d = expected_distances[ci - 1].unwrap_or(20.0); // default spacer if no link
            cum_dist[ci] = cum_dist[ci - 1] + prev_len + link_d;
        }

        // For each pair of found sites, check distance
        let mut n_pairs_ok = 0usize;

        // Collect (chain_index, all_positions) for found sites
        let found_chain_indices: Vec<usize> = chain_order.iter().enumerate()
            .filter(|(_, &si)| !positions[si].is_empty())
            .map(|(ci, _)| ci)
            .collect();

        for wi in 0..found_chain_indices.len() {
            for wj in (wi + 1)..found_chain_indices.len() {
                let ci_a = found_chain_indices[wi];
                let ci_b = found_chain_indices[wj];
                let si_a = chain_order[ci_a];
                let si_b = chain_order[ci_b];

                let _expected_d = cum_dist[ci_b] - cum_dist[ci_a] - site_seqs[si_a].len() as f64;
                // But we want spacer distance (not including site_a length)
                // Actually cum_dist already includes site lengths, so:
                let expected_total = cum_dist[ci_b] - cum_dist[ci_a]; // includes sites + spacers from A to B

                let tol = (expected_total * distance_tolerance).max(20.0);

                // Check if any occurrence pair matches
                let mut pair_ok = false;
                for &pa in &positions[si_a] {
                    for &pb in &positions[si_b] {
                        if pb > pa {
                            let obs = (pb - pa) as f64;
                            if (obs - expected_total).abs() <= tol {
                                pair_ok = true;
                                break;
                            }
                        }
                    }
                    if pair_ok { break; }
                }

                if pair_ok { n_pairs_ok += 1; }
            }
        }

        let n_links_confirmed = n_pairs_ok;

        // Chain compatible if: ≥3 sites found AND ≥1 pair at correct distance
        let chain_compatible = n_sites_present >= 3 && n_pairs_ok >= 1;

        let fields: Vec<&str> = name.split('_').collect();
        let period_ann: usize = fields.last().and_then(|s| s.parse().ok()).unwrap_or(0);
        let chr = if fields.len() >= 5 { fields[..fields.len()-4].join("_") } else { name.clone() };

        ArraySignature {
            name: name.clone(),
            length: seq.len(),
            period: period_ann,
            chr,
            best_strand: strand,
            n_sites_present,
            n_chain_links_confirmed: n_links_confirmed,
            chain_compatible,
        }
    }).collect();

    let n_compatible = enriched.iter().filter(|a| a.chain_compatible).count();
    let n_site_only = enriched.iter().filter(|a| a.n_sites_present >= 3 && !a.chain_compatible).count();
    eprintln!("  Chain-compatible arrays: {} / {}", n_compatible, all_arrays.len());
    eprintln!("  Rejected (sites present but wrong order/distances): {}", n_site_only);

    // Period distribution of compatible arrays
    let mut pcounts: HashMap<usize, usize> = HashMap::new();
    for a in enriched.iter().filter(|a| a.chain_compatible) {
        *pcounts.entry(a.period).or_insert(0) += 1;
    }
    let mut periods: Vec<(usize, usize)> = pcounts.into_iter().collect();
    periods.sort_by_key(|b| std::cmp::Reverse(b.1));
    eprintln!("  Periods: {}", periods.iter().take(8)
        .map(|(p, c)| format!("{}x{}", p, c)).collect::<Vec<_>>().join(", "));

    // Save results from this iteration
    let _n_candidates = candidates.len();
    let _ = links;
    final_sites = sites;
    final_chains = chains;
    final_n_compatible = n_compatible;

    // Update working_arrays: add chain-compatible arrays from ALL arrays
    let compatible_names: std::collections::HashSet<String> = enriched.iter()
        .filter(|a| a.chain_compatible)
        .map(|a| a.name.clone())
        .collect();
    let current_names: std::collections::HashSet<String> = working_arrays.iter()
        .map(|(n, _)| n.clone())
        .collect();

    let mut new_arrays: Vec<(String, Vec<u8>)> = Vec::new();
    for (name, seq) in &all_arrays {
        if compatible_names.contains(name) && !current_names.contains(name) {
            new_arrays.push((name.clone(), seq.clone()));
        }
    }

    if new_arrays.is_empty() {
        eprintln!("  No new arrays to add. Converged.");
        final_enriched = enriched;
        break;
    }

    eprintln!("  Adding {} new arrays for next iteration", new_arrays.len());
    working_arrays.extend(new_arrays);
    final_enriched = enriched;

    } // end EM loop

    // Claim compatible arrays for this SF
    let compatible_in_enrichment: std::collections::HashSet<String> = final_enriched.iter()
        .filter(|a| a.chain_compatible)
        .map(|a| a.name.clone())
        .collect();

    let newly_claimed: Vec<String> = compatible_in_enrichment.iter()
        .filter(|n| !claimed_names.contains(*n))
        .cloned()
        .collect();
    let n_new_claimed = newly_claimed.len();
    claimed_names.extend(compatible_in_enrichment);

    eprintln!("\n  SF{}: {} compatible arrays, {} newly claimed (total claimed: {})",
        sf_id, final_n_compatible, n_new_claimed, claimed_names.len());

    // Build motifs for this SF
    let sf_motifs: Vec<MotifCompat> = if !final_chains.is_empty() {
        final_chains[0].sites.iter().enumerate().map(|(i, &si)| {
            MotifCompat {
                name: format!("P{}SF{}M{}", target_period, sf_id, i + 1),
                sequence: final_sites[si].sequence.clone(),
                position_in_monomer: final_sites[si].position_mod_p,
                conservation: final_sites[si].family_support,
            }
        }).collect()
    } else {
        final_sites.iter().take(5).enumerate().map(|(i, s)| {
            MotifCompat {
                name: format!("P{}SF{}M{}", target_period, sf_id, i + 1),
                sequence: s.sequence.clone(),
                position_in_monomer: s.position_mod_p,
                conservation: s.family_support,
            }
        }).collect()
    };

    eprintln!("  SF{} motifs:", sf_id);
    for m in &sf_motifs {
        eprintln!("    {}: {} ({}bp)", m.name, m.sequence, m.sequence.len());
    }

    let family_id = format!("P{}SF{}", target_period, sf_id);
    all_sf_motifs.push((family_id.clone(), sf_motifs.clone()));

    // Record membership: only the arrays newly claimed by THIS SF.
    // Strip any "_rc" suffix so the name matches the original FASTA.
    let mut arrays_for_family: Vec<String> = newly_claimed
        .into_iter()
        .map(|n| n.strip_suffix("_rc").map(|s| s.to_string()).unwrap_or(n))
        .collect();
    arrays_for_family.sort();
    arrays_for_family.dedup();
    let n_family_arrays = arrays_for_family.len();
    sf_families.push(SfFamily {
        id: family_id,
        period: target_period,
        sf_index: sf_id,
        motifs: sf_motifs,
        arrays: arrays_for_family,
        n_arrays: n_family_arrays,
    });
    global_sf_id += 1;

    if n_new_claimed == 0 {
        eprintln!("  No new arrays claimed. Stopping SF discovery.");
        break;
    }

    } // end SF loop

    } // end PERIOD loop

    eprintln!("\n{}", "=".repeat(70));
    eprintln!("=== SUMMARY: {} families across {} periods, {} total arrays claimed / {} ===",
        all_sf_motifs.len(), periods_to_process.len(), claimed_names.len(), all_arrays.len());
    for (label, motifs) in &all_sf_motifs {
        eprintln!("  {}: {} motifs, first={}", label, motifs.len(),
            motifs.first().map(|m| m.sequence.as_str()).unwrap_or("?"));
    }

    // Output: use SF0 (dominant) as primary motif set
    let n_arrays = all_arrays.len();
    let n_enriched_families = all_sf_motifs.len();
    let n_compatible = claimed_names.len();
    let families: Vec<ProvisionalFamily> = Vec::new();

    // Collect all motifs from all SFs for output
    let motifs: Vec<MotifCompat> = all_sf_motifs.iter()
        .flat_map(|(_, m)| m.clone())
        .collect();

    eprintln!("\n=== MOTIFS (for motif_cut compatibility) ===");
    for m in &motifs {
        eprintln!("  {}: pos={}, seq={}, support={:.2}", m.name, m.position_in_monomer, m.sequence, m.conservation);
    }

    let report = ChainReport {
        target_period: args.period,
        n_arrays,
        n_candidate_anchors: 0,
        n_sites: motifs.len(),
        n_links: 0,
        n_chains: all_sf_motifs.len(),
        n_enriched_arrays: n_compatible,
        n_enriched_families,
        sites: Vec::new(),
        links: Vec::new(),
        chains: Vec::new(),
        motifs,
        families: sf_families,
        provisional_families: families,
    };

    std::fs::write(&args.output, serde_json::to_string_pretty(&report).unwrap()).unwrap();
    eprintln!("\nWritten to {}", args.output);
}

// ============ ANCHOR CANDIDATE ============

struct AnchorCandidate {
    hash: u64,
    n_arrays: usize,
    mod_p_positions: Vec<usize>,
    total_hits: usize,
    positional_concentration: f64,
    mode_position: usize,
    inter_occurrence_distances: Vec<usize>,
}

// Need default for the extra fields
impl AnchorCandidate {
    fn score(&self) -> f64 {
        (self.n_arrays as f64) * self.positional_concentration
    }
}

// ============ SITE GROWTH ============

fn grow_sites(sites: Vec<Site>, arrays: &[(String, Vec<u8>)], _n_arrays: usize, _min_support_frac: f64, period: usize) -> Vec<Site> {
    let max_extend = 30; // generous: entropy will tell us the real boundary

    sites.into_iter().map(|site| {
        let seed = site.sequence.as_bytes();
        let seed_len = seed.len();
        let context_len = seed_len + 2 * max_extend;

        // Collect flanking contexts from all occurrences (exact match)
        let contexts: Vec<Vec<u8>> = arrays.par_iter().flat_map(|(_, seq)| {
            let mut found: Vec<Vec<u8>> = Vec::new();
            for i in 0..seq.len().saturating_sub(seed_len) {
                if &seq[i..i + seed_len] == seed {
                    let ctx_start = i.saturating_sub(max_extend);
                    let ctx_end = (i + seed_len + max_extend).min(seq.len());
                    let offset = i - ctx_start;
                    let mut ctx = vec![b'N'; context_len];
                    let dst_start = max_extend - offset;
                    let src_len = ctx_end - ctx_start;
                    let copy_len = src_len.min(context_len.saturating_sub(dst_start));
                    ctx[dst_start..dst_start + copy_len]
                        .copy_from_slice(&seq[ctx_start..ctx_start + copy_len]);
                    found.push(ctx);
                }
            }
            found
        }).collect();

        if contexts.len() < 50 {
            eprintln!("  {} ({} bp) — {} contexts, keeping as-is", site.sequence, site.length, contexts.len());
            return site;
        }

        // Compute Shannon entropy per column across the full context window
        let entropy: Vec<f64> = (0..context_len).map(|pos| {
            let mut counts = [0u32; 4];
            let mut valid = 0u32;
            for ctx in &contexts {
                let idx = match ctx[pos] {
                    b'A' | b'a' => 0, b'C' | b'c' => 1,
                    b'G' | b'g' => 2, b'T' | b't' => 3, _ => continue,
                };
                counts[idx] += 1;
                valid += 1;
            }
            if valid == 0 { return 2.0; }
            let mut h = 0.0f64;
            for &c in &counts {
                if c > 0 {
                    let p = c as f64 / valid as f64;
                    h -= p * p.log2();
                }
            }
            h
        }).collect();

        // Seed occupies positions max_extend .. max_extend+seed_len in the context
        // Grow RIGHT: walk from max_extend+seed_len onward, stop at entropy jump
        let entropy_threshold = 1.0; // ~50% single-base = boundary

        let mut right_ext = 0;
        for ext in 1..=max_extend {
            let pos = max_extend + seed_len + ext - 1;
            if pos >= context_len { break; }
            if entropy[pos] > entropy_threshold { break; }
            right_ext = ext;
        }

        // Grow LEFT: walk from max_extend-1 backward
        let mut left_ext = 0;
        for ext in 1..=max_extend {
            let pos = max_extend - ext;
            if entropy[pos] > entropy_threshold { break; }
            left_ext = ext;
        }

        // Cap total site length at 11bp — if entropy allows more, trim to most conserved 11bp center
        let max_site_len = 11;
        let total_ext = left_ext + seed_len + right_ext;
        if total_ext > max_site_len {
            // Trim: keep the most conserved window of max_site_len
            let full_start = max_extend - left_ext;
            let full_end = max_extend + seed_len + right_ext;
            // Find lowest-entropy window of max_site_len within the grown region
            let mut best_start = full_start;
            let mut best_entropy = f64::MAX;
            for ws in full_start..=(full_end - max_site_len) {
                let h: f64 = entropy[ws..ws + max_site_len].iter().sum();
                if h < best_entropy {
                    best_entropy = h;
                    best_start = ws;
                }
            }
            left_ext = max_extend - best_start;
            right_ext = (best_start + max_site_len) - (max_extend + seed_len);
            // Clamp to avoid underflow
            if best_start > max_extend { left_ext = 0; }
            if best_start + max_site_len < max_extend + seed_len { right_ext = 0; }
        }

        // Build consensus for extended region
        let new_start = max_extend - left_ext;
        let new_end = max_extend + seed_len + right_ext;
        let new_len = new_end - new_start;

        let bases = [b'A', b'C', b'G', b'T'];
        let new_seq: Vec<u8> = (new_start..new_end).map(|pos| {
            let mut counts = [0u32; 4];
            for ctx in &contexts {
                if pos < ctx.len() {
                    let idx = match ctx[pos] {
                        b'A' | b'a' => 0, b'C' | b'c' => 1,
                        b'G' | b'g' => 2, b'T' | b't' => 3, _ => continue,
                    };
                    counts[idx] += 1;
                }
            }
            bases[counts.iter().enumerate().max_by_key(|(_, &c)| c).unwrap().0]
        }).collect();

        let new_sequence = String::from_utf8(new_seq).unwrap();
        // Growing left means the site's canonical start shifts backwards by left_ext.
        // Do this modulo period so the invariant (position_mod_p ∈ [0, period)) holds
        // even when the seed was near position 0.
        let p_i = period as i64;
        let new_position = ((site.position_mod_p as i64 - left_ext as i64).rem_euclid(p_i)) as usize;
        debug_assert!(new_position < period,
            "grow_sites produced position {} >= period {}", new_position, period);

        // Report with entropy at boundaries (positions just outside the grown region)
        let h_left = if left_ext > 0 && left_ext < max_extend {
            entropy[max_extend - left_ext - 1]
        } else { 2.0 };
        let h_right = {
            let p = max_extend + seed_len + right_ext;
            if right_ext > 0 && p < context_len { entropy[p] } else { 2.0 }
        };
        let h_core = entropy[max_extend..max_extend + seed_len].iter().sum::<f64>() / seed_len as f64;

        eprintln!("  {} → {} ({} → {} bp, L+{} R+{}) H_core={:.2} H_left_boundary={:.2} H_right_boundary={:.2}",
            site.sequence, new_sequence, site.length, new_len,
            left_ext, right_ext, h_core, h_left, h_right);

        Site {
            sequence: new_sequence,
            length: new_len,
            position_mod_p: new_position,
            ..site
        }
    }).collect()
}

fn enforce_spacers(mut sites: Vec<Site>, _period: usize) -> Vec<Site> {
    // Sort by position
    sites.sort_by(|a, b| a.position_mod_p.cmp(&b.position_mod_p).then_with(|| a.sequence.cmp(&b.sequence)));

    // 1. Remove exact duplicates (same sequence) — keep highest support
    let _before = sites.len();
    let mut seen_seqs: HashMap<String, usize> = HashMap::new(); // seq -> index of best
    let mut keep_dedup = vec![true; sites.len()];
    for (i, s) in sites.iter().enumerate() {
        if let Some(&prev_i) = seen_seqs.get(&s.sequence) {
            // Duplicate: keep higher support
            if s.family_support > sites[prev_i].family_support {
                keep_dedup[prev_i] = false;
                seen_seqs.insert(s.sequence.clone(), i);
            } else {
                keep_dedup[i] = false;
            }
        } else {
            seen_seqs.insert(s.sequence.clone(), i);
        }
    }
    let removed_exact = keep_dedup.iter().filter(|&&k| !k).count();
    sites = sites.into_iter().zip(keep_dedup).filter(|(_, k)| *k).map(|(s, _)| s).collect();
    if removed_exact > 0 {
        eprintln!("  Removed {} exact duplicates", removed_exact);
    }
    let after_dedup = sites.len();

    // 2. Remove sites whose sequence is a substring of another
    let seqs: Vec<String> = sites.iter().map(|s| s.sequence.clone()).collect();
    let mut keep = vec![true; sites.len()];
    for i in 0..sites.len() {
        for j in 0..sites.len() {
            if i == j || !keep[j] { continue; }
            if seqs[j].len() > seqs[i].len() && seqs[j].contains(&seqs[i]) {
                keep[i] = false; // i is substring of j, remove i
                break;
            }
        }
    }
    let removed_substr = keep.iter().filter(|&&k| !k).count();
    if removed_substr > 0 {
        eprintln!("  Removed {} substring sites", removed_substr);
    }
    let mut filtered: Vec<Site> = sites.into_iter().zip(keep).filter(|(_, k)| *k).map(|(s, _)| s).collect();

    // 3. Merge overlapping sites (position + length overlaps with next)
    filtered.sort_by(|a, b| a.position_mod_p.cmp(&b.position_mod_p).then_with(|| a.sequence.cmp(&b.sequence)));
    let mut merged: Vec<Site> = Vec::new();
    for site in filtered {
        if let Some(prev) = merged.last() {
            let prev_end = prev.position_mod_p + prev.length;
            if site.position_mod_p < prev_end {
                // Overlap: keep the one with higher support
                if site.family_support > prev.family_support {
                    *merged.last_mut().unwrap() = site;
                }
                continue;
            }
        }
        merged.push(site);
    }

    let after_merge = merged.len();
    if after_merge != after_dedup - removed_substr {
        eprintln!("  Merged {} overlapping sites", (after_dedup - removed_substr) - after_merge);
    }

    // Re-index
    for (i, s) in merged.iter_mut().enumerate() {
        s.id = i;
    }

    eprintln!("  {} sites after cleanup", merged.len());
    merged
}

// ============ LINK DISCOVERY ============

fn discover_links(
    arrays: &[(String, Vec<u8>)],
    sites: &[Site],
    site_hashes: &[u64],
    k: usize,
    p: usize,
    n_arrays: usize,
    max_sd: usize,
    min_support: f64,
) -> Vec<Link> {
    let n_sites = sites.len();
    if n_sites < 2 { return vec![]; }

    // For each array, find positions of each site using EXACT k-mer hash matching (fast)
    let array_site_positions: Vec<Vec<Vec<usize>>> = arrays.par_iter().map(|(_, seq)| {
        let mut positions = vec![Vec::new(); n_sites];
        for i in 0..seq.len().saturating_sub(k - 1) {
            if let Some(hash) = kmer_hash(&seq[i..i + k], k) {
                for (si, &sh) in site_hashes.iter().enumerate() {
                    if si >= n_sites { break; }
                    if hash == sh {
                        positions[si].push(i);
                    }
                }
            }
        }
        positions
    }).collect();

    // For each pair of sites, compute spacer distribution
    // Outer pair loop flattened to a single 0..n_sites*n_sites range and parallelized
    // via rayon. Each (si, sj) is independent: only reads from immutable refs, builds
    // a local spacers vec, returns Option<Link>. Final sort below makes order irrelevant.
    let mut links: Vec<Link> = (0..n_sites * n_sites)
        .into_par_iter()
        .filter_map(|idx| {
            let si = idx / n_sites;
            let sj = idx % n_sites;
            if si == sj { return None; }

            let expected_spacer = {
                let pi = sites[si].position_mod_p;
                let pj = sites[sj].position_mod_p;
                if pj > pi { pj - pi - k } else { p - pi + pj - k }
            };

            if expected_spacer > p { return None; }

            let mut spacers: Vec<i32> = Vec::new();
            let mut arrays_with_link = 0;

            for arr_positions in &array_site_positions {
                let pos_i = &arr_positions[si];
                let pos_j = &arr_positions[sj];

                if pos_i.is_empty() || pos_j.is_empty() { continue; }

                let mut found_in_array = false;
                // pos_j is sorted ascending by construction (Phase 1 fills it left-to-right).
                // For each pi, only pj in window [pi + k + expected_spacer ± max_sd*3] can
                // produce d = pj - pi - k with |d - expected_spacer| <= max_sd*3. Binary
                // search via partition_point reduces inner cost from O(K_i * K_j) to
                // O(K_i * log K_j + matches). Lower bound clamped to pi+1 to preserve the
                // pj > pi requirement (and avoid usize underflow when computing d).
                let half_window = max_sd * 3;
                for &pi in pos_i {
                    let target = pi + k + expected_spacer;
                    let lo = std::cmp::max(target.saturating_sub(half_window), pi + 1);
                    let hi = target.saturating_add(half_window);
                    let i_lo = pos_j.partition_point(|&x| x < lo);
                    let i_hi = pos_j.partition_point(|&x| x <= hi);
                    for &pj in &pos_j[i_lo..i_hi] {
                        let d = (pj - pi) as i32 - k as i32;
                        spacers.push(d);
                        found_in_array = true;
                    }
                }
                if found_in_array { arrays_with_link += 1; }
            }

            let support = arrays_with_link as f64 / n_arrays as f64;
            if support < min_support || spacers.len() < 5 { return None; }

            let mean = spacers.iter().sum::<i32>() as f64 / spacers.len() as f64;
            let variance = spacers.iter().map(|&s| (s as f64 - mean).powi(2)).sum::<f64>() / spacers.len() as f64;
            let sd = variance.sqrt();

            if sd <= max_sd as f64 {
                Some(Link {
                    from_site: si,
                    to_site: sj,
                    spacer_mean: mean,
                    spacer_sd: sd,
                    support,
                    n_arrays: arrays_with_link,
                })
            } else {
                None
            }
        })
        .collect();

    links.sort_by(|a, b| b.support.partial_cmp(&a.support).unwrap().then_with(|| (a.from_site, a.to_site).cmp(&(b.from_site, b.to_site))));
    links
}

// ============ CHAIN ASSEMBLY ============

fn assemble_chains(sites: &[Site], links: &[Link], period: usize) -> Vec<Chain> {
    if sites.is_empty() { return vec![]; }

    // Simple approach: sort all sites by position_mod_p → this IS the chain order.
    // Then find which consecutive pairs have supporting links.
    let mut ordered: Vec<usize> = (0..sites.len()).collect();
    ordered.sort_by(|&i, &j| sites[i].position_mod_p.cmp(&sites[j].position_mod_p).then_with(|| sites[i].sequence.cmp(&sites[j].sequence)).then_with(|| i.cmp(&j)));

    // Build link lookup
    let mut link_map: HashMap<(usize, usize), usize> = HashMap::new();
    for (li, link) in links.iter().enumerate() {
        link_map.insert((link.from_site, link.to_site), li);
    }

    // Build the full chain: all sites in positional order
    let mut chain_links: Vec<usize> = Vec::new();
    let mut linked_count = 0;
    for i in 0..ordered.len() {
        let from = ordered[i];
        let to = ordered[(i + 1) % ordered.len()];
        if let Some(&li) = link_map.get(&(from, to)) {
            chain_links.push(li);
            linked_count += 1;
        } else {
            chain_links.push(usize::MAX); // no link found
        }
    }

    let total_site_bp: usize = sites.iter().map(|s| s.length).sum();
    let coverage = total_site_bp as f64 / period as f64;

    let min_support = ordered.iter()
        .map(|&si| sites[si].family_support)
        .fold(f64::MAX, f64::min);

    eprintln!("  Full chain: {} sites in positional order, {}/{} consecutive links found",
        ordered.len(), linked_count, ordered.len());

    let full_chain = Chain {
        sites: ordered.clone(),
        links: chain_links,
        period_coverage: coverage.min(1.0),
        family_support: min_support,
        n_arrays: 0,
    };

    vec![full_chain]
}

// ============ SITE CLUSTERING ============

fn cluster_sites(mut sites: Vec<Site>, k: usize, p: usize) -> Vec<Site> {
    // Merge sites with overlapping mod-P positions (within k)
    sites.sort_by(|a, b| a.position_mod_p.cmp(&b.position_mod_p).then_with(|| a.sequence.cmp(&b.sequence)));

    let mut merged: Vec<Site> = Vec::new();
    for site in sites {
        let should_merge = merged.last().map(|prev: &Site| {
            let d = site.position_mod_p.abs_diff(prev.position_mod_p);
            let d = d.min(p - d);
            d < k
        }).unwrap_or(false);

        if should_merge {
            // Keep the one with higher support
            let last = merged.last_mut().unwrap();
            if site.family_support > last.family_support {
                *last = site;
            }
        } else {
            merged.push(site);
        }
    }

    // Re-index
    for (i, s) in merged.iter_mut().enumerate() {
        s.id = i;
    }
    merged
}

fn decode_kmer_sequences(arrays: &[(String, Vec<u8>)], k: usize, hashes: &[u64]) -> HashMap<u64, String> {
    let hash_set: std::collections::HashSet<u64> = hashes.iter().copied().collect();
    let mut result: HashMap<u64, String> = HashMap::new();

    // Just scan first array until we find all
    for (_, seq) in arrays {
        for i in 0..seq.len().saturating_sub(k - 1) {
            if let Some(hash) = kmer_hash(&seq[i..i + k], k) {
                if hash_set.contains(&hash) && !result.contains_key(&hash) {
                    result.insert(hash, String::from_utf8_lossy(&seq[i..i + k]).to_string());
                }
            }
        }
        if result.len() == hashes.len() { break; }
    }
    result
}

