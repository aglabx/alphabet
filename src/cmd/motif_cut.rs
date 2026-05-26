use std::collections::HashMap;
use clap::Parser;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use crate::monomer::revcomp;
use crate::io::read_fasta;

#[derive(Parser)]
#[command(name = "motif_cut", about = "Cut arrays into monomers at motif boundaries, classify by motif fingerprint + linker lengths")]
struct Args {
    /// Input: .10kb.fasta file
    input: String,

    /// Motifs JSON from discover_motifs (if not provided, uses built-in primate alpha satellite motifs)
    #[arg(short, long)]
    motifs: Option<String>,

    /// Output: monomers TSV
    #[arg(short, long, default_value = "monomers.tsv")]
    output: String,

    /// Output: families report JSON
    #[arg(long, default_value = "families.json")]
    report: String,

    /// Output: family consensus FASTA
    #[arg(long, default_value = "family_consensus.fa")]
    consensus_fa: String,

    /// Output directory for per-letter monomer dumps. One file per
    /// (family, letter); each row is `seq \t metadata...`. Default:
    /// `<dir-of-output-tsv>/letters/`.
    #[arg(long)]
    letters_dir: Option<String>,

    /// Max hamming distance for fuzzy motif match
    #[arg(long, default_value_t = 2)]
    max_mismatch: u32,

    /// Min distinct motifs to keep an array
    #[arg(long, default_value_t = 3)]
    min_motifs: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Deserialize)]
struct ChainsJson {
    #[serde(default)]
    motifs: Vec<MotifDef>,
    #[serde(default)]
    families: Vec<FamilyDef>,
}

#[derive(Debug, Deserialize)]
struct MotifDef {
    name: String,
    sequence: String,
    #[serde(default)]
    position_in_monomer: usize,
}

#[derive(Debug, Deserialize)]
struct FamilyDef {
    id: String,
    period: usize,
    #[serde(default)]
    sf_index: usize,
    motifs: Vec<MotifDef>,
    #[serde(default)]
    arrays: Vec<String>,
}

#[derive(Clone)]
struct AnchorSet {
    names: Vec<String>,
    sequences: Vec<Vec<u8>>,
    rc_sequences: Vec<Vec<u8>>,
    expected_spacers: Vec<usize>, // spacer to next motif (last = wrap)
    motif_len: usize,
    n_motifs: usize,
}

/// One SF family: a motif set plus its claimed arrays.
/// Populated either from chains.json `families` (new) or synthesized
/// from the legacy flat motif list / built-in defaults when no
/// per-family info is available.
#[derive(Clone)]
struct Family {
    id: String,
    period: usize,
    anchors: AnchorSet,
    /// Original array names claimed by this family. Empty means "match
    /// any array" (legacy single-family mode).
    array_whitelist: Option<std::collections::HashSet<String>>,
}

fn anchors_from_motifs(defs: &[MotifDef]) -> AnchorSet {
    let names: Vec<String> = defs.iter().map(|m| m.name.clone()).collect();
    let seqs: Vec<Vec<u8>> = defs.iter().map(|m| m.sequence.as_bytes().to_vec()).collect();
    let rcs: Vec<Vec<u8>> = seqs.iter().map(|s| revcomp(s)).collect();
    let positions: Vec<usize> = defs.iter().map(|m| m.position_in_monomer).collect();
    let mlen = seqs.first().map(|s| s.len()).unwrap_or(11);
    let n = seqs.len();
    let expected_spacers: Vec<usize> = (0..n).map(|i| {
        if i + 1 < n {
            let next_pos = positions[i + 1];
            let this_end = positions[i] + seqs[i].len();
            next_pos.saturating_sub(this_end)
        } else {
            20
        }
    }).collect();
    AnchorSet { names, sequences: seqs, rc_sequences: rcs, expected_spacers, motif_len: mlen, n_motifs: n }
}

/// Resolve motifs into one or more families.
///
/// Preference order:
///   1. chains.json with `families` → one Family per SF, each restricted
///      to its claimed arrays.
///   2. chains.json with flat `motifs` only (legacy) → single Family
///      matching any array.
///   3. no --motifs → built-in primate alpha, single Family matching any.
fn load_families(path: Option<&str>) -> Vec<Family> {
    match path {
        Some(p) => {
            let data: ChainsJson = serde_json::from_str(
                &std::fs::read_to_string(p).unwrap_or_else(|e| panic!("Cannot read {}: {}", p, e))
            ).unwrap();
            if !data.families.is_empty() {
                eprintln!("Loaded {} families from {}:", data.families.len(), p);
                let mut out = Vec::new();
                for fam in data.families {
                    let anchors = anchors_from_motifs(&fam.motifs);
                    eprintln!(
                        "  {}: period={} sf={} motifs={} arrays={}",
                        fam.id, fam.period, fam.sf_index, anchors.n_motifs, fam.arrays.len()
                    );
                    for (i, name) in anchors.names.iter().enumerate() {
                        eprintln!(
                            "    {}: {} ({}bp)",
                            name, String::from_utf8_lossy(&anchors.sequences[i]), anchors.sequences[i].len()
                        );
                    }
                    let whitelist: std::collections::HashSet<String> = fam.arrays.into_iter().collect();
                    out.push(Family {
                        id: fam.id,
                        period: fam.period,
                        anchors,
                        array_whitelist: Some(whitelist),
                    });
                }
                out
            } else {
                let anchors = anchors_from_motifs(&data.motifs);
                eprintln!(
                    "Loaded {} motifs from {} (legacy flat, no families) — running single-family mode",
                    anchors.n_motifs, p
                );
                vec![Family {
                    id: "default".into(),
                    period: 0,
                    anchors,
                    array_whitelist: None,
                }]
            }
        }
        None => {
            let defaults: Vec<(&str, &str)> = vec![
                ("M1", "ACATCACAAAG"),
                ("M2", "AGAATGCTTCT"),
                ("M3", "GAAGATATTTC"),
                ("M4", "TCCACTTGCAG"),
                ("M5", "AAAGAGTGTTT"),
            ];
            let defs: Vec<MotifDef> = defaults.iter().enumerate().map(|(i, (n, s))| MotifDef {
                name: n.to_string(),
                sequence: s.to_string(),
                position_in_monomer: [21, 41, 67, 112, 132][i],
            }).collect();
            let anchors = anchors_from_motifs(&defs);
            eprintln!("Using built-in primate alpha satellite motifs (5 x 11bp)");
            vec![Family {
                id: "default".into(),
                period: 171,
                anchors,
                array_whitelist: None,
            }]
        }
    }
}

// Cut at M1 (first motif) occurrences — monomer = M1...next_M1

#[derive(Debug, Clone)]
struct MotifHit {
    motif_idx: usize,
    position: usize,
    mismatches: u32,
    strand: char,
}

/// A monomer cut from the array
#[derive(Debug, Clone)]
struct Monomer {
    array_id: String,
    /// Family that cut this monomer (e.g. "P171SF0", or "default" in
    /// legacy single-family mode).
    family_id: String,
    monomer_idx: u32,
    start: usize,
    end: usize,
    length: usize,
    sequence: String,
    // Which chain sites are present
    site_present: Vec<bool>,
    // Distances between consecutive present sites
    distances: Vec<(usize, usize, i32)>, // (from_site, to_site, distance)
    // Site order string: actual order of found sites as letters (a,b,c,d,e,...)
    // e.g. "abcde" = full, "ab-de" = c missing, "abcdeabcde" = dimer
    site_order: String,
    // Site structure: sites with inter-site distances
    // e.g. "a18b34c43d18e" or "a100--d18e"
    site_structure: String,
    // Letter key: site presence + binned distances (structural identity)
    letter_key: String,
    // Subtype key: letter_key + site sequences (includes motif mutations)
    subtype_key: String,
}

/// Letter = structural monomer type (sites + distances)
#[derive(Debug, Clone, Serialize)]
struct Letter {
    key: String,
    name: String,        // A, B, C, ... AA, AB, ...
    size: usize,
    site_present: Vec<bool>,
    n_sites: usize,
    distance_pattern: String,
    mean_length: f64,
    consensus: String,
    fraction: f64,
    n_subtypes: usize,
}

#[derive(Debug, Serialize)]
struct CutReport {
    n_arrays: usize,
    n_arrays_cut: usize,
    n_monomers: usize,
    n_letters: usize,
    n_subtypes_total: usize,
    families: Vec<FamilyReport>,
}

#[derive(Debug, Serialize)]
struct FamilyReport {
    id: String,
    period: usize,
    n_arrays: usize,
    n_monomers: usize,
    letters: Vec<Letter>,
}

/// Output of running the scan → normalize → cut → classify pipeline on a
/// single family's arrays.
struct FamilyResult {
    monomers: Vec<Monomer>,
    letters: Vec<Letter>,
    key_to_letter: HashMap<String, String>,
    key_to_subtype: HashMap<String, String>,
    arrays_cut: usize,
    total_subtypes: usize,
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let families = load_families(args.motifs.as_deref());

    // --- Read ALL arrays ---
    eprintln!("\nReading {}...", args.input);
    let all_arrays = read_fasta(&args.input).unwrap_or_else(|e| panic!("{:?}", e));
    eprintln!("  {} total arrays in fasta", all_arrays.len());

    // --- Partition arrays by family ---
    // An array that's whitelisted by multiple families (shouldn't happen
    // for discover_chains output, since SFs claim disjointly, but we
    // still guard against it) goes to the first matching family.
    // Arrays not claimed by any family are skipped.
    let single_family_mode = families.iter().any(|f| f.array_whitelist.is_none());
    let mut buckets: Vec<Vec<(String, Vec<u8>)>> = vec![Vec::new(); families.len()];
    let mut unclaimed = 0usize;
    for (name, seq) in &all_arrays {
        if single_family_mode {
            buckets[0].push((name.clone(), seq.clone()));
            continue;
        }
        let mut placed = false;
        for (fi, fam) in families.iter().enumerate() {
            if let Some(wl) = &fam.array_whitelist {
                if wl.contains(name) {
                    buckets[fi].push((name.clone(), seq.clone()));
                    placed = true;
                    break;
                }
            }
        }
        if !placed { unclaimed += 1; }
    }
    if !single_family_mode {
        eprintln!(
            "  array → family assignment: {} claimed, {} unclaimed (skipped)",
            all_arrays.len() - unclaimed, unclaimed
        );
    }

    // --- Run the per-family pipeline ---
    let mut family_results: Vec<(Family, FamilyResult)> = Vec::new();
    for (fam, fam_arrays) in families.into_iter().zip(buckets.into_iter()) {
        if fam_arrays.is_empty() {
            eprintln!("\n===== Family {} : no arrays, skipping =====", fam.id);
            continue;
        }
        eprintln!("\n===== Family {} : {} arrays =====", fam.id, fam_arrays.len());
        let result = process_family(&fam, &fam_arrays, &args);
        family_results.push((fam, result));
    }

    // --- Aggregate stats for top-level report ---
    let total_monomers: usize = family_results.iter().map(|(_, r)| r.monomers.len()).sum();
    let arrays_cut: usize = family_results.iter().map(|(_, r)| r.arrays_cut).sum();
    let total_letters: usize = family_results.iter().map(|(_, r)| r.letters.len()).sum();
    let total_subtypes: usize = family_results.iter().map(|(_, r)| r.total_subtypes).sum();

    // --- Write outputs ---
    eprintln!("\nWriting outputs...");

    // Monomers TSV — self-describing manifest, then one row per monomer.
    // Contract:
    //   `family`  = SF family that cut this monomer (e.g. "P171SF0")
    //   `array_id` = original FASTA name (no _rc suffix)
    //   `strand`  = '+' if kept forward, '-' if motif_cut revcomped the
    //              array during strand normalization.
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.output).unwrap());
        writeln!(f, "#alphabet v{} / cut", env!("CARGO_PKG_VERSION")).unwrap();
        writeln!(f, "#input: {}", args.input).unwrap();
        writeln!(f, "#motifs_source: {}", args.motifs.as_deref().unwrap_or("<built-in primate alpha>")).unwrap();
        for (fam, _) in &family_results {
            let chain_sites_str: String = fam.anchors.names.iter().zip(fam.anchors.sequences.iter())
                .map(|(n, s)| format!("{}={}", n, String::from_utf8_lossy(s)))
                .collect::<Vec<_>>().join(" ");
            writeln!(f, "#family {} period={} sites: {}", fam.id, fam.period, chain_sites_str).unwrap();
        }
        writeln!(f, "#columns: family array_id strand monomer_idx start end length letter subtype site_order site_structure sites distances sequence").unwrap();
        writeln!(f, "family\tarray_id\tstrand\tmonomer_idx\tstart\tend\tlength\tletter\tsubtype\tsite_order\tsite_structure\tsites\tdistances\tsequence").unwrap();
        for (fam, res) in &family_results {
            for m in &res.monomers {
                let letter = res.key_to_letter.get(&m.letter_key).map(|s| s.as_str()).unwrap_or("?");
                let subtype = res.key_to_subtype.get(&m.subtype_key).map(|s| s.as_str()).unwrap_or("?");
                let sites_str: String = (0..fam.anchors.n_motifs).map(|i| if m.site_present[i] { '1' } else { '0' }).collect();
                let dist_str: String = m.distances.iter()
                    .map(|(_, _, d)| format!("{}", d))
                    .collect::<Vec<_>>().join(",");
                let (raw_id, strand) = match m.array_id.strip_suffix("_rc") {
                    Some(base) => (base, '-'),
                    None => (m.array_id.as_str(), '+'),
                };
                writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    m.family_id, raw_id, strand, m.monomer_idx, m.start, m.end, m.length,
                    letter, subtype, m.site_order, m.site_structure, sites_str, dist_str, m.sequence).unwrap();
            }
        }
    }

    // Consensus FASTA — prefix letter name with family_id so records
    // don't collide between families.
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(&args.consensus_fa).unwrap());
        for (fam, res) in &family_results {
            for lt in &res.letters {
                if lt.consensus.is_empty() { continue; }
                let sites_str: String = (0..fam.anchors.n_motifs).map(|i| if lt.site_present[i] { '1' } else { '0' }).collect();
                writeln!(f, ">{}:letter_{} n={} sites={} len={:.0} subtypes={}",
                    fam.id, lt.name, lt.size, sites_str, lt.mean_length, lt.n_subtypes).unwrap();
                writeln!(f, "{}", lt.consensus).unwrap();
            }
        }
    }

    // HOR strings per array (pipe-separated letters, per-family)
    {
        use std::io::Write;
        let hor_path = args.output.replace(".tsv", "_hor.tsv");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&hor_path).unwrap());
        writeln!(f, "family\tarray_id\tn_monomers\thor_string").unwrap();

        for (_fam, res) in &family_results {
            let mut array_monomers: HashMap<String, Vec<(u32, String)>> = HashMap::new();
            for m in &res.monomers {
                let letter = res.key_to_letter.get(&m.letter_key).map(|s| s.as_str()).unwrap_or("?");
                array_monomers.entry(m.array_id.clone()).or_default().push((m.monomer_idx, letter.to_string()));
            }
            let mut sorted_arrays: Vec<_> = array_monomers.iter().collect();
            sorted_arrays.sort_by(|a, b| b.1.len().cmp(&a.1.len()));
            for (array_id, monomers) in &sorted_arrays {
                let mut mono_sorted: Vec<_> = monomers.to_vec();
                mono_sorted.sort_by_key(|m| m.0);
                let hor_str: String = mono_sorted.iter().map(|(_, l)| l.as_str()).collect::<Vec<_>>().join("|");
                let raw_id = array_id.strip_suffix("_rc").unwrap_or(array_id);
                let fam_id = &res.monomers.first().map(|m| m.family_id.as_str()).unwrap_or("?");
                writeln!(f, "{}\t{}\t{}\t{}", fam_id, raw_id, mono_sorted.len(), hor_str).unwrap();
            }
        }

        eprintln!("  HOR strings: {}", hor_path);
    }

    // Per-letter monomer dumps: one file per (family, letter). Format is
    // one row per monomer, sequence first so downstream `cut`/grep/awk
    // stays easy. File layout: <letters_dir>/<family>/<letter>.tsv.
    {
        use std::io::Write;
        let letters_dir = match &args.letters_dir {
            Some(d) => std::path::PathBuf::from(d),
            None => {
                let out = std::path::PathBuf::from(&args.output);
                out.parent().unwrap_or(std::path::Path::new(".")).join("letters")
            }
        };
        std::fs::create_dir_all(&letters_dir).expect("create letters dir");
        let mut n_files = 0usize;
        for (fam, res) in &family_results {
            let fam_dir = letters_dir.join(&fam.id);
            std::fs::create_dir_all(&fam_dir).expect("create family letters dir");

            // Group monomers by letter name (A, B, C, ...).
            let mut per_letter: HashMap<String, Vec<&Monomer>> = HashMap::new();
            for m in &res.monomers {
                let letter = res.key_to_letter.get(&m.letter_key).cloned().unwrap_or_else(|| "?".into());
                per_letter.entry(letter).or_default().push(m);
            }
            for (letter_name, mons) in &per_letter {
                let path = fam_dir.join(format!("{}.tsv", letter_name));
                let mut f = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
                writeln!(f, "#alphabet v{} / cut letters", env!("CARGO_PKG_VERSION")).unwrap();
                writeln!(f, "#family: {} period={}", fam.id, fam.period).unwrap();
                writeln!(f, "#letter: {} (n={})", letter_name, mons.len()).unwrap();
                writeln!(f, "#columns: sequence array_id strand start end length monomer_idx subtype site_order site_structure").unwrap();
                for m in mons {
                    let subtype = res.key_to_subtype.get(&m.subtype_key).map(|s| s.as_str()).unwrap_or("?");
                    let (raw_id, strand) = match m.array_id.strip_suffix("_rc") {
                        Some(base) => (base, '-'),
                        None => (m.array_id.as_str(), '+'),
                    };
                    writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        m.sequence, raw_id, strand, m.start, m.end, m.length,
                        m.monomer_idx, subtype, m.site_order, m.site_structure).unwrap();
                }
                n_files += 1;
            }
        }
        eprintln!("  per-letter dumps: {} files under {}", n_files, letters_dir.display());
    }

    // Report JSON
    let report = CutReport {
        n_arrays: all_arrays.len(),
        n_arrays_cut: arrays_cut,
        n_monomers: total_monomers,
        n_letters: total_letters,
        n_subtypes_total: total_subtypes,
        families: family_results.iter().map(|(fam, res)| FamilyReport {
            id: fam.id.clone(),
            period: fam.period,
            n_arrays: res.arrays_cut,
            n_monomers: res.monomers.len(),
            letters: res.letters.clone(),
        }).collect(),
    };
    std::fs::write(&args.report, serde_json::to_string_pretty(&report).unwrap()).unwrap();

    eprintln!("Done.");
}

/// Run scan → normalize → cut → classify on one family's arrays.
fn process_family(fam: &Family, family_arrays: &[(String, Vec<u8>)], args: &Args) -> FamilyResult {
    let anchors = &fam.anchors;
    let motifs: Vec<(Vec<u8>, Vec<u8>)> = anchors.sequences.iter().zip(anchors.rc_sequences.iter())
        .map(|(f, r)| (f.clone(), r.clone())).collect();
    let n_motifs = anchors.n_motifs;
    let motif_len = anchors.motif_len;

    eprintln!("  scanning motifs (min {} distinct)...", args.min_motifs);
    let scan_results: Vec<(String, Vec<u8>, Vec<MotifHit>, usize)> = family_arrays.par_iter().map(|(name, seq)| {
        let hits = scan_motifs(seq, &motifs, args.max_mismatch);
        let distinct = {
            let mut seen = vec![false; n_motifs];
            for h in &hits { seen[h.motif_idx] = true; }
            seen.iter().filter(|&&b| b).count()
        };
        (name.clone(), seq.clone(), hits, distinct)
    }).collect();

    let arrays_passing = scan_results.iter().filter(|(_, _, _, dm)| *dm >= args.min_motifs).count();
    eprintln!("  arrays passing motif filter: {}", arrays_passing);

    eprintln!("  strand normalizing...");
    let normalized: Vec<(String, Vec<u8>, Vec<MotifHit>, usize)> = scan_results.into_par_iter()
        .filter(|(_, _, _, dm)| *dm >= args.min_motifs)
        .map(|(name, seq, hits, dm)| {
            let fwd_count = hits.iter().filter(|h| h.strand == '+').count();
            let rc_count = hits.iter().filter(|h| h.strand == '-').count();
            if rc_count > fwd_count {
                let rc_seq = revcomp(&seq);
                let new_hits = scan_motifs(&rc_seq, &motifs, args.max_mismatch);
                let new_dm = {
                    let mut seen = vec![false; n_motifs];
                    for h in &new_hits { seen[h.motif_idx] = true; }
                    seen.iter().filter(|&&b| b).count()
                };
                (format!("{}_rc", name), rc_seq, new_hits, new_dm)
            } else {
                (name, seq, hits, dm)
            }
        }).collect();

    let n_flipped = normalized.iter().filter(|(n, _, _, _)| n.ends_with("_rc")).count();
    eprintln!("  kept fwd: {}, flipped: {}", normalized.len() - n_flipped, n_flipped);

    eprintln!("  cutting monomers...");
    let ref_seqs: Vec<String> = anchors.sequences.iter().map(|s| String::from_utf8_lossy(s).to_string()).collect();
    let all_monomers: Vec<Vec<Monomer>> = normalized.par_iter()
        .map(|(name, seq, hits, _)| {
            cut_at_motifs(name, &fam.id, seq, hits, n_motifs, motif_len, &anchors.names, &anchors.expected_spacers, &ref_seqs)
        }).collect();
    let arrays_cut = all_monomers.iter().filter(|v| !v.is_empty()).count();
    let flat: Vec<Monomer> = all_monomers.into_iter().flatten().collect();
    eprintln!("  {} monomers from {} arrays", flat.len(), arrays_cut);

    // Classify within this family only.
    let n_total = flat.len();
    let mut letter_members: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, m) in flat.iter().enumerate() {
        letter_members.entry(m.letter_key.clone()).or_default().push(i);
    }
    let mut letter_list: Vec<(String, Vec<usize>)> = letter_members.into_iter().collect();
    letter_list.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    let mut letters: Vec<Letter> = Vec::new();
    let mut key_to_letter: HashMap<String, String> = HashMap::new();
    let mut key_to_subtype: HashMap<String, String> = HashMap::new();
    let mut total_subtypes = 0usize;
    for (rank, (key, members)) in letter_list.iter().enumerate() {
        let letter_name = make_letter_name(rank);
        key_to_letter.insert(key.clone(), letter_name.clone());
        let sample = &flat[members[0]];
        let mean_len = members.iter().map(|&i| flat[i].length as f64).sum::<f64>() / members.len() as f64;
        let dist_str = if !sample.distances.is_empty() {
            sample.distances.iter()
                .map(|(f, t, d)| format!("{}→{}:{}", anchors.names[*f], anchors.names[*t], d))
                .collect::<Vec<_>>().join(", ")
        } else {
            "no_sites".to_string()
        };
        let mut subtype_counts: HashMap<String, usize> = HashMap::new();
        for &mi in members {
            *subtype_counts.entry(flat[mi].subtype_key.clone()).or_insert(0) += 1;
        }
        let n_subtypes = subtype_counts.len();
        total_subtypes += n_subtypes;
        let mut subtype_list: Vec<(String, usize)> = subtype_counts.into_iter().collect();
        subtype_list.sort_by(|a, b| b.1.cmp(&a.1));
        for (si, (skey, _)) in subtype_list.iter().enumerate() {
            key_to_subtype.insert(skey.clone(), format!("{}{}", letter_name, si + 1));
        }
        let consensus = {
            let seqs: Vec<&[u8]> = members.iter()
                .take(5000)
                .map(|&i| flat[i].sequence.as_bytes())
                .filter(|s| s.len() == flat[members[0]].length)
                .collect();
            if !seqs.is_empty() && seqs[0].len() < 500 {
                build_consensus(&seqs, seqs[0].len())
            } else { String::new() }
        };
        letters.push(Letter {
            key: key.clone(),
            name: letter_name,
            size: members.len(),
            site_present: sample.site_present.clone(),
            n_sites: sample.site_present.iter().filter(|&&b| b).count(),
            distance_pattern: dist_str,
            mean_length: mean_len,
            consensus,
            fraction: if n_total > 0 { members.len() as f64 / n_total as f64 } else { 0.0 },
            n_subtypes,
        });
    }

    eprintln!("  {}: {} letters, {} subtypes, {} monomers", fam.id, letters.len(), total_subtypes, n_total);
    for lt in letters.iter().take(10) {
        let sites_str: String = (0..n_motifs).map(|i| if lt.site_present[i] { '●' } else { '·' }).collect();
        eprintln!("    {:<4} n={:>6} ({:>4.1}%) len={:>5.0} sites={} sub={:>3}",
            lt.name, lt.size, lt.fraction * 100.0, lt.mean_length, sites_str, lt.n_subtypes);
    }
    if letters.len() > 10 {
        eprintln!("    ... {} more letters", letters.len() - 10);
    }

    FamilyResult {
        monomers: flat,
        letters,
        key_to_letter,
        key_to_subtype,
        arrays_cut,
        total_subtypes,
    }
}

fn site_display_char(i: usize) -> char {
    match i {
        0..=25 => (b'a' + i as u8) as char,
        26..=51 => (b'A' + (i - 26) as u8) as char,
        _ => '?',
    }
}

fn make_letter_name(rank: usize) -> String {
    if rank < 26 {
        format!("{}", (b'A' + rank as u8) as char)
    } else {
        let first = (rank / 26 - 1).min(25);
        let second = rank % 26;
        format!("{}{}", (b'A' + first as u8) as char, (b'A' + second as u8) as char)
    }
}

// ============ CUTTING LOGIC ============

fn cut_at_motifs(array_id: &str, family_id: &str, seq: &[u8], hits: &[MotifHit], n_motifs: usize, motif_len: usize, anchor_names: &[String], expected_spacers: &[usize], anchors_ref: &[String]) -> Vec<Monomer> {
    if hits.is_empty() { return vec![]; }

    // Strategy: find M1 occurrences as monomer start points.
    // A monomer runs from one M1 to the next M1.
    // If no M1, fall back to the most frequent motif.

    // Find the best "start motif" — the one that appears most regularly
    let m1_hits: Vec<usize> = hits.iter().enumerate()
        .filter(|(_, h)| h.motif_idx == 0) // M1
        .map(|(i, _)| i)
        .collect();

    // If not enough M1 hits, use whatever motif is most frequent
    let (_start_motif, start_positions) = if m1_hits.len() >= 2 {
        (0usize, m1_hits)
    } else {
        // Find motif with most hits
        let mut counts = vec![0usize; n_motifs];
        for h in hits { counts[h.motif_idx] += 1; }
        let best = counts.iter().enumerate().max_by_key(|(_, &c)| c).unwrap().0;
        let positions: Vec<usize> = hits.iter().enumerate()
            .filter(|(_, h)| h.motif_idx == best)
            .map(|(i, _)| i)
            .collect();
        (best, positions)
    };

    // Cut: each monomer starts at a start_motif hit position
    let mut monomers: Vec<Monomer> = Vec::new();

    for wi in 0..start_positions.len() {
        let hit_idx = start_positions[wi];
        let start_pos = hits[hit_idx].position;
        let end_pos = if wi + 1 < start_positions.len() {
            hits[start_positions[wi + 1]].position
        } else {
            seq.len()
        };

        // Skip if too short or too long
        let length = end_pos - start_pos;
        if !(100..=400).contains(&length) { continue; }

        let mono_seq = &seq[start_pos..end_pos];

        // Find which sites are present + extract their actual sequences
        let mut site_present = vec![false; n_motifs];
        let mut site_sequences = vec![String::new(); n_motifs];
        let mut mono_hits: Vec<(usize, usize)> = Vec::new(); // (site_idx, relative_position)

        for h in hits {
            if h.position >= start_pos && h.position + motif_len <= end_pos {
                let si = h.motif_idx;
                if !site_present[si] { // take first hit per site
                    site_present[si] = true;
                    let site_start = h.position;
                    let ref_len = anchors_ref[si].len();
                    let ext_end = (site_start + ref_len).min(end_pos);
                    // Exact match found by scanner → sequence = reference (or very close)
                    site_sequences[si] = String::from_utf8_lossy(&seq[site_start..ext_end]).to_string();
                    // CIGAR: "=" for exact match sites (no indels by definition)
                }
                mono_hits.push((si, h.position - start_pos));
            }
        }
        mono_hits.sort_by_key(|&(_, p)| p);

        // === Distance-based presence inference ===
        // For pairs of found sites, check if missing sites between them are
        // present (distance matches expected) or truly deleted (distance collapsed)
        // Dedup mono_hits first for this analysis
        let mut deduped_hits: Vec<(usize, usize)> = Vec::new();
        {
            let mut seen = vec![false; n_motifs];
            for &(si, pos) in &mono_hits {
                if !seen[si] { seen[si] = true; deduped_hits.push((si, pos)); }
            }
        }

        if deduped_hits.len() >= 2 {
            // For each pair of found sites, infer intermediate presence by distance
            for wi in 0..deduped_hits.len() - 1 {
                let (si, pos_i) = deduped_hits[wi];
                let (sj, pos_j) = deduped_hits[wi + 1];

                // Sites between si and sj in chain order
                if sj <= si { continue; } // only forward chain order

                for mid in (si + 1)..sj {
                    if site_present[mid] { continue; } // already found exact

                    // Compute expected distance si→sj using real chain spacers
                    let total_expected = {
                        let mut d = 0usize;
                        for k in si..sj {
                            d += anchor_names[k].len(); // motif length
                            d += expected_spacers.get(k).copied().unwrap_or(10); // spacer
                        }
                        d
                    };
                    let total_observed = pos_j - pos_i;

                    // Tolerance: 5bp
                    let tolerance = 5u32;
                    let delta = total_observed as i32 - total_expected as i32;
                    if delta.unsigned_abs() <= tolerance {
                        site_present[mid] = true;
                        // Compute predicted position of mid using chain distances
                        let mut pred_offset = 0usize;
                        for k in si..mid {
                            pred_offset += anchor_names[k].len();
                            pred_offset += expected_spacers.get(k).copied().unwrap_or(10);
                        }
                        let abs_pos = pos_i + pred_offset;
                        let mid_ref = anchors_ref[mid].as_bytes();
                        let mid_len = mid_ref.len();
                        // Extract region: pad ±2bp to accommodate small indels
                        let pad = 2;
                        let ext_start = abs_pos.saturating_sub(pad);
                        let ext_end = (abs_pos + mid_len + pad).min(mono_seq.len());
                        if ext_end > ext_start {
                            let region = &mono_seq[ext_start..ext_end];
                            let (cigar, aligned_seq) = align_and_cigar(mid_ref, region);
                            site_sequences[mid] = format!("{}:{}", aligned_seq, cigar);
                        }
                    }
                }
            }
        }

        // Add inferred sites to mono_hits at their predicted positions
        for mid in 0..n_motifs {
            if site_present[mid] && !mono_hits.iter().any(|&(si, _)| si == mid) {
                // Inferred site: find its predicted position from nearest found neighbor
                if let Some(&(prev_si, prev_pos)) = deduped_hits.iter().rfind(|&&(si, _)| si < mid) {
                    let mut pred_offset = 0usize;
                    for k in prev_si..mid {
                        pred_offset += anchor_names[k].len();
                        pred_offset += expected_spacers.get(k).copied().unwrap_or(10);
                    }
                    mono_hits.push((mid, prev_pos + pred_offset));
                }
            }
        }
        mono_hits.sort_by_key(|&(_, p)| p);

        // Build site_order string: all sites (found + inferred) in positional order.
        // Single-char display: a..z (0..25), A..Z (26..51), '?' for overflow.
        // Real satellites have <<52 motifs, so the overflow path only triggers
        // for degenerate inputs (e.g. garbage-in arrays from misparsed periods).
        let site_letters: Vec<char> = (0..n_motifs).map(site_display_char).collect();
        let all_hits_order: String = mono_hits.iter()
            .map(|&(si, _)| if si < site_letters.len() { site_letters[si] } else { '?' })
            .collect();
        // Also build version with gaps: expected chain order, '-' for missing
        let chain_with_gaps: String = (0..n_motifs).map(|i| {
            if site_present[i] { site_letters[i] } else { '-' }
        }).collect();
        // Use the ALL hits order (shows duplications)
        let site_order = if all_hits_order.len() > n_motifs {
            all_hits_order.clone()
        } else {
            chain_with_gaps
        };

        // Build site_structure with gap_before:sites:gap_after
        let site_structure: String = {
            let mut parts: Vec<String> = Vec::new();
            let mono_len = mono_seq.len();

            // gap_before: distance from monomer start to first site
            let gap_before = if let Some(&(_, first_pos)) = mono_hits.first() {
                first_pos
            } else { mono_len };
            parts.push(format!("{}:", gap_before));

            let mut prev_end: Option<usize> = None;
            let mut prev_si: Option<usize> = None;
            for &(si, pos) in &mono_hits {
                if let (Some(pe), Some(ps)) = (prev_end, prev_si) {
                    let gap = pos as i32 - pe as i32;
                    let prev_chain_idx = (0..n_motifs).position(|i| i == ps).unwrap_or(0);
                    let curr_chain_idx = (0..n_motifs).position(|i| i == si).unwrap_or(0);
                    if curr_chain_idx > prev_chain_idx + 1 {
                        let n_missing = curr_chain_idx - prev_chain_idx - 1;
                        let dashes: String = "-".repeat(n_missing);
                        parts.push(format!("{}{}", gap, dashes));
                    } else {
                        parts.push(format!("{}", gap));
                    }
                }
                // Site letter + CIGAR if inferred with indel
                let letter = if si < site_letters.len() { site_letters[si] } else { '?' };
                let cigar_part = if site_sequences[si].contains(':') {
                    // Inferred site with CIGAR: "seq:CIGAR"
                    let cigar = site_sequences[si].split(':').next_back().unwrap_or("=");
                    if cigar != "=" { format!("{}[{}]", letter, cigar) } else { format!("{}", letter) }
                } else {
                    format!("{}", letter)
                };
                parts.push(cigar_part);
                prev_end = Some(pos + anchor_names[si].len());
                prev_si = Some(si);
            }
            // gap_after: distance from last site end to monomer end
            let gap_after = if let Some(pe) = prev_end {
                mono_len.saturating_sub(pe)
            } else { mono_len };
            parts.push(format!(":{}", gap_after));

            if parts.len() <= 2 {
                "-".to_string()
            } else {
                parts.join("")
            }
        };

        // Dedup for distances: keep only first hit per site
        let mut seen_sites = vec![false; n_motifs];
        mono_hits.retain(|&(si, _)| {
            if seen_sites[si] { false } else { seen_sites[si] = true; true }
        });

        // Compute distances between consecutive present sites
        let mut distances: Vec<(usize, usize, i32)> = Vec::new();
        for i in 0..mono_hits.len().saturating_sub(1) {
            let (from_s, from_p) = mono_hits[i];
            let (to_s, to_p) = mono_hits[i + 1];
            let from_len = anchor_names[from_s].len();
            let dist = (to_p as i32) - (from_p as i32) - (from_len as i32);
            distances.push((from_s, to_s, dist));
        }

        // Letter key: site presence + binned distances (structural identity)
        let presence_key: String = (0..n_motifs).map(|i| if site_present[i] { '1' } else { '0' }).collect();
        let distance_key: String = distances.iter()
            .map(|(_, _, d)| {
                let binned = ((*d + 10) / 20) * 20;
                format!("{}", binned)
            })
            .collect::<Vec<_>>().join(",");
        let letter_key = format!("{}:{}", presence_key, distance_key);

        // Subtype key: letter + actual site sequences (captures motif mutations)
        let sites_str: String = (0..n_motifs).map(|i| {
            if site_present[i] { site_sequences[i].as_str() } else { "-" }
        }).collect::<Vec<_>>().join("|");
        let subtype_key = format!("{}:{}", letter_key, sites_str);

        monomers.push(Monomer {
            array_id: array_id.to_string(),
            family_id: family_id.to_string(),
            monomer_idx: monomers.len() as u32,
            start: start_pos,
            end: end_pos,
            length,
            sequence: String::from_utf8_lossy(mono_seq).to_string(),
            site_present,
            distances,
            site_order,
            site_structure,
            letter_key,
            subtype_key,
        });
    }

    monomers
}

// ============ MOTIF SCANNING ============

fn scan_motifs(seq: &[u8], motifs: &[(Vec<u8>, Vec<u8>)], max_mismatch: u32) -> Vec<MotifHit> {
    let mut hits: Vec<MotifHit> = Vec::new();

    for (midx, (fwd, rc)) in motifs.iter().enumerate() {
        let mlen = fwd.len();
        if seq.len() < mlen { continue; }

        for i in 0..=(seq.len() - mlen) {
            let window = &seq[i..i + mlen];
            let fwd_mm = hamming_bytes(window, fwd);
            if fwd_mm <= max_mismatch {
                hits.push(MotifHit { motif_idx: midx, position: i, mismatches: fwd_mm, strand: '+' });
                continue;
            }
            let rc_mm = hamming_bytes(window, rc);
            if rc_mm <= max_mismatch {
                hits.push(MotifHit { motif_idx: midx, position: i, mismatches: rc_mm, strand: '-' });
            }
        }
    }

    hits.sort_by_key(|h| h.position);
    let ml = motifs.first().map(|(f, _)| f.len()).unwrap_or(11);
    dedup_hits(&mut hits, ml);
    hits
}

fn dedup_hits(hits: &mut Vec<MotifHit>, motif_len: usize) {
    if hits.len() < 2 { return; }
    let mut keep = vec![true; hits.len()];
    for i in 0..hits.len() {
        if !keep[i] { continue; }
        for j in (i + 1)..hits.len() {
            if !keep[j] { continue; }
            if hits[j].position < hits[i].position + motif_len {
                if hits[j].mismatches < hits[i].mismatches { keep[i] = false; break; }
                else { keep[j] = false; }
            } else { break; }
        }
    }
    let mut w = 0;
    for r in 0..hits.len() { if keep[r] { hits.swap(w, r); w += 1; } }
    hits.truncate(w);
}

fn hamming_bytes(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b.iter()).filter(|(x, y)| !x.eq_ignore_ascii_case(y)).count() as u32
}

/// Simple NW alignment of reference motif to extracted region, returns CIGAR string
fn align_and_cigar(reference: &[u8], query: &[u8]) -> (String, String) {
    let n = reference.len();
    let m = query.len();
    if n == 0 || m == 0 {
        return ("*".to_string(), String::new());
    }

    // Affine gap NW: gap_open expensive, gap_extend cheap → grouped indels
    let match_score: i32 = 2;
    let mismatch: i32 = -1;
    let gap_open: i32 = -6;
    let gap_extend: i32 = -1;

    // DP matrices: M (match/mismatch), X (gap in query/deletion), Y (gap in ref/insertion)
    let inf = i32::MIN / 2;
    let mut dm = vec![vec![inf; m + 1]; n + 1]; // best ending in M
    let mut dx = vec![vec![inf; m + 1]; n + 1]; // best ending in gap-in-query (D)
    let mut dy = vec![vec![inf; m + 1]; n + 1]; // best ending in gap-in-ref (I)

    dm[0][0] = 0;
    for i in 1..=n { dx[i][0] = gap_open + gap_extend * (i as i32 - 1); }
    for j in 1..=m { dy[0][j] = gap_open + gap_extend * (j as i32 - 1); }

    for i in 1..=n {
        for j in 1..=m {
            let s = if reference[i-1].eq_ignore_ascii_case(&query[j-1]) {
                match_score
            } else { mismatch };
            dm[i][j] = (dm[i-1][j-1].max(dx[i-1][j-1]).max(dy[i-1][j-1])) + s;
            dx[i][j] = (dm[i-1][j] + gap_open).max(dx[i-1][j] + gap_extend);
            dy[i][j] = (dm[i][j-1] + gap_open).max(dy[i][j-1] + gap_extend);
        }
    }

    // Combine for traceback
    let mut dp = vec![vec![0i32; m + 1]; n + 1];
    for i in 0..=n {
        for j in 0..=m {
            dp[i][j] = dm[i][j].max(dx[i][j]).max(dy[i][j]);
        }
    }

    // Traceback using affine matrices
    let mut cigar_ops: Vec<(char, u32)> = Vec::new();
    let mut aligned_query = Vec::new();
    let mut i = n;
    let mut j = m;
    // Determine which matrix we end in
    let mut state = if dm[n][m] >= dx[n][m] && dm[n][m] >= dy[n][m] { 'M' }
        else if dx[n][m] >= dy[n][m] { 'D' } else { 'I' };

    while i > 0 || j > 0 {
        match state {
            'M' if i > 0 && j > 0 => {
                let _s = if reference[i-1].eq_ignore_ascii_case(&query[j-1]) { match_score } else { mismatch };
                let op = if reference[i-1].eq_ignore_ascii_case(&query[j-1]) { 'M' } else { 'X' };
                aligned_query.push(query[j-1]);
                if let Some(last) = cigar_ops.last_mut() {
                    if last.0 == op { last.1 += 1; } else { cigar_ops.push((op, 1)); }
                } else { cigar_ops.push((op, 1)); }
                // Where did we come from?
                let prev = dm[i-1][j-1].max(dx[i-1][j-1]).max(dy[i-1][j-1]);
                state = if prev == dm[i-1][j-1] { 'M' } else if prev == dx[i-1][j-1] { 'D' } else { 'I' };
                i -= 1; j -= 1;
            }
            'D' if i > 0 => {
                if let Some(last) = cigar_ops.last_mut() {
                    if last.0 == 'D' { last.1 += 1; } else { cigar_ops.push(('D', 1)); }
                } else { cigar_ops.push(('D', 1)); }
                state = if dm[i-1][j] + gap_open >= dx[i-1][j] + gap_extend { 'M' } else { 'D' };
                i -= 1;
            }
            'I' if j > 0 => {
                aligned_query.push(query[j-1]);
                if let Some(last) = cigar_ops.last_mut() {
                    if last.0 == 'I' { last.1 += 1; } else { cigar_ops.push(('I', 1)); }
                } else { cigar_ops.push(('I', 1)); }
                state = if dm[i][j-1] + gap_open >= dy[i][j-1] + gap_extend { 'M' } else { 'I' };
                j -= 1;
            }
            _ => {
                // Fallback: consume remaining
                if i > 0 { cigar_ops.push(('D', i as u32)); i = 0; }
                if j > 0 { for k in (0..j).rev() { aligned_query.push(query[k]); } cigar_ops.push(('I', j as u32)); j = 0; }
            }
        }
    }

    cigar_ops.reverse();
    aligned_query.reverse();

    // Merge M and X into M, then compress consecutive same ops
    let mut merged: Vec<(char, u32)> = Vec::new();
    for &(op, len) in &cigar_ops {
        let c = if op == 'X' { 'M' } else { op };
        if let Some(last) = merged.last_mut() {
            if last.0 == c { last.1 += len; } else { merged.push((c, len)); }
        } else {
            merged.push((c, len));
        }
    }

    let has_indel = merged.iter().any(|(op, _)| *op == 'D' || *op == 'I');
    let final_cigar = if !has_indel {
        "=".to_string()
    } else {
        merged.iter().map(|(op, len)| format!("{}{}", len, op)).collect()
    };

    let aligned_str = String::from_utf8_lossy(&aligned_query).to_string();
    (final_cigar, aligned_str)
}

fn build_consensus(seqs: &[&[u8]], len: usize) -> String {
    let mut counts = vec![[0u32; 4]; len];
    for seq in seqs {
        for (pos, &b) in seq.iter().enumerate().take(len) {
            let idx = match b { b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => continue };
            counts[pos][idx] += 1;
        }
    }
    let bases = [b'A', b'C', b'G', b'T'];
    counts.iter().map(|c| bases[c.iter().enumerate().max_by_key(|(_, &v)| v).unwrap().0] as char).collect()
}
