//! `cut-v2` — multi-anchor position-aware monomer cutter for alpha satellite.
//!
//! Ported from the Python prototype at
//! `research/2026_03_30_alphasplitter/state/scripts/cut_v2_multi_anchor.py`.
//!
//! Differences from `cut` (motif_cut.rs):
//!   * Uses ALL panel anchors with weighted score, not just M1.
//!   * Tight length constraint via per-step ±tol_period window around `last + period`,
//!     not a global [100, 400] filter.
//!   * Bi-positional anchors handled (one anchor can vote for two canonical offsets).
//!   * Skipped monomers are interpolated rather than dropped, so arrays that lost
//!     the landmark in one or two consecutive monomers don't become outliers.
//!   * Explicit outlier classification (too_short / no_first_cut / only_one_cut).

use std::fs::File;
use std::io::{BufWriter, Write};

use clap::Parser;
use rayon::prelude::*;

use crate::io::read_fasta;

/// Data-derived 15-cluster anchor grid (May 2026 re-derivation on 1.18M
/// cut monomers from 7-haplotype T2T panel). Replaces the May 2025
/// empirical 11-entry panel (preserved as PANEL_LEGACY below for
/// comparison). Each entry uses the highest-DF 5-mer from its cluster
/// as the single representative.
pub const PANEL: &[(&[u8], &[usize], usize, f64)] = &[
    // (5-mer, &[canonical_positions], tolerance, weight)
    (b"GAAAC", &[57, 162], 3, 1.0),   // cluster_0  (66 members, DF 97.9%) — bi-pos
    (b"ACAGA", &[147, 22], 3, 1.0),   // cluster_1  (63 members, DF 97.0%) — bi-pos
    (b"TCTTT", &[64, 37],  3, 1.0),   // cluster_2  (61 members, DF 97.0%) — bi-pos
    (b"TTGGA", &[92, 114], 3, 1.0),   // cluster_7  (58 members, DF 95.8%) — bi-pos
    (b"CATTC", &[154, 12], 3, 1.0),   // cluster_3  (58 members, DF 95.6%) — bi-pos
    (b"TCTTC", &[129],     3, 0.9),   // cluster_4  (32 members, DF 93.0%)
    (b"GATGT", &[4],       3, 0.9),   // cluster_6  (26 members, DF 93.0%)
    (b"ATCTG", &[76],      3, 0.9),   // cluster_8  (23 members, DF 91.9%)
    (b"AGCAG", &[48],      3, 0.9),   // cluster_5  (30 members, DF 92.5%)
    (b"TAAAA", &[137],     3, 0.9),   // cluster_10 (35 members, DF 92.4%)
    (b"GTGGA", &[84],      3, 0.9),   // cluster_11 (25 members, DF 92.0%)
    (b"AGGCC", &[105],     3, 0.85),  // cluster_16 (26 members, DF 85.9%)
    (b"TGAAC", &[29],      3, 0.85),  // cluster_9  (20 members, DF 87.7%)
    (b"AAAGG", &[119],     3, 0.85),  // cluster_13 (12 members, DF 86.1%)
    (b"GCTTT", &[99],      3, 0.85),  // cluster_14 (17 members, DF 86.3%)
];

/// Old empirical 11-entry panel (May 2025). Preserved for benchmark
/// comparison; not used by cut-v2 default.
pub const PANEL_LEGACY: &[(&[u8], &[usize], usize, f64)] = &[
    (b"CTTTGTGATGT", &[0],         2,  2.0),
    (b"CAGAG",       &[25],        3,  1.0),
    (b"CTTTT",       &[40],        15, 0.5),
    (b"GTGGA",       &[86],        3,  1.0),
    (b"TGGAA",       &[117],       2,  1.0),
    (b"AGAAA",       &[162],       3,  1.0),
    (b"CATTC",       &[14, 155],   3,  0.8),
    (b"CAGAA",       &[149, 161],  3,  0.8),
    (b"GAAAC",       &[59, 163],   3,  0.7),
    (b"ACAGA",       &[24, 148],   3,  0.7),
    (b"AAACT",       &[141, 164],  3,  0.6),
];

#[derive(Parser)]
#[command(
    name = "cut-v2",
    about = "Multi-anchor position-aware monomer cutter (alpha satellite)"
)]
struct Args {
    /// Input FASTA file(s) — canonical-strand alpha satellite arrays.
    #[arg(long = "fasta", num_args = 1.., required = true, value_name = "PATH")]
    fasta: Vec<String>,

    /// Expected monomer period in bp.
    #[arg(long, default_value_t = 171)]
    period: usize,

    /// ± tolerance per cut step (window around last + period).
    #[arg(long = "tol-period", default_value_t = 30)]
    tol_period: usize,

    /// Minimum anchor score to accept a cut.
    #[arg(long = "min-score", default_value_t = 2.5)]
    min_score: f64,

    /// Maximum consecutive interpolated monomers before giving up on this array.
    #[arg(long = "max-skip", default_value_t = 10)]
    max_skip: usize,

    /// Rayon threads (0 = library default).
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,

    /// Retry the reverse complement when the forward pass finds no cut, and
    /// report which strand each array was cut on.
    ///
    /// The anchor table is canonical-strand, so an input written the other way
    /// round yields `no_first_cut` however much satellite it contains. That is
    /// correct for the documented input (assembly arrays in canonical
    /// orientation) and wrong for anything else: on PacBio HiFi reads, which
    /// arrive in random orientation, the forward-only pass refuses ~48% of
    /// alpha-bearing reads and the refused set is ~100% minus-strand.
    ///
    /// With this flag `start`/`end` stay in the INPUT frame while `mono_idx`
    /// follows the cut order along the strand that matched, so a `-` array has
    /// decreasing `start` as `mono_idx` rises. `sequence` is the monomer as
    /// cut, i.e. on the matching strand. Adds an `orient` column; without the
    /// flag the output is unchanged.
    #[arg(long = "auto-orient", default_value_t = false)]
    auto_orient: bool,

    /// Output: monomer rows TSV.
    #[arg(long = "out-monomers", default_value = "monomers.tsv")]
    out_monomers: String,

    /// Output: per-array summary TSV.
    #[arg(long = "out-summary", default_value = "summary.tsv")]
    out_summary: String,

    /// Output: outlier arrays TSV.
    #[arg(long = "out-outliers", default_value = "outliers.tsv")]
    out_outliers: String,

    /// Output: exceptions TSV — direct cuts whose gap_length deviates from
    /// `k_chosen * period` by more than `--exception-threshold` bp. These are
    /// inputs the cutter DID process but where the period assumption was
    /// measurably violated; downstream consumers should treat the surrounding
    /// monomers as suspect.
    #[arg(long = "out-exceptions", default_value = "exceptions.tsv")]
    out_exceptions: String,

    /// Residual threshold (bp) above which a gap counts as an exception.
    /// Default 20 — derived from empirical 98% containment of valid alpha
    /// monomers within ±20 bp of period=171 in the 7-hap T2T panel.
    #[arg(long = "exception-threshold", default_value_t = 20)]
    exception_threshold: usize,

    /// Bases of sequence context emitted around each exception (left and right).
    #[arg(long = "exception-context", default_value_t = 30)]
    exception_context: usize,
}

/// Cutter parameters — shared by CLI and tests.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    pub period: usize,
    pub tol_period: usize,
    pub min_score: f64,
    pub max_skip: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self { period: 171, tol_period: 30, min_score: 2.5, max_skip: 10 }
    }
}

/// All non-overlapping start positions of `sub` in `seq`.
///
/// Non-overlapping (advance by `sub.len()` after each hit) mirrors the Python
/// prototype — these panel anchors aren't self-overlapping at biologically
/// meaningful offsets, and overlap would double-count weight.
pub fn find_all_positions(seq: &[u8], sub: &[u8]) -> Vec<usize> {
    let mut out = Vec::new();
    let n = sub.len();
    if n == 0 || seq.len() < n {
        return out;
    }
    let mut i = 0usize;
    while i + n <= seq.len() {
        if &seq[i..i + n] == sub {
            out.push(i);
            i += n;
        } else {
            i += 1;
        }
    }
    out
}

/// `score[s]` = total panel weight contributed if a monomer were to start at `s`.
pub fn compute_score(seq: &[u8]) -> Vec<f64> {
    let l = seq.len();
    let mut score = vec![0.0f64; l];
    if l == 0 {
        return score;
    }
    let lmax = (l - 1) as isize;
    for (anchor, positions, tol, weight) in PANEL {
        let tol = *tol as isize;
        let w = *weight;
        for hit in find_all_positions(seq, anchor) {
            let hp = hit as isize;
            for &p in *positions {
                let p = p as isize;
                let s_lo_i = hp - p - tol;
                let s_hi_i = hp - p + tol;
                if s_hi_i < 0 {
                    continue;
                }
                let s_lo = s_lo_i.max(0) as usize;
                let s_hi = s_hi_i.min(lmax) as usize;
                for cell in &mut score[s_lo..=s_hi] {
                    *cell += w;
                }
            }
        }
    }
    score
}

/// One direct cut placement event during walk: (last_cut, next_cut, gap_length, k_chosen).
/// `k_chosen == 1` means no interpolation; `k_chosen >= 2` means `k_chosen - 1`
/// interpolated cuts were inserted between `last_cut` and `next_cut`.
#[derive(Debug, Clone, Copy)]
pub struct GapEvent {
    pub last_cut: usize,
    pub next_cut: usize,
    pub gap_length: usize,
    pub k_chosen: usize,
}

/// Walk left-to-right finding cuts at ~period spacing with score >= `min_score`.
/// Returns `(cuts, interpolated_flags, gap_events)`. Empty if no first cut found.
pub fn walk_cuts(score: &[f64], params: &Params) -> (Vec<usize>, Vec<bool>, Vec<GapEvent>) {
    let l = score.len();
    let mut cuts: Vec<usize> = Vec::new();
    let mut interp: Vec<bool> = Vec::new();
    let mut gap_events: Vec<GapEvent> = Vec::new();
    let eps = 1e-9_f64;
    let min_floor = params.min_score - 1e-3;

    // First cut: highest score in [0, period+tol_period], tie-break to smallest s.
    let first_hi = (params.period + params.tol_period).min(l);
    let mut best_s: Option<usize> = None;
    let mut best_score = min_floor;
    for (s, &v) in score.iter().enumerate().take(first_hi) {
        if v > best_score + eps {
            best_score = v;
            best_s = Some(s);
        }
    }
    if best_s.is_none() {
        // Fall back to scanning the first ~2 periods of the array.
        let scan_hi = (2 * params.period).min(l);
        for (s, &v) in score.iter().enumerate().take(scan_hi) {
            if v > best_score + eps {
                best_score = v;
                best_s = Some(s);
            }
        }
        if best_s.is_none() {
            return (cuts, interp, gap_events);
        }
    }
    cuts.push(best_s.unwrap());
    interp.push(false);

    // Subsequent cuts.
    loop {
        let last = *cuts.last().unwrap();
        if last + params.period < params.tol_period + 1 {
            // numerically impossible to advance — bail.
            break;
        }
        if last + params.period >= l + params.tol_period {
            break;
        }
        let mut advanced = false;
        for k in 1..=params.max_skip {
            let target = last + k * params.period;
            let window_lo = target.saturating_sub(params.tol_period);
            if window_lo >= l {
                break;
            }
            let window_hi = (target + params.tol_period).min(l - 1);
            if window_hi < window_lo {
                continue;
            }
            let mut bs: Option<usize> = None;
            let mut bscore = min_floor;
            for (offset, &v) in score[window_lo..=window_hi].iter().enumerate() {
                let s = window_lo + offset;
                if v > bscore + eps {
                    bscore = v;
                    bs = Some(s);
                } else if (v - bscore).abs() <= eps {
                    if let Some(cur) = bs {
                        let cur_d = (cur as isize - target as isize).abs();
                        let new_d = (s as isize - target as isize).abs();
                        if new_d < cur_d {
                            bs = Some(s);
                        }
                    }
                }
            }
            if let Some(b) = bs {
                let gap = b - last;
                gap_events.push(GapEvent {
                    last_cut: last,
                    next_cut: b,
                    gap_length: gap,
                    k_chosen: k,
                });
                for i in 1..k {
                    let frac = gap as f64 * i as f64 / k as f64;
                    let interp_pos = last + frac.round() as usize;
                    cuts.push(interp_pos);
                    interp.push(true);
                }
                cuts.push(b);
                interp.push(false);
                advanced = true;
                break;
            }
        }
        if !advanced {
            break;
        }
    }

    (cuts, interp, gap_events)
}

/// Cut result for a single array — used by tests and the CLI writer.
#[derive(Debug, Clone)]
pub struct ArrayCut {
    pub cuts: Vec<usize>,
    pub interpolated: Vec<bool>,
    pub score: Vec<f64>,
    pub gap_events: Vec<GapEvent>,
}

/// One-shot: compute score, walk cuts.
pub fn cut_array(seq: &[u8], params: &Params) -> ArrayCut {
    let score = compute_score(seq);
    let (cuts, interpolated, gap_events) = walk_cuts(&score, params);
    ArrayCut { cuts, interpolated, score, gap_events }
}

#[derive(Debug, Clone)]
struct ExceptionRow {
    last_cut: usize,
    next_cut: usize,
    gap_length: usize,
    k_chosen: usize,
    expected: usize,
    residual: i64,
    context_left: String,
    context_right: String,
}

#[derive(Debug, Clone)]
enum Outcome {
    Cut {
        array_id: String,
        seq_len: usize,
        rows: Vec<MonomerRow>,
        n_interp: usize,
        n_direct: usize,
        mean_score: f64,
        median_len: usize,
        mean_len: f64,
        exceptions: Vec<ExceptionRow>,
    },
    Outlier {
        array_id: String,
        seq_len: usize,
        reason: &'static str,
        exceptions: Vec<ExceptionRow>,
    },
}

#[derive(Debug, Clone)]
struct MonomerRow {
    mono_idx: usize,
    start: usize,
    end: usize,
    length: usize,
    interpolated: bool,
    score_left: f64,
    score_right: f64,
    sequence: String,
    /// '+' if cut on the input strand, '-' if on its reverse complement.
    /// Only written when --auto-orient is set.
    orient: char,
}

/// Cut `seq` on the strand it is given. `process_record` wraps this to add the
/// reverse-complement retry.
fn process_record_one_strand(
    name: &str,
    seq: &[u8],
    params: &Params,
    exc_threshold: usize,
    exc_context: usize,
) -> Outcome {
    let l = seq.len();
    if l < params.period {
        return Outcome::Outlier {
            array_id: name.to_string(),
            seq_len: l,
            reason: "too_short",
            exceptions: Vec::new(),
        };
    }
    let ArrayCut { cuts, interpolated, score, gap_events } = cut_array(seq, params);
    let exceptions = build_exceptions(seq, &gap_events, params.period, exc_threshold, exc_context);
    if cuts.is_empty() {
        return Outcome::Outlier {
            array_id: name.to_string(),
            seq_len: l,
            reason: "no_first_cut",
            exceptions,
        };
    }
    if cuts.len() < 2 {
        return Outcome::Outlier {
            array_id: name.to_string(),
            seq_len: l,
            reason: "only_one_cut",
            exceptions,
        };
    }

    let mut rows: Vec<MonomerRow> = Vec::with_capacity(cuts.len() - 1);
    let mut score_sum = 0.0;
    let mut score_n = 0usize;
    for &s in &cuts {
        score_sum += score[s];
        score_n += 1;
    }
    for i in 0..(cuts.len() - 1) {
        let s = cuts[i];
        let e = cuts[i + 1];
        let flag = interpolated[i] || interpolated[i + 1];
        let mono_seq = &seq[s..e];
        rows.push(MonomerRow {
            mono_idx: i,
            start: s,
            end: e,
            length: e - s,
            interpolated: flag,
            score_left: score[s],
            score_right: score[e],
            sequence: String::from_utf8_lossy(mono_seq).to_string(),
            orient: '+',
        });
    }

    let n_interp = rows.iter().filter(|r| r.interpolated).count();
    let n_direct = rows.len() - n_interp;
    let mean_score = if score_n > 0 { score_sum / score_n as f64 } else { 0.0 };
    let mut lens: Vec<usize> = rows.iter().map(|r| r.length).collect();
    lens.sort_unstable();
    let median_len = lens[lens.len() / 2];
    let mean_len = lens.iter().sum::<usize>() as f64 / lens.len() as f64;

    Outcome::Cut {
        array_id: name.to_string(),
        seq_len: l,
        rows,
        n_interp,
        n_direct,
        mean_score,
        median_len,
        mean_len,
        exceptions,
    }
}

/// Cut `seq`, optionally retrying its reverse complement.
///
/// The retry exists because the anchor table is canonical-strand: an array
/// written the other way round produces `no_first_cut` no matter how much
/// satellite it holds. Coordinates from a reverse hit are mapped back to the
/// input frame so the caller never has to flip anything; `mono_idx` keeps the
/// cut order along the matching strand, so a `-` array counts down in `start`.
fn process_record(
    name: &str,
    seq: &[u8],
    params: &Params,
    exc_threshold: usize,
    exc_context: usize,
    auto_orient: bool,
) -> Outcome {
    let fwd = process_record_one_strand(name, seq, params, exc_threshold, exc_context);
    if !auto_orient {
        return fwd;
    }
    if let Outcome::Cut { .. } = fwd {
        return fwd;
    }
    let l = seq.len();
    let rc = crate::monomer::revcomp(seq);
    match process_record_one_strand(name, &rc, params, exc_threshold, exc_context) {
        // both strands refused: say so rather than reporting the forward reason,
        // which would read as "this sequence has no monomer structure"
        Outcome::Outlier { array_id, seq_len, exceptions, .. } => Outcome::Outlier {
            array_id,
            seq_len,
            reason: "no_cut_either_strand",
            exceptions,
        },
        Outcome::Cut {
            array_id,
            seq_len,
            mut rows,
            n_interp,
            n_direct,
            mean_score,
            median_len,
            mean_len,
            exceptions,
        } => {
            for r in rows.iter_mut() {
                let (s, e) = (r.start, r.end);
                r.start = l - e;
                r.end = l - s;
                r.orient = '-';
            }
            Outcome::Cut {
                array_id,
                seq_len,
                rows,
                n_interp,
                n_direct,
                mean_score,
                median_len,
                mean_len,
                exceptions,
            }
        }
    }
}

fn build_exceptions(
    seq: &[u8],
    gap_events: &[GapEvent],
    period: usize,
    threshold: usize,
    context: usize,
) -> Vec<ExceptionRow> {
    let l = seq.len();
    let mut out = Vec::new();
    for ev in gap_events {
        let expected = ev.k_chosen * period;
        let residual = ev.gap_length as i64 - expected as i64;
        if residual.unsigned_abs() as usize <= threshold {
            continue;
        }
        let left_lo = ev.last_cut.saturating_sub(context);
        let left_hi = ev.last_cut.min(l);
        let right_lo = ev.next_cut.min(l);
        let right_hi = (ev.next_cut + context).min(l);
        let context_left = if left_hi > left_lo {
            String::from_utf8_lossy(&seq[left_lo..left_hi]).to_string()
        } else {
            String::new()
        };
        let context_right = if right_hi > right_lo {
            String::from_utf8_lossy(&seq[right_lo..right_hi]).to_string()
        } else {
            String::new()
        };
        out.push(ExceptionRow {
            last_cut: ev.last_cut,
            next_cut: ev.next_cut,
            gap_length: ev.gap_length,
            k_chosen: ev.k_chosen,
            expected,
            residual,
            context_left,
            context_right,
        });
    }
    out
}

pub fn run_from_args(argv: Vec<String>) {
    let args = Args::parse_from(&argv);

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let params = Params {
        period: args.period,
        tol_period: args.tol_period,
        min_score: args.min_score,
        max_skip: args.max_skip,
    };

    // Read all input FASTAs.
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for path in &args.fasta {
        eprintln!("Reading {}...", path);
        let r = read_fasta(path).unwrap_or_else(|e| panic!("{:?}", e));
        eprintln!("  {} records", r.len());
        records.extend(r);
    }
    eprintln!("Total records: {}", records.len());

    // Process in parallel.
    let exc_threshold = args.exception_threshold;
    let exc_context = args.exception_context;
    let outcomes: Vec<Outcome> = records
        .par_iter()
        .map(|(name, seq)| {
            process_record(name, seq, &params, exc_threshold, exc_context, args.auto_orient)
        })
        .collect();

    // Write outputs.
    let mut fout_m = BufWriter::new(File::create(&args.out_monomers).unwrap_or_else(|e| {
        panic!("Cannot create {}: {}", args.out_monomers, e)
    }));
    let mut fout_s = BufWriter::new(File::create(&args.out_summary).unwrap_or_else(|e| {
        panic!("Cannot create {}: {}", args.out_summary, e)
    }));
    let mut fout_o = BufWriter::new(File::create(&args.out_outliers).unwrap_or_else(|e| {
        panic!("Cannot create {}: {}", args.out_outliers, e)
    }));
    let mut fout_e = BufWriter::new(File::create(&args.out_exceptions).unwrap_or_else(|e| {
        panic!("Cannot create {}: {}", args.out_exceptions, e)
    }));

    writeln!(
        fout_m,
        "array_id\tmono_idx\tstart\tend\tlength\tinterpolated\tscore_left\tscore_right\tsequence{}",
        if args.auto_orient { "\torient" } else { "" }
    )
    .unwrap();
    writeln!(
        fout_s,
        "array_id\tlen_bp\tn_monomers\tn_interpolated\tn_direct\tmean_score\tmedian_length\tmean_length"
    )
    .unwrap();
    writeln!(fout_o, "array_id\tlen_bp\treason").unwrap();
    writeln!(
        fout_e,
        "array_id\tlast_cut\tnext_cut\tgap_length\tk_chosen\texpected_length\tresidual_bp\tabs_residual_bp\tcontext_left\tcontext_right"
    )
    .unwrap();

    let mut n_records = 0usize;
    let mut n_with = 0usize;
    let mut n_outliers = 0usize;
    let mut total_mono = 0usize;
    let mut total_direct = 0usize;
    let mut total_exceptions = 0usize;
    // too_short, no_first_cut, only_one_cut, no_cut_either_strand
    let mut breakdown = [0usize; 4];

    for o in &outcomes {
        n_records += 1;
        let (array_id_ref, exceptions_ref) = match o {
            Outcome::Cut { array_id, exceptions, .. } => (array_id.as_str(), exceptions),
            Outcome::Outlier { array_id, exceptions, .. } => (array_id.as_str(), exceptions),
        };
        for ex in exceptions_ref {
            writeln!(
                fout_e,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:+}\t{}\t{}\t{}",
                array_id_ref,
                ex.last_cut,
                ex.next_cut,
                ex.gap_length,
                ex.k_chosen,
                ex.expected,
                ex.residual,
                ex.residual.unsigned_abs(),
                ex.context_left,
                ex.context_right
            )
            .unwrap();
            total_exceptions += 1;
        }
        match o {
            Outcome::Cut {
                array_id, seq_len, rows, n_interp, n_direct, mean_score, median_len, mean_len,
                exceptions: _,
            } => {
                n_with += 1;
                total_mono += rows.len();
                total_direct += n_direct;
                for r in rows {
                    writeln!(
                        fout_m,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}{}",
                        array_id,
                        r.mono_idx,
                        r.start,
                        r.end,
                        r.length,
                        if r.interpolated { 1 } else { 0 },
                        r.score_left,
                        r.score_right,
                        r.sequence,
                        if args.auto_orient {
                            format!("\t{}", r.orient)
                        } else {
                            String::new()
                        }
                    )
                    .unwrap();
                }
                writeln!(
                    fout_s,
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.1}",
                    array_id,
                    seq_len,
                    rows.len(),
                    n_interp,
                    n_direct,
                    mean_score,
                    median_len,
                    mean_len
                )
                .unwrap();
            }
            Outcome::Outlier { array_id, seq_len, reason, exceptions: _ } => {
                n_outliers += 1;
                match *reason {
                    "too_short" => breakdown[0] += 1,
                    "no_first_cut" => breakdown[1] += 1,
                    "only_one_cut" => breakdown[2] += 1,
                    "no_cut_either_strand" => breakdown[3] += 1,
                    // an unclassified reason would vanish from the summary and
                    // the totals would stop adding up -- refuse instead
                    other => panic!("unhandled outlier reason {:?}", other),
                }
                writeln!(fout_o, "{}\t{}\t{}", array_id, seq_len, reason).unwrap();
            }
        }
    }

    eprintln!("records processed:     {}", n_records);
    eprintln!("records with monomers: {}", n_with);
    eprintln!(
        "outliers:              {}  (too_short={}, no_first_cut={}, only_one_cut={}, no_cut_either_strand={})",
        n_outliers, breakdown[0], breakdown[1], breakdown[2], breakdown[3]
    );
    eprintln!("total monomers:        {}", total_mono);
    let pct = if total_mono > 0 {
        100.0 * total_direct as f64 / total_mono as f64
    } else {
        0.0
    };
    eprintln!("direct (non-interp):   {} ({:.2}%)", total_direct, pct);
    eprintln!(
        "exceptions (|residual| > {} bp): {}",
        exc_threshold, total_exceptions
    );
}
