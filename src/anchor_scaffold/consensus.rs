//! Star consensus + EM driver for the anchor-scaffold MSA pipeline.
//!
//! Brief: `docs/anchor_scaffold_msa_brief.md` § Core modules / 4 (star_consensus)
//! and § Parameters (EM rounds / convergence).
//!
//! Pipeline:
//!   1. `process_monomer` — for one monomer + a current consensus, runs the
//!      detect → chain → extract_pieces → ond_align pipeline and returns a
//!      `MonomerOutcome` classified into STANDARD / OUTLIER / EXCEPTION with
//!      the brief's closed-enum cause.
//!   2. `update_tally` — accumulates a STANDARD outcome's edit-script
//!      contributions into a per-column `ColumnTally`. Outliers/exceptions
//!      are skipped (the partial alignment of outliers can be noisy and is
//!      emitted separately).
//!   3. `rebuild_consensus` — selects the highest-count base per column
//!      (ties: lexicographic A<C<G<T; gap is NOT a candidate base — single-
//!      base, byte-stable consensus per ASM).
//!   4. `run_em` — loops 1-3 up to `max_rounds` or until fewer than
//!      `convergence_threshold` fraction of columns change.
//!
//! Insertion handling: `Ins{col, base}` ops have no canonical column residency
//! (they sit BETWEEN columns), so they are recorded in the JSONL edit script
//! but do NOT contribute to the column-level tally. Inserts therefore can't
//! shift the consensus length. For the smoke fixture this is correct by
//! construction (counts/columns are defined in canonical-column space).

use std::collections::HashMap;

use serde::Serialize;

use crate::align::{ond_align, EditScript, Op};
use crate::anchor_detect::{detect_anchors, detect_anchors_chained, AnchorHit, PanelEntry};
use crate::anchor_scaffold::{extract_pieces, Chain, PieceKind};

/// Output stream a monomer is routed to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Stream {
    Standard,
    Outlier,
    Exception,
}

/// Closed enum of causes per the brief's I/O contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum Cause {
    Ok,
    TooFewAnchors,
    SegmentOverBand,
    GapUnalignable,
    NoncanonicalAnchorOrder,
}

/// JSON I/O record for one piece-level event. Only emitted for outliers.
#[derive(Debug, Clone, Serialize)]
pub struct PieceReport {
    pub span: (usize, usize),
    pub canonical_span: (usize, usize),
    /// `script.len()` when alignment succeeded within band; `None` when the
    /// piece exceeded the band (true D is unknown; `Over band` is the signal).
    #[serde(rename = "D", skip_serializing_if = "Option::is_none")]
    pub d: Option<usize>,
    pub status: PieceStatus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum PieceStatus {
    Ok,
    OverBand,
}

/// Per-monomer outcome — one JSONL line in standard/outliers/exceptions.
#[derive(Debug, Clone, Serialize)]
pub struct MonomerOutcome {
    pub id: String,
    pub stream: Stream,
    pub cause: Cause,
    pub n_present: usize,
    pub anchors: Vec<AnchorHit>,
    /// Empty for exceptions; partial for outliers (only succeeded pieces'
    /// edits); complete for standard. Columns are in CONSENSUS coordinates
    /// (already offset from the piece's canonical_span.0).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub edit_script: EditScript,
    /// Only populated for outliers, listing which pieces failed alignment.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub pieces: Vec<PieceReport>,
}

/// Process one monomer against the current consensus and panel. The result
/// classifies the monomer into a stream + cause and (for standard/outlier)
/// returns the edit script ready to feed `update_tally`.
pub fn process_monomer(
    id: String,
    monomer: &[u8],
    consensus: &[u8],
    panel: &[PanelEntry],
    min_anchors: usize,
    max_hd: usize,
    chained: bool,
) -> MonomerOutcome {
    let hits = if chained {
        detect_anchors_chained(panel, monomer, max_hd)
    } else {
        detect_anchors(panel, monomer, max_hd)
    };
    let chain = Chain::from_hits(hits);

    if !chain.has_enough_anchors(min_anchors) {
        return MonomerOutcome {
            id,
            stream: Stream::Exception,
            cause: Cause::TooFewAnchors,
            n_present: chain.n_primary,
            anchors: chain.hits,
            edit_script: Vec::new(),
            pieces: Vec::new(),
        };
    }

    if chain.noncanonical_order {
        return MonomerOutcome {
            id,
            stream: Stream::Outlier,
            cause: Cause::NoncanonicalAnchorOrder,
            n_present: chain.n_primary,
            anchors: chain.hits,
            edit_script: Vec::new(),
            pieces: Vec::new(),
        };
    }

    let pieces = extract_pieces(&chain, panel);
    let mut edit_script = EditScript::new();
    let mut failed: Vec<PieceReport> = Vec::new();
    let mut outlier_cause: Option<Cause> = None;

    for piece in &pieces {
        let cons_stretch = piece.consensus_seq(consensus);
        let mono_piece = piece.monomer_seq(monomer);
        let band = piece.band();
        match ond_align(cons_stretch, mono_piece, band) {
            None => {
                let cause = match piece.kind {
                    PieceKind::InterAnchor => Cause::SegmentOverBand,
                    PieceKind::MissingAnchorGap { .. } => Cause::GapUnalignable,
                };
                if outlier_cause.is_none() {
                    outlier_cause = Some(cause);
                }
                failed.push(PieceReport {
                    span: piece.monomer_span,
                    canonical_span: piece.canonical_span,
                    d: None,
                    status: PieceStatus::OverBand,
                });
            }
            Some(script) => {
                let offset = piece.canonical_span.0;
                for op in script {
                    let mapped = match op {
                        Op::Sub { col, to } => Op::Sub { col: col + offset, to },
                        Op::Del { col } => Op::Del { col: col + offset },
                        Op::Ins { col, base } => Op::Ins { col: col + offset, base },
                    };
                    edit_script.push(mapped);
                }
            }
        }
    }

    match outlier_cause {
        Some(cause) => MonomerOutcome {
            id,
            stream: Stream::Outlier,
            cause,
            n_present: chain.n_primary,
            anchors: chain.hits,
            edit_script,
            pieces: failed,
        },
        None => MonomerOutcome {
            id,
            stream: Stream::Standard,
            cause: Cause::Ok,
            n_present: chain.n_primary,
            anchors: chain.hits,
            edit_script,
            pieces: Vec::new(),
        },
    }
}

/// Per-canonical-column nucleotide + gap counts.
///
/// Indexing: `counts[col][0..=4]` = `[A, C, G, T, gap]`.
#[derive(Debug, Clone)]
pub struct ColumnTally {
    pub counts: Vec<[u32; 5]>,
}

const I_A: usize = 0;
const I_C: usize = 1;
const I_G: usize = 2;
const I_T: usize = 3;
const I_GAP: usize = 4;

impl ColumnTally {
    pub fn new(consensus_len: usize) -> Self {
        Self { counts: vec![[0u32; 5]; consensus_len] }
    }

    pub fn total_at(&self, col: usize) -> u32 {
        self.counts[col].iter().sum()
    }
}

fn base_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(I_A),
        b'C' => Some(I_C),
        b'G' => Some(I_G),
        b'T' => Some(I_T),
        _ => None,
    }
}

/// Canonical columns this monomer aligned and may vote on / be displayed at.
/// STANDARD = all columns; OUTLIER = the anchored extent
/// `[first_anchor ..= last_anchor]` minus the spans of pieces that failed the
/// band (and the flanks); EXCEPTION = none.
pub fn covered_mask(outcome: &MonomerOutcome, n: usize) -> Vec<bool> {
    match outcome.stream {
        Stream::Standard => vec![true; n],
        Stream::Exception => vec![false; n],
        Stream::Outlier => {
            let mut mask = vec![false; n];
            if let (Some(first), Some(last)) = (
                outcome.anchors.iter().map(|a| a.canonical_pos).min(),
                outcome.anchors.iter().map(|a| a.canonical_pos).max(),
            ) {
                let hi = last.min(n.saturating_sub(1));
                for c in mask.iter_mut().take(hi + 1).skip(first) {
                    *c = true;
                }
                for p in &outcome.pieces {
                    let (s, e) = p.canonical_span;
                    for c in mask.iter_mut().take(e.min(n)).skip(s) {
                        *c = false;
                    }
                }
            }
            mask
        }
    }
}

/// Reconstruct one monomer's aligned row over the canonical columns for a
/// human-viewable MSA: consensus base by default, `Sub` override, `-` for `Del`,
/// `.` for columns not covered by an aligned piece (failed pieces / flanks).
/// Insertions are between-column events and are omitted from the fixed grid.
pub fn aligned_row(outcome: &MonomerOutcome, consensus: &[u8]) -> Vec<u8> {
    let n = consensus.len();
    let mut row = consensus.to_vec();
    for op in &outcome.edit_script {
        match *op {
            Op::Sub { col, to } => {
                if col < n {
                    row[col] = to as u8;
                }
            }
            Op::Del { col } => {
                if col < n {
                    row[col] = b'-';
                }
            }
            Op::Ins { .. } => {}
        }
    }
    let covered = covered_mask(outcome, n);
    for (col, cell) in row.iter_mut().enumerate() {
        if !covered[col] {
            *cell = b'.';
        }
    }
    row
}

/// Accumulate a monomer's per-column contribution into `tally`, over the columns
/// it actually aligned (see `covered_mask`). Per column: default = `consensus[col]`,
/// `Sub{col,to}` overrides to `to`, `Del{col}` overrides to gap, `Ins` is ignored
/// (between-column). Exceptions contribute nothing.
pub fn update_tally(tally: &mut ColumnTally, outcome: &MonomerOutcome, consensus: &[u8]) {
    // Exceptions (no usable anchor chain) contribute nothing.
    if outcome.stream == Stream::Exception {
        return;
    }
    assert_eq!(
        tally.counts.len(),
        consensus.len(),
        "tally length must match consensus length"
    );

    // Bytes shorter than allocating a full per-column override map: use
    // u8 with sentinel 0xFF for "no override; use consensus[col]".
    const NONE: u8 = 0xFF;
    const GAP: u8 = b'-';
    let mut overrides = vec![NONE; consensus.len()];

    for op in &outcome.edit_script {
        match *op {
            Op::Sub { col, to } => {
                if col < overrides.len() {
                    overrides[col] = to as u8;
                }
            }
            Op::Del { col } => {
                if col < overrides.len() {
                    overrides[col] = GAP;
                }
            }
            Op::Ins { .. } => { /* between-column event; not tallied */ }
        }
    }

    // Partial-tally: a monomer votes only on the columns it actually aligned
    // (STANDARD = all; OUTLIER = anchored extent minus failed-piece spans).
    // Voting the consensus base on an unaligned column would fabricate data.
    // Realises "align between the present lattice anchors" — a divergent gap
    // costs only its own columns (BTN-scaffold-band-too-tight-starves-column-profile).
    let covered = covered_mask(outcome, consensus.len());

    for (col, &ovr) in overrides.iter().enumerate() {
        if !covered[col] {
            continue;
        }
        let b = if ovr == NONE { consensus[col] } else { ovr };
        if b == GAP {
            tally.counts[col][I_GAP] = tally.counts[col][I_GAP].saturating_add(1);
        } else if let Some(idx) = base_idx(b) {
            tally.counts[col][idx] = tally.counts[col][idx].saturating_add(1);
        }
        // Non-ACGT consensus bytes (e.g. N) are silently ignored.
    }
}

/// Build a new consensus from the tally. For each column the highest base
/// count wins; ties break by lexicographic A<C<G<T (i.e. by index, since
/// `counts` is stored in that order). A column with zero base counts (only
/// gaps, or empty) keeps the previous consensus byte — the consensus length
/// is byte-stable per the ASM constraint.
pub fn rebuild_consensus(tally: &ColumnTally, prev: &[u8]) -> Vec<u8> {
    let mut next = prev.to_vec();
    for (col, counts) in tally.counts.iter().enumerate() {
        // Only A/C/G/T are candidates for the consensus byte; `gap` is
        // excluded so the consensus stays single-base.
        let mut best_idx = 0usize;
        let mut best_count = counts[0];
        for (i, &c) in counts.iter().enumerate().take(4).skip(1) {
            if c > best_count {
                best_count = c;
                best_idx = i;
            }
        }
        if best_count == 0 {
            // No standard monomers voted at this column — keep prev.
            continue;
        }
        next[col] = match best_idx {
            I_A => b'A',
            I_C => b'C',
            I_G => b'G',
            I_T => b'T',
            _ => prev[col], // unreachable given the take(4) bound
        };
    }
    next
}

/// EM parameters; all defaults mirror the brief's ASM-msa-parameter-set.
#[derive(Debug, Clone, Copy)]
pub struct EmParams {
    /// Maximum EM rounds (brief: ≤3).
    pub max_rounds: usize,
    /// Convergence: stop when fewer than this fraction of columns change.
    /// Brief: < 0.1% (i.e. 0.001).
    pub convergence_threshold: f64,
    /// `min_anchors` for the EXCEPTION gate (brief: 6 of 10 primaries).
    pub min_anchors: usize,
    /// Global cap on per-slot HD threshold (brief: 0 for smoke).
    pub max_hd: usize,
    /// Use relative-chaining anchor detection (floats with indel drift) instead
    /// of fixed absolute-position windows. See `detect_anchors_chained` +
    /// BTN-scaffold-panel-fixed-windows-dont-transfer.
    pub chained: bool,
}

impl Default for EmParams {
    fn default() -> Self {
        Self { max_rounds: 3, convergence_threshold: 0.001, min_anchors: 6, max_hd: 0, chained: false }
    }
}

/// Telemetry for one EM round.
#[derive(Debug, Clone)]
pub struct EmRoundReport {
    pub round_idx: usize,
    pub n_columns_changed: usize,
    pub stream_counts: HashMap<Stream, usize>,
}

#[derive(Debug, Clone)]
pub struct EmResult {
    pub final_consensus: Vec<u8>,
    pub final_outcomes: Vec<MonomerOutcome>,
    pub rounds: Vec<EmRoundReport>,
}

/// Run the EM driver: process all monomers against the current consensus,
/// rebuild the consensus from STANDARD contributions, repeat until
/// convergence or `max_rounds`.
///
/// Returns the final consensus, the final round's per-monomer outcomes (so
/// the caller can emit JSONL), and a per-round telemetry log.
pub fn run_em(
    seed_consensus: &[u8],
    monomers: &[(String, Vec<u8>)],
    panel: &[PanelEntry],
    params: &EmParams,
) -> EmResult {
    let mut current = seed_consensus.to_vec();
    let mut rounds: Vec<EmRoundReport> = Vec::new();
    let mut last_outcomes: Vec<MonomerOutcome> = Vec::new();

    for round_idx in 0..params.max_rounds {
        let outcomes: Vec<MonomerOutcome> = monomers
            .iter()
            .map(|(id, seq)| {
                process_monomer(
                    id.clone(),
                    seq,
                    &current,
                    panel,
                    params.min_anchors,
                    params.max_hd,
                    params.chained,
                )
            })
            .collect();

        let mut tally = ColumnTally::new(current.len());
        for o in &outcomes {
            update_tally(&mut tally, o, &current);
        }

        let next_consensus = rebuild_consensus(&tally, &current);
        let n_changed = current
            .iter()
            .zip(&next_consensus)
            .filter(|(a, b)| a != b)
            .count();

        let mut stream_counts: HashMap<Stream, usize> = HashMap::new();
        for o in &outcomes {
            *stream_counts.entry(o.stream).or_insert(0) += 1;
        }

        rounds.push(EmRoundReport {
            round_idx,
            n_columns_changed: n_changed,
            stream_counts,
        });
        last_outcomes = outcomes;

        let change_rate = if current.is_empty() {
            0.0
        } else {
            n_changed as f64 / current.len() as f64
        };
        current = next_consensus;
        if change_rate < params.convergence_threshold {
            break;
        }
    }

    EmResult {
        final_consensus: current,
        final_outcomes: last_outcomes,
        rounds,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_detect::SCAFFOLD_PANEL;

    const CONS_171: &[u8] = b"CTTTGTGATGTCATCATTCCATGCCAGAGTGCATGCATGCCTTTTTGCATGCATGCATGGAAACATGCATGCATGCATGCATGCATGTGGACATGCATGCATGCATGCATGCATGCATGGAAGCATGCATGCATGCATGCAAAACTGCACAGATGCATTCACAGAAGCATG";

    #[test]
    fn clean_monomer_is_standard_with_empty_script() {
        let o = process_monomer(
            "M01_clean".into(),
            CONS_171,
            CONS_171,
            SCAFFOLD_PANEL,
            6,
            0,
            false,
        );
        assert_eq!(o.stream, Stream::Standard);
        assert_eq!(o.cause, Cause::Ok);
        assert!(o.edit_script.is_empty(), "{:?}", o.edit_script);
        assert!(o.pieces.is_empty());
        assert_eq!(o.n_present, 10);
    }

    #[test]
    fn substitution_monomer_yields_single_sub() {
        // Mirror M07_sub12 from the fixture: sub at col 12 (A -> T).
        let mut m = CONS_171.to_vec();
        m[12] = b'T';
        let o = process_monomer("sub12".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Standard);
        assert_eq!(o.edit_script, vec![Op::Sub { col: 12, to: 'T' }]);
    }

    #[test]
    fn deletion_monomer_yields_single_del() {
        let mut m = CONS_171.to_vec();
        m.remove(12);
        let o = process_monomer("del12".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Standard);
        // ond_align reports a Del at a canonical column near 12; the precise
        // column can land at 11 or 12 depending on traceback tie-break
        // against the surrounding filler. Either is correct as long as it's
        // a single Del op in the right neighbourhood.
        assert_eq!(o.edit_script.len(), 1, "{:?}", o.edit_script);
        match o.edit_script[0] {
            Op::Del { col } => assert!(
                (11..=14).contains(&col),
                "del col {col} out of expected range"
            ),
            ref op => panic!("expected Del, got {op:?}"),
        }
    }

    #[test]
    fn outlier_segment_yields_outlier_with_failed_piece() {
        // Randomize 91..117 deterministically (no need to use rng — same
        // pattern as the fixture's outlier monomers in spirit).
        let mut m = CONS_171.to_vec();
        for (i, p) in (91..117).enumerate() {
            m[p] = b"ACGT"[i % 3]; // intentionally skewed -> guaranteed D > band
        }
        let o = process_monomer("outlier".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Outlier);
        assert_eq!(o.cause, Cause::SegmentOverBand);
        assert_eq!(o.pieces.len(), 1, "{:?}", o.pieces);
        let p = &o.pieces[0];
        assert_eq!(p.canonical_span, (91, 117));
        assert!(matches!(p.status, PieceStatus::OverBand));
    }

    #[test]
    fn exception_monomer_short_circuits() {
        // Destroy 6 primaries (mirror exception_destroy_set logic).
        let mut m = CONS_171.to_vec();
        for (pos, k) in [(14_usize, 5_usize), (24, 5), (40, 5), (59, 5), (86, 5), (117, 5)] {
            for j in 0..k {
                m[pos + j] = b"ATGC"[(pos + j) % 4];
            }
        }
        let o = process_monomer("exc".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Exception);
        assert_eq!(o.cause, Cause::TooFewAnchors);
        assert!(o.edit_script.is_empty());
        assert_eq!(o.n_present, 4);
    }

    #[test]
    fn knockout_stays_standard_with_small_edit_script() {
        // Knock out slot4 GTGGA@86 — across-missing-anchor piece is the
        // (slot3 -> slot5) span. The 4 substitutions inside the filler
        // replacement are within band=14.
        let mut m = CONS_171.to_vec();
        for j in 0..5 {
            m[86 + j] = b"ATGC"[(86 + j) % 4];
        }
        let o = process_monomer("ko86".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Standard);
        assert_eq!(o.cause, Cause::Ok);
        assert_eq!(o.n_present, 9);
        // The edit_script's columns must all fall in the knocked-out anchor's
        // canonical span (86..91).
        for op in &o.edit_script {
            let col = match op {
                Op::Sub { col, .. } | Op::Del { col } | Op::Ins { col, .. } => *col,
            };
            assert!((86..91).contains(&col), "edit at col {col} outside KO zone");
        }
    }

    #[test]
    fn tally_increments_consensus_bases_for_clean() {
        let o = process_monomer("c".into(), CONS_171, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        let mut t = ColumnTally::new(CONS_171.len());
        update_tally(&mut t, &o, CONS_171);
        // Every column receives exactly one base vote — the consensus byte.
        for (col, c) in t.counts.iter().enumerate() {
            let total: u32 = c.iter().sum();
            assert_eq!(total, 1, "col {col}: {c:?}");
            // The voted base must be A/C/G/T (no gap from a clean monomer).
            assert_eq!(c[I_GAP], 0);
        }
    }

    #[test]
    fn tally_skips_outliers_and_exceptions() {
        let mut t = ColumnTally::new(CONS_171.len());
        let outlier = MonomerOutcome {
            id: "x".into(),
            stream: Stream::Outlier,
            cause: Cause::SegmentOverBand,
            n_present: 10,
            anchors: Vec::new(),
            edit_script: vec![Op::Sub { col: 0, to: 'A' }],
            pieces: Vec::new(),
        };
        let exception = MonomerOutcome {
            id: "y".into(),
            stream: Stream::Exception,
            cause: Cause::TooFewAnchors,
            n_present: 4,
            anchors: Vec::new(),
            edit_script: Vec::new(),
            pieces: Vec::new(),
        };
        update_tally(&mut t, &outlier, CONS_171);
        update_tally(&mut t, &exception, CONS_171);
        assert!(t.counts.iter().all(|c| c.iter().sum::<u32>() == 0));
    }

    #[test]
    fn rebuild_picks_majority_base_lex_tiebreak() {
        let mut t = ColumnTally::new(3);
        // col 0: A wins outright
        t.counts[0] = [5, 1, 0, 0, 0];
        // col 1: A and C tied — A wins by lex (index 0 < index 1)
        t.counts[1] = [3, 3, 0, 0, 0];
        // col 2: gap dominates but isn't a candidate -> falls back to G (the
        // only non-zero base)
        t.counts[2] = [0, 0, 1, 0, 10];
        let prev = b"TTT";
        let next = rebuild_consensus(&t, prev);
        assert_eq!(next, b"AAG");
    }

    #[test]
    fn rebuild_keeps_prev_when_no_votes() {
        let t = ColumnTally::new(4);
        let prev = b"ACGT";
        let next = rebuild_consensus(&t, prev);
        assert_eq!(next, prev);
    }

    #[test]
    fn em_converges_in_one_round_on_clean_data() {
        let mons = vec![("m1".into(), CONS_171.to_vec()), ("m2".into(), CONS_171.to_vec())];
        let r = run_em(CONS_171, &mons, SCAFFOLD_PANEL, &EmParams::default());
        assert_eq!(r.final_consensus, CONS_171);
        assert!(r.rounds.len() == 1, "expected 1 round, got {}", r.rounds.len());
        assert_eq!(r.rounds[0].n_columns_changed, 0);
        assert_eq!(r.rounds[0].stream_counts.get(&Stream::Standard), Some(&2));
    }

    #[test]
    fn outlier_partial_tally_covers_aligned_columns_only() {
        // Scramble the slot4..slot5 inter-anchor segment (cols 91..117) so that
        // piece fails the band → OUTLIER, but the rest of the lattice aligns.
        // Partial-tally must count the aligned columns and skip the failed span.
        let mut m = CONS_171.to_vec();
        for (i, b) in m.iter_mut().enumerate().take(117).skip(91) {
            *b = b"ACGT"[i % 4];
        }
        let o = process_monomer("out".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Outlier, "{:?}", o.cause);
        let mut tally = ColumnTally::new(CONS_171.len());
        update_tally(&mut tally, &o, CONS_171);
        assert!(tally.total_at(14) >= 1, "aligned anchor column 14 must be tallied");
        assert_eq!(tally.total_at(100), 0, "failed-piece column 100 must be skipped");
    }

    #[test]
    fn aligned_row_clean_equals_consensus_and_outlier_dots_failed_span() {
        // clean monomer -> row identical to consensus
        let o = process_monomer("c".into(), CONS_171, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(aligned_row(&o, CONS_171), CONS_171.to_vec());

        // outlier with scrambled slot4..slot5 segment -> '.' inside failed span,
        // real base at an aligned anchor column.
        let mut m = CONS_171.to_vec();
        for (i, b) in m.iter_mut().enumerate().take(117).skip(91) {
            *b = b"ACGT"[i % 4];
        }
        let o = process_monomer("o".into(), &m, CONS_171, SCAFFOLD_PANEL, 6, 0, false);
        assert_eq!(o.stream, Stream::Outlier);
        let row = aligned_row(&o, CONS_171);
        assert_eq!(row.len(), CONS_171.len());
        assert_eq!(row[100], b'.', "failed-piece column must be '.'");
        assert_ne!(row[14], b'.', "aligned column must carry a base");
    }
}
