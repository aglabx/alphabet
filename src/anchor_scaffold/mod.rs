//! Chain + piece extraction for the anchor-scaffold MSA pipeline.
//!
//! Submodules:
//!   * (this file) — Chain, Piece, extract_pieces
//!   * `consensus` — ColumnTally, process_monomer, run_em (star MSA + EM)
//!
//! Brief: `docs/anchor_scaffold_msa_brief.md` § Core modules / 1-2.
//!
//! Stages handled here:
//!   * `Chain::from_hits` — orders present anchors, counts primaries, flags
//!     non-canonical order. Caller checks `has_enough_anchors(min_anchors)`
//!     to route to the EXCEPTIONS stream.
//!   * `extract_pieces` — between every pair of consecutive present anchors
//!     emits one `Piece`. If the panel slots between them are not adjacent
//!     in the canonical order, the piece is tagged `MissingAnchorGap` with
//!     the indices of the skipped panel entries (later: cause-of-outlier
//!     when D > band on this longer piece).
//!
//! Pieces are pure metadata: they carry monomer/canonical spans and panel
//! anchor handles. Sequence slicing is done by helpers on the monomer/
//! consensus buffers the caller already owns.

pub mod consensus;
pub use consensus::*;

use crate::anchor_detect::{AnchorHit, PanelEntry, Role};

/// Ordered list of present anchors for one monomer + derived facts.
#[derive(Debug, Clone)]
pub struct Chain {
    /// Hits sorted by `found_pos` (== sorted by canonical column for any
    /// canonical-order chain).
    pub hits: Vec<AnchorHit>,
    /// Number of `Primary`-role hits. Compared against `min_anchors` to
    /// decide EXCEPTION routing.
    pub n_primary: usize,
    /// True if at least one consecutive pair has decreasing panel index
    /// (i.e. a later panel slot precedes an earlier one in monomer space).
    /// Routes the monomer to OUTLIERS with `cause = NONCANONICAL_ANCHOR_ORDER`.
    pub noncanonical_order: bool,
}

impl Chain {
    /// Build a chain from `detect_anchors` output.
    ///
    /// `hits` must already be sorted by `found_pos` (which is what
    /// `detect_anchors` returns).
    pub fn from_hits(hits: Vec<AnchorHit>) -> Self {
        let n_primary = hits.iter().filter(|h| h.role == Role::Primary).count();
        let noncanonical_order = hits
            .windows(2)
            .any(|w| w[0].panel_idx >= w[1].panel_idx);
        Self { hits, n_primary, noncanonical_order }
    }

    /// True when `n_primary >= min_anchors` — the gate between
    /// STANDARD/OUTLIER and EXCEPTION.
    pub fn has_enough_anchors(&self, min_anchors: usize) -> bool {
        self.n_primary >= min_anchors
    }
}

/// What kind of piece sits between two present anchors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PieceKind {
    /// The bounding anchors are adjacent panel slots — straight inter-anchor
    /// filler.
    InterAnchor,
    /// At least one canonical panel slot is missing between the bounding
    /// anchors. The skipped slots are reported by their panel index so the
    /// outlier emitter can name which columns were spanned.
    MissingAnchorGap { skipped_slots: Vec<usize> },
}

/// One alignable piece of a monomer: the region between two consecutive
/// present anchors, in both monomer and canonical-column coordinates.
#[derive(Debug, Clone)]
pub struct Piece {
    /// Half-open range in the monomer (`start..end`), excluding the bounding
    /// anchor k-mers themselves.
    pub monomer_span: (usize, usize),
    /// Half-open range in canonical columns. `canonical_span.1 -
    /// canonical_span.0` is the canonical "piece length" the brief refers to
    /// when computing `band = max(4, ceil(0.25 * piece_len))`.
    pub canonical_span: (usize, usize),
    /// Panel index of the left flanking anchor.
    pub left_anchor_idx: usize,
    /// Panel index of the right flanking anchor.
    pub right_anchor_idx: usize,
    pub kind: PieceKind,
}

impl Piece {
    /// Banded edit-distance budget per the brief:
    /// `band = max(4, ceil(0.25 * canonical_piece_len))`.
    pub fn band(&self) -> usize {
        let len = self.canonical_span.1.saturating_sub(self.canonical_span.0);
        4.max(((len as f64) * 0.25).ceil() as usize)
    }

    /// Borrow the piece's bytes from the monomer buffer. Returns an empty
    /// slice for degenerate spans (e.g. anchors that overlap after large
    /// indel shifts).
    pub fn monomer_seq<'m>(&self, monomer: &'m [u8]) -> &'m [u8] {
        let (s, e) = self.monomer_span;
        if e <= s || e > monomer.len() {
            return &[];
        }
        &monomer[s..e]
    }

    /// Borrow the piece's bytes from the consensus buffer. Same degenerate
    /// handling as `monomer_seq`.
    pub fn consensus_seq<'c>(&self, consensus: &'c [u8]) -> &'c [u8] {
        let (s, e) = self.canonical_span;
        if e <= s || e > consensus.len() {
            return &[];
        }
        &consensus[s..e]
    }
}

/// For each pair of consecutive present anchors in `chain`, emit one
/// `Piece`. The chain is assumed canonical (`!chain.noncanonical_order`); on
/// a non-canonical chain the caller should short-circuit to OUTLIER instead.
///
/// Returns an empty vector when the chain has fewer than two hits — there's
/// nothing to align between.
pub fn extract_pieces(chain: &Chain, panel: &[PanelEntry]) -> Vec<Piece> {
    let mut pieces = Vec::with_capacity(chain.hits.len().saturating_sub(1));
    for w in chain.hits.windows(2) {
        let left = &w[0];
        let right = &w[1];
        let left_kmer_len = panel[left.panel_idx].kmer.len();

        // Spans excluding the bounding anchor k-mers themselves. We clamp
        // monomer_end >= monomer_start so a degenerate (anchor-overlap) span
        // becomes empty, not negative.
        let monomer_start = left.found_pos + left_kmer_len;
        let monomer_end = right.found_pos.max(monomer_start);
        let canonical_start = left.canonical_pos + left_kmer_len;
        let canonical_end = right.canonical_pos.max(canonical_start);

        let skipped: Vec<usize> = (left.panel_idx + 1..right.panel_idx).collect();
        let kind = if skipped.is_empty() {
            PieceKind::InterAnchor
        } else {
            PieceKind::MissingAnchorGap { skipped_slots: skipped }
        };

        pieces.push(Piece {
            monomer_span: (monomer_start, monomer_end),
            canonical_span: (canonical_start, canonical_end),
            left_anchor_idx: left.panel_idx,
            right_anchor_idx: right.panel_idx,
            kind,
        });
    }
    pieces
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchor_detect::{detect_anchors, SCAFFOLD_PANEL};

    const CONS_171: &[u8] = b"CTTTGTGATGTCATCATTCCATGCCAGAGTGCATGCATGCCTTTTTGCATGCATGCATGGAAACATGCATGCATGCATGCATGCATGTGGACATGCATGCATGCATGCATGCATGCATGGAAGCATGCATGCATGCATGCAAAACTGCACAGATGCATTCACAGAAGCATG";

    fn build_chain(monomer: &[u8]) -> Chain {
        Chain::from_hits(detect_anchors(SCAFFOLD_PANEL, monomer, 0))
    }

    #[test]
    fn empty_chain_safe() {
        let c = Chain::from_hits(Vec::new());
        assert_eq!(c.n_primary, 0);
        assert!(!c.noncanonical_order);
        assert!(!c.has_enough_anchors(1));
    }

    #[test]
    fn canonical_consensus_chain() {
        let c = build_chain(CONS_171);
        assert_eq!(c.hits.len(), 11);
        assert_eq!(c.n_primary, 10);
        assert!(!c.noncanonical_order);
        assert!(c.has_enough_anchors(6));
        assert!(!c.has_enough_anchors(11));
    }

    #[test]
    fn canonical_consensus_yields_ten_inter_anchor_pieces() {
        let c = build_chain(CONS_171);
        let pieces = extract_pieces(&c, SCAFFOLD_PANEL);
        assert_eq!(pieces.len(), 10);
        for p in &pieces {
            assert!(matches!(p.kind, PieceKind::InterAnchor), "{p:?}");
        }
        // First piece: boundary (CTTTGTGATGT, len 11, pos 0) -> slot0 (CATTC@14).
        // monomer_span = canonical_span = (11, 14) — 3 bp filler.
        assert_eq!(pieces[0].monomer_span, (11, 14));
        assert_eq!(pieces[0].canonical_span, (11, 14));
        assert_eq!(pieces[0].band(), 4); // max(4, ceil(0.25*3)) = 4
        // Last piece: slot8 (CATTC@155, len 5) -> slot9 (CAGAA@161).
        // ends at 160 (= canonical_pos of slot9 == 161 minus... wait, end is 161).
        let last = pieces.last().unwrap();
        assert_eq!(last.monomer_span, (160, 161));
        assert_eq!(last.canonical_span, (160, 161));
        assert_eq!(last.band(), 4);
    }

    #[test]
    fn knockout_produces_missing_anchor_gap() {
        // Knockout slot4 (GTGGA@86) with the fixture's filler rule.
        let mut m = CONS_171.to_vec();
        for j in 0..5 {
            m[86 + j] = b"ATGC"[(86 + j) % 4];
        }
        let c = build_chain(&m);
        assert_eq!(c.n_primary, 9);
        assert!(c.has_enough_anchors(6));
        let pieces = extract_pieces(&c, SCAFFOLD_PANEL);
        // 10 hits -> 9 pieces.
        assert_eq!(pieces.len(), 9);

        let gaps: Vec<&Piece> = pieces
            .iter()
            .filter(|p| matches!(p.kind, PieceKind::MissingAnchorGap { .. }))
            .collect();
        assert_eq!(gaps.len(), 1, "{pieces:?}");

        // The gap spans slot3 (GAAAC@59, ends at 64) -> slot5 (TGGAA@117).
        // canonical_span = (64, 117), len 53, band = max(4, ceil(0.25*53)) = 14.
        let gap = gaps[0];
        assert_eq!(gap.canonical_span, (64, 117));
        assert_eq!(gap.band(), 14);
        match &gap.kind {
            PieceKind::MissingAnchorGap { skipped_slots } => {
                // panel_idx of slot4 in SCAFFOLD_PANEL is 5 (boundary=0, slot0=1, ..., slot4=5).
                assert_eq!(skipped_slots, &[5]);
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn outlier_segment_is_inter_anchor_with_band_seven() {
        // Outlier monomer: positions 91..117 randomized; anchors intact.
        // The slot4 -> slot5 piece is the targeted segment.
        let mut m = CONS_171.to_vec();
        for (i, p) in (91..117).enumerate() {
            m[p] = b"ACGT"[i % 4]; // arbitrary deterministic randomization
        }
        let c = build_chain(&m);
        assert_eq!(c.n_primary, 10);
        let pieces = extract_pieces(&c, SCAFFOLD_PANEL);
        assert_eq!(pieces.len(), 10);
        // slot4 (GTGGA@86, ends at 91) -> slot5 (TGGAA@117). canonical_span = (91, 117).
        let outlier_piece = pieces.iter().find(|p| p.canonical_span == (91, 117)).expect("found");
        assert!(matches!(outlier_piece.kind, PieceKind::InterAnchor));
        assert_eq!(outlier_piece.band(), 7); // brief: 26 bp -> band 7
    }

    #[test]
    fn exception_chain_gates_on_min_anchors() {
        // Knock out 6 of 10 primaries (the fixture's exception-set logic).
        let mut m = CONS_171.to_vec();
        for (pos, k) in [(14_usize, 5_usize), (24, 5), (40, 5), (59, 5), (86, 5), (117, 5)] {
            for j in 0..k {
                m[pos + j] = b"ATGC"[(pos + j) % 4];
            }
        }
        let c = build_chain(&m);
        assert_eq!(c.n_primary, 4);
        assert!(!c.has_enough_anchors(6), "must route to exceptions");
        // Pieces would still extract from the remaining 5 hits, but caller
        // should short-circuit before reaching extract_pieces.
        let pieces = extract_pieces(&c, SCAFFOLD_PANEL);
        assert_eq!(pieces.len(), 4); // 5 hits -> 4 windows
    }

    #[test]
    fn piece_band_thresholds() {
        let mk = |a: usize, b: usize| Piece {
            monomer_span: (0, 0),
            canonical_span: (a, b),
            left_anchor_idx: 0,
            right_anchor_idx: 1,
            kind: PieceKind::InterAnchor,
        };
        assert_eq!(mk(0, 0).band(), 4, "empty span -> floor 4");
        assert_eq!(mk(0, 15).band(), 4, "ceil(3.75)=4 == floor");
        assert_eq!(mk(0, 16).band(), 4, "ceil(4)=4 == floor");
        assert_eq!(mk(0, 17).band(), 5, "ceil(4.25)=5 > floor");
        assert_eq!(mk(0, 26).band(), 7, "brief: 26 bp -> 7");
        assert_eq!(mk(0, 53).band(), 14, "brief: 53 bp -> 14");
    }

    #[test]
    fn piece_seq_helpers_clamp_degenerate_span() {
        // Anchors-overlap pathology: monomer_end < monomer_start would be
        // clamped to == monomer_start by extract_pieces (yielding an empty
        // span). The helpers must return empty slices for such cases.
        let degenerate = Piece {
            monomer_span: (10, 10),
            canonical_span: (10, 10),
            left_anchor_idx: 0,
            right_anchor_idx: 1,
            kind: PieceKind::InterAnchor,
        };
        assert!(degenerate.monomer_seq(CONS_171).is_empty());
        assert!(degenerate.consensus_seq(CONS_171).is_empty());
    }

    #[test]
    fn piece_seq_helpers_borrow_correct_bytes() {
        // First piece in the canonical chain: monomer[11..14] == "CAT" — the
        // 3-bp filler between the boundary anchor (CTTTGTGATGT, ends at 11)
        // and slot0 (CATTC, starts at 14). Filler base i is `"ATGC"[i%4]`,
        // overridden by the boundary k-mer through index 10, then resuming:
        // 11='C' (ATGC[3]), 12='A' (ATGC[0]), 13='T' (ATGC[1]) -> "CAT".
        let p = Piece {
            monomer_span: (11, 14),
            canonical_span: (11, 14),
            left_anchor_idx: 0,
            right_anchor_idx: 1,
            kind: PieceKind::InterAnchor,
        };
        assert_eq!(p.monomer_seq(CONS_171), b"CAT");
        assert_eq!(p.consensus_seq(CONS_171), b"CAT");
    }

    #[test]
    fn noncanonical_order_detected() {
        // Hand-craft a chain where the second hit has lower panel_idx than
        // the first — simulating an anchor reordering.
        use crate::anchor_detect::AnchorHit;
        let hits = vec![
            AnchorHit {
                panel_idx: 3,
                label: "slot2",
                canonical_pos: 40,
                found_pos: 20,
                hd: 0,
                role: Role::Primary,
            },
            AnchorHit {
                panel_idx: 2,
                label: "slot1",
                canonical_pos: 24,
                found_pos: 50,
                hd: 0,
                role: Role::Primary,
            },
        ];
        let c = Chain::from_hits(hits);
        assert!(c.noncanonical_order);
    }
}
