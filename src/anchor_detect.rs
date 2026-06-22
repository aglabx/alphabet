//! Per-monomer anchor detection for the anchor-scaffold MSA pipeline.
//!
//! Brief: `docs/anchor_scaffold_msa_brief.md` § Core modules / 1 (chain) +
//! § Parameters. This module owns the canonical scaffold panel and the
//! per-slot windowed Hamming scan that every consumer of anchor hits should
//! call.
//!
//! Why this is separate from `cmd::cut_v2::PANEL`:
//!   * `cut_v2::PANEL` is the data-derived V5.1 grid used for ARRAY-level
//!     cut placement (a score-spreading algorithm across the whole array,
//!     looking for where a monomer's start could plausibly land).
//!   * `SCAFFOLD_PANEL` here lists canonical MSA-column landmarks for an
//!     ALREADY-CUT, canonical-rotation monomer (per-slot windowed match,
//!     fixed columns). Different positions, different algorithm.
//!
//! The brief's "REUSE the cut-v2 path; do not re-implement detection" is
//! satisfied by this single source of truth for the scaffold scan, not by
//! sharing the cut-v2 PANEL data (whose positions don't line up with the
//! canonical-anchor positions the smoke fixture plants).

use serde::{Deserialize, Serialize};

/// Role of a panel slot for the chain / min_anchors logic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Role {
    /// Frames the monomer; excluded from the `min_anchors` count.
    Boundary,
    /// Counted toward `min_anchors`.
    Primary,
    /// Presence-only; also excluded from the `min_anchors` count.
    Alternate,
}

/// One slot in a scaffold panel.
#[derive(Debug, Clone)]
pub struct PanelEntry {
    pub label: &'static str,
    pub kmer: &'static [u8],
    pub canonical_pos: usize,
    pub tolerance: usize,
    pub hd_threshold: usize,
    pub role: Role,
}

/// Canonical 11-slot scaffold panel: 1 boundary + 10 primaries.
///
/// Positions and k-mers are derived from
/// `tests/data/make_smoke_fixture.py`'s `ANCHORS` table — the fixture is the
/// ground truth for what the scaffolder is expected to find. `tolerance=3`
/// mirrors cut-v2's per-anchor window; `hd_threshold=0` is the smoke default
/// and will be raised via CLI for real data (see brief's `*(recalibrate)*`).
///
/// CATTC@14 (slot0) and CATTC@155 (slot8) are bi-positional: the same k-mer
/// fills two slots at distinct columns; detection is per-slot and scoped to
/// that slot's window, so they never cross-confuse.
pub const SCAFFOLD_PANEL: &[PanelEntry] = &[
    PanelEntry { label: "boundary", kmer: b"CTTTGTGATGT", canonical_pos:   0, tolerance: 3, hd_threshold: 0, role: Role::Boundary },
    PanelEntry { label: "slot0",    kmer: b"CATTC",       canonical_pos:  14, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot1",    kmer: b"CAGAG",       canonical_pos:  24, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot2",    kmer: b"CTTTT",       canonical_pos:  40, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot3",    kmer: b"GAAAC",       canonical_pos:  59, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot4",    kmer: b"GTGGA",       canonical_pos:  86, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot5",    kmer: b"TGGAA",       canonical_pos: 117, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot6",    kmer: b"AAACT",       canonical_pos: 141, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot7",    kmer: b"ACAGA",       canonical_pos: 148, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot8",    kmer: b"CATTC",       canonical_pos: 155, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
    PanelEntry { label: "slot9",    kmer: b"CAGAA",       canonical_pos: 161, tolerance: 3, hd_threshold: 0, role: Role::Primary  },
];

/// One anchor hit in a monomer.
///
/// JSON shape (brief I/O contract):
/// `{"slot":"slot0","canonical_pos":14,"pos":14,"hd":0,"role":"primary"}`
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnchorHit {
    #[serde(skip)]
    pub panel_idx: usize,
    #[serde(rename = "slot")]
    pub label: &'static str,
    pub canonical_pos: usize,
    #[serde(rename = "pos")]
    pub found_pos: usize,
    pub hd: usize,
    pub role: Role,
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Scan `monomer` for each panel slot independently and return one hit per
/// slot where a k-mer match within the slot's `hd_threshold` and `tolerance`
/// exists.
///
/// `max_hd` further caps each slot's HD threshold (so the caller can globally
/// tighten matching at runtime without rebuilding the panel).
///
/// Tie-break within a slot's window (per ASM-msa-parameter-set):
///   1. minimum Hamming distance,
///   2. minimum `|found_pos - canonical_pos|`,
///   3. minimum `found_pos` (deterministic leftmost when fully tied).
///
/// Cross-slot tie-break (primary > alternate, then slot index) is implicit:
/// each slot is scanned independently in its own window, so a single
/// `monomer` position can only contribute to slots whose windows include it.
/// In practice for `tolerance=3` and the canonical panel the slot windows are
/// disjoint, so cross-slot collisions don't arise.
///
/// Hits are returned sorted by `found_pos` so the chain stage can walk them
/// left-to-right.
pub fn detect_anchors(panel: &[PanelEntry], monomer: &[u8], max_hd: usize) -> Vec<AnchorHit> {
    let mut hits = Vec::with_capacity(panel.len());
    for (idx, entry) in panel.iter().enumerate() {
        let k = entry.kmer.len();
        if k == 0 || monomer.len() < k {
            continue;
        }
        let last_start = monomer.len() - k;
        let lo = entry.canonical_pos.saturating_sub(entry.tolerance);
        let hi = (entry.canonical_pos + entry.tolerance).min(last_start);
        if lo > hi {
            continue;
        }
        let hd_cap = entry.hd_threshold.min(max_hd);
        let mut best: Option<(usize, usize)> = None;
        for pos in lo..=hi {
            let hd = hamming(&monomer[pos..pos + k], entry.kmer);
            if hd > hd_cap {
                continue;
            }
            let better = match best {
                None => true,
                Some((bp, bh)) => match hd.cmp(&bh) {
                    std::cmp::Ordering::Less => true,
                    std::cmp::Ordering::Greater => false,
                    std::cmp::Ordering::Equal => {
                        let cur_d = pos.abs_diff(entry.canonical_pos);
                        let best_d = bp.abs_diff(entry.canonical_pos);
                        match cur_d.cmp(&best_d) {
                            std::cmp::Ordering::Less => true,
                            std::cmp::Ordering::Greater => false,
                            std::cmp::Ordering::Equal => pos < bp,
                        }
                    }
                },
            };
            if better {
                best = Some((pos, hd));
            }
        }
        if let Some((pos, hd)) = best {
            hits.push(AnchorHit {
                panel_idx: idx,
                label: entry.label,
                canonical_pos: entry.canonical_pos,
                found_pos: pos,
                hd,
                role: entry.role,
            });
        }
    }
    hits.sort_by_key(|h| h.found_pos);
    hits
}

/// Number of `Primary`-role hits (the quantity compared against `min_anchors`).
pub fn count_primary(hits: &[AnchorHit]) -> usize {
    hits.iter().filter(|h| h.role == Role::Primary).count()
}

/// Allowed deviation between an observed gap and the canonical gap between two
/// chained anchors. Base tolerance + 20% of the canonical distance, so longer
/// (or skipped-slot) gaps tolerate more accumulated indel drift. Frame-offset
/// invariant: only *differences* of positions enter, never absolute positions.
fn gap_tol(canon_gap: i64) -> i64 {
    3 + (canon_gap.abs() * 2) / 10
}

/// Relative-chaining anchor detection (BTN-scaffold-panel-fixed-windows-dont-transfer).
///
/// Unlike [`detect_anchors`] (fixed absolute ±tolerance windows from monomer
/// start — which collapses on real indel-bearing monomers), this collects every
/// candidate k-mer match per slot across the WHOLE monomer, then selects the
/// longest co-linear chain whose inter-anchor gaps match the canonical spacings
/// within [`gap_tol`]. Because only position *differences* are constrained, the
/// result floats with indel drift and is invariant to the monomer's rotation
/// frame offset. Ties broken by lower total Hamming distance.
///
/// O(c²) per monomer in the candidate count `c`; `c` is bounded by `max_hd`.
pub fn detect_anchors_chained(panel: &[PanelEntry], monomer: &[u8], max_hd: usize) -> Vec<AnchorHit> {
    // 1. all candidate matches per slot, genome-wide.
    struct Cand {
        panel_idx: usize,
        pos: usize,
        hd: usize,
    }
    let mut cands: Vec<Cand> = Vec::new();
    for (idx, entry) in panel.iter().enumerate() {
        let k = entry.kmer.len();
        if k == 0 || monomer.len() < k {
            continue;
        }
        let hd_cap = entry.hd_threshold.min(max_hd);
        let last = monomer.len() - k;
        for pos in 0..=last {
            let hd = hamming(&monomer[pos..pos + k], entry.kmer);
            if hd <= hd_cap {
                cands.push(Cand { panel_idx: idx, pos, hd });
            }
        }
    }
    if cands.is_empty() {
        return Vec::new();
    }
    // 2. order by position, then panel index.
    cands.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.panel_idx.cmp(&b.panel_idx)));
    let n = cands.len();

    // 3. DP for the longest gap-consistent co-linear chain (tie: min total hd).
    let mut dp = vec![1usize; n];
    let mut hdsum: Vec<i64> = cands.iter().map(|c| c.hd as i64).collect();
    let mut prev = vec![usize::MAX; n];
    let canon = |idx: usize| panel[idx].canonical_pos as i64;
    let mut best_i = 0usize;
    for i in 0..n {
        for j in 0..i {
            if cands[j].panel_idx >= cands[i].panel_idx || cands[j].pos >= cands[i].pos {
                continue;
            }
            let canon_gap = canon(cands[i].panel_idx) - canon(cands[j].panel_idx);
            let obs_gap = cands[i].pos as i64 - cands[j].pos as i64;
            if (obs_gap - canon_gap).abs() > gap_tol(canon_gap) {
                continue;
            }
            let cand_len = dp[j] + 1;
            let cand_hd = hdsum[j] + cands[i].hd as i64;
            if cand_len > dp[i] || (cand_len == dp[i] && cand_hd < hdsum[i]) {
                dp[i] = cand_len;
                hdsum[i] = cand_hd;
                prev[i] = j;
            }
        }
        if dp[i] > dp[best_i] || (dp[i] == dp[best_i] && hdsum[i] < hdsum[best_i]) {
            best_i = i;
        }
    }

    // 4. backtrack the winning chain.
    let mut chain_idx = Vec::new();
    let mut cur = best_i;
    while cur != usize::MAX {
        chain_idx.push(cur);
        cur = prev[cur];
    }
    chain_idx.reverse();

    // 5. materialize hits (already in increasing pos by construction).
    chain_idx
        .iter()
        .map(|&ci| {
            let c = &cands[ci];
            let e = &panel[c.panel_idx];
            AnchorHit {
                panel_idx: c.panel_idx,
                label: e.label,
                canonical_pos: e.canonical_pos,
                found_pos: c.pos,
                hd: c.hd,
                role: e.role,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The synthetic 171bp consensus produced by
    /// `tests/data/make_smoke_fixture.py`'s `build_consensus()`. Inlined so
    /// the unit test is self-contained (no file I/O).
    const CONS_171: &[u8] = b"CTTTGTGATGTCATCATTCCATGCCAGAGTGCATGCATGCCTTTTTGCATGCATGCATGGAAACATGCATGCATGCATGCATGCATGTGGACATGCATGCATGCATGCATGCATGCATGGAAGCATGCATGCATGCATGCAAAACTGCACAGATGCATTCACAGAAGCATG";

    #[test]
    fn consensus_is_171bp() {
        assert_eq!(CONS_171.len(), 171);
    }

    #[test]
    fn panel_has_canonical_layout() {
        assert_eq!(SCAFFOLD_PANEL.len(), 11);
        assert_eq!(SCAFFOLD_PANEL[0].role, Role::Boundary);
        assert_eq!(SCAFFOLD_PANEL[0].canonical_pos, 0);
        for entry in &SCAFFOLD_PANEL[1..] {
            assert_eq!(entry.role, Role::Primary, "{} not primary", entry.label);
        }
        // Positions are strictly increasing — order matches MSA-column order.
        for w in SCAFFOLD_PANEL.windows(2) {
            assert!(w[0].canonical_pos < w[1].canonical_pos);
        }
    }

    #[test]
    fn consensus_yields_all_eleven_anchors_at_canonical_pos() {
        let hits = detect_anchors(SCAFFOLD_PANEL, CONS_171, 0);
        assert_eq!(hits.len(), 11, "missing hits: {hits:?}");
        for h in &hits {
            assert_eq!(h.found_pos, h.canonical_pos, "{} not at canonical_pos", h.label);
            assert_eq!(h.hd, 0);
        }
        assert_eq!(count_primary(&hits), 10);
    }

    #[test]
    fn hits_sorted_by_found_pos() {
        let hits = detect_anchors(SCAFFOLD_PANEL, CONS_171, 0);
        for w in hits.windows(2) {
            assert!(w[0].found_pos < w[1].found_pos, "not sorted: {hits:?}");
        }
    }

    #[test]
    fn bipositional_cattc_disambiguated_by_window() {
        // CATTC fills both slot0 (canonical_pos 14) and slot8 (155). Per-slot
        // detection must report exactly these two, never cross-confuse them.
        let hits = detect_anchors(SCAFFOLD_PANEL, CONS_171, 0);
        let cattc: Vec<_> = hits
            .iter()
            .filter(|h| SCAFFOLD_PANEL[h.panel_idx].kmer == b"CATTC")
            .collect();
        assert_eq!(cattc.len(), 2);
        assert_eq!(cattc[0].label, "slot0");
        assert_eq!(cattc[0].found_pos, 14);
        assert_eq!(cattc[1].label, "slot8");
        assert_eq!(cattc[1].found_pos, 155);
    }

    #[test]
    fn knockout_loses_only_targeted_anchor() {
        // Knockout slot4 (GTGGA@86) using the fixture's filler rule
        // (`ATGC[(pos+j)%4]`). Expect 10 hits, slot4 missing.
        let mut m = CONS_171.to_vec();
        let ko = 86usize;
        for j in 0..5 {
            m[ko + j] = b"ATGC"[(ko + j) % 4];
        }
        let hits = detect_anchors(SCAFFOLD_PANEL, &m, 0);
        assert_eq!(hits.len(), 10, "{hits:?}");
        assert!(!hits.iter().any(|h| h.label == "slot4"));
    }

    #[test]
    fn exception_drops_six_primaries_yields_five_hits() {
        // Mirror the fixture's first exception monomer: knock out prim[0..6]
        // (slot0..slot5 in MSA-column order). Boundary + 4 remaining primaries
        // = 5 total hits.
        let mut m = CONS_171.to_vec();
        let knock = [(14_usize, 5_usize), (24, 5), (40, 5), (59, 5), (86, 5), (117, 5)];
        for (pos, k) in knock {
            for j in 0..k {
                m[pos + j] = b"ATGC"[(pos + j) % 4];
            }
        }
        let hits = detect_anchors(SCAFFOLD_PANEL, &m, 0);
        assert_eq!(hits.len(), 5, "{hits:?}");
        assert_eq!(count_primary(&hits), 4);
        assert!(hits.iter().any(|h| h.role == Role::Boundary));
    }

    #[test]
    fn empty_monomer_yields_no_hits() {
        assert!(detect_anchors(SCAFFOLD_PANEL, b"", 0).is_empty());
    }

    #[test]
    fn max_hd_zero_blocks_single_substitution() {
        // 1-base sub inside slot0 (CATTC@14): max_hd=0 must drop the slot.
        let mut m = CONS_171.to_vec();
        m[18] = b'G'; // CATTC[4] = C -> G
        let hits = detect_anchors(SCAFFOLD_PANEL, &m, 0);
        assert_eq!(hits.len(), 10);
        assert!(!hits.iter().any(|h| h.label == "slot0"));
    }

    #[test]
    fn relaxed_panel_admits_single_substitution() {
        // With per-slot hd_threshold=1 AND max_hd=1, slot0 should be found at
        // pos 14 with hd=1.
        let relaxed: &[PanelEntry] = &[PanelEntry {
            label: "slot0",
            kmer: b"CATTC",
            canonical_pos: 14,
            tolerance: 3,
            hd_threshold: 1,
            role: Role::Primary,
        }];
        let mut m = CONS_171.to_vec();
        m[18] = b'G';
        let hits = detect_anchors(relaxed, &m, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].found_pos, 14);
        assert_eq!(hits[0].hd, 1);
    }

    #[test]
    fn json_shape_matches_brief() {
        let hits = vec![AnchorHit {
            panel_idx: 1,
            label: "slot0",
            canonical_pos: 14,
            found_pos: 14,
            hd: 0,
            role: Role::Primary,
        }];
        let s = serde_json::to_string(&hits).unwrap();
        // panel_idx is intentionally not serialized (it's a runtime handle,
        // not a stable identifier outside the panel object).
        assert!(!s.contains("panel_idx"));
        assert!(s.contains(r#""slot":"slot0""#), "{s}");
        assert!(s.contains(r#""pos":14"#), "{s}");
        assert!(s.contains(r#""canonical_pos":14"#), "{s}");
        assert!(s.contains(r#""hd":0"#), "{s}");
        assert!(s.contains(r#""role":"primary""#), "{s}");
    }

    #[test]
    fn chained_clean_consensus_recovers_all_primaries() {
        let hits = detect_anchors_chained(SCAFFOLD_PANEL, CONS_171, 0);
        assert_eq!(count_primary(&hits), 10, "{hits:?}");
        // bi-positional CATTC still split across slot0 and slot8 by spacing.
        let cattc: Vec<_> = hits
            .iter()
            .filter(|h| SCAFFOLD_PANEL[h.panel_idx].kmer == b"CATTC")
            .collect();
        assert_eq!(cattc.len(), 2);
    }

    #[test]
    fn chained_beats_windowed_under_indel_drift() {
        // Insert 5 bases after slot0 so slot1..slot9 shift +5 — past the ±3
        // absolute windows the windowed detector uses, but the inter-anchor
        // GAPS (what chaining keys on) are preserved.
        let mut m = CONS_171[..16].to_vec();
        m.extend_from_slice(b"AAAAA");
        m.extend_from_slice(&CONS_171[16..]);
        let windowed = count_primary(&detect_anchors(SCAFFOLD_PANEL, &m, 1));
        let chained = count_primary(&detect_anchors_chained(SCAFFOLD_PANEL, &m, 1));
        assert!(
            chained > windowed,
            "chained {chained} should beat windowed {windowed} under drift"
        );
        assert!(chained >= 9, "chained should recover most anchors, got {chained}");
    }

    #[test]
    fn chained_empty_on_no_candidates() {
        assert!(detect_anchors_chained(SCAFFOLD_PANEL, b"NNNNNNNNNN", 0).is_empty());
    }
}
