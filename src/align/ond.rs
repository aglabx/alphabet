use serde::{Deserialize, Serialize};

/// One edit operation in CONSENSUS coordinates.
///
/// `col` indexes positions in the consensus stretch (`a` in `ond_align`).
/// - `Sub{col, to}`         — consensus base at `col` replaced with `to` in the monomer.
/// - `Del{col}`             — consensus base at `col` is absent in the monomer.
/// - `Ins{col, base}`       — extra base in the monomer inserted BEFORE consensus column `col`.
///
/// JSON shape matches the brief's I/O contract:
///   `{"op":"SUB","col":12,"to":"T"}`
///   `{"op":"DEL","col":70}`
///   `{"op":"INS","col":33,"base":"C"}`
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "UPPERCASE")]
pub enum Op {
    Sub { col: usize, to: char },
    Del { col: usize },
    Ins { col: usize, base: char },
}

pub type EditScript = Vec<Op>;

/// Banded Levenshtein edit distance + script.
///
/// `a` = consensus stretch, `b` = monomer piece. Returns `None` when the edit
/// distance exceeds `band` (the caller flags the piece, the monomer is routed
/// to the OUTLIERS stream).
///
/// Naming note: the brief calls this `ond_align` after Myers's O(ND) diff. The
/// golden cases in the brief require `D=1` for a single substitution, which is
/// Levenshtein semantics (sub cost = 1), not Myers's insert+delete-only D. The
/// algorithm here is banded Levenshtein DP; the name is kept for spec
/// traceability with EXP-anchor-scaffold-msa / ASM-msa-parameter-set.
///
/// Traceback tie-break order — diagonal > deletion > insertion. This is what
/// turns adjacent sub-cost-1 cells into a single `Sub` op and lets the brief's
/// golden ins/del cases land on the canonical column.
pub fn ond_align(a: &[u8], b: &[u8], band: usize) -> Option<EditScript> {
    let n = a.len();
    let m = b.len();

    // Minimum edit distance >= |n - m| (the length difference must be paid in
    // pure indels). Early-out keeps the hot OUTLIERS path cheap.
    if n.abs_diff(m) > band {
        return None;
    }

    const INF: usize = usize::MAX / 4;
    let mut dp = vec![vec![INF; m + 1]; n + 1];
    dp[0][0] = 0;
    for i in 1..=n.min(band) {
        dp[i][0] = i;
    }
    for j in 1..=m.min(band) {
        dp[0][j] = j;
    }

    for i in 1..=n {
        let j_min = i.saturating_sub(band).max(1);
        let j_max = (i + band).min(m);
        if j_min > j_max {
            continue;
        }
        for j in j_min..=j_max {
            let cost_sub = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let d_diag = dp[i - 1][j - 1].saturating_add(cost_sub);
            let d_del = dp[i - 1][j].saturating_add(1);
            let d_ins = dp[i][j - 1].saturating_add(1);
            dp[i][j] = d_diag.min(d_del).min(d_ins);
        }
    }

    let total = dp[n][m];
    if total > band {
        return None;
    }

    let mut script: EditScript = Vec::new();
    let mut i = n;
    let mut j = m;
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost_sub = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            if dp[i][j] == dp[i - 1][j - 1].saturating_add(cost_sub) {
                if cost_sub == 1 {
                    script.push(Op::Sub {
                        col: i - 1,
                        to: b[j - 1] as char,
                    });
                }
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i - 1][j].saturating_add(1) {
            script.push(Op::Del { col: i - 1 });
            i -= 1;
            continue;
        }
        if j > 0 && dp[i][j] == dp[i][j - 1].saturating_add(1) {
            script.push(Op::Ins {
                col: i,
                base: b[j - 1] as char,
            });
            j -= 1;
            continue;
        }
        // dp inconsistent — should be unreachable for a well-formed table.
        return None;
    }
    script.reverse();
    Some(script)
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Golden cases verbatim from docs/anchor_scaffold_msa_brief.md ---

    #[test]
    fn golden_exact_match() {
        let r = ond_align(b"ACGTACGT", b"ACGTACGT", 4).expect("D=0 <= band");
        assert!(r.is_empty(), "expected empty script, got {r:?}");
    }

    #[test]
    fn golden_single_sub() {
        let r = ond_align(b"ACGTACGT", b"ACTTACGT", 4).expect("D=1 <= band");
        assert_eq!(r, vec![Op::Sub { col: 2, to: 'T' }]);
    }

    #[test]
    fn golden_single_ins() {
        let r = ond_align(b"ACGTACGT", b"ACGTGACGT", 4).expect("D=1 <= band");
        assert_eq!(r, vec![Op::Ins { col: 4, base: 'G' }]);
    }

    #[test]
    fn golden_single_del() {
        let r = ond_align(b"ACGTACGT", b"ACGACGT", 4).expect("D=1 <= band");
        assert_eq!(r, vec![Op::Del { col: 3 }]);
    }

    #[test]
    fn golden_over_band_returns_none() {
        assert!(ond_align(b"AAAAAAAA", b"CCCCCCCC", 2).is_none());
    }

    // --- Boundary cases (not in the brief, but exercised by piece extraction
    // at monomer ends and across long missing-anchor gaps) ---

    #[test]
    fn empty_pair_zero_distance() {
        let r = ond_align(b"", b"", 0).expect("D=0");
        assert!(r.is_empty());
    }

    #[test]
    fn empty_consensus_inserts_all() {
        let r = ond_align(b"", b"ACGT", 4).expect("D=4 = band");
        assert_eq!(
            r,
            vec![
                Op::Ins { col: 0, base: 'A' },
                Op::Ins { col: 0, base: 'C' },
                Op::Ins { col: 0, base: 'G' },
                Op::Ins { col: 0, base: 'T' },
            ]
        );
    }

    #[test]
    fn empty_piece_deletes_all() {
        let r = ond_align(b"ACGT", b"", 4).expect("D=4 = band");
        assert_eq!(
            r,
            vec![
                Op::Del { col: 0 },
                Op::Del { col: 1 },
                Op::Del { col: 2 },
                Op::Del { col: 3 },
            ]
        );
    }

    #[test]
    fn length_diff_alone_over_band_short_circuits() {
        // |3 - 9| = 6 > band=3 — must return None before filling the DP.
        assert!(ond_align(b"ACG", b"ACGTACGTT", 3).is_none());
    }

    // --- JSON shape: the on-disk JSONL records have to match the brief ---

    #[test]
    fn json_shape_matches_brief() {
        let script = vec![
            Op::Sub { col: 12, to: 'T' },
            Op::Del { col: 70 },
            Op::Ins { col: 33, base: 'C' },
        ];
        let s = serde_json::to_string(&script).expect("serialize");
        assert!(
            s.contains(r#"{"op":"SUB","col":12,"to":"T"}"#),
            "Sub shape drifted: {s}"
        );
        assert!(
            s.contains(r#"{"op":"DEL","col":70}"#),
            "Del shape drifted: {s}"
        );
        assert!(
            s.contains(r#"{"op":"INS","col":33,"base":"C"}"#),
            "Ins shape drifted: {s}"
        );
    }
}
