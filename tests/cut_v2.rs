//! Smoke tests for `cut-v2` — ports the Python prototype's test cases.
//!
//! Each test builds a synthetic alpha satellite array, drives the cutter
//! directly via the library API (`alphabet::cmd::cut_v2::cut_array`), and
//! checks the same properties the Python tests check. Test 3 is the
//! load-bearing rotation-consistency check.

use alphabet::cmd::cut_v2::{cut_array, ArrayCut, Params};

/// Build a 171-bp synthetic alpha monomer.
///
/// All panel anchors are placed at their canonical positions, filler is `N`
/// (no panel anchor contains N, so no incidental hits). A label byte at
/// position 100 distinguishes letters A/B/V/… in cross-array tests.
fn build_monomer(label: u8, drop_landmark: bool, drop_islands: &[&[u8]]) -> Vec<u8> {
    let length = 171;
    let mut seq = vec![b'N'; length];
    seq[100] = label;

    if drop_landmark {
        for s in seq.iter_mut().take(11) {
            *s = b'G';
        }
    } else {
        for (i, c) in b"CTTTGTGATGT".iter().enumerate() {
            seq[i] = *c;
        }
    }

    let placements: &[(&[u8], usize)] = &[
        (b"CATTC", 14),
        (b"CAGAG", 25),
        (b"CTTTT", 40),
        (b"GTGGA", 86),
        (b"TGGAA", 117),
        (b"CAGAA", 149),
        (b"AGAAA", 162),
    ];
    for (anchor, pos) in placements {
        if drop_islands.iter().any(|d| d == anchor) {
            continue;
        }
        if pos + anchor.len() > length {
            continue;
        }
        for (i, c) in anchor.iter().enumerate() {
            seq[pos + i] = *c;
        }
    }
    seq
}

fn concat(parts: &[&[u8]]) -> Vec<u8> {
    let mut v = Vec::with_capacity(parts.iter().map(|p| p.len()).sum());
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

fn monomers(seq: &[u8], cut: &ArrayCut) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    if cut.cuts.len() < 2 {
        return out;
    }
    for i in 0..cut.cuts.len() - 1 {
        out.push(seq[cut.cuts[i]..cut.cuts[i + 1]].to_vec());
    }
    out
}

fn interp_flag_per_monomer(cut: &ArrayCut) -> Vec<bool> {
    let mut out = Vec::new();
    if cut.cuts.len() < 2 {
        return out;
    }
    for i in 0..cut.cuts.len() - 1 {
        out.push(cut.interpolated[i] || cut.interpolated[i + 1]);
    }
    out
}

#[test]
fn test_homo_aaaa() {
    let a = build_monomer(b'A', false, &[]);
    let arr = concat(&[&a, &a, &a, &a]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);

    assert_eq!(mons.len(), 3, "expected 3 monomers from 4-mer homo array");
    assert_eq!(mons[0].len(), 171, "monomer length");
    let set: std::collections::HashSet<&Vec<u8>> = mons.iter().collect();
    assert_eq!(set.len(), 1, "homo array → all monomers identical");
    assert!(
        mons[0].starts_with(b"CTTTGTGATGT"),
        "monomer should start at landmark"
    );
}

#[test]
fn test_hor_abav() {
    let a = build_monomer(b'A', false, &[]);
    let b = build_monomer(b'B', false, &[]);
    let v = build_monomer(b'V', false, &[]);
    let arr = concat(&[&a, &b, &a, &v]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);

    assert_eq!(mons.len(), 3);
    assert_eq!(mons[0], mons[2], "HOR: m0 == m2 (both A)");
    assert_ne!(mons[0], mons[1], "HOR: m0 != m1 (A vs B)");
}

#[test]
fn test_both_share_letter_a() {
    // THE CRITICAL ROTATION-CONSISTENCY TEST.
    // Cutting AAAA and ABAV in one run must yield byte-identical A monomers
    // across arrays. cut-v1 fails this; cut-v2 must pass.
    let a = build_monomer(b'A', false, &[]);
    let b = build_monomer(b'B', false, &[]);
    let v = build_monomer(b'V', false, &[]);
    let homo = concat(&[&a, &a, &a, &a]);
    let hor = concat(&[&a, &b, &a, &v]);

    let params = Params::default();
    let homo_mons = monomers(&homo, &cut_array(&homo, &params));
    let hor_mons = monomers(&hor, &cut_array(&hor, &params));

    let homo_set: std::collections::HashSet<&Vec<u8>> = homo_mons.iter().collect();
    assert_eq!(homo_set.len(), 1, "homo uniform");
    assert_eq!(hor_mons[0], homo_mons[0], "hor-A pos 0 == homo-A");
    assert_eq!(hor_mons[2], homo_mons[0], "hor-A pos 2 == homo-A");
    assert_ne!(hor_mons[1], homo_mons[0], "hor-B != A");
}

#[test]
fn test_divergent_no_landmark() {
    // Whole array has the 11-bp landmark mutated (paralog-divergent); only
    // island anchors remain. cut-v1 dumps these as outliers; cut-v2 cuts them.
    let a_div = build_monomer(b'A', true, &[]);
    let arr = concat(&[&a_div, &a_div, &a_div, &a_div]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);

    assert_eq!(mons.len(), 3, "divergent array still cuts to 3 monomers");
    assert_eq!(mons[0].len(), 171);
    let set: std::collections::HashSet<&Vec<u8>> = mons.iter().collect();
    assert_eq!(set.len(), 1, "all identical");
    assert!(
        mons[0].starts_with(b"GGGGGGGGGGG"),
        "monomer starts where the (deleted) landmark used to live"
    );
}

#[test]
fn test_partial_anchor_loss() {
    // Mixed: m0 full, m1 lost landmark, m2 lost landmark AND CATTC, m3 full.
    // Islands alone (CAGAG+CTTTT+GTGGA+TGGAA+CAGAA+AGAAA) score >> min_score,
    // so all cuts should be direct (no interpolation needed).
    let a = build_monomer(b'A', false, &[]);
    let a_no_lm = build_monomer(b'A', true, &[]);
    let a_no_lm_no_cattc = build_monomer(b'A', true, &[b"CATTC"]);
    let arr = concat(&[&a, &a_no_lm, &a_no_lm_no_cattc, &a]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);
    let flags = interp_flag_per_monomer(&cut);

    assert_eq!(mons.len(), 3);
    assert!(
        flags.iter().all(|&f| !f),
        "all 3 monomers should be direct (no interpolation): {:?}",
        flags
    );
}

#[test]
fn test_anchor_free_middle_interpolated() {
    // Middle monomer has no anchors at all (just T's). Cutter must skip it
    // with k=2, interpolate exactly one cut, and not collapse neighbours.
    let a = build_monomer(b'A', false, &[]);
    let noisy = vec![b'T'; 171];
    let arr = concat(&[&a, &noisy, &a, &a]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);
    let flags = interp_flag_per_monomer(&cut);

    assert_eq!(mons.len(), 3, "expected 3 monomers");
    assert!(
        flags.iter().any(|&f| f),
        "at least one monomer should be flagged interpolated: {:?}",
        flags
    );
}

#[test]
fn test_outlier_recovery() {
    // The 390 v1-outlier case: array has 0 landmarks anywhere but each monomer
    // carries the full island panel. v1 flags it as outlier; v2 recovers all
    // 4 monomers (5 templates → 4 cut→cut spans).
    let a_div = build_monomer(b'A', true, &[]);
    let b_div = build_monomer(b'B', true, &[]);
    let arr = concat(&[&a_div, &b_div, &a_div, &b_div, &a_div]);
    let cut = cut_array(&arr, &Params::default());
    let mons = monomers(&arr, &cut);

    assert_eq!(mons.len(), 4, "outlier-class array should recover 4 monomers");
}
