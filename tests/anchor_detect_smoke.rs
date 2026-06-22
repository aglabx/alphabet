//! Integration test: run `detect_anchors` over the smoke fixture and assert
//! per-monomer hit counts match the brief's expected streams.
//!
//! This is the first end-to-end gate connecting the canonical scaffold panel
//! to the ground-truth fixture. If it stays green, downstream chain/MSA
//! stages can trust that detection is wired up correctly.

use std::collections::HashMap;
use std::path::PathBuf;

use alphabet::anchor_detect::{count_primary, detect_anchors, Role, SCAFFOLD_PANEL};
use alphabet::io::read_fasta;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("fixtures")
        .join(name)
}

#[test]
fn smoke_fixture_anchor_counts_match_expected_stream() {
    let records = read_fasta(fixture("smoke.fasta").to_str().unwrap()).expect("read smoke.fasta");
    let inv_bytes = std::fs::read(fixture("smoke_invariants.json")).expect("read invariants");
    let inv: serde_json::Value = serde_json::from_slice(&inv_bytes).expect("parse invariants");
    let monomers = inv["monomers"].as_array().expect("monomers array");

    assert_eq!(records.len(), 24, "fixture must have 24 monomers");
    assert_eq!(monomers.len(), 24);

    let mut per_stream: HashMap<&str, usize> = HashMap::new();
    let mut primary_histogram: HashMap<(String, usize), usize> = HashMap::new();

    for ((id, seq), inv_entry) in records.iter().zip(monomers) {
        let expected_id = inv_entry["id"].as_str().expect("id");
        let expected_stream = inv_entry["expected_stream"].as_str().expect("stream");
        assert_eq!(id, expected_id, "fixture order mismatch");

        let hits = detect_anchors(SCAFFOLD_PANEL, seq, 0);
        let n_primary = count_primary(&hits);

        // Boundary CTTTGTGATGT@0 must be present in every monomer —
        // the fixture never alters position 0.
        assert!(
            hits.iter().any(|h| h.role == Role::Boundary),
            "{id} ({expected_stream}): missing boundary anchor"
        );

        // Per-stream invariants per the brief's fixture design table.
        match expected_stream {
            "standard" => {
                let expected = if id.contains("knockout") { 9 } else { 10 };
                assert_eq!(
                    n_primary, expected,
                    "{id}: standard stream expected {expected} primaries, got {n_primary} ({hits:?})"
                );
            }
            "outlier" => {
                // 10/10 anchors intact; D>band lives in the segment, not here.
                assert_eq!(
                    n_primary, 10,
                    "{id}: outlier must keep all 10 primary anchors, got {n_primary}"
                );
            }
            "exception" => {
                // The fixture knocks out 6 of 10; 4 remain (< min_anchors=6).
                assert_eq!(
                    n_primary, 4,
                    "{id}: exception expected 4 primaries, got {n_primary}"
                );
                assert!(
                    n_primary < 6,
                    "{id}: must be < min_anchors (6) to be exception"
                );
            }
            other => panic!("{id}: unknown stream {other:?}"),
        }

        *per_stream.entry(expected_stream).or_insert(0) += 1;
        *primary_histogram
            .entry((expected_stream.to_string(), n_primary))
            .or_insert(0) += 1;
    }

    // Mirror the fixture's stream_counts.
    assert_eq!(per_stream.get("standard"), Some(&18));
    assert_eq!(per_stream.get("outlier"), Some(&3));
    assert_eq!(per_stream.get("exception"), Some(&3));

    // Detailed n_primary histogram (catches regressions where, e.g., a
    // tolerance change moves anchors out of window).
    assert_eq!(
        primary_histogram.get(&("standard".into(), 10)),
        Some(&15),
        "standard non-knockout count drifted"
    );
    assert_eq!(
        primary_histogram.get(&("standard".into(), 9)),
        Some(&3),
        "knockout (n_primary=9) count drifted"
    );
    assert_eq!(
        primary_histogram.get(&("outlier".into(), 10)),
        Some(&3),
        "outlier count drifted"
    );
    assert_eq!(
        primary_histogram.get(&("exception".into(), 4)),
        Some(&3),
        "exception count drifted"
    );
}
