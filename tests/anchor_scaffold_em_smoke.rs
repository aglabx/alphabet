//! End-to-end smoke test for the EM driver: run `run_em` on the canonical
//! 24-monomer fixture against the synthetic consensus seed, assert
//! brief-design invariants (stream counts, consensus stability, single-round
//! convergence) and byte-identical re-run.

use std::collections::HashMap;
use std::path::PathBuf;

use alphabet::anchor_detect::SCAFFOLD_PANEL;
use alphabet::anchor_scaffold::{run_em, EmParams, Stream};
use alphabet::io::read_fasta;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("fixtures")
        .join(name)
}

fn load_monomers() -> Vec<(String, Vec<u8>)> {
    read_fasta(fixture("smoke.fasta").to_str().unwrap()).expect("read smoke.fasta")
}

fn load_seed_consensus() -> Vec<u8> {
    let recs = read_fasta(fixture("smoke_consensus.fa").to_str().unwrap()).expect("read seed");
    assert_eq!(recs.len(), 1, "seed consensus must be a single record");
    recs.into_iter().next().unwrap().1
}

#[test]
fn em_run_matches_brief_design_table_and_is_deterministic() {
    let seed = load_seed_consensus();
    assert_eq!(seed.len(), 171, "seed must be 171 bp");

    let monomers = load_monomers();
    assert_eq!(monomers.len(), 24);

    let params = EmParams::default();
    let result = run_em(&seed, &monomers, SCAFFOLD_PANEL, &params);

    // Convergence: 0 columns change in the first round on the clean fixture,
    // so EM stops after a single round (the smoke data isn't a moving target).
    assert_eq!(result.rounds.len(), 1, "{:?}", result.rounds);
    assert_eq!(result.rounds[0].n_columns_changed, 0);

    // Consensus stays byte-stable.
    assert_eq!(
        result.final_consensus, seed,
        "smoke fixture must not drift the consensus"
    );

    // Stream counts mirror the fixture's invariants.
    let counts: HashMap<Stream, usize> = result
        .final_outcomes
        .iter()
        .fold(HashMap::new(), |mut h, o| {
            *h.entry(o.stream).or_insert(0) += 1;
            h
        });
    assert_eq!(counts.get(&Stream::Standard), Some(&18));
    assert_eq!(counts.get(&Stream::Outlier), Some(&3));
    assert_eq!(counts.get(&Stream::Exception), Some(&3));

    // Determinism: same inputs -> identical outcomes serialized as JSONL.
    let result2 = run_em(&seed, &monomers, SCAFFOLD_PANEL, &params);
    let j1 = serde_json::to_string(&result.final_outcomes).unwrap();
    let j2 = serde_json::to_string(&result2.final_outcomes).unwrap();
    assert_eq!(j1, j2, "two runs must produce byte-identical JSONL");

    // Spot-check specific monomers from the fixture by id.
    let by_id: HashMap<&str, &alphabet::anchor_scaffold::MonomerOutcome> = result
        .final_outcomes
        .iter()
        .map(|o| (o.id.as_str(), o))
        .collect();

    // M01_clean: standard, empty edit script.
    let m01 = by_id["M01_clean"];
    assert_eq!(m01.stream, Stream::Standard);
    assert!(m01.edit_script.is_empty(), "{:?}", m01.edit_script);

    // M07_sub12: substitution at col 12 (A -> T).
    let m07 = by_id["M07_sub12"];
    assert_eq!(m07.stream, Stream::Standard);
    assert_eq!(m07.edit_script.len(), 1);

    // Outliers: cause is SegmentOverBand (slot4 -> slot5 piece exceeds band 7).
    for id in ["M19_outlier", "M20_outlier", "M21_outlier"] {
        let o = by_id[id];
        assert_eq!(o.stream, Stream::Outlier, "{id}");
        assert_eq!(
            o.cause,
            alphabet::anchor_scaffold::Cause::SegmentOverBand,
            "{id}"
        );
        assert_eq!(o.pieces.len(), 1, "{id}: {:?}", o.pieces);
        assert_eq!(o.pieces[0].canonical_span, (91, 117));
    }

    // Exceptions: cause is TooFewAnchors, n_present == 4.
    for id in ["M22_exception", "M23_exception", "M24_exception"] {
        let o = by_id[id];
        assert_eq!(o.stream, Stream::Exception, "{id}");
        assert_eq!(
            o.cause,
            alphabet::anchor_scaffold::Cause::TooFewAnchors,
            "{id}"
        );
        assert_eq!(o.n_present, 4, "{id}");
        assert!(o.edit_script.is_empty(), "{id}");
    }

    // Knockouts stay STANDARD: their edit scripts live inside the knocked
    // span (canonical cols 86..91 / 117..122 / 141..146 respectively).
    for (id, expected_zone) in [
        ("M16_knockout86", 86..91),
        ("M17_knockout117", 117..122),
        ("M18_knockout141", 141..146),
    ] {
        let o = by_id[id];
        assert_eq!(o.stream, Stream::Standard, "{id}");
        assert!(!o.edit_script.is_empty(), "{id}: expected ≥1 edit");
        for op in &o.edit_script {
            let col = match op {
                alphabet::align::Op::Sub { col, .. }
                | alphabet::align::Op::Del { col }
                | alphabet::align::Op::Ins { col, .. } => *col,
            };
            assert!(
                expected_zone.contains(&col),
                "{id}: edit at col {col} outside KO zone {expected_zone:?}"
            );
        }
    }
}
