//! Integration test: chain + piece extraction over every monomer in the
//! smoke fixture. Each monomer is checked against the brief's fixture design
//! table (clean/sub/del/ins/knockout/outlier/exception).

use std::collections::HashMap;
use std::path::PathBuf;

use alphabet::anchor_detect::{detect_anchors, SCAFFOLD_PANEL};
use alphabet::anchor_scaffold::{extract_pieces, Chain, PieceKind};
use alphabet::io::read_fasta;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("fixtures")
        .join(name)
}

#[test]
fn smoke_fixture_pieces_match_brief_design() {
    let records = read_fasta(fixture("smoke.fasta").to_str().unwrap()).expect("read smoke");
    let inv_bytes = std::fs::read(fixture("smoke_invariants.json")).expect("read invariants");
    let inv: serde_json::Value = serde_json::from_slice(&inv_bytes).expect("parse");
    let monomers = inv["monomers"].as_array().expect("monomers");
    assert_eq!(records.len(), 24);

    let mut totals = HashMap::<&str, (usize, usize, usize)>::new(); // stream -> (n_monomers, n_inter, n_gap)

    for ((id, seq), inv_entry) in records.iter().zip(monomers) {
        let stream = inv_entry["expected_stream"].as_str().expect("stream");
        let chain = Chain::from_hits(detect_anchors(SCAFFOLD_PANEL, seq, 0));
        assert!(
            !chain.noncanonical_order,
            "{id}: fixture is not supposed to plant noncanonical order"
        );

        if stream == "exception" {
            assert!(
                !chain.has_enough_anchors(6),
                "{id}: must fail min_anchors gate"
            );
            // Skip piece extraction — caller routes whole monomer to exceptions.
            continue;
        }

        assert!(
            chain.has_enough_anchors(6),
            "{id}: stream={stream} but chain failed min_anchors gate"
        );
        let pieces = extract_pieces(&chain, SCAFFOLD_PANEL);
        let n_inter = pieces
            .iter()
            .filter(|p| matches!(p.kind, PieceKind::InterAnchor))
            .count();
        let n_gap = pieces
            .iter()
            .filter(|p| matches!(p.kind, PieceKind::MissingAnchorGap { .. }))
            .count();

        // Expected piece counts per stream:
        //   standard non-knockout:  10 inter, 0 gap (n_hits=11 -> 10 windows)
        //   standard knockout:       8 inter, 1 gap (n_hits=10 -> 9 windows; 1 of them spans the knocked slot)
        //   outlier:                10 inter, 0 gap (all 10 anchors present)
        match stream {
            "standard" => {
                if id.contains("knockout") {
                    assert_eq!(n_inter + n_gap, 9, "{id}: total pieces");
                    assert_eq!(n_gap, 1, "{id}: knockout must produce one MissingAnchorGap");
                    // Verify gap canonical span lines up with the knocked
                    // anchor's canonical position. Knockout positions and
                    // their expected gap spans (from the brief):
                    //   slot4 KO @ 86  -> gap (64, 117), band=14
                    //   slot5 KO @ 117 -> gap (91, 141), band=13
                    //   slot6 KO @ 141 -> gap (122, 148), band=7
                    let gap = pieces
                        .iter()
                        .find(|p| matches!(p.kind, PieceKind::MissingAnchorGap { .. }))
                        .unwrap();
                    let expected = if id.contains("knockout86") {
                        ((64usize, 117usize), 14usize)
                    } else if id.contains("knockout117") {
                        ((91, 141), 13)
                    } else if id.contains("knockout141") {
                        ((122, 148), 7)
                    } else {
                        panic!("{id}: unexpected knockout variant")
                    };
                    assert_eq!(gap.canonical_span, expected.0, "{id}: gap span");
                    assert_eq!(gap.band(), expected.1, "{id}: gap band");
                } else {
                    assert_eq!(n_inter + n_gap, 10, "{id}: total pieces");
                    assert_eq!(n_gap, 0, "{id}: non-knockout standard must be all inter");
                }
            }
            "outlier" => {
                assert_eq!(n_inter + n_gap, 10, "{id}: total pieces");
                assert_eq!(n_gap, 0, "{id}: outlier keeps all anchors");
                // The randomized span 91..117 lives in the slot4 -> slot5
                // inter-anchor piece. Verify the piece exists with band 7
                // (brief: 26 bp -> band 7).
                let outlier_piece = pieces
                    .iter()
                    .find(|p| p.canonical_span == (91, 117))
                    .unwrap_or_else(|| panic!("{id}: expected slot4->slot5 piece"));
                assert_eq!(outlier_piece.band(), 7);
                assert!(matches!(outlier_piece.kind, PieceKind::InterAnchor));
            }
            other => panic!("{id}: unknown stream {other:?}"),
        }

        let entry = totals.entry(stream).or_insert((0, 0, 0));
        entry.0 += 1;
        entry.1 += n_inter;
        entry.2 += n_gap;
    }

    // Aggregate consistency against the fixture design:
    //   18 standard total (15 non-knockout * 10 + 3 knockout * 9 = 177 pieces, 3 gaps)
    //   3 outlier total (3 * 10 = 30 pieces, 0 gaps)
    let (n_std, std_inter, std_gap) = totals["standard"];
    assert_eq!(n_std, 18);
    assert_eq!(std_gap, 3, "exactly 3 knockout gaps across the fixture");
    assert_eq!(std_inter, 15 * 10 + 3 * 8, "{std_inter}"); // 174

    let (n_out, out_inter, out_gap) = totals["outlier"];
    assert_eq!(n_out, 3);
    assert_eq!(out_inter, 30);
    assert_eq!(out_gap, 0);
}
