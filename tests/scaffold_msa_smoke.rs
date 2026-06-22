//! End-to-end smoke test for the `scaffold-msa` subcommand.
//!
//! Drives the built binary against the smoke fixture, verifies the six
//! output artifacts the brief requires, and re-runs to confirm byte-
//! identical output (deterministic).

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_alphabet")
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("fixtures")
        .join(name)
}

fn tmp_outdir(name: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("test-tmp")
        .join(name);
    if dir.exists() {
        std::fs::remove_dir_all(&dir).ok();
    }
    std::fs::create_dir_all(&dir).expect("create outdir");
    dir
}

fn run_scaffold(outdir: &Path) {
    let status = Command::new(bin())
        .args([
            "scaffold-msa",
            "--monomers",
            fixture("smoke.fasta").to_str().unwrap(),
            "--seed",
            fixture("smoke_consensus.fa").to_str().unwrap(),
            "-o",
            outdir.to_str().unwrap(),
        ])
        .status()
        .expect("spawn alphabet scaffold-msa");
    assert!(status.success(), "scaffold-msa exited non-zero");
}

#[test]
fn scaffold_msa_writes_all_six_artifacts() {
    let outdir = tmp_outdir("scaffold_msa_smoke");
    run_scaffold(&outdir);

    for name in [
        "standard.jsonl",
        "outliers.jsonl",
        "exceptions.jsonl",
        "consensus.fa",
        "columns.tsv",
        "summary.tsv",
    ] {
        let path = outdir.join(name);
        let meta = std::fs::metadata(&path).unwrap_or_else(|e| panic!("missing {name}: {e}"));
        assert!(meta.len() > 0, "{name} is empty");
    }
}

#[test]
fn scaffold_msa_jsonl_line_counts_match_streams() {
    let outdir = tmp_outdir("scaffold_msa_jsonl_counts");
    run_scaffold(&outdir);

    for (name, expected) in [
        ("standard.jsonl", 18usize),
        ("outliers.jsonl", 3),
        ("exceptions.jsonl", 3),
    ] {
        let txt = std::fs::read_to_string(outdir.join(name)).expect("read");
        let n = txt.lines().filter(|l| !l.is_empty()).count();
        assert_eq!(n, expected, "{name}: got {n} lines, expected {expected}");
    }
}

#[test]
fn scaffold_msa_summary_matches_brief() {
    let outdir = tmp_outdir("scaffold_msa_summary");
    run_scaffold(&outdir);

    let txt = std::fs::read_to_string(outdir.join("summary.tsv")).expect("read summary");
    // Stream block.
    assert!(txt.contains("standard\t18"), "summary: {txt}");
    assert!(txt.contains("outlier\t3"), "summary: {txt}");
    assert!(txt.contains("exception\t3"), "summary: {txt}");
    // Cause block.
    assert!(txt.contains("OK\t18"), "summary: {txt}");
    assert!(txt.contains("SEGMENT_OVER_BAND\t3"), "summary: {txt}");
    assert!(txt.contains("TOO_FEW_ANCHORS\t3"), "summary: {txt}");
    // None of these should ever fire on the smoke fixture.
    assert!(!txt.contains("GAP_UNALIGNABLE"), "no gap-unalignable expected");
    assert!(!txt.contains("NONCANONICAL_ANCHOR_ORDER"), "no noncanonical expected");
}

#[test]
fn scaffold_msa_consensus_is_byte_stable_against_seed() {
    let outdir = tmp_outdir("scaffold_msa_consensus");
    run_scaffold(&outdir);

    let seed_text = std::fs::read_to_string(fixture("smoke_consensus.fa")).expect("read seed");
    let out_text = std::fs::read_to_string(outdir.join("consensus.fa")).expect("read out");

    // Strip headers, compare sequence bodies.
    let seed_seq: String = seed_text
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect();
    let out_seq: String = out_text
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect();
    assert_eq!(
        seed_seq, out_seq,
        "refined consensus must match seed on the clean smoke fixture"
    );
}

#[test]
fn scaffold_msa_is_deterministic_across_runs() {
    let a = tmp_outdir("scaffold_msa_det_a");
    let b = tmp_outdir("scaffold_msa_det_b");
    run_scaffold(&a);
    run_scaffold(&b);
    for name in [
        "standard.jsonl",
        "outliers.jsonl",
        "exceptions.jsonl",
        "consensus.fa",
        "columns.tsv",
        "summary.tsv",
    ] {
        let aa = std::fs::read(a.join(name)).expect("a");
        let bb = std::fs::read(b.join(name)).expect("b");
        assert_eq!(aa, bb, "{name} differs between two runs");
    }
}
