//! End-to-end smoke tests for `discover` and `cut` subcommands.
//!
//! Drives the built binary against the in-repo `tests/data/test_chm13.fasta`
//! fixture and checks output format + non-empty results. Tests run with
//! `-t 1` and minimal EM (`--max-sf 1 --max-iter 1`) to keep CI fast.
//!
//! Outputs land under `target/test-tmp/<test-name>/` (gitignored,
//! ephemeral). Following the project rule: nothing in `/tmp`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_alphabet")
}

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("test_chm13.fasta")
}

fn tmp_outdir(name: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("test-tmp")
        .join(name);
    if dir.exists() {
        std::fs::remove_dir_all(&dir).ok();
    }
    std::fs::create_dir_all(&dir).expect("create tmp outdir");
    dir
}

fn assert_nonempty(path: &Path) {
    let meta = std::fs::metadata(path).unwrap_or_else(|e| panic!("missing {}: {}", path.display(), e));
    assert!(meta.len() > 0, "{} is empty", path.display());
}

#[test]
fn discover_smoke() {
    let out = tmp_outdir("discover_smoke");
    let chains = out.join("chains.json");

    let status = Command::new(bin())
        .args([
            "discover",
            fixture().to_str().unwrap(),
            "-o", chains.to_str().unwrap(),
            "-t", "1",
            "--max-sf", "1",
            "--max-iter", "1",
        ])
        .status()
        .expect("spawn alphabet discover");
    assert!(status.success(), "discover exited non-zero");

    assert_nonempty(&chains);

    // Sanity-parse the JSON and require non-zero arrays + non-empty chains
    // or families. (n_chains==0 + n_arrays==0 would mean discover silently
    // produced an empty pipeline, which is what this test guards against.)
    let bytes = std::fs::read(&chains).expect("read chains.json");
    let v: serde_json::Value = serde_json::from_slice(&bytes).expect("parse chains.json");
    let n_arrays = v.get("n_arrays").and_then(|x| x.as_u64()).unwrap_or(0);
    assert!(n_arrays > 0, "n_arrays must be > 0, got {}", n_arrays);
    let n_chains = v.get("n_chains").and_then(|x| x.as_u64()).unwrap_or(0);
    let n_fams = v.get("n_enriched_families").and_then(|x| x.as_u64()).unwrap_or(0);
    assert!(
        n_chains > 0 || n_fams > 0,
        "chains.json must have ≥1 chain or family, got n_chains={n_chains}, n_enriched_families={n_fams}"
    );
}

#[test]
fn cut_smoke() {
    // cut needs a chains.json — produce one inline so this test stands alone.
    let out = tmp_outdir("cut_smoke");
    let chains = out.join("chains.json");
    let monomers = out.join("monomers.tsv");
    let report = out.join("families.json");
    let consensus = out.join("family_consensus.fa");
    let letters = out.join("letters");

    let status = Command::new(bin())
        .args([
            "discover",
            fixture().to_str().unwrap(),
            "-o", chains.to_str().unwrap(),
            "-t", "1",
            "--max-sf", "1",
            "--max-iter", "1",
        ])
        .status()
        .expect("spawn discover");
    assert!(status.success(), "discover exited non-zero");
    assert_nonempty(&chains);

    let status = Command::new(bin())
        .args([
            "cut",
            fixture().to_str().unwrap(),
            "-m", chains.to_str().unwrap(),
            "-o", monomers.to_str().unwrap(),
            "--report", report.to_str().unwrap(),
            "--consensus-fa", consensus.to_str().unwrap(),
            "--letters-dir", letters.to_str().unwrap(),
            "-t", "1",
        ])
        .status()
        .expect("spawn cut");
    assert!(status.success(), "cut exited non-zero");

    assert_nonempty(&monomers);

    // Verify header + at least one data row. The TSV starts with `#`
    // comment lines (provenance: tool version, input, motifs source);
    // the first non-comment line is the column header.
    let text = std::fs::read_to_string(&monomers).expect("read monomers.tsv");
    let mut non_comment = text.lines().filter(|l| !l.starts_with('#'));
    let header = non_comment.next().expect("column header line");
    assert!(header.contains("array_id"), "header missing array_id: {header:?}");
    assert!(header.contains("sequence"), "header missing sequence: {header:?}");
    let n_rows = non_comment.count();
    assert!(n_rows > 0, "monomers.tsv has only a header, no data rows");
}
