# Changelog

## v0.2.0 — 2026-05-27

### Features
- **cut-v2**: replace empirical PANEL with data-derived 15-cluster grid
  (preserve `PANEL_LEGACY` for comparison). Re-derivation on 1.18M cut
  monomers from a 7-haplotype T2T panel; independently rediscovered all 5
  known bi-positional anchors (CATTC, GAAAC, ACAGA, TCTTT, TTGGA) from raw
  enumeration. Max score ~14.0 (was 6.7); `--min-score` default unchanged.
- **cut-v2**: add `--out-exceptions` flag — emits rows for direct cuts
  where `|gap_length − k*period|` exceeds the threshold, with bp context
  strings.
- **cut-v2** subcommand itself: multi-anchor position-aware monomer cutter
  for alpha satellite, ported from the Python prototype.

### Maturity / infrastructure
- **CI**: GitHub Actions workflow runs build + test + clippy with
  `-D warnings` on push and PR.
- **Lints**: zero clippy warnings under `cargo clippy --all-targets
  -- -D warnings`. Three stylistic lints with broad disagreement
  (`needless_range_loop`, `too_many_arguments`, `type_complexity`) are
  silenced crate-wide with rationale in `src/lib.rs`.
- **Tests**: 9 integration tests passing. Added `tests/discover_cut_smoke.rs`
  exercising `discover` and `cut` end-to-end against the in-repo
  `tests/data/test_chm13.fasta` fixture.
- **Errors**: `src/io.rs` now returns `anyhow::Result<T>` with `Context`
  attached at every fallible step (file path, line number, parse cause).
  Callers `unwrap_or_else(|e| panic!("{:?}", e))` to surface the full
  chain instead of bare unwrap traces.
- **Naming**: binary/crate renamed `alphasplitter` → `alphabet`.

## v0.1.0

Initial release as **AlphaSplitter** at
[github.com/aglabx/alphasplitter](https://github.com/aglabx/alphasplitter).
