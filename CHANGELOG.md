# Changelog

## Unreleased

### Features
- **scaffold-msa** subcommand: anchor-scaffolded MSA over canonical-rotation
  monomers. See `docs/anchor_scaffold_msa_brief.md`. Pipeline: detect
  anchors → chain → extract pieces (inter-anchor + missing-anchor gap) →
  banded Levenshtein per piece → star-consensus column tally → rebuild →
  EM (≤ 3 rounds, < 0.1% column-change convergence). Three output streams
  (`standard.jsonl`, `outliers.jsonl`, `exceptions.jsonl`) plus
  `consensus.fa`, `columns.tsv`, `summary.tsv`. Closed-enum cause
  (OK | TOO_FEW_ANCHORS | SEGMENT_OVER_BAND | GAP_UNALIGNABLE |
  NONCANONICAL_ANCHOR_ORDER). Byte-identical across runs (tie-break fixed:
  HD → |Δpos| → leftmost; consensus tie-break lex A<C<G<T).
- **Smoke fixture**: `tests/data/fixtures/{smoke.fasta, smoke_consensus.fa,
  smoke_invariants.json}` — 24 synthetic monomers exercising each branch
  (clean/sub/del/ins/anchor-knockout/outlier/exception). Generator at
  `tests/data/make_smoke_fixture.py`, `make smoke-fixture` regenerates.
- **Library modules**: `align::ond` (banded Levenshtein + edit script),
  `anchor_detect` (per-monomer SCAFFOLD_PANEL scan), `anchor_scaffold`
  (chain + piece extraction + tally + rebuild + EM driver).

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
