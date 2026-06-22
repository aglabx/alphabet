# Programmer brief — anchor-scaffolded MSA (Rust, Myers O(ND))

> Source: research_compiler artifact `ART-programmer-brief-anchor-scaffold-msa`.
> Copied here verbatim so spec and implementation live together in this repo.
> Companion artifacts (in research_compiler tree, read-only from here):
> `EXP-anchor-scaffold-msa`, `ASM-msa-parameter-set`.

Implement an anchor-scaffolded MSA in the alphasplitter repo. Rust (tens of
millions of monomers). Logic/architecture below.

## Inputs
- Monomer corpus (canonical-rotation FASTA / packed) — ART-alpha-canonical-dataset.
- Anchor detections per monomer — REUSE the cut-v2 path (RES-cut-v2-rust-port);
  do not re-implement detection. Each monomer -> ordered (slot_id, found_pos, hd).
- Parameters from ASM-msa-parameter-set (CLI flags; no hardcoded magic numbers).

## Core modules
1. chain — order present anchors; if n_present < min_anchors -> exceptions stream.
2. segment extractor — between consecutive present anchors emit a piece; across a
   missing anchor emit the longer gap piece (tag which canonical columns it spans).
3. ond_align — banded Myers O(ND) of a piece vs the current per-column consensus
   stretch; returns edit script + D. Reuse an existing O(ND) impl if present.
   D > band -> mark piece, monomer -> outliers.
4. star_consensus — accumulate per-column base/gap counts across all monomers;
   rebuild consensus. EM driver: seed -> align all -> rebuild -> realign, <=3
   rounds or until < 0.1% columns change.
5. emit — three streams: standard / outliers / exceptions + summary. Nothing
   dropped silently.

## Pipeline (ASCII)

 monomer corpus (canonical rotation, FASTA/packed)
        |
        v
 [ detect ]   reuse cut-v2 detector -> per monomer: ordered [(slot_id,pos,hd)]
        |
        v
 [ chain ] --- n_present < min_anchors ? --- yes ---> EXCEPTIONS
        |
        v  (present anchors fix their canonical columns)
 [ emit pieces between fixed columns ]
   - inter-anchor segment   (consecutive present anchors)
   - missing-anchor gap     (skip >=1 absent anchor, longer piece)
        |
        v  for each piece
 [ ond_align ]  banded Myers O(ND) vs consensus stretch
        |        D > band ? -- yes --> mark piece, monomer -> OUTLIERS
        v  edit scripts (small D)
 [ star_consensus ]  accumulate per-column counts; rebuild consensus
        ^___ EM: realign all (<=3 rounds / <0.1% cols change) ___
        |
        v
 [ emit ]  STANDARD: consensus.fa + *.jsonl + columns.tsv
           OUTLIERS / EXCEPTIONS + summary.tsv (counts + enum cause)

## Module layout (recommendation; final repo layout is yours)
- `src/align/ond.rs` — LIBRARY: pure banded Myers O(ND).
  `fn ond_align(a: &[u8], b: &[u8], band: usize) -> Option<EditScript>`
  `a` = consensus stretch, `b` = monomer piece; returns `None` when D > band.
  Unit-tested standalone, no I/O.
- `src/anchor_scaffold/mod.rs` — LIBRARY: chaining, piece extraction (segment +
  missing-anchor gap), star consensus, EM driver. Testable without I/O.
- `src/bin/anchor_scaffold_msa.rs` (or a subcommand) — ORCHESTRATOR: CLI flags
  (defaults from ASM), reads corpus, calls the libs, serializes the streams.

## I/O contract between modules + on disk
`EditScript = ordered Vec<Op>` in CONSENSUS coordinates, where
`Op in { Sub{col,to}, Del{col}, Ins{col,base} }` (`Ins{col}` = base in the monomer
before consensus column `col`). This is the on-disk diff and the unit-test target.

Per-monomer records are JSONL (uniform across streams); aggregates TSV/FASTA:
- `standard.jsonl` — one object per line:
  ```
  {"id","stream":"standard","n_present":<int>,"anchors":[{"slot","pos","hd"}],
   "edit_script":[{"op":"SUB","col":12,"to":"T"},{"op":"DEL","col":70},{"op":"INS","col":33,"base":"C"}],
   "cause":"OK"}
  ```
- `outliers.jsonl` — same shape + per-piece detail:
  ```
  ... "edit_script":[...partial...],"pieces":[{"span":[90,117],"D":19,"status":"OVER_BAND"}],"cause":"SEGMENT_OVER_BAND"
  ```
- `exceptions.jsonl` — `{"id","stream":"exception","n_present":4,"anchors":[...],"cause":"TOO_FEW_ANCHORS"}` (no edit_script)
- `consensus.fa` — refined canonical-column consensus (single record).
- `columns.tsv` — `canonical_col \t A \t C \t G \t T \t gap \t entropy \t indel_rate`
- `summary.tsv` — two blocks: `stream\tcount`, then `cause\tcount`.

`cause` is a closed ENUM (not free-text):
`OK | TOO_FEW_ANCHORS | SEGMENT_OVER_BAND | GAP_UNALIGNABLE | NONCANONICAL_ANCHOR_ORDER`

## ond_align unit test — golden cases (self-contained; no CI dep)
`a` = consensus, `b` = piece:
```
a=ACGTACGT  b=ACGTACGT   band=4  -> D=0, empty script
a=ACGTACGT  b=ACTTACGT   band=4  -> D=1, [Sub{col:2,to:'T'}]
a=ACGTACGT  b=ACGTGACGT  band=4  -> D=1, [Ins{col:4,base:'G'}]   (G between T,A: unambiguous)
a=ACGTACGT  b=ACGACGT    band=4  -> D=1, [Del{col:3}]
a=AAAAAAAA  b=CCCCCCCC   band=2  -> None (D=8 > band)
```

## Parameters (ASM-msa-parameter-set) — concrete defaults; *(recalibrate)* on CHM13 pilot
- scaffold panel: 15-slot V5.1 (`anchor_panel.tsv`); 10 primary slots + boundary X.
- anchor positions/tolerances/HD: exactly per `anchor_panel.tsv` (per-slot
  canonical_pos, tolerance_bp, hd_threshold). Same windows as cut-v2.
- min_anchors present: 6 of 10 primary (boundary X excluded) *(recalibrate)*.
- O(ND) band per piece: `band = max(4, ceil(0.25 * piece_len))`; D > band => piece
  flagged *(recalibrate)*. (Smoke: 26 bp outlier seg -> band 7, D~19>7 => outlier;
  ~53 bp knockout gap -> band 14, D~4<=14 => standard.)
- anchor collision tie-break: (1) min HD to canonical k-mer -> (2) min |pos-canonical_pos|
  -> (3) primary role over alternate -> (4) min slot_idx.
- bi-positional matching: detection is per-slot, scoped to that slot's position
  window ONLY; same k-mer string fills two slots at distinct positions
  (`CATTC` = slot0@14 AND slot8@155; `ACAGA`/`GAAAC`/`AAACT` also repeat). `CATTC` near 14
  maps to slot0, never slot8.
- outlier vs exception: `n_present < min_anchors` => EXCEPTION (whole monomer).
  Else if >=1 piece D > band (or gap unalignable within band) => OUTLIER, but the
  partial alignment is kept/emitted. Non-canonical anchor order => OUTLIER.
- EM rounds/convergence: <=3 rounds; stop when < 0.1% of columns change.
- segment-consensus tie-break: highest base count; ties lexicographic A<C<G<T
  (no IUPAC; single-base, byte-stable consensus).
- seed consensus: existing canonical alpha consensus.
- monomer rotation: canonical (AnchorGrid Stage-2 frame) — upstream dependency.

## Smoke fixture (ground-truth; you run, <30 s, no production data/threads)
Deterministic synthetic fixture; each monomer's expected stream + column is fixed
by construction. Generator -> `fixtures/{smoke.fasta, smoke_consensus.fa, smoke_invariants.json}`.
Synthetic: 171 bp consensus = ATGC filler (no CpG; only homopolymer is the `CTTTT`
anchor/filler junction, <=5, not at any test site), 15-slot anchors planted at
canonical positions. CODE test only, not biology.

24 monomers, expected `{standard:18, outlier:3, exception:3}`:
```
clean         x6  trivial align              D=0, 10/10 primary anchors
substitution  x4  inter-anchor segment       D=1 at col 12/33/70/130, exact base
deletion      x3  segment 1-bp gap           len 170; gap unique col 12/70/130
insertion     x2  segment 1-bp insert        len 172; col 33/130, base != neighbours
anchor-knockout x3 across-missing-anchor O(ND) 9/10 present; gap aligns (D<=band) -> standard
outlier       x3  segment D>band             10/10 present; seg 91-116 random -> outliers
exception     x3  < min_anchors              6 anchors -> filler, 4 present -> exceptions
```

fixtures md5 (verify your regen):
```
smoke.fasta            d8b31801c8a57cf4e14a95b7a71ab9e4
smoke_consensus.fa     1a67ab0484f53c79621323393833eb5f
smoke_invariants.json  6f37c976216d1d1ea776582799ed2586
```

Smoke assertions:
- `n_standard + n_outliers + n_exceptions == 24`; stream counts == `{18,3,3}`.
- consensus length == 171.
- planted substitutions appear at the expected column with the expected base.
- each deletion/insertion produces exactly one gap/insert at the expected column.
- the three knockout monomers are in STANDARD (proves the missing-anchor path), not exceptions.
- two runs of the tool are byte-identical.
- unit test on `ond_align`: known pair -> known D + script (golden cases above).

## Acceptance criteria — done when:
- [ ] Smoke fixture passes ALL assertions (counts, columns, knockout->standard, determinism).
- [ ] `ond_align` unit test passes; matches the golden banded edit distances.
- [ ] Three streams emitted; counts sum to `n_input`; outliers/exceptions carry enum `cause`.
- [ ] Output byte-identical across two runs (tie-break fixed per ASM).
- [ ] All imposed params are CLI flags with defaults from ASM — no hardcoded magic numbers.
- [ ] `--shuffle-input` hook (or accepts a pre-shuffled file) for the human's shuffle null on aglab0.
- [ ] Wall-clock scales ~linearly on smoke->CHM13 ramp (report monomers/s).
- [ ] No writes to `research/.../state/**`; production run left to the human.

## Out of scope
- `research/.../state/**` (research-compiler owns it).
- Production run on the full corpus (the human runs it on aglab0).
