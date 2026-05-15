# Alphabet

Ab initio discovery of satellite DNA monomer alphabets from tandem repeat arrays.

Alphabet identifies periodic conserved anchor sites in satellite DNA, assembles them into chain templates, and classifies monomers by structural indels — without reference HMMs, curated databases, or prior knowledge of the satellite family.

## Quick Start

Alphabet is a single `alphabet` CLI with subcommands. Build once, then invoke.

```bash
# Build
cargo build --release          # or: make

# End-to-end pipeline (discover → cut → annotate)
./target/release/alphabet run arrays.10kb.fasta -o results/ -t 96

# Or stage-by-stage
./target/release/alphabet discover arrays.10kb.fasta -o chains.json -t 96
./target/release/alphabet cut      arrays.10kb.fasta -m chains.json -o monomers.tsv -t 96
./target/release/alphabet annotate monomers.tsv annotated.tsv
```

Run `alphabet --help` for the full subcommand list; every subcommand takes `--help` as well. `make run INPUT=arrays.10kb.fasta OUTDIR=results THREADS=96` is the Makefile shortcut for the end-to-end path.

## What it does

**Input:** FASTA of satellite DNA arrays (e.g., from [Satellome](https://github.com/ad3002/satellome) or TRASH).

**Output (under `-o <outdir>`):**
- `chains.json` — per-SF motif bundles with claimed array membership
- `monomers.tsv` — one row per monomer: `family`, `array_id`, coords, `letter`, `subtype`, `site_order`, `site_structure`, sequence
- `monomers_hor.tsv` — HOR strings per array
- `annotated.tsv` — `monomers.tsv` + CENP-B columns
- `families.json` — per-family alphabet report
- `family_consensus.fa` — one consensus per (family, letter)
- `letters/<family>/<letter>.tsv` — all monomers of each letter, one row per monomer (sequence first, then metadata) — ready to feed into muscle/hmmbuild

## Pipeline

```
satellome .10kb.fasta
  │
  ▼
alphabet discover      (per-period, iterative SF discovery)
  ├── Scan 8-mers → filter by periodicity + support + complexity
  ├── Grow by Shannon entropy (cap 11bp)
  ├── Discover links (exact co-occurrence, stable spacing)
  ├── Assemble chain → enrichment (pairwise distances + occupancy ≥30%)
  └── → chains.json (families[]: motifs + claimed array whitelist)
  │
  ▼
alphabet cut -m chains.json
  ├── Partition arrays by family (each array gets cut with only its
  │     family's motifs — no cross-SF contamination)
  ├── Scan both strands, strand-normalize
  ├── Cut at first-site positions → monomers
  ├── Classify per family: letter (indel-defined) + subtype (site variants)
  ├── site_order: abcde, ab-de, abdebcde (dimers)
  ├── site_structure: 0:a18b24c43d18e:57
  └── → monomers.tsv, monomers_hor.tsv, letters/<family>/<letter>.tsv,
       family_consensus.fa, families.json
  │
  ▼
alphabet annotate
  ├── CENP-B box scoring (nTTCGnnnnAnnCGGGn, both strands)
  ├── 9-char and 17-char CIGAR notation
  └── → annotated.tsv
```

## Core Model

A satellite monomer is a chain of conserved anchor sites separated by variable spacers:

```
[--site_A--]---spacer---[--site_B--]---spacer---[--site_C--]...
  conserved    LENGTH      conserved    LENGTH     conserved
  sequence     matters!    sequence     matters!   sequence
```

**Letter** = defined by indels only (which sites are present/deleted + spacer length changes), NOT by substitutions. Letters are assigned locally within each SF family: "A" in `P171SF0` is not the same object as "A" in `P171SF1`; the `family` column disambiguates.

**Site presence detection:** a site is present if the distance between its neighbors matches expected — no fuzzy sequence matching needed.

## Subcommands

| Subcommand | Purpose |
|---|---|
| `run` | End-to-end: discover → cut → annotate |
| `discover` | Chain-first anchor discovery + enrichment + multi-SF/period |
| `cut` | Array → monomers, letter/subtype classification, per-letter dumps |
| `annotate` | Add CENP-B box columns to a monomers TSV |
| `find-box` | Search for any degenerate box pattern (CENP-B, TIGD4, etc.) |
| `spacing` | CENP-B box spacing analysis at single-base resolution |
| `find-periodic-boxes` | De novo periodic box discovery from paired k-mers |
| `reads {build-hmms,scan,classify,extract,alphabet}` | ONT read processing (HPC-HMM workflow) |
| `dev {find-phase,motif-graph}` | Hidden research subcommands |

## Validated on

| Species | Satellite | Period | Sites | Arrays | Monomers | Letters |
|---------|-----------|--------|-------|--------|----------|---------|
| Human (CHM13) | Alpha satellite | 171bp | 12 | 522/1223 | 411K | 585 |
| Human (CHM13) | Human sat I | 2420bp | 87 | 70 | — | — |
| Mouse (C57BL/6J) | Minor satellite | 120bp | 8 | 17/28 | 44K | 105 |
| Zebra finch | 191bp satellite | 191bp | 11 | 262 | 104K | 279 |
| Zebra finch | 717bp satellite | 717bp | 33 | 70 | 82K | 127 |
| Zebra finch | 3772bp satellite | 3772bp | 131 | 108 | — | — |
| 7 primate species | Alpha satellite | 171bp | 5 (v1) | — | 8M | 26 |

Alexandrov validation: 38/39 SF types recovered, 96.3% of monomers assigned.

## Key Findings

- **585 structural letter types** in human alpha satellite (15× more than HumAS-HMMER's 39)
- **CENP-B variable positions** are chromosome-specific suprachromosomal family markers
- **CENP-B spacing within CDR** = chromosome-specific nucleosome positioning signature
- **Mouse centromeric satellite saturated** with CENP-B boxes (89% of monomers)
- **Zebra finch has two alternative box systems** (BoxA-D at 355bp, TIGD4 at 715bp) on different chromosomes
- **Convergent architecture**: human (340bp), ZF (355bp), mouse (120bp) → 1 box ≈ 1 nucleosome repeat

## Output Format (annotated TSV)

Manifest header (one `#` block from `cut`, one from `annotate`):

```
#alphabet v0.2.0 / cut
#input: arrays.10kb.fasta
#motifs_source: chains.json
#family P171SF0 period=171 sites: M1=ACATCACAAAG M2=AGAATGCTTCT ...
#columns: family array_id strand monomer_idx start end length letter subtype
#          site_order site_structure sites distances sequence
#alphabet v0.2.0 / annotate
#CENP-B_box_ref: nTTCGnnnnAnnCGGGn (17bp, searched on BOTH strands)
#cenpb_labels: B+ (score==9), B? (score 7-8), B- (score<=6)
```

Columns:

| Column | Description |
|--------|-------------|
| family | SF family that cut this monomer (e.g. `P171SF0`) |
| array_id | Original FASTA name (no `_rc` suffix) |
| strand | `+` kept forward, `-` array was revcomped during strand normalization |
| start, end | Coordinates within the array |
| letter | Structural type (indel-defined), local to family |
| subtype | Letter + site sequence variant |
| site_order | Chain sites: `abcde`, `ab-de`, `abdebcde` |
| site_structure | `gap_before:sites_with_distances:gap_after` |
| cenpb | `B+` (9/9), `B?` (7-8/9), `B-` (≤6/9) |
| cenpb_cigar17 | Full 17-char CIGAR |
| sequence | Monomer DNA sequence |

## Requirements

- Rust 1.70+
- Input: satellite arrays in FASTA (recommend ≥10kb arrays)

## Citation

Komissarov A. (2026). Alphabet: Ab initio discovery of satellite DNA monomer alphabets from tandem repeat arrays. *In preparation.*

## License

MIT
