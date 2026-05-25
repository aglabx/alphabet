#!/bin/bash
# End-to-end test for Alphabet pipeline
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$DIR/.." && pwd)"
BIN="$ROOT/target/release"
DATA="$DIR/data"
OUT="$DIR/output"

mkdir -p "$OUT"

echo "=== Building ==="
cd "$ROOT"
cargo build --release 2>&1 | tail -3

echo ""
echo "=== Step 1: discover_chains (CHM13 alpha satellite) ==="
$BIN/discover_chains "$DATA/test_chm13.fasta" --period 171 -o "$OUT/chains.json" -t 4 2>&1 | tail -10
test -f "$OUT/chains.json" && echo "OK: chains.json created" || { echo "FAIL: no chains.json"; exit 1; }

echo ""
echo "=== Step 2: motif_cut ==="
$BIN/motif_cut "$DATA/test_chm13.fasta" -m "$OUT/chains.json" -o "$OUT/monomers.tsv" --report "$OUT/alphabet.json" -t 4 2>&1 | tail -10
NMON=$(wc -l < "$OUT/monomers.tsv")
echo "Monomers: $NMON lines"
test "$NMON" -gt 100 && echo "OK: monomers generated" || { echo "FAIL: too few monomers ($NMON)"; exit 1; }

echo ""
echo "=== Step 3: annotate_cenpb ==="
$BIN/annotate_cenpb "$OUT/monomers.tsv" "$OUT/annotated.tsv"
NANN=$(wc -l < "$OUT/annotated.tsv")
echo "Annotated: $NANN lines"

# Verify CENP-B columns added
NCOL=$(head -1 "$OUT/annotated.tsv" | awk -F'\t' '{print NF}')
echo "Columns: $NCOL"
test "$NCOL" -ge 18 && echo "OK: CENP-B columns present" || { echo "FAIL: missing columns ($NCOL)"; exit 1; }

# Verify B+ hits exist
BPLUS=$(grep -c "B+" "$OUT/annotated.tsv" || true)
echo "B+ monomers: $BPLUS"
test "$BPLUS" -gt 0 && echo "OK: B+ hits found" || { echo "FAIL: no B+ hits"; exit 1; }

echo ""
echo "=== Step 4: cenpb_spacing ==="
$BIN/cenpb_spacing "$DATA/test_chm13.fasta" 7 4 2>&1 | head -10

echo ""
echo "=== Step 5: find_box (CENP-B pattern) ==="
$BIN/find_box "$DATA/test_chm13.fasta" --pattern NTTCGNNNNANNCGGGN --min-score 7 2>&1 | head -10

echo ""
echo "=== ALL TESTS PASSED ==="
echo "Output in: $OUT/"
ls -lh "$OUT/"
