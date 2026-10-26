#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/chain"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== 1. Running chain on circular reference (primary alignments) ==="
$BIN -r data/ref_circular.fa -q data/reads.fa -k 15 -w 10 -c -p -o alignments.paf

echo "Generated alignments.paf:"
cat alignments.paf

echo -e "\n=== 2. Plotting Synteny Dot-Plot (SVG) ==="
python3 scripts/plot_chains.py alignments.paf dotplot.svg

rm -f alignments.paf dotplot.svg
echo -e "\nExample run completed successfully!"
