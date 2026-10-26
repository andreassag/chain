#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/chain"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== Benchmark: Banded DP vs. Fenwick RMQ Chaining ==="

echo "1. Banded DP mode (-k 15 -w 10):"
/usr/bin/time -v $BIN -r data/ref_circular.fa -q data/reads.fa -k 15 -w 10 -c -o /dev/null 2>&1 | grep -E "Elapsed|Maximum resident"

echo "2. Fenwick RMQ mode (--fenwick):"
/usr/bin/time -v $BIN -r data/ref_circular.fa -q data/reads.fa -k 15 -w 10 -c --fenwick -o /dev/null 2>&1 | grep -E "Elapsed|Maximum resident"

echo "Benchmark completed."
