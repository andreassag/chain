> [!NOTE]
> **Archive Notice**: This project is an educational and reference implementation of minimizer seed-and-chaining in C++20. It is finished and archived; no active maintenance or feature development is planned.

# chain: Minimizer Seed-and-Chain Engine

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`chain` is a high-performance C++20 co-linear seed-and-chain engine inspired by `minimap2`. It extracts sliding-window $(k, w)$-minimizers in $O(N)$ time via monotonic queues, identifies exact seed anchors across forward and reverse complement strands, supports native circular chromosome origin-wrapping, and computes optimal chaining dynamic programming using both banded DP ($O(N \cdot B)$) and coordinate-compressed Fenwick Range Maximum Query trees ($O(N \log N)$).

---

## Key Features

1. **Sliding-Window Canonical $(k, w)$-Minimizers**:
   - Monotonic min-queue maintaining the minimum canonical hash in $O(1)$ amortized time per base.
   - Strand-aware anchor pairing with automatic reverse complement coordination.
2. **Native Circular Chromosome Origin Wrapping (`--circular`)**:
   - Maps long reads or syntenic contigs that span across the circular bacterial origin ($L_{\text{ref}} - 1 \leftrightarrow 0$).
   - Seamlessly stitches anchors in $[0, 2L_{\text{ref}})$ and normalizes coordinates modulo chromosome length, tagged with `circ:Z:true`.
3. **Dual Chaining Solvers**:
   - **Banded DP ($O(N \cdot B)$)**: Fast exact chaining over bounded lookback windows ($B = 64$).
   - **Fenwick Tree RMQ ($O(N \log N)$)**: Global co-linear chaining with gap cost optimization.
4. **Standard 12-Column PAF & Standalone SVG Dot-Plot**:
   - Emits standard PAF records (`tp:A:P/S`, `cm:i:`, `s1:i:`, `circ:Z:`).
   - Bundled Python script (`scripts/plot_chains.py`) turns PAF files into interactive SVG synteny dot-plots.

---

## Architectural Overview

```
Reference & Query FASTAs
          |
  [ extract_minimizers(k, w) ]  <--- NtHash + Monotonic Min-Queue
          |
  [ Reference Seed Index ] (Hash Table with Circular Origin Duplication)
          |
  [ Seed Anchor Collection ] (Forward & Reverse-Complement Matching)
          |
          +-----------------------+
          |                       |
  [ Banded DP O(N*B) ]   [ Fenwick RMQ O(N log N) ]
          |                       |
          +-----------+-----------+
                      |
           [ Optimal Chain Backtrace ]
                      |
              +-------+-------+
              |               |
         PAF Formatter    SVG Dot-Plot
```

---

## Build & Test

### Building via Makefile
```bash
make all      # Builds bin/chain
make test     # Runs comprehensive unit tests
```

### Building via CMake
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
ctest --output-on-failure
```

---

## CLI Usage

```bash
# Basic alignment
./bin/chain -r reference.fa -q reads.fa -k 15 -w 10 -o alignments.paf

# Circular chromosome alignment with primary-only filter
./bin/chain -r bacterial_chr.fa -q pacbio_reads.fa -k 15 -w 10 -c -p -o circular.paf

# Fast Fenwick RMQ mode for ultra-long synteny
./bin/chain -r ref.fa -q assembly.fa --fenwick -o synteny.paf

# Generate SVG dot-plot
python3 scripts/plot_chains.py circular.paf dotplot.svg
```

---

## Verification & Benchmarks

Run automated test suite and benchmarks:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Performance on 50 kb Circular Genome:
- **Minimizer Extraction**: ~2 ms for 50 kb reference.
- **Banded DP Chaining**: ~10 ms elapsed time for 16 kb total read sequences.
- **Fenwick RMQ Chaining**: ~10 ms elapsed time.
- **Memory Footprint**: Peak RSS ~ 6 MB.

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
