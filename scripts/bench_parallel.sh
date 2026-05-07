#!/usr/bin/env bash
set -euo pipefail
# Phase 3a parallel speedup measurement on a 10k-read lambda corpus.

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

cargo build --release --bin rusttie >/dev/null 2>&1

if [ ! -f /tmp/big_r1.fq ]; then
    wgsim -e 0.005 -r 0 -R 0 -N 10000 -1 50 -2 50 -S 42 \
      validate/fixtures/lambda_virus.fa /tmp/big_r1.fq /tmp/big_r2.fq >/dev/null 2>&1
fi

echo "=== single-threaded (-p 1) ==="
time target/release/rusttie -p 1 -x validate/fixtures/lambda_virus \
    -U /tmp/big_r1.fq -S /tmp/rt_p1.sam --no-head

echo "=== 4 threads (-p 4) ==="
time target/release/rusttie -p 4 -x validate/fixtures/lambda_virus \
    -U /tmp/big_r1.fq -S /tmp/rt_p4.sam --no-head

echo "=== 8 threads (-p 8) ==="
time target/release/rusttie -p 8 -x validate/fixtures/lambda_virus \
    -U /tmp/big_r1.fq -S /tmp/rt_p8.sam --no-head

echo "=== bowtie2 default ==="
time bowtie2 -x validate/fixtures/lambda_virus \
    -U /tmp/big_r1.fq -S /tmp/bt2.sam --no-head 2>/dev/null

echo "=== diff: rusttie -p1 vs -p8 ==="
if diff -q /tmp/rt_p1.sam /tmp/rt_p8.sam >/dev/null; then
    echo "identical (parallel preserves byte-equivalence)"
else
    diff /tmp/rt_p1.sam /tmp/rt_p8.sam | head -5
    exit 1
fi
