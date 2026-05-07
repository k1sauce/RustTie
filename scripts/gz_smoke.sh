#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

cargo build --bin rusttie 2>&1 | tail -2

# Need a reads file; regenerate if missing.
if [ ! -f /tmp/rusttie_self_reads.fq ]; then
    wgsim -e 0.005 -r 0 -R 0 -N 200 -1 50 -2 50 -S 42 \
        validate/fixtures/lambda_virus.fa \
        /tmp/rusttie_self_reads.fq /tmp/rusttie_self_reads.r2.fq >/dev/null 2>&1
fi
gzip -k -f /tmp/rusttie_self_reads.fq

target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq.gz -S /tmp/rt_gz.sam --no-head
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_plain.sam --no-head

if diff -q /tmp/rt_gz.sam /tmp/rt_plain.sam >/dev/null; then
    echo "gzipped and plain produce identical SAM output"
else
    diff /tmp/rt_gz.sam /tmp/rt_plain.sam | head -5
    exit 1
fi
