#!/usr/bin/env bash
# Debug a few chr22 reads with RUSTTIE_DEBUG set.
set -uo pipefail
cd "$(dirname "$0")/.."

# First three reads from R1 (BT2 maps these all successfully)
head -12 /tmp/rusttie_chr22/reads_R1.fq > /tmp/dbg_reads.fq
echo "=== reads under test ==="
awk 'NR % 4 == 1' /tmp/dbg_reads.fq

echo ""
echo "=== rusttie with RUSTTIE_DEBUG=1 ==="
RUSTTIE_DEBUG=1 target/release/rusttie -p 1 \
    -x /tmp/rusttie_chr22/rt_chr22 \
    -U /tmp/dbg_reads.fq -S /tmp/dbg.sam --no-head 2>&1 | head -50

echo ""
echo "=== resulting SAM (POS column) ==="
awk -F'\t' '{print $1, "flag="$2, "pos="$4}' /tmp/dbg.sam
