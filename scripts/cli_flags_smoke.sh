#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if [ ! -f /tmp/rusttie_self_reads.fq ]; then
    wgsim -e 0.005 -r 0 -R 0 -N 200 -1 50 -2 50 -S 42 \
        validate/fixtures/lambda_virus.fa \
        /tmp/rusttie_self_reads.fq /tmp/rusttie_self_reads.r2.fq >/dev/null 2>&1
fi

echo "=== --no-unal: skip unmapped records ==="
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_unal.sam --no-head
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_no_unal.sam --no-head --no-unal
total=$(wc -l < /tmp/rt_unal.sam)
filt=$(wc -l < /tmp/rt_no_unal.sam)
unal=$(awk -F'\t' '$2 == 4 {n++} END{print n+0}' /tmp/rt_unal.sam)
echo "  full: $total records, $unal unmapped; --no-unal: $filt records (expected: $total - $unal)"

echo "=== --mp 4,4: every mismatch costs exactly 4 (instead of 2-6) ==="
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_mp.sam --no-head --mp 4,4
echo "  AS distribution:"
awk -F'\t' '$2 != 4 {for(i=12;i<=NF;i++) if($i ~ /^AS:i:/) print $i}' /tmp/rt_mp.sam | sort | uniq -c | head -5

echo "=== --score-min L,-100,0: very loose threshold ==="
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_loose.sam --no-head --score-min "L,-100,0"
loose_mapped=$(awk -F'\t' '$2 != 4 {n++} END{print n+0}' /tmp/rt_loose.sam)
default_mapped=$(awk -F'\t' '$2 != 4 {n++} END{print n+0}' /tmp/rt_unal.sam)
echo "  loose mapped: $loose_mapped, default mapped: $default_mapped"

echo "=== --very-sensitive: accepted as preset (no-op for now) ==="
target/debug/rusttie -x validate/fixtures/lambda_virus \
    -U /tmp/rusttie_self_reads.fq -S /tmp/rt_vs.sam --no-head --very-sensitive
diff -q /tmp/rt_vs.sam /tmp/rt_unal.sam >/dev/null && echo "  SAM identical to default (preset is a no-op)"
