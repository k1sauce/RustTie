#!/usr/bin/env bash
# Can BT2 itself load our rusttie-build index? Ultimate drop-in test.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

target/debug/rusttie-build validate/fixtures/lambda_virus.fa /tmp/rt_lambda

if [ ! -f /tmp/rusttie_self_reads.fq ]; then
    wgsim -e 0.005 -r 0 -R 0 -N 200 -1 50 -2 50 -S 42 \
        validate/fixtures/lambda_virus.fa \
        /tmp/rusttie_self_reads.fq /tmp/rusttie_self_reads.r2.fq >/dev/null 2>&1
fi

# rusttie-build doesn't yet produce the reverse index; copy from BT2's.
cp validate/fixtures/lambda_virus.rev.1.bt2 /tmp/rt_lambda.rev.1.bt2
cp validate/fixtures/lambda_virus.rev.2.bt2 /tmp/rt_lambda.rev.2.bt2

echo "=== bowtie2 attempting to load our self-built index ==="
bowtie2 -x /tmp/rt_lambda -U /tmp/rusttie_self_reads.fq -S /tmp/bt2_on_ours.sam 2>&1 | tail -10

echo
echo "=== diff (bt2 on our index) vs (bt2 on BT2's index) ==="
bowtie2 -x validate/fixtures/lambda_virus -U /tmp/rusttie_self_reads.fq \
    -S /tmp/bt2_on_bt2.sam 2>/dev/null
diff <(grep -v '^@' /tmp/bt2_on_ours.sam) <(grep -v '^@' /tmp/bt2_on_bt2.sam) | head -10
echo "(empty diff means BT2 produces identical output on both indexes)"
