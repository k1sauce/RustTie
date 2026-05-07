#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."

cargo build --release --bin rusttie --bin rusttie-build 2>&1 | tail -2

# Build the multi_n_long fixture
target/release/rusttie-build validate/fixtures/multi_n_long.fa /tmp/mn_idx >/dev/null

# Probe with the chr1[43..] stretch content
cat > /tmp/mn_probe.fq << 'EOF'
@chr1_43_in_stretch_1
CGGAATAGCATGCATGCATGCATGCATGCATGCATG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@chr1_0_in_stretch_0
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@chr2_2_in_stretch_3
CGGTACGTACGTACGTACGTACGTACGTACGTACG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

echo "=== rusttie with RUSTTIE_DEBUG=1 ==="
RUSTTIE_DEBUG=1 target/release/rusttie -p 1 -x /tmp/mn_idx -U /tmp/mn_probe.fq -S /tmp/mn_out.sam --no-head 2>&1
echo ""
echo "=== resulting SAM ==="
cat /tmp/mn_out.sam
