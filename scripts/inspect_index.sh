#!/usr/bin/env bash
set -euo pipefail
echo "=== file sizes ==="
ls -la /tmp/rusttie_self_lambda* validate/fixtures/lambda_virus.{1,2,3,4}.bt2 | awk '{printf "%10d %s\n", $5, $9}'
echo
echo "=== our .1.bt2 header (first 48 bytes) ==="
xxd /tmp/rusttie_self_lambda.1.bt2 | head -3
echo
echo "=== bt2 .1.bt2 header (first 48 bytes) ==="
xxd validate/fixtures/lambda_virus.1.bt2 | head -3
echo
echo "=== our .2.bt2 header ==="
xxd /tmp/rusttie_self_lambda.2.bt2 | head -2
echo
echo "=== bt2 .2.bt2 header ==="
xxd validate/fixtures/lambda_virus.2.bt2 | head -2
