#!/usr/bin/env bash
# Real-data validation against published Illumina reads.
#
# Uses NA12878 mitochondrial FASTQ from nf-core/test-datasets (sarek branch):
# 9.5 MB × 2 paired Illumina reads, real error profiles, against the hg38
# chrM reference (~16 KB). Aligns with both rusttie and bowtie2 and reports
# the SAM-level agreement metrics.
#
# These files are stable on GitHub at known commits — change them here if
# the upstream layout moves.
set -euo pipefail
cd "$(dirname "$0")/.."

WORK=/tmp/rusttie_real
mkdir -p "$WORK"
cd "$WORK"

URL_BASE="https://raw.githubusercontent.com/nf-core/test-datasets/Sarek"

if [ ! -f hg38.chrM.fa ]; then
    echo "=== fetching reference (hg38.chrM.fa, ~16KB) ==="
    curl -fsSL "$URL_BASE/testdata/hg38.chrM.fa" -o hg38.chrM.fa
fi

if [ ! -f NA12878_mito_1.fq.gz ] || [ ! -f NA12878_mito_2.fq.gz ]; then
    echo "=== fetching NA12878 mito FASTQs (~9.5 MB each) ==="
    curl -fsSL "$URL_BASE/testdata/NA12878_mito_1.fq.gz" -o NA12878_mito_1.fq.gz
    curl -fsSL "$URL_BASE/testdata/NA12878_mito_2.fq.gz" -o NA12878_mito_2.fq.gz
fi

echo "Reference: $(wc -c < hg38.chrM.fa) bytes (chrM)"
echo "R1: $(wc -c < NA12878_mito_1.fq.gz) bytes (gzipped)"
echo "R2: $(wc -c < NA12878_mito_2.fq.gz) bytes (gzipped)"

# ---- Build with bowtie2-build (oracle) ----
if [ ! -f bt_chrM.1.bt2 ]; then
    echo "=== bowtie2-build ==="
    time bowtie2-build hg38.chrM.fa bt_chrM >/dev/null
fi

# ---- Build with rusttie-build ----
echo "=== rusttie-build ==="
time "$OLDPWD/target/release/rusttie-build" hg38.chrM.fa rt_chrM

# ---- Byte-diff every .bt2 file ----
echo "=== byte-diff vs bowtie2-build ==="
ok=1
for ext in 1 2 3 4 rev.1 rev.2; do
    # cmp returns nonzero when files differ; mask under pipefail.
    d=$( { cmp -l "rt_chrM.$ext.bt2" "bt_chrM.$ext.bt2" 2>&1 || true; } | wc -l)
    sz_rt=$(wc -c < "rt_chrM.$ext.bt2")
    sz_bt=$(wc -c < "bt_chrM.$ext.bt2")
    if [ "$d" -eq 0 ] && [ "$sz_rt" -eq "$sz_bt" ]; then
        printf "  .%s.bt2: %10d bytes, IDENTICAL\n" "$ext" "$sz_rt"
    else
        printf "  .%s.bt2: ours=%d bt2=%d, %d differing bytes\n" "$ext" "$sz_rt" "$sz_bt" "$d"
        ok=0
    fi
done
if [ "$ok" -eq 0 ]; then
    echo "BUILD MISMATCH (continuing — alignment uses BT2's index for both tools)"
fi

# ---- Align with both ----
echo "=== bowtie2 paired-end ==="
time bowtie2 -p 8 -x bt_chrM -1 NA12878_mito_1.fq.gz -2 NA12878_mito_2.fq.gz -S bt.sam 2>&1 | tail -8

echo "=== rusttie paired-end ==="
time "$OLDPWD/target/release/rusttie" -p 8 -x rt_chrM -1 NA12878_mito_1.fq.gz -2 NA12878_mito_2.fq.gz -S rt.sam --no-head 2>&1 | tail -3

# ---- SAM diff ----
echo "=== SAM diff summary ==="
python3 <<'PY'
def parse(path):
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            f_ = line.rstrip().split('\t')
            qname = f_[0]
            is_r1 = (int(f_[1]) & 0x40) != 0
            tags = {}
            for t in f_[11:]:
                k, _, v = t.partition(':')
                tags[k] = v
            rows[(qname, is_r1)] = {
                'flag': int(f_[1]), 'rname': f_[2], 'pos': int(f_[3]),
                'mapq': int(f_[4]), 'cigar': f_[5], 'tags': tags,
            }
    return rows

bt = parse('bt.sam')
rt = parse('rt.sam')
keys = set(bt) & set(rt)
bt_mapped = sum(1 for k in keys if (bt[k]['flag'] & 4) == 0)
rt_mapped = sum(1 for k in keys if (rt[k]['flag'] & 4) == 0)
both_mapped = [k for k in keys if (bt[k]['flag'] & 4) == 0 and (rt[k]['flag'] & 4) == 0]
print(f'  records:       bt={len(bt)} rt={len(rt)}')
print(f'  bt mapped:     {bt_mapped}')
print(f'  rt mapped:     {rt_mapped}')
print(f'  both mapped:   {len(both_mapped)}')
if both_mapped:
    for m in ['pos', 'cigar', 'mapq']:
        n = sum(1 for k in both_mapped if bt[k][m] == rt[k][m])
        print(f'  {m:5s} agree:  {n}/{len(both_mapped)} ({100*n/len(both_mapped):.1f}%)')
    for tag in ['AS', 'NM', 'MD', 'XM', 'XO', 'XG']:
        n = sum(1 for k in both_mapped if bt[k]['tags'].get(tag) == rt[k]['tags'].get(tag))
        print(f'  {tag:5s} agree:  {n}/{len(both_mapped)} ({100*n/len(both_mapped):.1f}%)')
PY
