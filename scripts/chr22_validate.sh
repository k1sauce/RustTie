#!/usr/bin/env bash
# chr22 validation: download GRCh38 chr22, build with both bowtie2-build and
# rusttie-build, byte-diff every .bt2 file, then align ~10k synthetic reads
# with both tools and SAM-diff.
set -euo pipefail
cd "$(dirname "$0")/.."

WORK=/tmp/rusttie_chr22
mkdir -p "$WORK"
cd "$WORK"

FA_GZ=Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
FA=Homo_sapiens.GRCh38.dna.chromosome.22.fa
URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/$FA_GZ"

if [ ! -f "$FA" ]; then
    if [ ! -f "$FA_GZ" ]; then
        echo "=== downloading chr22 FASTA from Ensembl ==="
        curl -fsSL -o "$FA_GZ" "$URL"
    fi
    gunzip -k "$FA_GZ"
fi
echo "chr22 FASTA: $(wc -c < "$FA") bytes"

# ---- Build with bowtie2-build (oracle) ----
if [ ! -f bt_chr22.1.bt2 ]; then
    echo "=== bowtie2-build ==="
    time bowtie2-build "$FA" bt_chr22 >/dev/null
fi

# ---- Build with rusttie-build ----
echo "=== rusttie-build ==="
time "$OLDPWD/target/release/rusttie-build" "$FA" rt_chr22

# ---- Byte-diff every .bt2 file ----
echo "=== byte-diff vs bowtie2-build ==="
ok=1
for ext in 1 2 3 4 rev.1 rev.2; do
    d=$(cmp -l "rt_chr22.$ext.bt2" "bt_chr22.$ext.bt2" 2>&1 | wc -l)
    sz_rt=$(wc -c < "rt_chr22.$ext.bt2")
    sz_bt=$(wc -c < "bt_chr22.$ext.bt2")
    if [ "$d" -eq 0 ] && [ "$sz_rt" -eq "$sz_bt" ]; then
        printf "  .%s.bt2: %10d bytes, IDENTICAL\n" "$ext" "$sz_rt"
    else
        printf "  .%s.bt2: ours=%d bt2=%d, %d differing bytes\n" "$ext" "$sz_rt" "$sz_bt" "$d"
        ok=0
    fi
done
[ "$ok" -eq 0 ] && { echo "BUILD MISMATCH — skipping alignment diff"; exit 1; }

# ---- Generate synthetic reads with wgsim ----
if [ ! -f reads_R1.fq ]; then
    echo "=== wgsim 10000 paired reads ==="
    wgsim -e 0.005 -r 0 -R 0 -N 10000 -1 100 -2 100 -d 350 -s 30 -S 42 \
        "$FA" reads_R1.fq reads_R2.fq >/dev/null
fi

# ---- Align with both ----
echo "=== bowtie2 paired-end ==="
time bowtie2 -p 8 -x bt_chr22 -1 reads_R1.fq -2 reads_R2.fq -S bt.sam 2>&1 | tail -8

echo "=== rusttie paired-end ==="
time "$OLDPWD/target/release/rusttie" -p 8 -x rt_chr22 -1 reads_R1.fq -2 reads_R2.fq -S rt.sam --no-head 2>&1 | tail -3

# ---- SAM diff (excluding headers, since rusttie omits them) ----
echo "=== SAM diff summary ==="
python3 <<'PY'
import collections, sys

def parse(path):
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith('@') or not line.strip():
                continue
            f_ = line.rstrip('\n').split('\t')
            qname = f_[0]
            is_r1 = (int(f_[1]) & 0x40) != 0
            tags = {}
            for t in f_[11:]:
                k, _, v = t.partition(':')
                tags[k] = v
            rows[(qname, is_r1)] = {
                'flag': int(f_[1]),
                'rname': f_[2],
                'pos': int(f_[3]),
                'mapq': int(f_[4]),
                'cigar': f_[5],
                'rnext': f_[6],
                'pnext': int(f_[7]),
                'tlen': int(f_[8]),
                'tags': tags,
            }
    return rows

bt = parse('bt.sam')
rt = parse('rt.sam')
print(f'  total records  bt2={len(bt)} rt={len(rt)}')

# Treat all keys; intersect
keys = set(bt) | set(rt)
only_bt = sum(1 for k in keys if k in bt and k not in rt)
only_rt = sum(1 for k in keys if k in rt and k not in bt)
print(f'  only in bt2:   {only_bt}')
print(f'  only in rt:    {only_rt}')

both = [k for k in keys if k in bt and k in rt]
both_mapped = [k for k in both if (bt[k]["flag"] & 4) == 0 and (rt[k]["flag"] & 4) == 0]
bt_only_mapped = [k for k in both if (bt[k]["flag"] & 4) == 0 and (rt[k]["flag"] & 4) != 0]
rt_only_mapped = [k for k in both if (bt[k]["flag"] & 4) != 0 and (rt[k]["flag"] & 4) == 0]
both_unmapped = [k for k in both if (bt[k]["flag"] & 4) != 0 and (rt[k]["flag"] & 4) != 0]
print(f'  both mapped:   {len(both_mapped)}')
print(f'  bt2 mapped, rt unmapped: {len(bt_only_mapped)}')
print(f'  rt mapped, bt2 unmapped: {len(rt_only_mapped)}')
print(f'  both unmapped: {len(both_unmapped)}')

# Field-level agreement on both-mapped
metrics = ['pos', 'rname', 'cigar', 'tlen']
for m in metrics:
    n = sum(1 for k in both_mapped if bt[k][m] == rt[k][m])
    pct = 100*n/len(both_mapped) if both_mapped else 0
    print(f'  {m:10s} agree: {n}/{len(both_mapped)} ({pct:.1f}%)')

# Tag-level
tag_keys = ['AS', 'NM', 'MD', 'XM', 'XO', 'XG']
for tag in tag_keys:
    n = sum(1 for k in both_mapped if bt[k]['tags'].get(tag) == rt[k]['tags'].get(tag))
    pct = 100*n/len(both_mapped) if both_mapped else 0
    print(f'  {tag:10s} agree: {n}/{len(both_mapped)} ({pct:.1f}%)')

# MAPQ
n = sum(1 for k in both_mapped if bt[k]['mapq'] == rt[k]['mapq'])
pct = 100*n/len(both_mapped) if both_mapped else 0
print(f'  MAPQ       agree: {n}/{len(both_mapped)} ({pct:.1f}%)')

# Sample a few divergences
divergent = [k for k in both_mapped if bt[k]['pos'] != rt[k]['pos'] or bt[k]['cigar'] != rt[k]['cigar']]
print(f'\n  divergent both-mapped: {len(divergent)}')
for k in divergent[:3]:
    print(f'    {k[0]} (r1={k[1]}):')
    print(f'      bt2: pos={bt[k]["pos"]} cigar={bt[k]["cigar"]} mapq={bt[k]["mapq"]} as={bt[k]["tags"].get("AS")}')
    print(f'      rt:  pos={rt[k]["pos"]} cigar={rt[k]["cigar"]} mapq={rt[k]["mapq"]} as={rt[k]["tags"].get("AS")}')
PY

echo
echo "=== MAPQ disagreement breakdown ==="
python3 "$OLDPWD/scripts/mapq_diff.py" bt.sam rt.sam --read-len 100
