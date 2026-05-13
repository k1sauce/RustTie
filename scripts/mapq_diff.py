#!/usr/bin/env python3
# MAPQ disagreement diagnostic.
#
# Reads two paired-end SAMs aligned from the same reads (BT2 oracle + RustTie)
# and characterizes WHERE MAPQ disagreement comes from. Split by whether AS
# agrees: AS-agree + MAPQ-disagree is pure alternate-set drift (secbest); AS
# disagree means the chosen alignment itself differs, which is a different bug.
#
# Usage:
#   python3 scripts/mapq_diff.py bt.sam rt.sam            # chr22 paths
#   python3 scripts/mapq_diff.py --read-len 100 a.sam b.sam
#
# The "both_low" bucket the README calls out is rows where both MAPQ < 30.
# Paired-end MAPQ in BT2's BowtieMapq2 (monotone=true, --score-min L,-0.6,-0.6,
# 100bp pairs) has pair_diff = 120 (= -2 * mate_score_min), pair_best =
# AS_R1 + AS_R2, pair_bestOver = pair_best + 120. From a given MAPQ bin we
# can back out the secbest range BT2 must have observed; comparing to our
# bin tells us how the alternate-set sizes differ.

from __future__ import annotations

import argparse
import collections
import sys
from typing import Iterable


def parse_sam(path: str) -> dict[tuple[str, bool], dict]:
    rows: dict[tuple[str, bool], dict] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            qname = cols[0]
            flag = int(cols[1])
            if flag & 4:
                continue  # unmapped — no MAPQ to compare
            is_r1 = (flag & 0x40) != 0
            tags: dict[str, str] = {}
            for t in cols[11:]:
                k, _, v = t.partition(":")
                tags[k] = v
            rows[(qname, is_r1)] = {
                "flag": flag,
                "rname": cols[2],
                "pos": int(cols[3]),
                "mapq": int(cols[4]),
                "cigar": cols[5],
                "tlen": int(cols[8]),
                "as": int(tags["AS"].split(":")[-1]) if "AS" in tags else None,
            }
    return rows


# BT2 bin table for paired-end end-to-end mode. Each entry:
# (bestdiff_min_frac, bestover_min_frac, mapq).
# bestdiff bins are descending; within each bestdiff bin, bestover bins are
# descending. None for bestover_min_frac means "no secbest" (this column unused
# when we have both AS values). Top-of-bin uses bestOver == diff in BT2 source,
# which equals bestOver >= diff in end-to-end (scores <= 0 cap bestOver at diff).
#
# Format: (bestdiff_lo, bins_within), where bins_within is
# [(bestover_lo, mapq), ...] sorted bestover_lo descending. 1.0 is the
# "== diff" top-of-bin sentinel.
END_TO_END_BINS_PAIRED = [
    (0.9, [(1.0, 39), (0.0, 33)]),
    (0.8, [(1.0, 38), (0.0, 27)]),
    (0.7, [(1.0, 37), (0.0, 26)]),
    (0.6, [(1.0, 36), (0.0, 22)]),
    (0.5, [(1.0, 35), (0.84, 25), (0.68, 16), (0.0, 5)]),
    (0.4, [(1.0, 34), (0.84, 21), (0.68, 14), (0.0, 4)]),
    (0.3, [(1.0, 32), (0.88, 18), (0.67, 15), (0.0, 3)]),
    (0.2, [(1.0, 31), (0.88, 17), (0.67, 11), (0.0, 0)]),
    (0.1, [(1.0, 30), (0.88, 12), (0.67, 7), (0.0, 0)]),
    # bestdiff in (0, 0.1*diff): only two leaves.
    (0.0001, [(0.67, 6), (0.0, 2)]),
    # bestdiff == 0
    (0.0, [(0.67, 1), (0.0, 0)]),
]

NO_SECBEST_BINS_PAIRED = [
    (0.8, 42), (0.7, 40), (0.6, 24), (0.5, 23), (0.4, 8), (0.3, 3), (0.0, 0),
]


def bin_for_mapq(mapq: int) -> list[tuple[float, float | None, float | None]]:
    """Return list of (bestdiff_lo_frac, bestover_lo_frac, bestover_hi_frac)
    bins of `diff` that BT2 would have to be in to emit this MAPQ. Used to
    back-solve the implied secbest range.

    Each returned tuple is *inclusive lo, exclusive hi* in fractions of diff.
    A None bestover_hi means "<= 1.0" (top of bin). Multiple bins are returned
    because the same MAPQ can come from different bestdiff buckets in rare
    cases (it doesn't, actually, but we keep the API general).
    """
    out = []
    # Compute upper bound for each bestdiff bin from the entry above it.
    prev_diff_lo = 1.0  # top sentinel (bestdiff cannot exceed diff)
    for bd_lo, leaves in END_TO_END_BINS_PAIRED:
        # leaves sorted bestover_lo descending; compute hi from the entry above
        prev_bo_lo = 1.0001  # exclusive sentinel above 1.0
        for bo_lo, m in leaves:
            if m == mapq:
                out.append((bd_lo, prev_diff_lo, bo_lo, prev_bo_lo))
            prev_bo_lo = bo_lo
        prev_diff_lo = bd_lo
    return out


def implied_secbest_score(
    pair_best: int, pair_smin: int, mapq: int
) -> tuple[int, int] | None:
    """Given a paired pair_best and pair_smin, return the (lo, hi) range of
    pair_secbest scores (inclusive lo, exclusive hi) that would produce this
    MAPQ. None if MAPQ unreachable. Returns the union if multiple bins match
    (first match — they don't overlap in practice).
    """
    diff = max(1, 0 - pair_smin)
    best_over = pair_best - pair_smin
    best_over_frac = best_over / diff if diff > 0 else 0.0
    bins = bin_for_mapq(mapq)
    for bd_lo_f, bd_hi_f, bo_lo_f, bo_hi_f in bins:
        # We need bestover_frac in [bo_lo_f, bo_hi_f). If not, skip.
        if not (bo_lo_f <= best_over_frac < bo_hi_f):
            continue
        # bestdiff range maps to secbest range:
        #   bestdiff = best - secbest
        #   secbest = best - bestdiff
        # bestdiff in [bd_lo_f*diff, bd_hi_f*diff) → secbest in
        #   (best - bd_hi_f*diff, best - bd_lo_f*diff]
        sec_hi = pair_best - bd_lo_f * diff  # inclusive
        sec_lo = pair_best - bd_hi_f * diff  # exclusive
        return (int(sec_lo) + 1, int(sec_hi) + 1)  # convert to int-inclusive lo, exclusive hi
    return None


def fmt_pct(n: int, d: int) -> str:
    return f"{n}/{d} ({100*n/d:.1f}%)" if d else "0/0 (n/a)"


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bt_sam", help="Bowtie 2 oracle SAM")
    ap.add_argument("rt_sam", help="RustTie SAM")
    ap.add_argument(
        "--read-len",
        type=int,
        default=100,
        help="Read length in bp (paired symmetric). Used to compute pair_smin "
        "from BT2 default --score-min L,-0.6,-0.6.",
    )
    ap.add_argument(
        "--samples",
        type=int,
        default=5,
        help="qnames to print per bucket for spot-checking.",
    )
    args = ap.parse_args(argv[1:])

    # BT2 default --score-min L,-0.6,-0.6 → per-mate smin = floor(-0.6 + -0.6 * rdlen)
    # truncates toward zero in BT2's TAlScore cast. For 100bp: -60.6 → -60.
    per_mate_smin = int(-0.6 + -0.6 * args.read_len)
    pair_smin = per_mate_smin * 2

    print(f"== MAPQ diff: bt={args.bt_sam} rt={args.rt_sam}")
    print(f"   per-mate smin={per_mate_smin}, pair_smin={pair_smin}, pair_diff={max(1, -pair_smin)}")

    bt = parse_sam(args.bt_sam)
    rt = parse_sam(args.rt_sam)

    # Build paired pair_best (AS sum) per qname for back-solving secbest.
    def pair_score(rows: dict, qname: str) -> int | None:
        a = rows.get((qname, True))
        b = rows.get((qname, False))
        if a is None or b is None or a["as"] is None or b["as"] is None:
            return None
        return a["as"] + b["as"]

    # Only consider records mapped in BOTH bt and rt.
    keys = set(bt) & set(rt)
    print(f"   records mapped in both: {len(keys)}")

    agree = 0
    disagree_same_as = []  # qname, is_r1
    disagree_diff_as = []
    for k in keys:
        if bt[k]["mapq"] == rt[k]["mapq"]:
            agree += 1
            continue
        if bt[k]["as"] == rt[k]["as"]:
            disagree_same_as.append(k)
        else:
            disagree_diff_as.append(k)
    total_disagree = len(disagree_same_as) + len(disagree_diff_as)
    print(f"   MAPQ agree:                 {fmt_pct(agree, len(keys))}")
    print(f"   MAPQ disagree (total):      {total_disagree}")
    print(f"     ├─ AS agrees (secbest):   {fmt_pct(len(disagree_same_as), total_disagree)}")
    print(f"     └─ AS disagrees (chosen): {fmt_pct(len(disagree_diff_as), total_disagree)}")

    # The README hint: "both_low" — both MAPQ < 30 (the multi-mapper regime).
    both_low = [k for k in disagree_same_as if bt[k]["mapq"] < 30 and rt[k]["mapq"] < 30]
    both_high = [k for k in disagree_same_as if bt[k]["mapq"] >= 30 and rt[k]["mapq"] >= 30]
    mixed = [k for k in disagree_same_as if k not in both_low and k not in both_high]
    print(f"\n   Within AS-agree disagreements:")
    print(f"     both_low  (<30 each): {fmt_pct(len(both_low), len(disagree_same_as))}")
    print(f"     both_high (>=30):     {fmt_pct(len(both_high), len(disagree_same_as))}")
    print(f"     mixed (one each):     {fmt_pct(len(mixed), len(disagree_same_as))}")

    # Confusion matrix on AS-agree disagreements: rt_mapq -> bt_mapq counts.
    cm = collections.Counter((rt[k]["mapq"], bt[k]["mapq"]) for k in disagree_same_as)
    print(f"\n   Top (rt_mapq, bt_mapq) bins on AS-agree disagreements:")
    for (rt_m, bt_m), n in cm.most_common(15):
        print(f"     rt={rt_m:>3d}  bt={bt_m:>3d}  n={n}")

    # Spot-check sample qnames for both_low, with implied secbest ranges.
    print(f"\n   both_low samples (AS-agree, both MAPQ<30) — qname, AS_pair, implied secbest ranges:")
    for k in both_low[: args.samples]:
        qname, is_r1 = k
        pair_best = pair_score(bt, qname)
        if pair_best is None:
            continue
        bt_range = implied_secbest_score(pair_best, pair_smin, bt[k]["mapq"])
        rt_range = implied_secbest_score(pair_best, pair_smin, rt[k]["mapq"])
        print(
            f"     {qname:30s} r1={is_r1!s:5s} pair_AS={pair_best:>4d}  "
            f"rt_mapq={rt[k]['mapq']:>3d} (implies pair_sec ∈ {rt_range})  "
            f"bt_mapq={bt[k]['mapq']:>3d} (implies pair_sec ∈ {bt_range})"
        )

    # Diff-direction: are we systematically too high, too low, or noisy?
    delta = collections.Counter()
    for k in disagree_same_as:
        delta[rt[k]["mapq"] - bt[k]["mapq"]] += 1
    rt_higher = sum(n for d, n in delta.items() if d > 0)
    rt_lower = sum(n for d, n in delta.items() if d < 0)
    print(
        f"\n   Direction (AS-agree):  rt > bt: {rt_higher},  rt < bt: {rt_lower}"
        f"  → {'rt confident-er' if rt_higher > rt_lower else 'rt conservative-er'}"
    )

    # AS-disagree disagreements deserve to be flagged separately — they are
    # NOT the README's "both_low" bucket; they indicate the chosen alignment
    # differs, which is upstream of MAPQ.
    if disagree_diff_as:
        print(f"\n   AS-disagree sample (different alignment chosen):")
        for k in disagree_diff_as[: args.samples]:
            qname, is_r1 = k
            print(
                f"     {qname:30s} r1={is_r1!s:5s}  "
                f"rt AS={rt[k]['as']:>4d} cigar={rt[k]['cigar']:<6s} mapq={rt[k]['mapq']:>3d}  "
                f"bt AS={bt[k]['as']:>4d} cigar={bt[k]['cigar']:<6s} mapq={bt[k]['mapq']:>3d}"
            )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
