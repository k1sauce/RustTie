//! Joint paired-mode descent — a first port of BT2's `extendSeedsPaired`
//! (`vendor/bowtie2/aligner_sw_driver.cpp:1582`).
//!
//! Difference from `align::align_read_with_descent` used on each mate
//! independently: this function takes BOTH mates' reads + qualities,
//! collects seed anchors for *both* mates and merges them into a single
//! priority queue, and for each successfully-extended anchor it
//! immediately attempts mate-rescue of the OTHER mate in the concordance
//! window. Each `(anchor, rescued)` success becomes one entry in the
//! pair-candidate pool — matching the `(rs1_[i], rs2_[i])` parallel-list
//! structure BT2 populates via `AlnSinkWrap::report(rs1, rs2)`
//! (`vendor/bowtie2/aln_sink.cpp:1413`).
//!
//! This addresses the residual ~7.7% MAPQ gap on chr22: BT2's pair pool
//! contains close-alternate pair candidates that emerge from joint
//! extension and don't appear in BT2's output SAM either, so we cannot
//! recover them by post-processing per-mate alignment lists. They must be
//! generated *during* descent.
//!
//! First-pass scope: ungapped extension only, no SW fallback, no
//! re-seeding. Seed-hit cap + `-D` failure budget still apply. Gapped /
//! reseeding will be added incrementally in follow-up sessions.
//!
//! Wire-in: gated by `RUSTTIE_JOINT_DESCENT=1` in [`rusttie-cli`] so we
//! can A/B against the existing independent-then-rescue path without
//! risking regressions.

use std::collections::HashSet;

use rusttie_index::{Bt2Index, BitPairReference};

use crate::align::{
    DEFAULT_SEED_LEN, PrioritizedCandidate, REPETITIVE_HITS_THRESHOLD, Scoring, Strand,
    collect_prioritized, mate_rescue, score_candidate_gapped, score_candidate_ungapped,
    seed_interval,
};
use crate::bt2_random::{RandomSource, ascii_to_bt2_base, gen_rand_seed};
use crate::paired::{FRAG_LEN_MAX, FRAG_LEN_MIN, PairCandidate};
use crate::revcomp::reverse_complement;
use crate::align::Alignment;

/// Cap on pair-pool size — matches BT2's default `mhits + 1 = 51`
/// (`vendor/bowtie2/bt2_search.cpp:343`).
pub const PAIR_POOL_CAP: usize = 51;

/// Output of joint-descent: pair pool (for paired MAPQ), per-mate lists
/// (for unpaired / discordant fallback + secbest threading), and per-mate
/// secbests (for the single-end fallback MAPQ path).
#[derive(Debug, Default)]
pub struct JointDescentResult {
    pub pair_pool: Vec<PairCandidate>,
    pub r1_alns: Vec<Alignment>,
    pub r2_alns: Vec<Alignment>,
    pub r1_secbest: Option<i32>,
    pub r2_secbest: Option<i32>,
}

/// Identifies which mate an anchor candidate belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AnchorMate {
    R1,
    R2,
}

/// Anchor candidate + which mate produced it. Sorted across both mates by
/// `sa_range_size` ascending so the priority queue matches BT2's
/// least-repetitive-first ordering.
#[derive(Debug, Clone, Copy)]
struct JointCandidate {
    mate: AnchorMate,
    cand: PrioritizedCandidate,
}

/// Joint paired-mode descent. See module docs.
#[allow(clippy::too_many_arguments)]
pub fn align_pair_jointly(
    idx: &Bt2Index,
    refs: &BitPairReference,
    r1_seq: &[u8],
    r1_qual: &[u8],
    r1_name: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    r2_name: &[u8],
    scoring: &Scoring,
    seed_hit_cap: u32,
    descent_budget: u32,
    descent_reseed: u32,
) -> JointDescentResult {
    // Precompute reverse-complements + reversed qualities for both mates
    // (used for reverse-strand candidate extensions + mate-rescue).
    let r1_rc = reverse_complement(r1_seq);
    let r2_rc = reverse_complement(r2_seq);
    let r1_qual_rev: Vec<u8> = r1_qual.iter().rev().copied().collect();
    let r2_qual_rev: Vec<u8> = r2_qual.iter().rev().copied().collect();

    let r1_smin = scoring.score_min(r1_seq.len() as u32);
    let r2_smin = scoring.score_min(r2_seq.len() as u32);

    // BT2-faithful per-read PRNG seed (paired: XOR of both mates' seeds).
    // The PRNG is used by the optional bt2_descent path inside
    // `collect_prioritized` to sample large SA ranges in BT2's order.
    // Qname is not threaded down here yet — use empty name; the seed
    // still depends on sequence + quality which is the dominant input.
    let r1_bt2_seq: Vec<u8> = r1_seq.iter().map(|&b| ascii_to_bt2_base(b)).collect();
    let r2_bt2_seq: Vec<u8> = r2_seq.iter().map(|&b| ascii_to_bt2_base(b)).collect();
    let r1_qual_raw: Vec<u8> = r1_qual.iter().map(|&q| q.saturating_sub(33)).collect();
    let r2_qual_raw: Vec<u8> = r2_qual.iter().map(|&q| q.saturating_sub(33)).collect();
    let seed_r1 = gen_rand_seed(&r1_bt2_seq, &r1_qual_raw, r1_name, 0);
    let seed_r2 = gen_rand_seed(&r2_bt2_seq, &r2_qual_raw, r2_name, 0);
    let mut rnd = RandomSource::new(seed_r1 ^ seed_r2);

    let r1_strands: [(Strand, &[u8]); 2] = [
        (Strand::Forward, r1_seq),
        (Strand::Reverse, r1_rc.as_slice()),
    ];
    let r2_strands: [(Strand, &[u8]); 2] = [
        (Strand::Forward, r2_seq),
        (Strand::Reverse, r2_rc.as_slice()),
    ];

    let mut out = JointDescentResult::default();

    // Two-tier dedup: `cand_seen_*` tracks which seed-inferred anchor
    // positions we've already tried to extend (avoids redundant SW work
    // across passes). `aln_seen_*` tracks which final alignment positions
    // landed in the per-mate `*_alns` lists (so the same alignment found
    // by two anchors isn't double-counted). Pair-pool dedup is its own
    // set on the full pair key.
    let mut cand_seen_r1: HashSet<(u32, u32, Strand)> = HashSet::new();
    let mut cand_seen_r2: HashSet<(u32, u32, Strand)> = HashSet::new();
    let mut aln_seen_r1: HashSet<(u32, u32, Strand)> = HashSet::new();
    let mut aln_seen_r2: HashSet<(u32, u32, Strand)> = HashSet::new();
    let mut seen_pairs: HashSet<(u32, u32, Strand, u32, Strand)> = HashSet::new();

    // Mate-rescue redundancy filter: per-(ref_id, anchor_strand) sorted list
    // of anchor ref_offs that already had a mate-rescue attempt. The
    // rescue window for an FR anchor spans `frag_max` bp; two anchors
    // within `frag_max` cover overlapping windows and one rescue is
    // sufficient. This is the structural cost driver — at high
    // `--seed-hit-cap` / `-D` we extend hundreds of anchors per read, and
    // running SW mate-rescue on each was driving the 47s wall time.
    // Skipping redundant windows is the analog of BT2's `redMate1_` /
    // `redMate2_` per-cell dedup applied at the *rescue-window* layer
    // (where it actually saves work for ungapped anchors).
    // Hard cap on total mate-rescue attempts per read (across both mates,
    // across passes). Mirrors BT2's `maxDp` knob in `extendSeedsPaired`.
    // Without it, hi-cap settings spend ~50× the wall time of the
    // default-cap path doing SW on every successfully-extended anchor
    // that lands at a distinct genomic locus.
    //
    // 50 is the knee from the chr22 sweep — MAPQ plateaus at this value
    // (matches K=200's 92.6%) while wall time keeps growing past it.
    // Override at runtime with `RUSTTIE_MAX_RESCUE=<n>` for tuning.
    let max_rescue_attempts: usize = std::env::var("RUSTTIE_MAX_RESCUE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(50);
    let mut rescue_attempts: usize = 0;

    let mut r1_best: i32 = i32::MIN;
    let mut r1_secbest: i32 = i32::MIN;
    let mut r2_best: i32 = i32::MIN;
    let mut r2_secbest: i32 = i32::MIN;

    // Accumulator for cross-pass Phase 2 (gapped SW) fallback. Mirrors
    // `align_read_with_descent`'s pattern: only iterate gapped if Phase 1
    // ungapped found nothing for both mates.
    let mut all_cands: Vec<JointCandidate> = Vec::new();

    // BT2 `-R` re-seeding: each pass shifts seed offsets so a fresh set
    // of seeds covers the read. Mirrors `align_read_with_descent` for
    // parity. Per-pass-`-D`-budget terminates within a pass; re-seed gate
    // (BT2's `repetitive` flag) decides whether to launch the next pass.
    let n_passes = descent_reseed + 1;
    let interval = seed_interval((r1_seq.len() as u32).min(r2_seq.len() as u32));

    'pass_loop: for pass in 0..n_passes {
        let shift = if descent_reseed == 0 {
            0
        } else {
            pass * interval / n_passes
        };

        let mut pass_cands: Vec<JointCandidate> = Vec::new();
        let mut cap_fired = false;
        let mut pass_total_hits: u64 = 0;
        let mut pass_aligned_seeds: u32 = 0;

        let mut buf: Vec<PrioritizedCandidate> = Vec::new();
        for (strand, query) in &r1_strands {
            buf.clear();
            collect_prioritized(
                idx,
                query,
                DEFAULT_SEED_LEN,
                shift,
                seed_hit_cap,
                *strand,
                &mut buf,
                &mut pass_total_hits,
                &mut pass_aligned_seeds,
                &mut cap_fired,
                Some(&mut rnd),
            );
            for c in buf.drain(..) {
                pass_cands.push(JointCandidate {
                    mate: AnchorMate::R1,
                    cand: c,
                });
            }
        }
        for (strand, query) in &r2_strands {
            buf.clear();
            collect_prioritized(
                idx,
                query,
                DEFAULT_SEED_LEN,
                shift,
                seed_hit_cap,
                *strand,
                &mut buf,
                &mut pass_total_hits,
                &mut pass_aligned_seeds,
                &mut cap_fired,
                Some(&mut rnd),
            );
            for c in buf.drain(..) {
                pass_cands.push(JointCandidate {
                    mate: AnchorMate::R2,
                    cand: c,
                });
            }
        }

        pass_cands.sort_by_key(|jc| (jc.cand.sa_range_size, jc.cand.ref_id, jc.cand.ref_off));

        let mut consecutive_failures: u32 = 0;
        for jc in &pass_cands {
            if consecutive_failures >= descent_budget {
                break;
            }
            if out.pair_pool.len() >= PAIR_POOL_CAP {
                break 'pass_loop;
            }

            let (anchor_seq, anchor_qual, anchor_smin) = match (jc.mate, jc.cand.strand) {
                (AnchorMate::R1, Strand::Forward) => (r1_seq, r1_qual, r1_smin),
                (AnchorMate::R1, Strand::Reverse) => {
                    (r1_rc.as_slice(), r1_qual_rev.as_slice(), r1_smin)
                }
                (AnchorMate::R2, Strand::Forward) => (r2_seq, r2_qual, r2_smin),
                (AnchorMate::R2, Strand::Reverse) => {
                    (r2_rc.as_slice(), r2_qual_rev.as_slice(), r2_smin)
                }
            };

            // Skip if we've already tried to extend this anchor candidate
            // (across passes — re-seeding shifts can rediscover the same
            // seed-inferred position).
            let cand_key = (jc.cand.ref_id, jc.cand.ref_off, jc.cand.strand);
            let cand_seen = match jc.mate {
                AnchorMate::R1 => &mut cand_seen_r1,
                AnchorMate::R2 => &mut cand_seen_r2,
            };
            if !cand_seen.insert(cand_key) {
                continue;
            }

            let anchor_aln = score_candidate_ungapped(
                refs,
                jc.cand.ref_id,
                jc.cand.ref_off,
                anchor_seq,
                anchor_qual,
                jc.cand.strand,
                scoring,
                anchor_smin,
            );
            let Some(anchor_aln) = anchor_aln else {
                consecutive_failures += 1;
                continue;
            };

            match jc.mate {
                AnchorMate::R1 => {
                    update_score_window(anchor_aln.score, &mut r1_best, &mut r1_secbest)
                }
                AnchorMate::R2 => {
                    update_score_window(anchor_aln.score, &mut r2_best, &mut r2_secbest)
                }
            };

            // Append to per-mate list (deduped by final alignment key).
            let aln_key = (anchor_aln.ref_id, anchor_aln.ref_off, anchor_aln.strand);
            match jc.mate {
                AnchorMate::R1 => {
                    if aln_seen_r1.insert(aln_key) {
                        out.r1_alns.push(anchor_aln.clone());
                    }
                }
                AnchorMate::R2 => {
                    if aln_seen_r2.insert(aln_key) {
                        out.r2_alns.push(anchor_aln.clone());
                    }
                }
            }

            // Mate-rescue the OTHER mate (BT2's `extendSeedsPaired`
            // inner mate-find step). Bounded by `MAX_RESCUE_ATTEMPTS` so
            // hi-cap settings don't pay ~50× the wall time running SW on
            // every successfully-extended anchor at a distinct repetitive
            // locus. A window-overlap dedup was tried but cut into
            // legitimately-distinct pair candidates — anchors A and B
            // within `frag_max` of each other can still produce two
            // distinct `(A, rescued)` and `(B, rescued)` pair candidates
            // even when their rescue *windows* overlap.
            if rescue_attempts >= max_rescue_attempts {
                continue;
            }
            rescue_attempts += 1;

            let (other_seq, other_qual, other_rc, other_qual_rev) = match jc.mate {
                AnchorMate::R1 => (r2_seq, r2_qual, r2_rc.as_slice(), r2_qual_rev.as_slice()),
                AnchorMate::R2 => (r1_seq, r1_qual, r1_rc.as_slice(), r1_qual_rev.as_slice()),
            };
            let rescued = mate_rescue(
                refs,
                &anchor_aln,
                other_seq,
                other_qual,
                other_rc,
                other_qual_rev,
                scoring,
                FRAG_LEN_MIN,
                FRAG_LEN_MAX,
            );
            let Some(rescued) = rescued else {
                // Anchor extension succeeded but no mate found in the
                // concordance window — partial progress. Don't tick the
                // failure budget (the anchor alignment will still
                // participate in the discordant/unpaired fallback).
                continue;
            };

            // Append rescued to OTHER mate's list (deduped).
            match jc.mate {
                AnchorMate::R1 => {
                    let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
                    if aln_seen_r2.insert(key) {
                        out.r2_alns.push(rescued.clone());
                    }
                    update_score_window(rescued.score, &mut r2_best, &mut r2_secbest);
                }
                AnchorMate::R2 => {
                    let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
                    if aln_seen_r1.insert(key) {
                        out.r1_alns.push(rescued.clone());
                    }
                    update_score_window(rescued.score, &mut r1_best, &mut r1_secbest);
                }
            }

            let (r1_aln, r2_aln) = match jc.mate {
                AnchorMate::R1 => (anchor_aln, rescued),
                AnchorMate::R2 => (rescued, anchor_aln),
            };
            let pair_key = (
                r1_aln.ref_id,
                r1_aln.ref_off,
                r1_aln.strand,
                r2_aln.ref_off,
                r2_aln.strand,
            );
            if seen_pairs.insert(pair_key)
                && let Some(cand) = PairCandidate::try_new(r1_aln, r2_aln)
            {
                out.pair_pool.push(cand);
                consecutive_failures = 0;
            }
            // Don't tick `consecutive_failures` for a duplicate pair —
            // it's a re-discovery, not a missed extension.
        }

        // Save this pass's candidates so the Phase 2 gapped fallback can
        // re-iterate them with SW DP if Phase 1 found nothing.
        all_cands.extend(pass_cands);

        // BT2's re-seed gate: continue iterating passes only if we don't
        // already have a non-repetitive set of paired hits. Use the same
        // threshold as `align_read_with_descent` for consistency.
        let avg_hits = if pass_aligned_seeds > 0 {
            pass_total_hits / pass_aligned_seeds as u64
        } else {
            0
        };
        let repetitive = avg_hits > REPETITIVE_HITS_THRESHOLD as u64;
        if !out.pair_pool.is_empty() && !repetitive {
            break;
        }
    }

    // Phase 2 (SW rescue): if Phase 1 found nothing for either mate, try
    // gapped extension on every collected candidate. Handles indel reads
    // where ungapped Hamming can't score. Mirrors
    // `align_read_with_descent`'s Phase 2 logic.
    if out.r1_alns.is_empty() && out.r2_alns.is_empty() {
        for jc in &all_cands {
            let (anchor_seq, anchor_qual, anchor_smin) = match (jc.mate, jc.cand.strand) {
                (AnchorMate::R1, Strand::Forward) => (r1_seq, r1_qual, r1_smin),
                (AnchorMate::R1, Strand::Reverse) => {
                    (r1_rc.as_slice(), r1_qual_rev.as_slice(), r1_smin)
                }
                (AnchorMate::R2, Strand::Forward) => (r2_seq, r2_qual, r2_smin),
                (AnchorMate::R2, Strand::Reverse) => {
                    (r2_rc.as_slice(), r2_qual_rev.as_slice(), r2_smin)
                }
            };
            let Some(anchor_aln) = score_candidate_gapped(
                refs,
                jc.cand.ref_id,
                jc.cand.ref_off,
                anchor_seq,
                anchor_qual,
                jc.cand.strand,
                scoring,
            ) else {
                continue;
            };
            if anchor_aln.score < anchor_smin {
                continue;
            }

            let aln_key = (anchor_aln.ref_id, anchor_aln.ref_off, anchor_aln.strand);
            match jc.mate {
                AnchorMate::R1 => {
                    update_score_window(anchor_aln.score, &mut r1_best, &mut r1_secbest);
                    if aln_seen_r1.insert(aln_key) {
                        out.r1_alns.push(anchor_aln.clone());
                    }
                }
                AnchorMate::R2 => {
                    update_score_window(anchor_aln.score, &mut r2_best, &mut r2_secbest);
                    if aln_seen_r2.insert(aln_key) {
                        out.r2_alns.push(anchor_aln.clone());
                    }
                }
            }

            // Also attempt mate-rescue from this gapped anchor — gives a
            // pair candidate for indel-only reads. Reuses the rescue cap.
            if rescue_attempts >= max_rescue_attempts {
                continue;
            }
            rescue_attempts += 1;
            let (other_seq, other_qual, other_rc, other_qual_rev) = match jc.mate {
                AnchorMate::R1 => (r2_seq, r2_qual, r2_rc.as_slice(), r2_qual_rev.as_slice()),
                AnchorMate::R2 => (r1_seq, r1_qual, r1_rc.as_slice(), r1_qual_rev.as_slice()),
            };
            let Some(rescued) = mate_rescue(
                refs,
                &anchor_aln,
                other_seq,
                other_qual,
                other_rc,
                other_qual_rev,
                scoring,
                FRAG_LEN_MIN,
                FRAG_LEN_MAX,
            ) else {
                continue;
            };
            match jc.mate {
                AnchorMate::R1 => {
                    let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
                    if aln_seen_r2.insert(key) {
                        out.r2_alns.push(rescued.clone());
                    }
                    update_score_window(rescued.score, &mut r2_best, &mut r2_secbest);
                }
                AnchorMate::R2 => {
                    let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
                    if aln_seen_r1.insert(key) {
                        out.r1_alns.push(rescued.clone());
                    }
                    update_score_window(rescued.score, &mut r1_best, &mut r1_secbest);
                }
            }
            let (r1_aln, r2_aln) = match jc.mate {
                AnchorMate::R1 => (anchor_aln, rescued),
                AnchorMate::R2 => (rescued, anchor_aln),
            };
            let pair_key = (
                r1_aln.ref_id,
                r1_aln.ref_off,
                r1_aln.strand,
                r2_aln.ref_off,
                r2_aln.strand,
            );
            if seen_pairs.insert(pair_key)
                && let Some(cand) = PairCandidate::try_new(r1_aln, r2_aln)
            {
                out.pair_pool.push(cand);
            }
        }
    }

    // Sort per-mate lists score-descending, leftmost-first — matches the
    // ordering invariant `classify_pair_set`'s Cartesian fallback assumes.
    let sort_alns = |alns: &mut Vec<Alignment>| {
        alns.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.ref_id.cmp(&b.ref_id))
                .then(a.ref_off.cmp(&b.ref_off))
        });
    };
    sort_alns(&mut out.r1_alns);
    sort_alns(&mut out.r2_alns);

    out.r1_secbest = if r1_secbest > i32::MIN { Some(r1_secbest) } else { None };
    out.r2_secbest = if r2_secbest > i32::MIN { Some(r2_secbest) } else { None };

    // Sort + cap the pair pool — BT2's mhits+1 limit on `rs1_`/`rs2_`.
    out.pair_pool.sort_by(|a, b| b.score_sum.cmp(&a.score_sum));
    out.pair_pool.truncate(PAIR_POOL_CAP);

    out
}

/// Fold a new score `s` into running `(best, secbest)`. Same semantics as
/// `align::update_score_window` — duplicated here because that helper is
/// private and we need the same logic on the joint-descent tracking.
#[inline]
fn update_score_window(s: i32, best: &mut i32, secbest: &mut i32) {
    if s > *best {
        if *best != i32::MIN {
            *secbest = *best;
        }
        *best = s;
    } else if (s == *best && *secbest < *best) || (s < *best && s > *secbest) {
        *secbest = s;
    }
}
