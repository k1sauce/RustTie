//! BT2-faithful seed prioritization, the next building block toward
//! closing the chr22 MAPQ gap (#3, see `rusttie.md` Phase 0).
//!
//! Ports the algorithmic shape of `SwDriver::prioritizeSATupsRands` from
//! `vendor/bowtie2/aligner_sw_driver.cpp:492-738`. This is the function
//! that decides — given a set of seed hits — *which* SA-range rows to
//! enumerate and in what order before the descent loop's `-D` budget
//! fires.
//!
//! The algorithm has two phases:
//!
//! 1. **Smalls** — SA ranges with `size ≤ nsm` (BT2 default 5) are
//!    processed exhaustively. All rows resolved.
//! 2. **Non-smalls** — large SA ranges are sampled one row at a time
//!    via a weighted [`RowSampler`] (chooses *which* range to sample
//!    next, weighted toward smaller / more-specific ranges) plus a
//!    per-range [`Random1toN`] (chooses *which* row in that range).
//!
//! This shape diverges fundamentally from our current
//! `align::collect_prioritized` which uses a `seed_hit_cap` cutoff and
//! drops large ranges entirely. The cutoff is empirically tuned for
//! performance but produces a different alignment-candidate population
//! than BT2, which is the root cause of the residual MAPQ gap.
//!
//! This module provides the algorithm in isolation; integration into
//! the wider descent pipeline is gated on `#3` Phase 1 completion.

use crate::bt2_random::{Random1toN, RandomSource, RowSampler};

/// Input to [`rank_seed_hits`]: per-seed-offset SA range sizes for each
/// strand. `None` (or `Some(0)`) means "no hits at this offset/strand".
/// Length of both vectors must equal `n_offs`.
#[derive(Debug, Clone)]
pub struct SeedRankInput {
    pub n_offs: usize,
    /// `fw_sizes[i]` = SA range size for forward seed at offset index `i`.
    pub fw_sizes: Vec<Option<u32>>,
    /// `rc_sizes[i]` = SA range size for revcomp seed at offset index `i`.
    pub rc_sizes: Vec<Option<u32>>,
}

/// Port of `SeedResults::rankSeedHits` (`aligner_seed.h:1019-1080`).
/// Returns `(offset_index, fw_strand)` tuples in BT2's exact iteration
/// order: smallest-SA-range first, with PRNG-driven tie-breaking via
/// `nextBool` (strand preference) and `nextU32` (wrap-around scan start).
///
/// Critically advances the PRNG state in the exact same call sequence
/// as BT2 — every rank step does 1 nextBool + 2 nextU32. Downstream
/// `Random1toN` consumption then aligns with BT2's by construction.
pub fn rank_seed_hits(
    input: &SeedRankInput,
    rnd: &mut RandomSource,
) -> Vec<(usize, bool)> {
    let n = input.n_offs;
    assert_eq!(input.fw_sizes.len(), n);
    assert_eq!(input.rc_sizes.len(), n);
    let mut sorted_fw = vec![false; n];
    let mut sorted_rc = vec![false; n];

    let nonempty = |s: &Option<u32>| s.map(|x| x > 0).unwrap_or(false);
    let total: usize = input
        .fw_sizes
        .iter()
        .filter(|s| nonempty(s))
        .count()
        + input.rc_sizes.iter().filter(|s| nonempty(s)).count();

    let mut out: Vec<(usize, bool)> = Vec::with_capacity(total);
    if n == 0 || total == 0 {
        return out;
    }

    let dump = std::env::var_os("RUSTTIE_DUMP_RANK").is_some();
    while out.len() < total {
        let mut min_sz = u32::MAX;
        let mut min_idx: usize = 0;
        let mut min_fw = true;
        let pre_last = rnd.last();
        let rb = rnd.next_bool();
        let post_bool_last = rnd.last();
        if dump {
            eprintln!(
                "[rt-rank step={} pre_last={} rb={} post_bool_last={}]",
                out.len(),
                pre_last,
                rb as u32,
                post_bool_last
            );
        }
        for fwi in 0..2 {
            let fw = fwi == (if rb { 1 } else { 0 });
            let sizes = if fw { &input.fw_sizes } else { &input.rc_sizes };
            let sorted = if fw { &sorted_fw } else { &sorted_rc };
            let pre_u32_last = rnd.last();
            let nu = rnd.next_u32();
            let post_u32_last = rnd.last();
            let i_start = (nu as usize) % n;
            if dump {
                eprintln!(
                    "[rt-rank step={} fwi={} fw={} pre_u32_last={} nu={} post_u32_last={} n={} i_start={}]",
                    out.len(),
                    fwi,
                    fw,
                    pre_u32_last,
                    nu,
                    post_u32_last,
                    n,
                    i_start
                );
            }
            let mut i = i_start;
            for _ in 0..n {
                if let Some(sz) = sizes[i]
                    && sz > 0
                    && !sorted[i]
                    && sz < min_sz
                {
                    min_sz = sz;
                    min_idx = i;
                    min_fw = fw;
                }
                i = if i + 1 == n { 0 } else { i + 1 };
            }
        }
        debug_assert!(min_sz != u32::MAX, "should always find a remaining seed");
        if dump {
            eprintln!(
                "[rt-rank step={} picked={{idx={},fw={}}}]",
                out.len(),
                min_idx,
                min_fw
            );
        }
        if min_fw {
            sorted_fw[min_idx] = true;
        } else {
            sorted_rc[min_idx] = true;
        }
        out.push((min_idx, min_fw));
    }
    out
}

/// A single seed hit in the suffix array: an SA range `[sa_lo, sa_hi)`
/// with metadata about which seed produced it. Equivalent to BT2's
/// `SATupleAndPos` (`aligner_sw_driver.h` near `SATupleAndPos`).
#[derive(Debug, Clone)]
pub struct SeedHit {
    /// Inclusive start of the SA range.
    pub sa_lo: u32,
    /// Exclusive end of the SA range.
    pub sa_hi: u32,
    /// Seed offset in the read (0-based, from the read's 5' end).
    pub rdoff: u32,
    /// Seed length (typically 22 for BT2 default).
    pub seedlen: u32,
    /// Whether the seed matched the forward (fw=true) or reverse-complement
    /// (fw=false) of the read.
    pub fw: bool,
    /// Per-seed extension info — number of additional matching bases to
    /// the left/right of the seed. Used as the weight numerator
    /// (`(nlex + nrex + 1)^2`) in [`RowSampler`]. Set to 0 if no
    /// extension is performed.
    pub nlex: u32,
    pub nrex: u32,
}

impl SeedHit {
    pub fn size(&self) -> u32 {
        self.sa_hi - self.sa_lo
    }

    /// Extension length used by [`RowSampler`] weight: `nlex + nrex + 1`.
    pub fn extended_len(&self) -> u32 {
        self.nlex + self.nrex + 1
    }
}

/// A single resolved sample from [`prioritize_sa_tups_rands`] — a
/// specific row in a specific SA range, along with the originating seed
/// metadata. Caller resolves `sa_row` to a text position via
/// `resolve_text_pos` and extends from there.
#[derive(Debug, Clone, Copy)]
pub struct PrioritizedRow {
    /// BWT row to resolve to a text offset.
    pub sa_row: u32,
    /// Seed offset in the read (passed through from `SeedHit::rdoff`).
    pub rdoff: u32,
    pub seedlen: u32,
    pub fw: bool,
    /// Size of the SA range this row was sampled from (after dedup of
    /// covered seeds). Useful for downstream priority decisions.
    pub sa_range_size: u32,
}

/// Port of BT2's `prioritizeSATupsRands` (`aligner_sw_driver.cpp:492-738`).
/// Returns a list of rows to process in BT2-equivalent order: all rows
/// from "small" ranges first (in sorted-by-size order), then weighted
/// random sampling from large ranges until `maxelt` rows have been
/// emitted or all ranges are exhausted.
///
/// **Parameters match BT2's defaults:**
/// - `nsm = 5`: SA ranges with size ≤ 5 are "small" and processed whole.
/// - `lensq = true, szsq = true`: square both terms in the weight.
/// - `maxelt`: cap on total rows emitted.
///
/// This function does *not* do seed extension or covered-seed filtering
/// (BT2's `extend()` + `ExtendRange` logic). Those would be added when
/// integrating with our actual descent loop; the algorithmic shape here
/// captures the row-sampling order which is the part that drives the
/// MAPQ gap.
pub fn prioritize_sa_tups_rands(
    seeds: Vec<SeedHit>,
    nsm: u32,
    lensq: bool,
    szsq: bool,
    maxelt: usize,
    rnd: &mut RandomSource,
) -> Vec<PrioritizedRow> {
    let mut out: Vec<PrioritizedRow> = Vec::new();
    if seeds.is_empty() || maxelt == 0 {
        return out;
    }

    // BT2's `satpos.sort()` (`aligner_sw_driver.cpp:612`) is a stable
    // size-only sort. Mirroring that exactly regressed the legacy path
    // from 94.4% → 94.3% MAPQ (-8 chosen-disagree); the rnd-tagged sort
    // empirically samples tie ordering in a way that better matches
    // BT2's net pool composition once the legacy mate-rescue loop runs.
    // Keeping the rnd-tagged version pending the row-sampling port that
    // would make the byte-exact variant net-positive.
    let mut tagged: Vec<(SeedHit, u32)> =
        seeds.into_iter().map(|s| (s, rnd.next_u32())).collect();
    tagged.sort_by_key(|(s, tag)| (s.size(), *tag));
    let sorted: Vec<SeedHit> = tagged.into_iter().map(|(s, _)| s).collect();

    let nsmall = sorted.iter().take_while(|s| s.size() <= nsm).count();

    // Phase 1 — smalls: emit all rows from each small range. BT2 uses
    // Random1toN even within small ranges (`aligner_sw_driver.cpp:1859`
    // — `rands_[i].next(rnd)`) but consumes that PRNG only during the
    // extension loop in `extendSeedsPaired`, NOT inside `prioritizeSATupsRands`.
    // To match BT2's PRNG state at the next call site (the next
    // `rankSeedHits`), we must NOT consume PRNG for these small-range
    // picks here — the caller is responsible for doing so per anchor
    // extension if it wants byte-exact match.
    //
    // RUSTTIE_RT_SMALLS_RND=1 reinstates the old "consume PRNG inline"
    // behavior for A/B regression testing against the prior 94.4% baseline.
    let smalls_rand = std::env::var_os("RUSTTIE_RT_SMALLS_RND").is_some();
    for seed in sorted.iter().take(nsmall) {
        if out.len() >= maxelt {
            return out;
        }
        let remaining = maxelt - out.len();
        let sz = seed.size() as usize;
        let take = sz.min(remaining);
        if smalls_rand {
            let mut sampler = Random1toN::new(sz);
            for _ in 0..take {
                let r = sampler.next(rnd);
                out.push(PrioritizedRow {
                    sa_row: seed.sa_lo + r as u32,
                    rdoff: seed.rdoff,
                    seedlen: seed.seedlen,
                    fw: seed.fw,
                    sa_range_size: seed.size(),
                });
            }
        } else {
            // Deterministic BWT-row order. No PRNG consumed.
            for r in 0..take {
                out.push(PrioritizedRow {
                    sa_row: seed.sa_lo + r as u32,
                    rdoff: seed.rdoff,
                    seedlen: seed.seedlen,
                    fw: seed.fw,
                    sa_range_size: seed.size(),
                });
            }
        }
    }
    if out.len() >= maxelt || nsmall == sorted.len() {
        return out;
    }

    // Phase 2 — non-smalls: weighted-random sampling. The RowSampler
    // picks WHICH range to draw from; a per-range Random1toN picks the
    // specific row.
    let large_bins: Vec<(usize, usize)> = sorted
        .iter()
        .skip(nsmall)
        .map(|s| (s.extended_len() as usize, s.size() as usize))
        .collect();
    let mut row_samp = RowSampler::new(large_bins, lensq, szsq);
    let mut row_choosers: Vec<Option<Random1toN>> = vec![None; sorted.len() - nsmall];

    while out.len() < maxelt {
        // Pick a non-small range.
        if row_samp.total_mass() <= 0.0 {
            // All ranges exhausted.
            break;
        }
        let bin = row_samp.next(rnd);
        // Float drift can leave `total_mass > 0` while every bin is
        // eliminated; the sampler then returns its sentinel. Treat as
        // exhausted.
        if bin == usize::MAX {
            break;
        }
        // Lazy-initialize the per-range row chooser.
        if row_choosers[bin].is_none() {
            row_choosers[bin] = Some(Random1toN::new(sorted[nsmall + bin].size() as usize));
        }
        let chooser = row_choosers[bin].as_mut().unwrap();
        debug_assert!(!chooser.done());
        let r = chooser.next(rnd);
        if chooser.done() {
            row_samp.finished_range(bin);
        }
        let seed = &sorted[nsmall + bin];
        out.push(PrioritizedRow {
            sa_row: seed.sa_lo + r as u32,
            rdoff: seed.rdoff,
            seedlen: seed.seedlen,
            fw: seed.fw,
            sa_range_size: seed.size(),
        });
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_seed(sa_lo: u32, sa_hi: u32, rdoff: u32) -> SeedHit {
        SeedHit {
            sa_lo,
            sa_hi,
            rdoff,
            seedlen: 22,
            fw: true,
            nlex: 0,
            nrex: 0,
        }
    }

    /// All-small-range case: every row gets emitted exactly once, and
    /// ranges are processed in size-ascending order (size 1 first, then
    /// size 3, then size 5). Within each range the order is PRNG-driven
    /// so we only check membership, not specific row positions.
    #[test]
    fn smalls_only_emits_all_rows() {
        let mut rnd = RandomSource::new(42);
        let seeds = vec![
            mk_seed(100, 103, 0),  // 3 rows
            mk_seed(200, 205, 22), // 5 rows
            mk_seed(300, 301, 44), // 1 row
        ];
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 100, &mut rnd);
        assert_eq!(out.len(), 9);
        // First entry is the only row from the size-1 range.
        assert_eq!(out[0].sa_row, 300);
        // Next three are some permutation of {100, 101, 102}.
        let mut got_3: Vec<u32> = out[1..4].iter().map(|r| r.sa_row).collect();
        got_3.sort_unstable();
        assert_eq!(got_3, vec![100, 101, 102]);
        // Final five are some permutation of {200..205}.
        let mut got_5: Vec<u32> = out[4..9].iter().map(|r| r.sa_row).collect();
        got_5.sort_unstable();
        assert_eq!(got_5, vec![200, 201, 202, 203, 204]);
    }

    /// `maxelt` truncates output.
    #[test]
    fn maxelt_truncates() {
        let mut rnd = RandomSource::new(42);
        let seeds = vec![mk_seed(0, 10, 0)]; // 10 rows, all small (sz>nsm=5)
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 4, &mut rnd);
        // sz=10 > nsm=5 so this becomes a "non-small". Only `maxelt` rows
        // get emitted, sampled randomly.
        assert_eq!(out.len(), 4);
    }

    /// Large-range case: weighted sampling produces output of size `maxelt`,
    /// each row unique, all from the same range.
    #[test]
    fn large_range_unique_rows() {
        let mut rnd = RandomSource::new(12345);
        let seeds = vec![mk_seed(0, 50, 0)]; // 50 rows, sz > nsm so large
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 10, &mut rnd);
        assert_eq!(out.len(), 10);
        let mut seen: Vec<u32> = out.iter().map(|r| r.sa_row).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 10, "rows must be unique");
        // All in [0, 50).
        assert!(seen.iter().all(|&r| r < 50));
    }

    /// Mixed small + large: small range processed first (exhaustively
    /// but in random per-range order), then large range sampled.
    #[test]
    fn mixed_smalls_then_large() {
        let mut rnd = RandomSource::new(42);
        let seeds = vec![
            mk_seed(0, 100, 0),   // 100 rows, LARGE
            mk_seed(1000, 1003, 22), // 3 rows, SMALL
        ];
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 50, &mut rnd);
        assert_eq!(out.len(), 50);
        // First 3 should be the small range, in some PRNG order.
        let mut got_small: Vec<u32> = out[0..3].iter().map(|r| r.sa_row).collect();
        got_small.sort_unstable();
        assert_eq!(got_small, vec![1000, 1001, 1002]);
        // Remaining 47 from the large range, all in [0, 100).
        for row in &out[3..] {
            assert!(row.sa_row < 100);
        }
    }

    /// Determinism: same seeds + same PRNG state = same output.
    #[test]
    fn deterministic() {
        let collect = |seed: u32| {
            let mut rnd = RandomSource::new(seed);
            let seeds = vec![mk_seed(0, 50, 0), mk_seed(100, 110, 22)];
            prioritize_sa_tups_rands(seeds, 5, true, true, 30, &mut rnd)
                .into_iter()
                .map(|r| r.sa_row)
                .collect::<Vec<_>>()
        };
        assert_eq!(collect(42), collect(42));
    }

    /// rank_seed_hits: every non-empty seed appears exactly once, in
    /// size-ascending order with PRNG-driven ties.
    #[test]
    fn rank_covers_all_nonempty() {
        let mut rnd = RandomSource::new(42);
        let input = SeedRankInput {
            n_offs: 5,
            fw_sizes: vec![Some(10), Some(0), Some(3), None, Some(7)],
            rc_sizes: vec![Some(2), Some(5), None, Some(15), Some(1)],
        };
        let ranks = rank_seed_hits(&input, &mut rnd);
        // fw non-empty (size>0): offsets 0, 2, 4 → 3.
        // rc non-empty (size>0): offsets 0, 1, 3, 4 → 4.
        // Total = 7.
        assert_eq!(ranks.len(), 7);
        // Verify uniqueness.
        let mut seen: Vec<(usize, bool)> = ranks.clone();
        seen.sort();
        seen.dedup();
        assert_eq!(seen.len(), 7);
    }

    /// Ranks are sorted by size ascending (modulo PRNG tie-breaking).
    #[test]
    fn rank_sorts_by_size() {
        let mut rnd = RandomSource::new(42);
        let input = SeedRankInput {
            n_offs: 4,
            fw_sizes: vec![Some(100), Some(2), Some(50), Some(20)],
            rc_sizes: vec![None, None, None, None],
        };
        let ranks = rank_seed_hits(&input, &mut rnd);
        // All forward; sizes 2, 20, 50, 100 → expected order: offset 1, 3, 2, 0.
        let sizes: Vec<u32> = ranks
            .iter()
            .map(|(i, _)| input.fw_sizes[*i].unwrap())
            .collect();
        assert_eq!(sizes, vec![2, 20, 50, 100]);
    }

    /// Determinism: same seed → same rank order.
    #[test]
    fn rank_deterministic() {
        let run = |seed: u32| {
            let mut rnd = RandomSource::new(seed);
            let input = SeedRankInput {
                n_offs: 5,
                fw_sizes: vec![Some(10), Some(10), Some(10), Some(5), Some(2)],
                rc_sizes: vec![Some(7), Some(7), Some(7), Some(3), Some(1)],
            };
            rank_seed_hits(&input, &mut rnd)
        };
        assert_eq!(run(42), run(42));
        assert_ne!(run(42), run(43)); // different seed → likely different order on ties
    }

    /// Exhaustion: with `maxelt > total elements`, all rows are emitted
    /// once.
    #[test]
    fn exhaustion() {
        let mut rnd = RandomSource::new(42);
        let seeds = vec![mk_seed(0, 50, 0)]; // 50 rows
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 1000, &mut rnd);
        assert_eq!(out.len(), 50);
        let mut seen: Vec<u32> = out.iter().map(|r| r.sa_row).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 50);
    }
}
