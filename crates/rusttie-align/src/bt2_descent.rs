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

    // BT2 sorts SA ranges by size ascending. Stable sort so per-mate
    // input order is preserved among same-size ranges (matches BT2's
    // `EList::sort` which is `std::stable_sort`).
    let mut sorted: Vec<SeedHit> = seeds;
    sorted.sort_by_key(|s| s.size());

    let nsmall = sorted.iter().take_while(|s| s.size() <= nsm).count();

    // Phase 1 — smalls: emit all rows from each small range in order.
    for seed in sorted.iter().take(nsmall) {
        if out.len() >= maxelt {
            return out;
        }
        let remaining = maxelt - out.len();
        let sz = seed.size() as usize;
        let take = sz.min(remaining);
        for row in seed.sa_lo..(seed.sa_lo + take as u32) {
            out.push(PrioritizedRow {
                sa_row: row,
                rdoff: seed.rdoff,
                seedlen: seed.seedlen,
                fw: seed.fw,
                sa_range_size: seed.size(),
            });
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

    /// All-small-range case: every row gets emitted in sort-by-size order.
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
        // Should be sorted by SA range size: size 1 (300), then 3 (100s),
        // then 5 (200s).
        assert_eq!(out[0].sa_row, 300);
        assert_eq!(out[1].sa_row, 100);
        assert_eq!(out[2].sa_row, 101);
        assert_eq!(out[3].sa_row, 102);
        // 200s range — 5 rows.
        assert_eq!(out[4].sa_row, 200);
        assert_eq!(out[8].sa_row, 204);
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

    /// Mixed small + large: smalls emitted exhaustively first, then large
    /// is sampled.
    #[test]
    fn mixed_smalls_then_large() {
        let mut rnd = RandomSource::new(42);
        let seeds = vec![
            mk_seed(0, 100, 0),   // 100 rows, LARGE
            mk_seed(1000, 1003, 22), // 3 rows, SMALL
        ];
        let out = prioritize_sa_tups_rands(seeds, 5, true, true, 50, &mut rnd);
        assert_eq!(out.len(), 50);
        // First 3 should be the small range, in order.
        assert_eq!(out[0].sa_row, 1000);
        assert_eq!(out[1].sa_row, 1001);
        assert_eq!(out[2].sa_row, 1002);
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
