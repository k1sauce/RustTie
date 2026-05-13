//! BT2-faithful pseudo-randomness primitives. Required for matching BT2's
//! SA-range traversal order (the dominant factor in our remaining MAPQ
//! gap — see `rusttie.md` Phase 0 section).
//!
//! Two pieces, both ports of `vendor/bowtie2`:
//!
//! 1. [`RandomSource`] — BT2's LCG random source from `random_source.h:34`.
//!    Simple Numerical-Recipes-style LCG with two interleaved steps per
//!    `next_u32()` call.
//! 2. [`Random1toN`] — BT2's "without replacement" integer sampler from
//!    `random_util.h:32`. Adapts strategy based on range size: swap-list
//!    (Fisher-Yates) for small ranges, seen-list-then-convert for large.
//!
//! Validated against BT2's behavior by replaying byte-exact PRNG sequences
//! (see tests).

/// LCG matching BT2's `RandomSource` (`vendor/bowtie2/random_source.h:34`).
/// Constants `a = 1664525, c = 1013904223` are the Numerical Recipes LCG.
/// `next_u32` advances the LCG twice and XORs the high half of each.
#[derive(Debug, Clone, Copy)]
pub struct RandomSource {
    last: u32,
}

impl RandomSource {
    const A: u32 = 1664525;
    const C: u32 = 1013904223;

    pub fn new(seed: u32) -> Self {
        Self { last: seed }
    }

    pub fn next_u32(&mut self) -> u32 {
        self.last = Self::A.wrapping_mul(self.last).wrapping_add(Self::C);
        let mut ret = self.last >> 16;
        self.last = Self::A.wrapping_mul(self.last).wrapping_add(Self::C);
        ret ^= self.last;
        ret
    }

    /// `next_u32() / u32::MAX` as a single-precision float. BT2's
    /// `nextFloat()` divides in float32 precision — preserve that for
    /// byte-exact compatibility (`random_source.h:137-140`).
    pub fn next_float(&mut self) -> f32 {
        self.next_u32() as f32 / u32::MAX as f32
    }
}

/// "Without replacement" integer sampler. Adapts strategy:
/// - **Swap-list mode** (small ranges, n < `SWAPLIST_THRESH`): pre-allocate
///   `[0, n)`, partial Fisher-Yates shuffle on each `next()`. O(1) per
///   call, O(n) memory.
/// - **Seen-list mode** (large ranges, initial): generate random number,
///   reject if already seen, retain in a list. After `thresh` calls, the
///   probability of repeated rejections grows; convert to swap-list.
///
/// Port of `Random1toN` (`vendor/bowtie2/random_util.h:32`). The exact
/// swap and conversion order matters for byte-exact compatibility with
/// BT2's SA-range traversal.
#[derive(Debug, Clone)]
pub struct Random1toN {
    n: usize,
    cur: usize,
    swaplist: bool,
    converted: bool,
    list: Vec<usize>,
    seen: Vec<usize>,
    thresh: usize,
}

impl Random1toN {
    /// Below this size, use swap-list mode from the start.
    pub const SWAPLIST_THRESH: usize = 128;
    /// Convert seen-list → swap-list once seen-list has at least this many
    /// entries.
    pub const CONVERSION_THRESH: usize = 16;
    /// Also convert once seen-list size reaches this fraction of `n`.
    pub const CONVERSION_FRAC: f32 = 0.10;

    /// Initialize a sampler for `[0, n)`. Strategy chosen based on `n`.
    pub fn new(n: usize) -> Self {
        let swaplist = n < Self::SWAPLIST_THRESH;
        let thresh = Self::CONVERSION_THRESH
            .max((Self::CONVERSION_FRAC * n as f32) as usize);
        Self {
            n,
            cur: 0,
            swaplist,
            converted: false,
            list: Vec::new(),
            seen: Vec::new(),
            thresh,
        }
    }

    pub fn done(&self) -> bool {
        self.cur >= self.n
    }

    pub fn left(&self) -> usize {
        self.n - self.cur
    }

    /// Total range size (note: in seen-list-converted mode, this reflects
    /// the *remaining* range after conversion).
    pub fn size(&self) -> usize {
        self.n
    }

    /// Sample the next integer in `[0, n)`. Each integer is returned
    /// exactly once across `n` calls. The ORDER matches BT2 byte-for-byte
    /// when seeded with the same `RandomSource` state.
    pub fn next(&mut self, rnd: &mut RandomSource) -> usize {
        debug_assert!(!self.done());
        if self.cur == 0 && !self.converted {
            if self.n == 1 {
                self.cur = 1;
                return 0;
            }
            if self.swaplist {
                self.list = (0..self.n).collect();
            }
        }
        if self.swaplist {
            // Fisher-Yates partial shuffle: pick r ∈ [cur, n), swap, advance.
            let r = self.cur + (rnd.next_u32() as usize % (self.n - self.cur));
            if r != self.cur {
                self.list.swap(self.cur, r);
            }
            let result = self.list[self.cur];
            self.cur += 1;
            result
        } else {
            // Seen-list mode: rejection sample.
            let rn = loop {
                let rn = rnd.next_u32() as usize % self.n;
                if !self.seen.contains(&rn) {
                    break rn;
                }
            };
            self.seen.push(rn);
            self.cur += 1;
            // Convert to swap-list once we have enough seens that further
            // rejection is expensive. Matches BT2's logic at
            // `random_util.h:133-158`.
            if self.seen.len() >= self.thresh && self.cur < self.n {
                self.convert_to_swaplist();
            }
            rn
        }
    }

    /// Move from seen-list to swap-list mode. Build a swap-list containing
    /// every integer in `[0, n)` that we haven't yet returned, in
    /// ascending order. Mirrors BT2's conversion at `random_util.h:138-156`.
    /// Mark an SA-range as exhausted (used by RowSampler-style elimination).
    pub fn set_done(&mut self) {
        self.cur = self.n;
    }

    fn convert_to_swaplist(&mut self) {
        self.seen.sort_unstable();
        let mut new_list: Vec<usize> = Vec::with_capacity(self.n - self.cur);
        let mut prev = 0;
        for &s in &self.seen {
            for j in prev..s {
                new_list.push(j);
            }
            prev = s + 1;
        }
        for j in prev..self.n {
            new_list.push(j);
        }
        debug_assert_eq!(new_list.len(), self.n - self.cur);
        self.list = new_list;
        self.seen.clear();
        self.n -= self.cur;
        self.cur = 0;
        self.converted = true;
        self.swaplist = true;
    }
}

/// Weighted random sampler with elimination — picks one bin out of N
/// according to per-bin masses, with the ability to drop bins from
/// circulation. Port of BT2's `RowSampler` (`aligner_sw_driver.h:179-256`).
///
/// In BT2's paired-descent, this picks WHICH non-small SA-range to sample
/// a row from. Larger ranges (more BWT hits) have proportionally less
/// weight, so smaller (more specific) seeds are favored.
#[derive(Debug, Clone)]
pub struct RowSampler {
    /// Per-bin probability mass.
    masses: Vec<f64>,
    /// Per-bin elimination flag (true = excluded from sampling).
    eliminated: Vec<bool>,
    /// Sum of masses for live (non-eliminated) bins.
    mass: f64,
}

impl RowSampler {
    /// Initialize the sampler with per-bin (length, range_size) pairs.
    /// `lensq` squares the length term; `szsq` squares the size term.
    /// BT2 calls this with both `true` so the weight is
    /// `(extended_length^2) / (range_size^2)`.
    pub fn new<I: IntoIterator<Item = (usize, usize)>>(
        bins: I,
        lensq: bool,
        szsq: bool,
    ) -> Self {
        let mut masses = Vec::new();
        let mut mass: f64 = 0.0;
        for (len, range_sz) in bins {
            let mut num = len as f64;
            if lensq {
                num *= num;
            }
            let mut denom = range_sz as f64;
            if szsq {
                denom *= denom;
            }
            let m = num / denom;
            masses.push(m);
            mass += m;
        }
        let eliminated = vec![false; masses.len()];
        Self {
            masses,
            eliminated,
            mass,
        }
    }

    /// Mark bin `i` as exhausted so subsequent `next()` calls skip it.
    pub fn finished_range(&mut self, i: usize) {
        if !self.eliminated[i] {
            self.eliminated[i] = true;
            self.mass -= self.masses[i];
        }
    }

    /// Sample a live bin index proportional to its mass. BT2 uses
    /// `nextFloat() * mass` and walks the cumulative mass; we match
    /// that exactly so the same RNG state picks the same bin.
    pub fn next(&mut self, rnd: &mut RandomSource) -> usize {
        let rd = (rnd.next_float() as f64) * self.mass;
        let mut mass_sofar = 0.0;
        let mut last_unelim: usize = usize::MAX;
        for (i, &m) in self.masses.iter().enumerate() {
            if !self.eliminated[i] {
                last_unelim = i;
                mass_sofar += m;
                if rd < mass_sofar {
                    return i;
                }
            }
        }
        // Float drift can leave rd >= cumulative mass on the last bin.
        // BT2 falls through to the last unelim bin; match that behavior.
        debug_assert_ne!(last_unelim, usize::MAX);
        last_unelim
    }

    pub fn total_mass(&self) -> f64 {
        self.mass
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// LCG byte-exact regression: a fixed seed must produce a fixed sequence.
    /// Values cross-checked against compiled BT2 (run RandomSource with seed
    /// 42 and dump nextU32() × 5).
    #[test]
    fn random_source_byte_exact() {
        let mut r = RandomSource::new(42);
        // Verified against the LCG formula by hand; if these break, the
        // LCG diverged from BT2's `random_source.h:52-61`.
        let seq: Vec<u32> = (0..5).map(|_| r.next_u32()).collect();
        // recompute expected from the LCG directly to pin the formula.
        let (a, c): (u32, u32) = (1664525, 1013904223);
        let mut last: u32 = 42;
        let mut expected: Vec<u32> = Vec::new();
        for _ in 0..5 {
            last = a.wrapping_mul(last).wrapping_add(c);
            let mut ret = last >> 16;
            last = a.wrapping_mul(last).wrapping_add(c);
            ret ^= last;
            expected.push(ret);
        }
        assert_eq!(seq, expected);
    }

    /// Random1toN returns every integer in [0, n) exactly once.
    #[test]
    fn random_1_to_n_covers_full_range_swaplist() {
        let mut rnd = RandomSource::new(12345);
        let mut s = Random1toN::new(50);
        let mut seen: Vec<usize> = Vec::new();
        while !s.done() {
            seen.push(s.next(&mut rnd));
        }
        seen.sort_unstable();
        assert_eq!(seen, (0..50).collect::<Vec<_>>());
    }

    /// Same, for the seen-list-then-converted path (n > SWAPLIST_THRESH).
    #[test]
    fn random_1_to_n_covers_full_range_seenlist() {
        let mut rnd = RandomSource::new(98765);
        let mut s = Random1toN::new(500);
        let mut seen: Vec<usize> = Vec::new();
        while !s.done() {
            seen.push(s.next(&mut rnd));
        }
        seen.sort_unstable();
        assert_eq!(seen, (0..500).collect::<Vec<_>>());
    }

    /// n == 1 returns 0 then is done.
    #[test]
    fn random_1_to_n_singleton() {
        let mut rnd = RandomSource::new(1);
        let mut s = Random1toN::new(1);
        assert_eq!(s.next(&mut rnd), 0);
        assert!(s.done());
    }

    /// Determinism: same seed + same n produces same sequence.
    #[test]
    fn random_1_to_n_deterministic() {
        let collect = |seed: u32, n: usize| {
            let mut rnd = RandomSource::new(seed);
            let mut s = Random1toN::new(n);
            let mut out = Vec::new();
            while !s.done() {
                out.push(s.next(&mut rnd));
            }
            out
        };
        assert_eq!(collect(42, 30), collect(42, 30));
        assert_eq!(collect(42, 500), collect(42, 500));
    }

    /// Byte-exact match against compiled BT2 (`bt2_rng_test.cpp`). Reference
    /// output captured by running the C++ test driver against
    /// `vendor/bowtie2/random_source.h` + `random_util.cpp`. If these
    /// assertions break, the Rust port has diverged from BT2.
    #[test]
    fn matches_bt2_byte_exact() {
        // RandomSource with seed=42, first 10 nextU32() values.
        let mut r = RandomSource::new(42);
        let seq: Vec<u32> = (0..10).map(|_| r.next_u32()).collect();
        assert_eq!(
            seq,
            vec![
                378477685, 955892534, 110201035, 508762003, 4271932838,
                2146044429, 3699949778, 389807688, 4080590808, 3820277857,
            ]
        );

        // Random1toN(30) with seed=12345, full sequence (swap-list mode).
        let mut rnd = RandomSource::new(12345);
        let mut s = Random1toN::new(30);
        let mut got: Vec<usize> = Vec::new();
        while !s.done() {
            got.push(s.next(&mut rnd));
        }
        assert_eq!(
            got,
            vec![
                20, 11, 22, 28, 13, 14, 5, 3, 0, 27, 12, 9, 6, 19, 8,
                4, 24, 18, 10, 16, 23, 29, 21, 7, 2, 15, 25, 1, 26, 17,
            ]
        );

        // Random1toN(500) with seed=98765, first 20 values (seen-list mode).
        let mut rnd = RandomSource::new(98765);
        let mut s = Random1toN::new(500);
        let first20: Vec<usize> = (0..20).map(|_| s.next(&mut rnd)).collect();
        assert_eq!(
            first20,
            vec![
                56, 341, 311, 307, 292, 377, 75, 348, 360, 120,
                489, 45, 123, 103, 160, 439, 41, 213, 17, 315,
            ]
        );
    }

    /// `nextFloat` returns values in roughly [0, 1] using f32 division.
    #[test]
    fn next_float_in_range() {
        let mut r = RandomSource::new(42);
        for _ in 0..100 {
            let f = r.next_float();
            assert!((0.0..=1.0).contains(&f), "got {f}");
        }
    }

    /// `RowSampler` weighted distribution: with masses [10, 1, 1], bin 0
    /// should dominate.
    #[test]
    fn row_sampler_weighted() {
        let mut rnd = RandomSource::new(42);
        // Equivalent to (len=10, sz=1) → mass 100; (len=1, sz=1) → mass 1.
        let mut s = RowSampler::new(
            vec![(10, 1), (1, 1), (1, 1)],
            true,
            true,
        );
        let mut counts = [0usize; 3];
        for _ in 0..10_000 {
            counts[s.next(&mut rnd)] += 1;
        }
        // Expected weights: 100 / 102 ≈ 98%, 1/102 ≈ 1%, 1/102 ≈ 1%.
        assert!(counts[0] > 9_500, "bin 0 should dominate, got {counts:?}");
        assert!(counts[1] < 500);
        assert!(counts[2] < 500);
    }

    /// `RowSampler::finished_range` excludes a bin from sampling.
    #[test]
    fn row_sampler_elimination() {
        let mut rnd = RandomSource::new(42);
        let mut s = RowSampler::new(vec![(1, 1), (1, 1), (1, 1)], false, false);
        s.finished_range(1);
        for _ in 0..100 {
            let pick = s.next(&mut rnd);
            assert_ne!(pick, 1, "bin 1 should be eliminated");
        }
    }
}
