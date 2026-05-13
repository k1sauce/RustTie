//! Seed-and-extend alignment.
//!
//! Phase 2b scope: multi-seed exact FM-index search → candidate positions →
//! reference-window comparison → mismatches-only scoring (no indels yet).
//! BT2's default seed length is 22; we mirror that. Seed interval defaults
//! to `1 + 1.15 * sqrt(L)` rounded down (BT2's `-S S,1,1.15`).
//!
//! What's deferred to 2c:
//! - Allowing mismatches inside seeds (`-N 1`).
//! - Indels in the extension (banded SW via `block-aligner`).
//! - Quality-scaled mismatch penalties (currently constant 2 = BT2 Q40).

use rusttie_index::search::{backward_search, joined_to_ref, resolve_text_pos};
use rusttie_index::{BitPairReference, Bt2Index, exact_hits};

use crate::extend::{extend, try_ungapped};
use crate::revcomp::reverse_complement;

/// Process-global counters of how many candidates resolve via the ungapped
/// fast path vs the SW rescue. Bumped from `score_candidate_*`; printed by
/// the CLI on `RUSTTIE_PROFILE=1`. Relaxed-ordering atomics, ~1 ns/op,
/// kept in for perf debugging since the cost is below the noise floor.
#[doc(hidden)]
pub mod profile {
    use std::sync::atomic::{AtomicU64, Ordering};
    pub static N_UNGAPPED_OK: AtomicU64 = AtomicU64::new(0);
    pub static N_SW_FALLBACK: AtomicU64 = AtomicU64::new(0);
    pub static N_READS_CAP_FIRED: AtomicU64 = AtomicU64::new(0);
    pub static N_READS_TOTAL: AtomicU64 = AtomicU64::new(0);
    pub fn print() {
        eprintln!(
            "[profile] ungapped_ok={} sw_fallback={} reads_total={} reads_cap_fired={}",
            N_UNGAPPED_OK.load(Ordering::Relaxed),
            N_SW_FALLBACK.load(Ordering::Relaxed),
            N_READS_TOTAL.load(Ordering::Relaxed),
            N_READS_CAP_FIRED.load(Ordering::Relaxed),
        );
    }
}

/// Default maximum SA-range size for a seed before we treat it as too
/// repetitive and skip it entirely. Repetitive seeds (e.g., low-complexity
/// AT-rich 22-mers) can have thousands of hits in a 50M-bp reference;
/// extending every one is the dominant runtime cost. BT2 uses a similar
/// heuristic in its descent driver.
///
/// Cap value vs chr22 perf/recall (10k paired 100bp reads, 0.5% error):
///   ∞  → 95s, 100.0% recall (every candidate extended)
///   300→ 15.7s, 99.9% recall
///   100→ 6.3s,  99.6% recall
///   50 → 3.8s,  99.5% recall  ← default
///   30 → 2.7s,  99.4% recall
///
/// 50 is the knee — 25× faster than uncapped at 0.5% recall cost. Tunable
/// via `--seed-hit-cap` on the CLI.
pub const PER_SEED_HIT_CAP: u32 = 50;

/// BT2 `-D` default: consecutive seed-extension failures before giving up
/// on the current seed set. An "extension fails" when its alignment
/// neither improves best nor secbest.
pub const DESCENT_D_DEFAULT: u32 = 15;

/// BT2 `-R` default: maximum re-seedings if the current seed set is
/// "repetitive" (avg hits per aligned seed > 300, per BT2 manual).
pub const DESCENT_R_DEFAULT: u32 = 2;

/// BT2's repetitiveness threshold: when avg(total_hits / n_aligned_seeds)
/// exceeds this, the read is considered to have repetitive seeds and a
/// re-seeding pass is triggered.
pub const REPETITIVE_HITS_THRESHOLD: u32 = 300;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub ref_id: u32,
    /// 0-based leftmost reference position (SAM convention).
    pub ref_off: u32,
    pub strand: Strand,
    pub read_len: u32,
    /// Mismatches (CIGAR M with unequal bases). Excludes indels.
    pub mismatches: u32,
    /// Number of distinct gap runs (`XO:i`).
    pub gap_opens: u32,
    /// Number of gap-extending bases (`XG:i` total = ins + del bases).
    pub gap_extends: u32,
    /// `AS:i` value. Negative under BT2's end-to-end scoring.
    pub score: i32,
    /// CIGAR string built from the SW alignment ops.
    pub cigar: String,
    /// `MD:Z` tag value (without the `MD:Z:` prefix).
    pub md: String,
}

/// BT2's default seed length (`-L 22`).
pub const DEFAULT_SEED_LEN: u32 = 22;

/// BT2's default mismatch penalty for backwards-compat constant.
pub const MM_PENALTY_MAX: i32 = 6;
pub const MM_PENALTY_MIN: i32 = 2;
pub const MM_PENALTY_MAX_Q: u8 = 40;

/// Reference-window slack on each side of a seed-inferred candidate position.
/// Allows for indels that shift the alignment relative to the seed match.
pub const EXTEND_SLACK: u32 = 15;

/// Runtime-configurable scoring parameters (BT2 `--mp`, `--rdg`, `--rfg`,
/// `--score-min`, `--np`). Defaults match BT2's defaults exactly, so omitting
/// every flag produces the same SAM as before this struct existed.
#[derive(Debug, Clone, Copy)]
pub struct Scoring {
    /// Max mismatch penalty (Q≥40), `--mp` first arg.
    pub mp_max: i32,
    /// Min mismatch penalty (Q=0), `--mp` second arg.
    pub mp_min: i32,
    /// Quality at which mismatch penalty saturates at `mp_max`. Always 40.
    pub mp_max_q: u8,
    /// Read gap open penalty, `--rdg` first arg.
    pub rdg_open: i32,
    /// Read gap extend penalty, `--rdg` second arg.
    pub rdg_extend: i32,
    /// Reference gap open penalty, `--rfg` first arg.
    pub rfg_open: i32,
    /// Reference gap extend penalty, `--rfg` second arg.
    pub rfg_extend: i32,
    /// `--score-min L,A,B`: minimum acceptable score is `A + B * read_len`
    /// for `L` (linear) function. We only support `L` for now.
    pub score_min_const: f64,
    pub score_min_coeff: f64,
}

impl Default for Scoring {
    fn default() -> Self {
        // BT2 end-to-end defaults: --mp 6,2 --rdg 5,3 --rfg 5,3 --score-min L,-0.6,-0.6
        Self {
            mp_max: 6,
            mp_min: 2,
            mp_max_q: 40,
            rdg_open: 5,
            rdg_extend: 3,
            rfg_open: 5,
            rfg_extend: 3,
            score_min_const: -0.6,
            score_min_coeff: -0.6,
        }
    }
}

impl Scoring {
    /// Mismatch penalty for a base of Phred quality `q` (clamped to [0, mp_max_q]).
    #[inline]
    pub fn mm_penalty(&self, q: u8) -> i32 {
        let q = q.min(self.mp_max_q) as i32;
        self.mp_min + (self.mp_max - self.mp_min) * q / self.mp_max_q as i32
    }

    /// Minimum acceptable alignment score for a read of length `read_len`.
    pub fn score_min(&self, read_len: u32) -> i32 {
        // BT2 truncates toward zero (C-style `(T)ret` cast in
        // `simple_func.h:110`), not nearest-even round. `-60.6` → `-60`,
        // not `-61`. The off-by-one per mate compounds to off-by-two on
        // paired-end pair_smin and shifts the `bestdiff/diff` MAPQ ratio
        // into the next bin down, costing ~6% MAPQ agreement on chr22.
        (self.score_min_const + self.score_min_coeff * read_len as f64) as i32
    }
}

/// Convert ASCII Phred+33 quality byte to Phred Q.
#[inline]
pub fn phred33_to_q(b: u8) -> u8 {
    b.saturating_sub(33)
}

/// Backwards-compat shims for the now-deprecated free functions. New code
/// should use [`Scoring::mm_penalty`] / [`Scoring::score_min`] directly.
#[inline]
pub fn mm_penalty(q: u8) -> i32 {
    Scoring::default().mm_penalty(q)
}

pub fn score_min(read_len: u32) -> i32 {
    Scoring::default().score_min(read_len)
}

/// BT2's default seed interval `S,1,1.15` for a read of length `L`:
/// `1 + 1.15 * sqrt(L)`, floored. Always at least 1.
pub fn seed_interval(read_len: u32) -> u32 {
    let v = 1.0 + 1.15 * (read_len as f64).sqrt();
    (v.floor() as u32).max(1)
}

/// Seed offsets within a read of length `read_len` for the given seed length.
pub fn seed_offsets(read_len: u32, seed_len: u32) -> Vec<u32> {
    seed_offsets_shifted(read_len, seed_len, 0)
}

/// Seed offsets shifted by `shift` chars from the read's left end. Used by
/// the descent re-seeding logic (`-R`): each retry pass uses a different
/// shift so the read's true alignment site has a fresh chance to be hit by
/// a non-repetitive seed.
pub fn seed_offsets_shifted(read_len: u32, seed_len: u32, shift: u32) -> Vec<u32> {
    if read_len < seed_len {
        return Vec::new();
    }
    let interval = seed_interval(read_len);
    let last = read_len - seed_len;
    let start = shift % interval.max(1);
    let mut out = Vec::new();
    let mut o = start;
    while o <= last {
        out.push(o);
        o += interval;
    }
    // Always include the final offset so the right end of the read is covered.
    if out.last().copied() != Some(last) {
        out.push(last);
    }
    out
}

/// Per-read alignment result: the chosen `best` alignment, the next-best
/// score (used for *single-end* MAPQ via `mapq_v2`), and the full set of
/// valid alignments that passed `score_min`. The full set is what
/// paired-end MAPQ needs — BT2 computes paired MAPQ from the second-best
/// **concordant pair**, not from per-mate secbests, so we must Cartesian
/// the two mates' alignment sets to find it.
#[derive(Debug, Clone)]
pub struct AlignResult {
    pub best: Alignment,
    /// Score of the next-best **per-mate** alignment, or `None` if no other
    /// candidate passes `score_min`. May equal `best.score` on ties.
    pub secbest_score: Option<i32>,
    /// All valid alignments for this read (deduped by `(ref_id, ref_off,
    /// strand)`), sorted score-descending then leftmost-first. `all[0]` is
    /// always equal to `best`. Ordering is deterministic.
    pub all: Vec<Alignment>,
}

/// Try to align a read on either strand. Returns the best alignment by
/// score (ties broken by leftmost ref position) plus the second-best score
/// (for MAPQ), or `None` if no candidate passes `score_min`.
///
/// `qual` is the read's Phred+33 ASCII quality string; if empty, every base
/// is treated as Q40 (BT2 default for unknown qualities).
pub fn align_read(
    idx: &Bt2Index,
    refs: &BitPairReference,
    read: &[u8],
    qual: &[u8],
    scoring: &Scoring,
) -> Option<AlignResult> {
    align_read_with_descent(
        idx,
        refs,
        read,
        qual,
        scoring,
        PER_SEED_HIT_CAP,
        DESCENT_D_DEFAULT,
        DESCENT_R_DEFAULT,
    )
}

/// Backwards-compat: `--seed-hit-cap` plumbing path. Uses BT2 default `-D`/`-R`.
pub fn align_read_with_cap(
    idx: &Bt2Index,
    refs: &BitPairReference,
    read: &[u8],
    qual: &[u8],
    scoring: &Scoring,
    seed_hit_cap: u32,
) -> Option<AlignResult> {
    align_read_with_descent(
        idx,
        refs,
        read,
        qual,
        scoring,
        seed_hit_cap,
        DESCENT_D_DEFAULT,
        DESCENT_R_DEFAULT,
    )
}

/// BT2-faithful descent driver: prioritizes least-repetitive seeds, breaks
/// after `d_budget` consecutive non-improving extensions, and re-seeds
/// (shifted offsets) up to `r_reseed` times if the current seed set is
/// repetitive (avg hits per aligned seed > [`REPETITIVE_HITS_THRESHOLD`]).
/// Mirrors the BT2 manual's description of `-D` / `-R`.
#[allow(clippy::too_many_arguments)]
pub fn align_read_with_descent(
    idx: &Bt2Index,
    refs: &BitPairReference,
    read: &[u8],
    qual: &[u8],
    scoring: &Scoring,
    seed_hit_cap: u32,
    d_budget: u32,
    r_reseed: u32,
) -> Option<AlignResult> {
    let read_len = read.len() as u32;
    let smin = scoring.score_min(read_len);
    let dbg = std::env::var_os("RUSTTIE_DEBUG").is_some();
    if dbg {
        eprintln!("[align] read_len={read_len} smin={smin} D={d_budget} R={r_reseed}",);
    }

    let rc_query = reverse_complement(read);
    let rc_qual: Vec<u8> = qual.iter().rev().copied().collect();
    let strands: [(Strand, &[u8], &[u8]); 2] = [
        (Strand::Forward, read, qual),
        (Strand::Reverse, rc_query.as_slice(), rc_qual.as_slice()),
    ];

    // Accumulate every passing alignment so paired-end MAPQ can enumerate
    // concordant pairs (BT2 needs the second-best concordant *pair* score,
    // which we can't recover from per-mate best/secbest alone).
    // - `cand_seen_ungapped`: dedups Phase 1 candidate extensions across
    //   re-seed passes (don't try the same seed-inferred position twice).
    // - `aln_seen`: dedups final alignments by their resolved (ref_id,
    //   ref_off, strand). Phase 2 SW rescue uses *only* this set so it
    //   isn't blocked by Phase 1's per-candidate dedup — different
    //   bookkeeping levels.
    // `current_best` / `current_secbest` track BT2's `-D` "improves best or
    // secbest" rule on the fly.
    let mut all_alns: Vec<Alignment> = Vec::new();
    let mut current_best: i32 = i32::MIN;
    let mut current_secbest: i32 = i32::MIN;
    let mut cand_seen_ungapped = std::collections::HashSet::new();
    let mut aln_seen = std::collections::HashSet::new();
    let mut all_cands: Vec<PrioritizedCandidate> = Vec::new();
    let mut cap_fired = false;

    // Pass 0 uses the default offsets; passes 1..=R shift them so that a
    // fresh set of seeds covers the read. R+1 distinct shifts in [0, interval).
    let interval = seed_interval(read_len);
    let n_passes = r_reseed + 1;

    for pass in 0..n_passes {
        let shift = if r_reseed == 0 {
            0
        } else {
            pass * interval / n_passes
        };

        let mut pass_cands: Vec<PrioritizedCandidate> = Vec::new();
        let mut pass_total_hits: u64 = 0;
        let mut pass_aligned_seeds: u32 = 0;
        for (strand, query, _q) in &strands {
            collect_prioritized(
                idx,
                query,
                DEFAULT_SEED_LEN,
                shift,
                seed_hit_cap,
                *strand,
                &mut pass_cands,
                &mut pass_total_hits,
                &mut pass_aligned_seeds,
                &mut cap_fired,
            );
        }
        pass_cands.sort_by_key(|c| (c.sa_range_size, c.ref_id, c.ref_off));

        if dbg {
            eprintln!(
                "[align] pass={pass} shift={shift} cands={} total_hits={pass_total_hits} aligned_seeds={pass_aligned_seeds}",
                pass_cands.len(),
            );
        }

        let mut consecutive_failures: u32 = 0;
        for cand in &pass_cands {
            if consecutive_failures >= d_budget {
                break;
            }
            let cand_key = (cand.ref_id, cand.ref_off, cand.strand);
            if !cand_seen_ungapped.insert(cand_key) {
                continue;
            }
            let (_, query, q) = strands[strand_idx(cand.strand)];
            let aln_opt = score_candidate_ungapped(
                refs,
                cand.ref_id,
                cand.ref_off,
                query,
                q,
                cand.strand,
                scoring,
                smin,
            );
            match aln_opt {
                Some(aln) => {
                    let aln_key = (aln.ref_id, aln.ref_off, aln.strand);
                    let s = aln.score;
                    let improved = update_score_window(s, &mut current_best, &mut current_secbest);
                    if aln_seen.insert(aln_key) {
                        all_alns.push(aln);
                    }
                    if improved {
                        consecutive_failures = 0;
                    } else {
                        consecutive_failures += 1;
                    }
                }
                None => {
                    consecutive_failures += 1;
                }
            }
        }

        all_cands.extend(pass_cands);

        // BT2's re-seed gate: avg_hits/aligned_seeds > REPETITIVE_THRESHOLD.
        // If we already have a hit and seeds are normal, no need to re-seed.
        let avg_hits = if pass_aligned_seeds > 0 {
            pass_total_hits / pass_aligned_seeds as u64
        } else {
            0
        };
        let repetitive = avg_hits > REPETITIVE_HITS_THRESHOLD as u64;
        if !all_alns.is_empty() && !repetitive {
            break;
        }
    }

    // Phase 2 (SW rescue) only if Phase 1 found nothing. Iterate the same
    // candidates Phase 1 saw but use SW DP — handles indel reads that
    // ungapped Hamming can't score.
    if all_alns.is_empty() {
        for cand in &all_cands {
            let (_, query, q) = strands[strand_idx(cand.strand)];
            if let Some(aln) = score_candidate_gapped(
                refs,
                cand.ref_id,
                cand.ref_off,
                query,
                q,
                cand.strand,
                scoring,
            ) {
                if aln.score < smin {
                    continue;
                }
                let aln_key = (aln.ref_id, aln.ref_off, aln.strand);
                if aln_seen.insert(aln_key) {
                    all_alns.push(aln);
                }
            }
        }
    }

    profile::N_READS_TOTAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    if cap_fired {
        profile::N_READS_CAP_FIRED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }

    if all_alns.is_empty() {
        return None;
    }

    // Sort score-descending, leftmost-first so `all[0]` is the displayed
    // best and downstream paired-end secbest selection is deterministic.
    all_alns.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.ref_id.cmp(&b.ref_id))
            .then(a.ref_off.cmp(&b.ref_off))
    });
    let best = all_alns[0].clone();
    let secbest_score = all_alns.get(1).map(|a| a.score);
    Some(AlignResult {
        best,
        secbest_score,
        all: all_alns,
    })
}

/// Fold a new score `s` into running `(current_best, current_secbest)`.
/// Returns `true` iff `s` improves either bound — i.e. counts as a
/// non-failed extension under BT2's `-D` rule. A tie with `current_best`
/// counts as a new secbest (secbest goes from `<best` to `==best`).
#[inline]
fn update_score_window(s: i32, current_best: &mut i32, current_secbest: &mut i32) -> bool {
    if s > *current_best {
        if *current_best != i32::MIN {
            *current_secbest = *current_best;
        }
        *current_best = s;
        true
    } else if (s == *current_best && *current_secbest < *current_best)
        || (s < *current_best && s > *current_secbest)
    {
        *current_secbest = s;
        true
    } else {
        false
    }
}

#[inline]
pub(crate) fn strand_idx(s: Strand) -> usize {
    match s {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

/// Candidate enriched with the source seed's SA-range size and strand. The
/// descent driver sorts by `sa_range_size` ascending so least-repetitive
/// seeds extend first — same priority as BT2 uses.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PrioritizedCandidate {
    pub(crate) ref_id: u32,
    pub(crate) ref_off: u32,
    pub(crate) strand: Strand,
    pub(crate) sa_range_size: u32,
}

/// Collect prioritized candidates from `query` seeded at offsets shifted by
/// `shift` (for re-seeding passes). Pushes into `out`; counts total
/// SA-range hits and the number of seeds that produced ≥1 hit so the
/// caller can compute BT2's repetitiveness criterion. Sets `cap_fired`
/// true if any seed was skipped due to `seed_hit_cap` — caller uses this
/// to pessimize MAPQ (unverified candidates may include a real alternate).
#[allow(clippy::too_many_arguments)]
pub(crate) fn collect_prioritized(
    idx: &Bt2Index,
    query: &[u8],
    seed_len: u32,
    shift: u32,
    seed_hit_cap: u32,
    strand: Strand,
    out: &mut Vec<PrioritizedCandidate>,
    total_hits: &mut u64,
    aligned_seeds: &mut u32,
    cap_fired: &mut bool,
) {
    let read_len = query.len() as u32;
    let dbg = std::env::var_os("RUSTTIE_DEBUG").is_some();
    for off in seed_offsets_shifted(read_len, seed_len, shift) {
        let seed = &query[off as usize..(off + seed_len) as usize];
        let Some(range) = backward_search(idx, seed) else {
            continue;
        };
        if range.is_empty() {
            continue;
        }
        let hits = range.len();
        *total_hits += hits as u64;
        *aligned_seeds += 1;
        if hits > seed_hit_cap {
            *cap_fired = true;
            if dbg {
                eprintln!(
                    "[seed] strand={strand:?} off={off} skipped — too repetitive ({hits} > cap {seed_hit_cap})",
                );
            }
            continue;
        }
        for row in range.lo..range.hi {
            let pos = resolve_text_pos(idx, row);
            let hit = joined_to_ref(idx, pos);
            if hit.ref_off < off {
                continue;
            }
            out.push(PrioritizedCandidate {
                ref_id: hit.ref_id,
                ref_off: hit.ref_off - off,
                strand,
                sa_range_size: hits,
            });
        }
    }
}

/// Backwards-compat: legacy export of `exact_hits` is still used by the
/// `rusttie_aligner_handles_all_stretches_in_multi_contig` test path.
#[allow(dead_code)]
fn _suppress_unused_exact_hits_warning(idx: &Bt2Index, seed: &[u8]) {
    let _ = exact_hits(idx, seed);
}

/// Verify a candidate by extracting a slack reference window around the
/// seed-inferred position and running semiglobal SW. Handles indels.
///
/// The window is clamped to the unambiguous stretch containing
/// `seed_ref_off` so we never read across N gaps. If the read can't fit
/// within the containing stretch starting at `seed_ref_off`, the candidate
/// is rejected (its alignment would span Ns, which BT2 also rejects).
/// Stretch + window math common to both score paths. Returns the
/// (window, ungapped_lo, win_start) bundle, or None if the candidate spans
/// an N gap or falls outside the index.
fn locate_window(
    refs: &BitPairReference,
    ref_id: u32,
    seed_ref_off: u32,
    read_len: u32,
) -> Option<(Vec<u8>, u32, u32)> {
    let (_joined, rec_idx) = refs.locate(ref_id, seed_ref_off)?;
    let stretch_start = refs.ref_offsets[rec_idx];
    let stretch_end = stretch_start + refs.records[rec_idx].len;
    if seed_ref_off + read_len > stretch_end {
        return None;
    }
    let actual_left = (seed_ref_off - stretch_start).min(EXTEND_SLACK);
    let actual_right = (stretch_end - seed_ref_off - read_len).min(EXTEND_SLACK);
    let win_start = seed_ref_off - actual_left;
    let win_len = read_len + actual_left + actual_right;
    let window = refs.extract(ref_id, win_start, win_len)?;
    Some((window, actual_left, win_start))
}

/// Ungapped (Hamming) score of a single candidate. Returns `Some` only if
/// the alignment passes smin — partial / spurious candidates return None
/// without paying SW cost.
#[allow(clippy::too_many_arguments)]
pub(crate) fn score_candidate_ungapped(
    refs: &BitPairReference,
    ref_id: u32,
    seed_ref_off: u32,
    query: &[u8],
    qual: &[u8],
    strand: Strand,
    scoring: &Scoring,
    smin: i32,
) -> Option<Alignment> {
    let read_len = query.len() as u32;
    let (window, actual_left, _win_start) = locate_window(refs, ref_id, seed_ref_off, read_len)?;
    let lo = actual_left as usize;
    let hi = lo + read_len as usize;
    let ext = try_ungapped(query, qual, &window[lo..hi], seed_ref_off, scoring, smin)?;
    profile::N_UNGAPPED_OK.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    Some(Alignment {
        ref_id,
        ref_off: ext.ref_off,
        strand,
        read_len,
        mismatches: ext.mismatches,
        gap_opens: 0,
        gap_extends: 0,
        score: ext.score,
        cigar: ext.cigar,
        md: ext.md,
    })
}

/// SW DP rescue for candidates that ungapped didn't resolve. Only invoked
/// when Phase 1 finds zero alignments for the whole read, so a false
/// positive here at least doesn't waste SW work behind a known-good
/// ungapped hit.
fn score_candidate_gapped(
    refs: &BitPairReference,
    ref_id: u32,
    seed_ref_off: u32,
    query: &[u8],
    qual: &[u8],
    strand: Strand,
    scoring: &Scoring,
) -> Option<Alignment> {
    let read_len = query.len() as u32;
    let (window, _actual_left, win_start) = locate_window(refs, ref_id, seed_ref_off, read_len)?;
    profile::N_SW_FALLBACK.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let ext = extend(query, qual, &window, win_start, scoring)?;
    Some(Alignment {
        ref_id,
        ref_off: ext.ref_off,
        strand,
        read_len,
        mismatches: ext.mismatches,
        gap_opens: ext.gap_opens,
        gap_extends: ext.n_ins + ext.n_del,
        score: ext.score,
        cigar: ext.cigar,
        md: ext.md,
    })
}

/// FR mate-rescue: given an anchor mate alignment, search the
/// concordance-allowed window for the OTHER mate via SW (with an ungapped
/// fast path). Mirrors BT2's `extendSeedsPaired` inner mate-find step
/// (`vendor/bowtie2/aligner_sw_driver.cpp:2226-2311` + `pe.cpp:161-354`).
///
/// `other_fwd` / `other_qual` is the other mate's FASTQ-orientation read.
/// We RC internally based on the anchor's strand (FR pair: opposite strand).
/// Returns one alignment if it scores ≥ smin, else `None`.
///
/// Window math (FR, see `pe.cpp:206-350`):
/// - Anchor on forward strand → other on reverse strand to the RIGHT.
///   Other's leftmost ref position ∈ `[anchor.ref_off + frag_min - olen,
///   anchor.ref_off + frag_max - olen]`. SW reference span =
///   `[anchor.ref_off + frag_min - olen, anchor.ref_off + frag_max - 1]`.
/// - Anchor on reverse strand → other on forward strand to the LEFT.
///   Other's leftmost ref position ∈ `[anchor.ref_off + alen - frag_max,
///   anchor.ref_off + alen - frag_min]`. SW reference span =
///   `[anchor.ref_off + alen - frag_max, anchor.ref_off + alen - frag_min
///   + olen - 1]`.
#[allow(clippy::too_many_arguments)]
pub fn mate_rescue(
    refs: &BitPairReference,
    anchor: &Alignment,
    other_fwd: &[u8],
    other_qual: &[u8],
    other_rc: &[u8],
    other_qual_rev: &[u8],
    scoring: &Scoring,
    frag_min: u32,
    frag_max: u32,
) -> Option<Alignment> {
    let alen = anchor.read_len as i64;
    let olen = other_fwd.len() as u64;
    let olen_u32 = olen as u32;
    let smin = scoring.score_min(olen_u32);

    // Compute the [win_lo, win_hi) reference window and the strand the
    // OTHER mate sits on (and the corresponding query orientation).
    let (other_strand, query, q, win_lo_signed, win_hi_signed) = match anchor.strand {
        Strand::Forward => {
            let lo = anchor.ref_off as i64 + frag_min as i64 - olen as i64;
            let hi = anchor.ref_off as i64 + frag_max as i64; // exclusive
            (Strand::Reverse, other_rc, other_qual_rev, lo, hi)
        }
        Strand::Reverse => {
            let lo = anchor.ref_off as i64 + alen - frag_max as i64;
            let hi = anchor.ref_off as i64 + alen - frag_min as i64 + olen as i64; // exclusive
            (Strand::Forward, other_fwd, other_qual, lo, hi)
        }
    };

    let win_lo = win_lo_signed.max(0) as u32;
    if win_hi_signed <= win_lo as i64 {
        return None;
    }
    let win_hi = win_hi_signed as u32;
    if win_hi.saturating_sub(win_lo) < olen_u32 {
        return None;
    }

    // Clamp to the unambiguous stretch containing win_lo (avoid spanning Ns,
    // mirrors `score_candidate_*`).
    let (_joined, rec_idx) = refs.locate(anchor.ref_id, win_lo)?;
    let stretch_start = refs.ref_offsets[rec_idx];
    let stretch_end = stretch_start + refs.records[rec_idx].len;
    let win_start = win_lo.max(stretch_start);
    let win_end = win_hi.min(stretch_end);
    if win_end <= win_start || win_end - win_start < olen_u32 {
        return None;
    }
    let win_len = win_end - win_start;
    let window = refs.extract(anchor.ref_id, win_start, win_len)?;

    // Fast path: ungapped sliding-Hamming over the window. Pick the best
    // start-position by score. Most chr22 reads finish here.
    let mut best_off: u32 = 0;
    let mut best_score: i32 = i32::MIN;
    let mut best_md: String = String::new();
    let mut best_mm: u32 = 0;
    let last_offset = win_len - olen_u32;
    for off in 0..=last_offset {
        let lo = off as usize;
        let hi = lo + olen_u32 as usize;
        // Inline Hamming with early bailout if score drops too low to beat
        // the running best (helps when window is wide and most positions
        // are far from the read).
        let mut score: i32 = 0;
        let mut mismatches: u32 = 0;
        let mut md = String::new();
        let mut md_run: u32 = 0;
        let bailout_threshold = best_score; // any score <= best_score is no improvement
        let slice = &window[lo..hi];
        let mut bailed = false;
        for i in 0..olen as usize {
            let r = query[i];
            let g = slice[i];
            let g_is_acgt = matches!(g, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
            if g_is_acgt && r.eq_ignore_ascii_case(&g) {
                md_run += 1;
            } else {
                use std::fmt::Write as _;
                let _ = write!(&mut md, "{md_run}");
                md_run = 0;
                md.push(g as char);
                mismatches += 1;
                let qv = q
                    .get(i)
                    .copied()
                    .map(phred33_to_q)
                    .unwrap_or(scoring.mp_max_q);
                score -= scoring.mm_penalty(qv);
                if score <= bailout_threshold || score < smin {
                    bailed = true;
                    break;
                }
            }
        }
        if bailed {
            continue;
        }
        use std::fmt::Write as _;
        let _ = write!(&mut md, "{md_run}");
        if score > best_score {
            best_score = score;
            best_off = off;
            best_md = md;
            best_mm = mismatches;
        }
    }

    if best_score >= smin {
        let mut cigar = String::new();
        use std::fmt::Write as _;
        let _ = write!(&mut cigar, "{}M", olen_u32);
        return Some(Alignment {
            ref_id: anchor.ref_id,
            ref_off: win_start + best_off,
            strand: other_strand,
            read_len: olen_u32,
            mismatches: best_mm,
            gap_opens: 0,
            gap_extends: 0,
            score: best_score,
            cigar,
            md: best_md,
        });
    }

    // Gapped fallback: full SW over the window (handles indels). This is
    // the heavy path; only fires when ungapped didn't clear smin.
    profile::N_SW_FALLBACK.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let ext = extend(query, q, &window, win_start, scoring)?;
    if ext.score < smin {
        return None;
    }
    Some(Alignment {
        ref_id: anchor.ref_id,
        ref_off: ext.ref_off,
        strand: other_strand,
        read_len: olen_u32,
        mismatches: ext.mismatches,
        gap_opens: ext.gap_opens,
        gap_extends: ext.n_ins + ext.n_del,
        score: ext.score,
        cigar: ext.cigar,
        md: ext.md,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_base() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
    }

    fn read_lambda() -> Vec<u8> {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../validate/fixtures/lambda_virus.fa");
        let text = std::fs::read_to_string(&path).unwrap();
        let mut s = Vec::new();
        for l in text.lines() {
            if !l.starts_with('>') {
                s.extend(l.trim().as_bytes());
            }
        }
        s
    }

    #[test]
    fn seed_offsets_50bp() {
        // L=50: interval = floor(1 + 1.15 * sqrt(50)) = floor(1 + 8.13) = 9.
        // Last valid start = 50-22 = 28. So 0, 9, 18, 27, plus 28.
        assert_eq!(seed_offsets(50, 22), vec![0, 9, 18, 27, 28]);
    }

    #[test]
    fn perfect_read_aligns() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let refs = BitPairReference::open(fixture_base()).unwrap();
        let seq = read_lambda();
        let read = &seq[10_000..10_050];
        let q = vec![b'I'; 50]; // Q40 throughout
        let a = align_read(&idx, &refs, read, &q, &Scoring::default())
            .expect("aligned")
            .best;
        assert_eq!(a.ref_off, 10_000);
        assert_eq!(a.strand, Strand::Forward);
        assert_eq!(a.mismatches, 0);
        assert_eq!(a.score, 0);
        assert_eq!(a.md, "50");
    }

    #[test]
    fn one_mismatch_high_q_aligns() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let refs = BitPairReference::open(fixture_base()).unwrap();
        let seq = read_lambda();
        let mut read = seq[20_000..20_050].to_vec();
        let orig = read[30];
        read[30] = if orig == b'A' { b'C' } else { b'A' };
        let q = vec![b'I'; 50]; // Q40
        let a = align_read(&idx, &refs, &read, &q, &Scoring::default())
            .expect("aligned")
            .best;
        assert_eq!(a.mismatches, 1);
        assert_eq!(a.score, -6); // Q40 → max penalty
        assert_eq!(a.md, format!("30{}19", orig as char));
    }

    #[test]
    fn one_mismatch_low_q_uses_low_penalty() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let refs = BitPairReference::open(fixture_base()).unwrap();
        let seq = read_lambda();
        let mut read = seq[20_000..20_050].to_vec();
        let orig = read[30];
        read[30] = if orig == b'A' { b'C' } else { b'A' };
        let mut q = vec![b'I'; 50];
        q[30] = b'!'; // Phred 0
        let a = align_read(&idx, &refs, &read, &q, &Scoring::default())
            .expect("aligned")
            .best;
        assert_eq!(a.mismatches, 1);
        assert_eq!(a.score, -2); // Q0 → min penalty
    }

    #[test]
    fn rc_read_aligns() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let refs = BitPairReference::open(fixture_base()).unwrap();
        let seq = read_lambda();
        let original = &seq[30_000..30_050];
        let read = reverse_complement(original);
        let q = vec![b'I'; 50];
        let a = align_read(&idx, &refs, &read, &q, &Scoring::default())
            .expect("aligned")
            .best;
        assert_eq!(a.ref_off, 30_000);
        assert_eq!(a.strand, Strand::Reverse);
        assert_eq!(a.mismatches, 0);
    }

    #[test]
    fn mm_penalty_endpoints_and_middle() {
        assert_eq!(mm_penalty(0), 2);
        assert_eq!(mm_penalty(40), 6);
        assert_eq!(mm_penalty(20), 4);
        // Saturation above 40.
        assert_eq!(mm_penalty(60), 6);
    }
}
