//! Banded extension via `bio::alignment::pairwise` (semiglobal SW), plus
//! BT2-faithful AS recomputation, CIGAR, and MD tag construction.
//!
//! Strategy: bio's aligner doesn't support per-position quality, so we run
//! it with a fixed mismatch penalty (Q40 = -6) to find the alignment
//! structure, then walk the ops to recompute AS with Q-scaling and to build
//! CIGAR / MD / NM / XM / XO / XG.
//!
//! [`try_ungapped`] is a fast path that skips SW entirely when the read
//! aligns to the seed-inferred position with no indels — that's >99% of
//! reads on Illumina-like data. SW is reserved for the indel rescue path.

use bio::alignment::AlignmentOperation as Op;
use bio::alignment::pairwise::{Aligner, MIN_SCORE};
use std::fmt::Write;

use crate::align::{Scoring, phred33_to_q};

/// Result of extending a candidate to a full alignment.
pub struct Extended {
    /// 0-based reference position where the alignment starts.
    pub ref_off: u32,
    /// Mismatches across the alignment (excludes indel ops).
    pub mismatches: u32,
    /// Insertions in the read relative to reference (CIGAR `I` count).
    pub n_ins: u32,
    /// Deletions in the read relative to reference (CIGAR `D` count).
    pub n_del: u32,
    /// Number of distinct gap runs (CIGAR `XO` value).
    pub gap_opens: u32,
    /// Q-scaled alignment score (`AS:i`).
    pub score: i32,
    /// CIGAR string (`50M`, `25M2I23M`, etc.).
    pub cigar: String,
    /// MD tag value (without the `MD:Z:` prefix).
    pub md: String,
}

/// Ungapped (Hamming) extension of `read` against `ref_window` of equal
/// length, starting at `ref_off`. Returns the alignment if its Q-scaled
/// score is `>= smin`, else `None` (caller should fall back to SW DP).
/// Bails early as soon as `score < smin` so a hopeless alignment costs
/// only a few comparisons.
///
/// This is the dominant path on real data: with seed-and-extend, the
/// seed-inferred ref position is exactly correct unless an indel shifts
/// the read between the seed and the read's start — uncommon outside of
/// long-read or pathological data. For wgsim with `-R 0` it's 100%.
pub fn try_ungapped(
    read: &[u8],
    qual: &[u8],
    ref_window: &[u8],
    ref_off: u32,
    sc: &Scoring,
    smin: i32,
) -> Option<Extended> {
    debug_assert_eq!(ref_window.len(), read.len());

    let mut score: i32 = 0;
    let mut mismatches: u32 = 0;
    let mut md = String::new();
    let mut md_run: u32 = 0;

    for i in 0..read.len() {
        let r = read[i];
        let g = ref_window[i];
        let g_is_acgt = matches!(g, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
        let is_match = g_is_acgt && r.eq_ignore_ascii_case(&g);
        if is_match {
            md_run += 1;
        } else {
            let _ = write!(&mut md, "{md_run}");
            md_run = 0;
            md.push(g as char);
            mismatches += 1;
            let q = qual
                .get(i)
                .copied()
                .map(phred33_to_q)
                .unwrap_or(sc.mp_max_q);
            score -= sc.mm_penalty(q);
            // Hopeless: every remaining op only loses more score.
            if score < smin {
                return None;
            }
        }
    }
    let _ = write!(&mut md, "{md_run}");

    let mut cigar = String::new();
    let _ = write!(&mut cigar, "{}M", read.len());

    Some(Extended {
        ref_off,
        mismatches,
        n_ins: 0,
        n_del: 0,
        gap_opens: 0,
        score,
        cigar,
        md,
    })
}

/// Run semiglobal SW (read fully aligned, ref free at ends) of `read`
/// against the reference window starting at `win_off`. Returns `None` if
/// the window is too small or the alignment can't be recovered.
pub fn extend(
    read: &[u8],
    qual: &[u8],
    ref_window: &[u8],
    win_off: u32,
    sc: &Scoring,
) -> Option<Extended> {
    if ref_window.is_empty() || ref_window.len() < read.len() {
        return None;
    }

    // Fixed -mp_max mismatch for the SW search to find the alignment
    // structure; we rescore with Q-scaling below to get the per-position
    // penalty BT2 actually computes.
    let mp_max = sc.mp_max;
    let match_fn = move |a: u8, b: u8| -> i32 {
        if a.eq_ignore_ascii_case(&b) {
            0
        } else {
            -mp_max
        }
    };
    // Pick gap penalties; if open differs across read/ref directions, the
    // bio::pairwise aligner doesn't support asymmetric gaps natively. Use
    // the read-gap (D in CIGAR — read missing bases) penalty here. For
    // BT2's default these are identical.
    let mut aligner = Aligner::with_capacity(
        read.len(),
        ref_window.len(),
        -sc.rdg_open,
        -sc.rdg_extend,
        match_fn,
    );
    let aln = aligner.semiglobal(read, ref_window);

    // Walk ops to build CIGAR / MD / counts and rescore with Q.
    let mut cigar = String::new();
    let mut md = String::new();
    let mut md_run: u32 = 0;
    let mut cigar_run: u32 = 0;
    let mut cigar_op: u8 = 0; // 'M', 'I', 'D'
    let mut mismatches: u32 = 0;
    let mut n_ins: u32 = 0;
    let mut n_del: u32 = 0;
    let mut gap_opens: u32 = 0;
    let mut score: i32 = 0;
    let mut prev_was_gap: bool = false;
    let mut read_pos = aln.xstart;
    let mut ref_pos = aln.ystart;
    // We're tracking deletion runs in MD: emit `^<bases>` on the run end.
    let mut del_run: Vec<u8> = Vec::new();
    let flush_del_md = |del_run: &mut Vec<u8>, md: &mut String, md_run: &mut u32| {
        if !del_run.is_empty() {
            md.push_str(&md_run.to_string());
            *md_run = 0;
            md.push('^');
            for b in del_run.drain(..) {
                md.push(b as char);
            }
        }
    };

    for op in &aln.operations {
        match op {
            Op::Match | Op::Subst => {
                flush_del_md(&mut del_run, &mut md, &mut md_run);
                let r = read[read_pos];
                let g = ref_window[ref_pos];
                let is_match = r.eq_ignore_ascii_case(&g)
                    && matches!(g, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
                if is_match {
                    md_run += 1;
                } else {
                    md.push_str(&md_run.to_string());
                    md_run = 0;
                    md.push(g as char);
                    mismatches += 1;
                    let q = qual
                        .get(read_pos)
                        .copied()
                        .map(phred33_to_q)
                        .unwrap_or(sc.mp_max_q);
                    score -= sc.mm_penalty(q);
                }
                push_cigar(&mut cigar, &mut cigar_op, &mut cigar_run, b'M');
                read_pos += 1;
                ref_pos += 1;
                prev_was_gap = false;
            }
            Op::Ins => {
                // Insertion in read = "ref gap" in BT2 terminology
                // (read has bases not in reference) → use rfg penalties.
                if !prev_was_gap {
                    gap_opens += 1;
                    score -= sc.rfg_open;
                }
                score -= sc.rfg_extend;
                n_ins += 1;
                push_cigar(&mut cigar, &mut cigar_op, &mut cigar_run, b'I');
                read_pos += 1;
                prev_was_gap = true;
            }
            Op::Del => {
                // Deletion from read = "read gap" in BT2 (read missing
                // bases the reference has) → use rdg penalties.
                if !prev_was_gap {
                    gap_opens += 1;
                    score -= sc.rdg_open;
                }
                score -= sc.rdg_extend;
                n_del += 1;
                del_run.push(ref_window[ref_pos]);
                push_cigar(&mut cigar, &mut cigar_op, &mut cigar_run, b'D');
                ref_pos += 1;
                prev_was_gap = true;
            }
            Op::Xclip(_) | Op::Yclip(_) => {
                // Semiglobal clip ops are filtered by bio's aligner already.
            }
        }
    }
    flush_del_md(&mut del_run, &mut md, &mut md_run);
    md.push_str(&md_run.to_string());
    flush_cigar(&mut cigar, cigar_op, cigar_run);

    let ref_off = win_off + aln.ystart as u32;
    Some(Extended {
        ref_off,
        mismatches,
        n_ins,
        n_del,
        gap_opens,
        score,
        cigar,
        md,
    })
}

fn push_cigar(cigar: &mut String, cur_op: &mut u8, cur_run: &mut u32, op: u8) {
    if *cur_op == op {
        *cur_run += 1;
    } else {
        flush_cigar(cigar, *cur_op, *cur_run);
        *cur_op = op;
        *cur_run = 1;
    }
}

fn flush_cigar(cigar: &mut String, op: u8, run: u32) {
    if op != 0 && run > 0 {
        cigar.push_str(&run.to_string());
        cigar.push(op as char);
    }
}

#[allow(dead_code)]
const _: i32 = MIN_SCORE; // keep import used in case bio is updated
