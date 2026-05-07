//! Paired-end concordance check + paired-alignment summary.
//!
//! BT2-faithful (Phase 3l): the per-mate alignment sets are combined via
//! Cartesian product to enumerate concordant pairs, sorted by `r1.score +
//! r2.score` descending. The best concordant pair is selected as the
//! displayed alignment and the second-best concordant *pair* score is
//! threaded into MAPQ as the `secbest` input — this matches BT2's
//! `selectByScore` + `BowtieMapq2::mapq` for paired reads (see
//! `vendor/bowtie2/aln_sink.cpp:1580` and `vendor/bowtie2/unique.h:218-235`).
//! Per-mate secbest is **not** used for paired MAPQ — only the pair sum is.

use crate::align::{Alignment, Strand};

/// BT2 default minimum fragment length (`-I 0`).
pub const FRAG_LEN_MIN: u32 = 0;
/// BT2 default maximum fragment length (`-X 500`).
pub const FRAG_LEN_MAX: u32 = 500;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairType {
    /// Concordant pair: FR orientation, fragment in [I, X], same chrom.
    Concordant,
    /// Both mates aligned but not concordantly.
    Discordant,
    /// At least one mate failed to align — pair is "unpaired".
    Unpaired,
}

/// Outcome bundle for a read pair.
#[derive(Debug, Clone)]
pub struct PairOutcome {
    pub r1: Option<Alignment>,
    pub r2: Option<Alignment>,
    /// Per-mate secbest scores — used only for **single-end** fallback paths
    /// (one mate mapped, mate's own MAPQ). Concordant pairs use
    /// `concordant_pair_secbest` instead.
    pub r1_secbest: Option<i32>,
    pub r2_secbest: Option<i32>,
    pub pair_type: PairType,
    /// Fragment length (always positive). 0 if not concordant or one unmapped.
    pub frag_len: u32,
    /// For concordant pairs: sum-of-mate-scores for the second-best
    /// concordant pair found (or `None` if only one concordant pair existed).
    /// `BowtieMapq2` uses this — paired MAPQ depends on the next-best
    /// **pair** score, not on per-mate alternates.
    pub concordant_pair_secbest: Option<i32>,
    /// Additional concordant pairs after the primary (`r1`/`r2`), sorted by
    /// pair score descending. Populated only for `-k <int>` / `-a`
    /// reporting. Each entry is `(r1, r2, frag_len)`.
    pub additional_concordant: Vec<(Alignment, Alignment, u32)>,
}

/// Test whether two mate alignments form a concordant `--fr` pair.
/// Returns `(true, frag_len)` if concordant under BT2's defaults
/// (one strand each, FR orientation, fragment in `[FRAG_LEN_MIN, FRAG_LEN_MAX]`,
/// same reference). The frag_len is always positive.
pub fn is_concordant(a1: &Alignment, a2: &Alignment) -> Option<u32> {
    if a1.ref_id != a2.ref_id || a1.strand == a2.strand {
        return None;
    }
    let (fwd, rev) = if a1.strand == Strand::Forward {
        (a1, a2)
    } else {
        (a2, a1)
    };
    if fwd.ref_off > rev.ref_off + rev.read_len.saturating_sub(1) {
        return None;
    }
    let frag_len = (rev.ref_off + rev.read_len) - fwd.ref_off;
    if !(FRAG_LEN_MIN..=FRAG_LEN_MAX).contains(&frag_len) {
        return None;
    }
    Some(frag_len)
}

/// Enumerate concordant pairs from per-mate alignment sets, pick the
/// best-scoring pair (ties broken by leftmost r1 position then leftmost r2
/// position) and return the second-best pair's score for paired MAPQ. When
/// `report_k > 1`, additional concordant pairs (in score-descending order)
/// are returned in `additional_concordant` for `-k`/`-a` reporting.
///
/// Falls back to discordant (best-scoring r1 + best-scoring r2 with no
/// concordance constraint) or unpaired if either or both mates failed to
/// align — matches BT2's reporting decision tree.
pub fn classify_pair_set(
    r1_alns: &[Alignment],
    r2_alns: &[Alignment],
    r1_secbest: Option<i32>,
    r2_secbest: Option<i32>,
    report_k: u32,
) -> PairOutcome {
    if r1_alns.is_empty() || r2_alns.is_empty() {
        let r1 = r1_alns.first().cloned();
        let r2 = r2_alns.first().cloned();
        return PairOutcome {
            r1,
            r2,
            r1_secbest,
            r2_secbest,
            pair_type: PairType::Unpaired,
            frag_len: 0,
            concordant_pair_secbest: None,
            additional_concordant: Vec::new(),
        };
    }

    // Cartesian product, filter by concordance, sort by pair score.
    let mut concordant: Vec<(usize, usize, u32, i32)> = Vec::new();
    for (i, a1) in r1_alns.iter().enumerate() {
        for (j, a2) in r2_alns.iter().enumerate() {
            if let Some(frag_len) = is_concordant(a1, a2) {
                let score_sum = a1.score + a2.score;
                concordant.push((i, j, frag_len, score_sum));
            }
        }
    }
    if !concordant.is_empty() {
        concordant.sort_by(|a, b| {
            b.3.cmp(&a.3)
                .then(r1_alns[a.0].ref_id.cmp(&r1_alns[b.0].ref_id))
                .then(r1_alns[a.0].ref_off.cmp(&r1_alns[b.0].ref_off))
                .then(r2_alns[a.1].ref_off.cmp(&r2_alns[b.1].ref_off))
        });
        let (i, j, frag_len, _best_score) = concordant[0];
        let concordant_pair_secbest = concordant.get(1).map(|c| c.3);
        // Collect up to (report_k - 1) additional concordant pairs after
        // the primary, deduped against the primary itself.
        let mut additional = Vec::new();
        if report_k > 1 {
            let want = (report_k as usize).saturating_sub(1);
            for c in concordant.iter().skip(1).take(want) {
                additional.push((r1_alns[c.0].clone(), r2_alns[c.1].clone(), c.2));
            }
        }
        return PairOutcome {
            r1: Some(r1_alns[i].clone()),
            r2: Some(r2_alns[j].clone()),
            r1_secbest,
            r2_secbest,
            pair_type: PairType::Concordant,
            frag_len,
            concordant_pair_secbest,
            additional_concordant: additional,
        };
    }

    // No concordant pair — report best per-mate.
    PairOutcome {
        r1: Some(r1_alns[0].clone()),
        r2: Some(r2_alns[0].clone()),
        r1_secbest,
        r2_secbest,
        pair_type: PairType::Discordant,
        frag_len: 0,
        concordant_pair_secbest: None,
        additional_concordant: Vec::new(),
    }
}

/// Backwards-compat shim: single-best-per-mate, no `-k` reporting.
pub fn classify_pair(
    r1: Option<Alignment>,
    r2: Option<Alignment>,
    r1_secbest: Option<i32>,
    r2_secbest: Option<i32>,
) -> PairOutcome {
    let r1_vec: Vec<Alignment> = r1.into_iter().collect();
    let r2_vec: Vec<Alignment> = r2.into_iter().collect();
    classify_pair_set(&r1_vec, &r2_vec, r1_secbest, r2_secbest, 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(ref_off: u32, strand: Strand, len: u32) -> Alignment {
        Alignment {
            ref_id: 0,
            ref_off,
            strand,
            read_len: len,
            mismatches: 0,
            gap_opens: 0,
            gap_extends: 0,
            score: 0,
            cigar: format!("{len}M"),
            md: format!("{len}"),
        }
    }

    #[test]
    fn concordant_fr() {
        // R1 forward at 100, R2 reverse at 350. Both 50bp. Fragment = 350+50-100 = 300.
        let r1 = aln(100, Strand::Forward, 50);
        let r2 = aln(350, Strand::Reverse, 50);
        let p = classify_pair(Some(r1), Some(r2), None, None);
        assert_eq!(p.pair_type, PairType::Concordant);
        assert_eq!(p.frag_len, 300);
    }

    #[test]
    fn discordant_same_strand() {
        let r1 = aln(100, Strand::Forward, 50);
        let r2 = aln(350, Strand::Forward, 50);
        assert_eq!(
            classify_pair(Some(r1), Some(r2), None, None).pair_type,
            PairType::Discordant
        );
    }

    #[test]
    fn discordant_too_far() {
        let r1 = aln(100, Strand::Forward, 50);
        let r2 = aln(700, Strand::Reverse, 50); // 700+50-100 = 650 > 500
        assert_eq!(
            classify_pair(Some(r1), Some(r2), None, None).pair_type,
            PairType::Discordant
        );
    }

    #[test]
    fn discordant_wrong_orientation() {
        // R1 reverse to the LEFT of R2 forward — RF, not FR.
        let r1 = aln(100, Strand::Reverse, 50);
        let r2 = aln(350, Strand::Forward, 50);
        assert_eq!(
            classify_pair(Some(r1), Some(r2), None, None).pair_type,
            PairType::Discordant
        );
    }

    #[test]
    fn unpaired_one_unmapped() {
        let r1 = aln(100, Strand::Forward, 50);
        assert_eq!(
            classify_pair(Some(r1), None, None, None).pair_type,
            PairType::Unpaired
        );
    }
}
