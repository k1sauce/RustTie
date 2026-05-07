//! BowTie 2's V2 MAPQ formula (`BowtieMapq2` in `vendor/bowtie2/unique.h`).
//!
//! End-to-end mode only. Inputs are alignment scores in BT2's signed
//! convention (0 = perfect, negative = penalized). `score_min` is the floor
//! from `--score-min L,-0.6,-0.6` (i.e., -0.6 * read_len for L=50 → -30).

/// Compute MAPQ for a single-end alignment in end-to-end mode.
///
/// `best` is the AS score of the chosen alignment.
/// `secbest` is the AS score of the next-best alignment, or `None` if there
/// is no other valid alignment.
/// `score_min` is the floor (`scoreMin.f<TAlScore>(rdlen)`).
///
/// Direct port of `BowtieMapq2::mapq` for `sc_.monotone = true`
/// (`vendor/bowtie2/unique.h:223-332`). `scPer = 0` in end-to-end.
pub fn mapq_v2(best: i32, secbest: Option<i32>, score_min: i32) -> u8 {
    let sc_per: i32 = 0;
    let diff = (sc_per - score_min).max(1) as f64;
    let best_over = (best - score_min) as f64;

    if let Some(sec) = secbest {
        let bestdiff = (best - sec).abs() as f64;
        if bestdiff >= diff * 0.9 {
            if best_over >= diff { 39 } else { 33 }
        } else if bestdiff >= diff * 0.8 {
            if best_over >= diff { 38 } else { 27 }
        } else if bestdiff >= diff * 0.7 {
            if best_over >= diff { 37 } else { 26 }
        } else if bestdiff >= diff * 0.6 {
            if best_over >= diff { 36 } else { 22 }
        } else if bestdiff >= diff * 0.5 {
            if best_over >= diff {
                35
            } else if best_over >= diff * 0.84 {
                25
            } else if best_over >= diff * 0.68 {
                16
            } else {
                5
            }
        } else if bestdiff >= diff * 0.4 {
            if best_over >= diff {
                34
            } else if best_over >= diff * 0.84 {
                21
            } else if best_over >= diff * 0.68 {
                14
            } else {
                4
            }
        } else if bestdiff >= diff * 0.3 {
            if best_over >= diff {
                32
            } else if best_over >= diff * 0.88 {
                18
            } else if best_over >= diff * 0.67 {
                15
            } else {
                3
            }
        } else if bestdiff >= diff * 0.2 {
            if best_over >= diff {
                31
            } else if best_over >= diff * 0.88 {
                17
            } else if best_over >= diff * 0.67 {
                11
            } else {
                0
            }
        } else if bestdiff >= diff * 0.1 {
            if best_over >= diff {
                30
            } else if best_over >= diff * 0.88 {
                12
            } else if best_over >= diff * 0.67 {
                7
            } else {
                0
            }
        } else if bestdiff > 0.0 {
            if best_over >= diff * 0.67 { 6 } else { 2 }
        } else {
            // bestdiff == 0
            if best_over >= diff * 0.67 { 1 } else { 0 }
        }
    } else {
        // No second-best.
        if best_over >= diff * 0.8 {
            42
        } else if best_over >= diff * 0.7 {
            40
        } else if best_over >= diff * 0.6 {
            24
        } else if best_over >= diff * 0.5 {
            23
        } else if best_over >= diff * 0.4 {
            8
        } else if best_over >= diff * 0.3 {
            3
        } else {
            0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// L=50, score_min=-30, diff=30. Perfect alignment, no secbest → 42.
    #[test]
    fn unique_perfect() {
        assert_eq!(mapq_v2(0, None, -30), 42);
    }

    /// 1 high-Q mismatch (best=-6), no secbest. bestOver=24, ratio=0.8 → 42.
    #[test]
    fn unique_one_mm() {
        assert_eq!(mapq_v2(-6, None, -30), 42);
    }

    /// 1 indel (best=-8). bestOver=22, ratio=0.733 → 40 (between 0.7 and 0.8).
    #[test]
    fn unique_one_indel() {
        assert_eq!(mapq_v2(-8, None, -30), 40);
    }

    /// Multi-mapping perfect (best=0, secbest=0, bestdiff=0). → 1 (top branch).
    #[test]
    fn multi_perfect() {
        // bestdiff = 0, best_over = 30 = diff, so best_over >= diff * 0.67 → 1.
        assert_eq!(mapq_v2(0, Some(0), -30), 1);
    }

    /// Best=0, secbest=-6, bestdiff=6, ratio=0.2.
    /// best_over=diff → 31 (per the table at bestdiff bin 0.2).
    #[test]
    fn perfect_with_close_secbest() {
        assert_eq!(mapq_v2(0, Some(-6), -30), 31);
    }
}
