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
///
/// IMPORTANT: BT2 writes the bin thresholds as `(double)0.9f` etc. — a
/// **single-precision** float literal cast to double. Using `0.9_f64`
/// directly diverges at boundary cases: for `bestdiff=12, diff=120`, BT2's
/// `120 * (double)0.1f = 12.0000001788...` makes `12 >= …` FALSE (bin
/// <0.1, MAPQ 6), while `120 * 0.1_f64 = 12.0` makes `12 >= 12.0` TRUE
/// (bin 0.1, MAPQ 12). On chr22 this off-by-one explained ~120 reads in
/// the `rt=12 bt=6` disagreement bin.
const C_09: f64 = 0.9_f32 as f64;
const C_08: f64 = 0.8_f32 as f64;
const C_07: f64 = 0.7_f32 as f64;
const C_06: f64 = 0.6_f32 as f64;
const C_05: f64 = 0.5_f32 as f64;
const C_04: f64 = 0.4_f32 as f64;
const C_03: f64 = 0.3_f32 as f64;
const C_02: f64 = 0.2_f32 as f64;
const C_01: f64 = 0.1_f32 as f64;
const C_088: f64 = 0.88_f32 as f64;
const C_084: f64 = 0.84_f32 as f64;
const C_068: f64 = 0.68_f32 as f64;
const C_067: f64 = 0.67_f32 as f64;

pub fn mapq_v2(best: i32, secbest: Option<i32>, score_min: i32) -> u8 {
    let sc_per: i32 = 0;
    let diff = (sc_per - score_min).max(1) as f64;
    let best_over = (best - score_min) as f64;

    if let Some(sec) = secbest {
        let bestdiff = (best - sec).abs() as f64;
        if bestdiff >= diff * C_09 {
            if best_over >= diff { 39 } else { 33 }
        } else if bestdiff >= diff * C_08 {
            if best_over >= diff { 38 } else { 27 }
        } else if bestdiff >= diff * C_07 {
            if best_over >= diff { 37 } else { 26 }
        } else if bestdiff >= diff * C_06 {
            if best_over >= diff { 36 } else { 22 }
        } else if bestdiff >= diff * C_05 {
            if best_over >= diff {
                35
            } else if best_over >= diff * C_084 {
                25
            } else if best_over >= diff * C_068 {
                16
            } else {
                5
            }
        } else if bestdiff >= diff * C_04 {
            if best_over >= diff {
                34
            } else if best_over >= diff * C_084 {
                21
            } else if best_over >= diff * C_068 {
                14
            } else {
                4
            }
        } else if bestdiff >= diff * C_03 {
            if best_over >= diff {
                32
            } else if best_over >= diff * C_088 {
                18
            } else if best_over >= diff * C_067 {
                15
            } else {
                3
            }
        } else if bestdiff >= diff * C_02 {
            if best_over >= diff {
                31
            } else if best_over >= diff * C_088 {
                17
            } else if best_over >= diff * C_067 {
                11
            } else {
                0
            }
        } else if bestdiff >= diff * C_01 {
            if best_over >= diff {
                30
            } else if best_over >= diff * C_088 {
                12
            } else if best_over >= diff * C_067 {
                7
            } else {
                0
            }
        } else if bestdiff > 0.0 {
            if best_over >= diff * C_067 { 6 } else { 2 }
        } else {
            // bestdiff == 0
            if best_over >= diff * C_067 { 1 } else { 0 }
        }
    } else {
        // No second-best.
        if best_over >= diff * C_08 {
            42
        } else if best_over >= diff * C_07 {
            40
        } else if best_over >= diff * C_06 {
            24
        } else if best_over >= diff * C_05 {
            23
        } else if best_over >= diff * C_04 {
            8
        } else if best_over >= diff * C_03 {
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

    /// 1 high-Q mismatch (best=-6), no secbest. bestOver=24, ratio=0.8.
    /// BT2's `(double)0.8f` evaluates 30 * 0.8f → 24.000000357..., so the
    /// `24 >= 24.0000003...` check is FALSE and we fall through to the
    /// 0.7-bin: MAPQ 40. (Pre-fix, with `0.8_f64`, this returned 42.)
    #[test]
    fn unique_one_mm() {
        assert_eq!(mapq_v2(-6, None, -30), 40);
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

    /// Best=0, secbest=-6, bestdiff=6, ratio=0.2 in true arithmetic.
    /// BT2's `(double)0.2f = 0.2000000029...` makes `30 * 0.2f =
    /// 6.000000089...`, so `6 >= 6.0000000...` is FALSE and we fall
    /// through to the 0.1-bin: best_over=diff → MAPQ 30. (Pre-fix returned
    /// 31 because `0.2_f64 * 30 = 6.0` exactly.)
    #[test]
    fn perfect_with_close_secbest() {
        assert_eq!(mapq_v2(0, Some(-6), -30), 30);
    }
}
