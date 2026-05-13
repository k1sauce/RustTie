//! Alignment: seed-and-extend with mismatches + indels (Phase 2d).
//! SIMD'd extension via `block-aligner` lands in Phase 3.

pub mod align;
pub mod bt2_descent;
pub mod bt2_random;
pub mod extend;
pub mod mapq;
pub mod paired;
pub mod paired_descent;
pub mod revcomp;

#[doc(hidden)]
pub use align::profile;

pub use align::{
    AlignResult, Alignment, DEFAULT_SEED_LEN, DESCENT_D_DEFAULT, DESCENT_R_DEFAULT, EXTEND_SLACK,
    MM_PENALTY_MAX, MM_PENALTY_MAX_Q, MM_PENALTY_MIN, PER_SEED_HIT_CAP, REPETITIVE_HITS_THRESHOLD,
    Scoring, Strand, align_read, align_read_with_cap, align_read_with_descent, mate_rescue,
    mm_penalty, phred33_to_q, score_min, seed_interval, seed_offsets,
};
pub use mapq::mapq_v2;
pub use paired::{
    FRAG_LEN_MAX, FRAG_LEN_MIN, PairCandidate, PairOutcome, PairType, classify_pair,
    classify_pair_set, is_concordant,
};
pub use paired_descent::{JointDescentResult, PAIR_POOL_CAP, align_pair_jointly};
pub use revcomp::reverse_complement;
