//! FM-index backward search and SA → reference-position resolution.

use crate::bwt::{DollarPos, lf};
use crate::format::ascii_to_2bit;
use crate::reader::Bt2Index;

/// Suffix-array range `[lo, hi)` produced by backward search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SaRange {
    pub lo: u32,
    pub hi: u32,
}

impl SaRange {
    pub fn len(&self) -> u32 {
        self.hi.saturating_sub(self.lo)
    }

    pub fn is_empty(&self) -> bool {
        self.hi <= self.lo
    }
}

/// Resolved hit: (reference id, 0-based offset within that reference).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefHit {
    pub ref_id: u32,
    pub ref_off: u32,
}

/// Exact-match backward search over `query` (ASCII A/C/G/T, any case).
/// Returns `None` if any non-ACGT char is in the query (BT2 treats N's as
/// non-matching for exact seed search; mismatch handling is for later phases).
pub fn backward_search(idx: &Bt2Index, query: &[u8]) -> Option<SaRange> {
    if query.is_empty() {
        return Some(SaRange {
            lo: 0,
            hi: idx.params.bwt_len,
        });
    }
    let dollar = DollarPos::new(idx.z_off, &idx.params);

    // Initial range: full BWT.
    let mut lo: u32 = 0;
    let mut hi: u32 = idx.params.bwt_len;

    // Walk query right-to-left, narrowing [lo, hi).
    for &b in query.iter().rev() {
        let c = ascii_to_2bit(b)?;
        lo = lf(idx, dollar, c, lo);
        hi = lf(idx, dollar, c, hi);
        if hi <= lo {
            return Some(SaRange { lo, hi: lo });
        }
    }
    Some(SaRange { lo, hi })
}

/// Resolve a single SA-range row to its text position (joined-text offset).
///
/// Walks `LF` from `row` until we hit either `z_off` (text position 0) or a
/// row whose offset is sampled in `idx.offs`.
pub fn resolve_text_pos(idx: &Bt2Index, mut row: u32) -> u32 {
    let dollar = DollarPos::new(idx.z_off, &idx.params);
    let mut steps: u32 = 0;
    loop {
        if row == idx.z_off {
            return steps;
        }
        // Sampled if low `off_rate` bits of row are 0.
        let off_mask = (1u32 << idx.params.off_rate) - 1;
        if row & off_mask == 0 {
            let sampled = idx.offs[(row >> idx.params.off_rate) as usize];
            return sampled + steps;
        }
        // Walk LF: read BWT[row] as c, then row = LF(c, row).
        let c = crate::bwt::bwt_char(&idx.ebwt, &idx.params, row);
        // BWT[z_off] is encoded as 'A' but represents $; we should not be
        // walking from there because it's handled above.
        row = lf(idx, dollar, c, row);
        steps += 1;
        debug_assert!(steps < idx.params.bwt_len, "LF walk did not terminate");
    }
}

/// Map a joined-text position to (ref_id, ref_off) using the `rstarts` table.
///
/// For the lambda spike there's a single fragment so the answer is trivial,
/// but the binary-search version handles arbitrary multi-fragment indexes.
pub fn joined_to_ref(idx: &Bt2Index, joined_pos: u32) -> RefHit {
    // rstarts is sorted by joined_off ascending. Find the largest entry
    // with joined_off <= joined_pos.
    let starts = &idx.rstarts;
    let mut lo = 0usize;
    let mut hi = starts.len();
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if starts[mid].joined_off <= joined_pos {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let frag = starts[lo];
    RefHit {
        ref_id: frag.ref_id,
        ref_off: frag.ref_off + (joined_pos - frag.joined_off),
    }
}

/// Convenience: enumerate all reference hits for an exact-match query.
pub fn exact_hits(idx: &Bt2Index, query: &[u8]) -> Vec<RefHit> {
    let Some(range) = backward_search(idx, query) else {
        return Vec::new();
    };
    if range.is_empty() {
        return Vec::new();
    }
    (range.lo..range.hi)
        .map(|row| {
            let pos = resolve_text_pos(idx, row);
            joined_to_ref(idx, pos)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_base() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
    }

    /// Trivial case: empty query → full BWT range.
    #[test]
    fn empty_query_full_range() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let r = backward_search(&idx, b"").unwrap();
        assert_eq!(r.lo, 0);
        assert_eq!(r.hi, idx.params.bwt_len);
    }

    /// Single-char query: hi - lo = count of that char in BWT = fchr[c+1] - fchr[c].
    #[test]
    fn single_char_count_matches_fchr() {
        let idx = Bt2Index::open(fixture_base()).unwrap();
        for (c, ch) in b"ACGT".iter().enumerate() {
            let r = backward_search(&idx, &[*ch]).unwrap();
            let expected = idx.fchr[c + 1] - idx.fchr[c];
            assert_eq!(
                r.len(),
                expected,
                "char {} ({c}): hi-lo={}, expected {expected}",
                *ch as char,
                r.len()
            );
        }
    }
}
