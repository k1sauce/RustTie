//! BWT character access and `LF(c, r) = C[c] + rank(c, r)` lookup.
//!
//! The BWT is laid out in fixed-size "sides" of `side_sz` bytes. Within a
//! side, the first `side_bwt_sz` bytes are 2-bit packed BWT data (LSB-first
//! within each byte for forward sides), followed by 4 × `OFF_SIZE` bytes of
//! cumulative occurrence checkpoints for A/C/G/T (count from BWT[0..side_start)).
//!
//! The `$` symbol is encoded as 'A' at row `z_off` but excluded from the
//! checkpoints (bt2_idx.h:2955-2963: `count = false` when `saElt == 0`).
//! So we adjust the in-side count, not the checkpoint.

use byteorder::{ByteOrder, LittleEndian};

use crate::format::{EbwtParams, OFF_SIZE};
use crate::reader::Bt2Index;

/// Position of the `$` row, decomposed into byte/bp coords.
#[derive(Debug, Clone, Copy)]
pub struct DollarPos {
    /// Absolute byte offset into `ebwt` where `z_off` lives.
    pub byte_off: u32,
    /// Bit-pair within that byte (0..=3).
    pub bp_off: u32,
}

impl DollarPos {
    pub fn new(z_off: u32, params: &EbwtParams) -> Self {
        let side_num = z_off / params.side_bwt_len;
        let side_char_off = z_off % params.side_bwt_len;
        let byte_off = side_num * params.side_sz + (side_char_off >> 2);
        let bp_off = side_char_off & 3;
        Self { byte_off, bp_off }
    }
}

/// Read the 2-bit char at BWT row `r` (LSB-first packing within the byte).
/// Returns the encoded character. Caller must check whether `r == z_off` to
/// decide if the value represents `$` rather than 'A'.
#[inline]
pub fn bwt_char(ebwt: &[u8], params: &EbwtParams, r: u32) -> u8 {
    let side_num = r / params.side_bwt_len;
    let off_in_side = r % params.side_bwt_len;
    let byte_idx = (side_num * params.side_sz + (off_in_side >> 2)) as usize;
    let bp = (off_in_side & 3) as u8;
    (ebwt[byte_idx] >> (bp * 2)) & 0b11
}

/// `LF(c, r) = C[c] + rank(c, r)` — the row in the BWT corresponding to the
/// character preceding the suffix at row `r` (when that character is `c`).
///
/// `r` may be in `[0, bwt_len]`; the upper bound is used as the exclusive
/// end of an SA range.
pub fn lf(idx: &Bt2Index, dollar: DollarPos, c: u8, r: u32) -> u32 {
    let p = &idx.params;
    debug_assert!(c < 4, "c must be 0..=3");
    debug_assert!(r <= p.bwt_len, "r out of range");

    // Decompose r into (side, byte-in-side, bp-in-byte).
    let side_num = r / p.side_bwt_len;
    let off_in_side = r % p.side_bwt_len;
    let byte_in_side = off_in_side >> 2; // full bytes to scan
    let bp_in_byte = off_in_side & 3; // partial chars in next byte

    let side_start = (side_num * p.side_sz) as usize;
    let bwt_region = &idx.ebwt[side_start..side_start + p.side_bwt_sz as usize];

    // Count `c` in full bytes [0, byte_in_side). Hot inner loop on every
    // FM-index step; SWAR-popcount processes 8 bytes (32 packed chars) per
    // iteration with no branches.
    let mut c_count: u32 = count_in_bytes_swar(&bwt_region[..byte_in_side as usize], c);
    // Partial byte: count first `bp_in_byte` chars of byte `byte_in_side`.
    if bp_in_byte > 0 && (byte_in_side as usize) < bwt_region.len() {
        c_count += count_in_byte_partial(bwt_region[byte_in_side as usize], c, bp_in_byte);
    }

    // Adjust for the `$` masquerading as 'A': subtract 1 if we passed z_off.
    if c == 0 {
        let cur_byte_abs = side_start as u32 + byte_in_side;
        let z_byte = dollar.byte_off;
        let z_bp = dollar.bp_off;
        let passed = cur_byte_abs > z_byte || (cur_byte_abs == z_byte && bp_in_byte > z_bp);
        if passed && side_start as u32 <= z_byte {
            c_count -= 1;
        }
    }

    // Add the per-side checkpoint: u32 at [side_start + side_bwt_sz + c*4].
    let occ_base = side_start + p.side_bwt_sz as usize;
    let checkpoint = LittleEndian::read_u32(
        &idx.ebwt[occ_base + (c as usize) * OFF_SIZE..occ_base + (c as usize + 1) * OFF_SIZE],
    );

    // C[c] (fchr) + rank(c, r).
    idx.fchr[c as usize] + checkpoint + c_count
}

/// Count occurrences of 2-bit code `c` in a packed byte (4 chars).
#[inline]
fn count_in_byte(byte: u8, c: u8) -> u32 {
    let mut n = 0;
    for bp in 0..4 {
        if (byte >> (bp * 2)) & 0b11 == c {
            n += 1;
        }
    }
    n
}

/// Count occurrences of 2-bit code `c` across `bytes`, processing 8 bytes
/// at a time via SWAR. Each u64 chunk is XORed with a replicated pattern of
/// `c`, then we mark each 2-bit field as 0/1 (1 = matched) and popcount.
/// 5 ops per 32 chars — much tighter than a per-bp branch loop.
#[inline]
fn count_in_bytes_swar(bytes: &[u8], c: u8) -> u32 {
    debug_assert!(c < 4);
    // Replicate c into all 32 2-bit positions: e.g. c=2 → 0xAAAA_AAAA_AAAA_AAAA.
    let one_byte: u8 = c | (c << 2) | (c << 4) | (c << 6);
    let pattern: u64 = u64::from(one_byte) * 0x0101_0101_0101_0101;
    const LO: u64 = 0x5555_5555_5555_5555;

    let mut count: u32 = 0;
    let chunks = bytes.chunks_exact(8);
    let rem = chunks.remainder();
    for chunk in chunks {
        // Safety: chunks_exact yields slices of length 8.
        let u = u64::from_le_bytes(chunk.try_into().unwrap());
        let x = u ^ pattern;
        // After XOR, each 2-bit field is 0 iff it matched c. Mark non-zero
        // fields with 1 in the low bit, then invert so matches → 1.
        let nonzero = (x | (x >> 1)) & LO;
        let matches = nonzero ^ LO;
        count += matches.count_ones();
    }
    for &b in rem {
        count += count_in_byte(b, c);
    }
    count
}

/// Count occurrences of `c` in the first `bp` characters of a packed byte
/// (the lowest `bp * 2` bits, LSB-first).
#[inline]
fn count_in_byte_partial(byte: u8, c: u8, bp: u32) -> u32 {
    let mut n = 0;
    for i in 0..bp {
        if (byte >> (i * 2)) & 0b11 == c {
            n += 1;
        }
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_base() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
    }

    #[test]
    fn lf_at_zero_is_fchr() {
        // LF(c, 0) = fchr[c] + 0 (rank at row 0 is 0 for all c).
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let dollar = DollarPos::new(idx.z_off, &idx.params);
        for c in 0..4 {
            assert_eq!(lf(&idx, dollar, c, 0), idx.fchr[c as usize]);
        }
    }

    #[test]
    fn lf_at_bwt_len_is_fchr_next() {
        // LF(c, bwt_len) = fchr[c+1]: total count of chars <= c.
        let idx = Bt2Index::open(fixture_base()).unwrap();
        let dollar = DollarPos::new(idx.z_off, &idx.params);
        for c in 0..4 {
            assert_eq!(
                lf(&idx, dollar, c, idx.params.bwt_len),
                idx.fchr[c as usize + 1],
                "c={c}",
            );
        }
    }

    #[test]
    fn count_in_byte_basic() {
        // Byte 0b11_10_01_00 (LSB-first: A, C, G, T) → one of each.
        let b = 0b11_10_01_00;
        assert_eq!(count_in_byte(b, 0), 1); // A
        assert_eq!(count_in_byte(b, 1), 1); // C
        assert_eq!(count_in_byte(b, 2), 1); // G
        assert_eq!(count_in_byte(b, 3), 1); // T

        // Partial: first 2 chars are A, C.
        assert_eq!(count_in_byte_partial(b, 0, 2), 1);
        assert_eq!(count_in_byte_partial(b, 1, 2), 1);
        assert_eq!(count_in_byte_partial(b, 2, 2), 0);
        assert_eq!(count_in_byte_partial(b, 3, 2), 0);
    }
}
