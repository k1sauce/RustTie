//! `.bt2` binary file format constants and derived parameters.
//!
//! Authoritative source: `vendor/bowtie2/bt2_idx.h` and `bt2_io.cpp`.
//! Spec: `docs/bt2-format.md`.
//!
//! Spike scope: 32-bit "small" index built by `bowtie2-build-s`.
//! 64-bit ("large") indexes use `OFF_SIZE = 8`; not yet supported.

/// Size of `TIndexOffU` on disk for the small index.
pub const OFF_SIZE: usize = 4;

/// Endianness sentinel: u32 == 1 → native LE; == 0x01000000 → swap.
pub const ENDIAN_NATIVE: u32 = 1;
pub const ENDIAN_SWAPPED: u32 = 0x0100_0000;

/// Flag bits in the header (bt2_idx.h:820-823).
pub const FLAG_VALID: i32 = 1;
pub const FLAG_COLOR: i32 = 2;
pub const FLAG_ENTIRE_REV: i32 = 4;

/// 2-bit DNA encoding (bt2_idx.h packing).
pub const A: u8 = 0;
pub const C: u8 = 1;
pub const G: u8 = 2;
pub const T: u8 = 3;

/// Convert ASCII nucleotide → 2-bit code. `None` for ambiguous.
pub fn ascii_to_2bit(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(A),
        b'C' | b'c' => Some(C),
        b'G' | b'g' => Some(G),
        b'T' | b't' => Some(T),
        _ => None,
    }
}

/// Geometry of the FM-index, derived from the on-disk header.
/// Mirrors `EbwtParams` (bt2_idx.h:112-254).
#[derive(Debug, Clone)]
pub struct EbwtParams {
    /// Text length (excluding `$`). Header field.
    pub len: u32,
    /// BWT length (= `len + 1`, includes `$`).
    pub bwt_len: u32,
    /// `log2(lineSz)`. Header field.
    pub line_rate: u32,
    /// Bytes per BWT side block.
    pub line_sz: u32,
    /// `log2(SA sample stride)`. Header field.
    pub off_rate: u32,
    /// Stride between sampled SA entries.
    pub off_stride: u32,
    /// Characters in the ftab key. Header field.
    pub ftab_chars: u32,
    /// Number of ftab entries: `(1 << (ftab_chars * 2)) + 1`.
    pub ftab_len: u32,
    /// Eftab length: `ftab_chars * 2`.
    pub eftab_len: u32,
    /// Number of sampled SA entries.
    pub offs_len: u32,
    /// Bytes per side (= `line_sz` since `linesPerSide` is deprecated to 1).
    pub side_sz: u32,
    /// BWT bytes within a side: `side_sz - 4*OFF_SIZE`.
    pub side_bwt_sz: u32,
    /// BWT chars within a side (4 chars per byte): `side_bwt_sz * 4`.
    pub side_bwt_len: u32,
    /// Total number of sides.
    pub num_sides: u32,
    /// Total BWT+checkpoint storage: `num_sides * side_sz`.
    pub ebwt_tot_len: u32,
    /// Colorspace index (we don't support this in the spike).
    pub color: bool,
    /// Reverse-entire flag.
    pub entire_reverse: bool,
}

impl EbwtParams {
    /// Compute derived params from the four header fields + flags.
    /// Mirrors `Ebwt::init()` (bt2_idx.h:133-167).
    pub fn from_header(
        len: u32,
        line_rate: u32,
        off_rate: u32,
        ftab_chars: u32,
        flags: i32,
    ) -> Self {
        let bwt_len = len + 1;
        let line_sz = 1u32 << line_rate;
        let off_stride = 1u32 << off_rate;
        let ftab_len = (1u32 << (ftab_chars * 2)) + 1;
        let eftab_len = ftab_chars * 2;
        let offs_len = (bwt_len + off_stride - 1) >> off_rate;

        let side_sz = line_sz;
        let side_bwt_sz = side_sz - 4 * OFF_SIZE as u32;
        let side_bwt_len = side_bwt_sz * 4;

        // bwt_sz here is the count of BWT chars (= bwt_len) needing storage.
        let num_sides = bwt_len.div_ceil(side_bwt_len);
        let ebwt_tot_len = num_sides * side_sz;

        let abs_flags = -flags;
        Self {
            len,
            bwt_len,
            line_rate,
            line_sz,
            off_rate,
            off_stride,
            ftab_chars,
            ftab_len,
            eftab_len,
            offs_len,
            side_sz,
            side_bwt_sz,
            side_bwt_len,
            num_sides,
            ebwt_tot_len,
            color: (abs_flags & FLAG_COLOR) != 0,
            entire_reverse: (abs_flags & FLAG_ENTIRE_REV) != 0,
        }
    }
}
