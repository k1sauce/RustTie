//! Reference sequence access from `.3.bt2` (records) + `.4.bt2` (packed
//! 2-bit DNA). This is what BT2 itself uses at align time, and lets us avoid
//! requiring a separate FASTA on the command line.
//!
//! Format (spec doc + reference.cpp:103-170, 178-225):
//! - `.3.bt2`: 4-byte endian sentinel, `num_records: u32`, then
//!   `num_records × (off: u32, len: u32, first: u8)` packed.
//! - `.4.bt2`: raw 2-bit packed bytes; total chars = sum of all record lens,
//!   padded up to a multiple of 4.

use std::io::{Cursor, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use byteorder::{LittleEndian, ReadBytesExt};

use crate::format::ENDIAN_NATIVE;

/// One unambiguous-stretch record from `.3.bt2`.
#[derive(Debug, Clone, Copy)]
pub struct RefRecord {
    /// Offset (chars) within the containing reference sequence.
    pub off: u32,
    /// Length (chars) of this stretch.
    pub len: u32,
    /// True iff this record starts a new reference sequence.
    pub first: bool,
}

/// 2-bit packed reference plus the records describing where each stretch
/// lies within its containing reference sequence.
pub struct BitPairReference {
    pub records: Vec<RefRecord>,
    /// 2-bit packed concatenation of all unambiguous stretches.
    /// LSB-first within each byte (matches `.4.bt2` packing).
    pub packed: Vec<u8>,
    /// Total unambiguous chars (sum of `records[i].len`).
    pub total_chars: u32,
    /// For each record, its starting offset in `packed` (in chars, not bytes).
    pub joined_offsets: Vec<u32>,
    /// For each record, the absolute starting offset within the containing
    /// reference sequence (= cumulative `off + len` of preceding records of
    /// the same reference, plus this record's `off`). Computed at load time
    /// so `locate` is O(records); BT2 stores `off` as a per-stretch run-length
    /// of preceding Ns rather than an absolute offset.
    pub ref_offsets: Vec<u32>,
    /// For each record, the reference index it belongs to. Derived from the
    /// `first` bits during load (each `first=true` starts a new reference).
    pub ref_ids: Vec<u32>,
    /// For each `ref_id`, `(lo, hi)` such that `records[lo..hi]` are the
    /// stretches belonging to that reference. Lets `locate` be O(log K)
    /// in records-per-ref instead of O(N) over all records.
    pub ref_record_ranges: Vec<(usize, usize)>,
}

impl BitPairReference {
    pub fn open(base: impl AsRef<Path>) -> Result<Self> {
        let base = base.as_ref();
        let p3 = path_with_ext(base, "3.bt2");
        let p4 = path_with_ext(base, "4.bt2");
        let mut buf3 = Vec::new();
        std::fs::File::open(&p3)
            .with_context(|| format!("opening {}", p3.display()))?
            .read_to_end(&mut buf3)?;
        let mut buf4 = Vec::new();
        std::fs::File::open(&p4)
            .with_context(|| format!("opening {}", p4.display()))?
            .read_to_end(&mut buf4)?;
        Self::from_bytes(&buf3, &buf4)
    }

    pub fn from_bytes(file3: &[u8], file4: &[u8]) -> Result<Self> {
        let mut c = Cursor::new(file3);
        let endian = c.read_u32::<LittleEndian>()?;
        if endian != ENDIAN_NATIVE {
            bail!("byte-swapped .3.bt2 not supported");
        }
        let n = c.read_u32::<LittleEndian>()?;
        if n == 0 {
            bail!("empty .3.bt2");
        }
        let mut records = Vec::with_capacity(n as usize);
        let mut joined_offsets = Vec::with_capacity(n as usize);
        let mut ref_offsets = Vec::with_capacity(n as usize);
        let mut ref_ids = Vec::with_capacity(n as usize);
        let mut total: u32 = 0;
        let mut cur_ref: i64 = -1;
        let mut ref_pos_acc: u32 = 0;
        for _ in 0..n {
            let off = c.read_u32::<LittleEndian>()?;
            let len = c.read_u32::<LittleEndian>()?;
            let mut first_buf = [0u8; 1];
            c.read_exact(&mut first_buf)?;
            let first = first_buf[0] != 0;
            if first {
                cur_ref += 1;
                ref_pos_acc = 0;
            }
            // Per BT2 (`bt2_idx.h:2721-2745`, `ref_read.cpp`): `off` is the
            // count of leading Ns BEFORE this stretch (within the current
            // reference sequence), `len` is the unambiguous-stretch length.
            // Absolute offset where this stretch's ACGT starts within the
            // reference = `ref_pos_acc + off`.
            let stretch_ref_off = ref_pos_acc + off;
            ref_offsets.push(stretch_ref_off);
            ref_ids.push(cur_ref as u32);
            joined_offsets.push(total);
            total += len;
            ref_pos_acc = stretch_ref_off + len;
            records.push(RefRecord { off, len, first });
        }

        // .4.bt2 should have ceil(total / 4) bytes.
        let expected_bytes = total.div_ceil(4) as usize;
        if file4.len() < expected_bytes {
            bail!(
                ".4.bt2 too short: have {} bytes, need at least {expected_bytes}",
                file4.len()
            );
        }

        // Records are emitted in ref order with `first=true` marking each
        // new reference, so per-ref ranges are simply the runs between
        // those markers.
        let mut ref_record_ranges: Vec<(usize, usize)> = Vec::new();
        let mut run_start: usize = 0;
        for (i, r) in records.iter().enumerate() {
            if r.first && i > 0 {
                ref_record_ranges.push((run_start, i));
                run_start = i;
            }
        }
        ref_record_ranges.push((run_start, records.len()));

        Ok(Self {
            records,
            packed: file4[..expected_bytes].to_vec(),
            total_chars: total,
            joined_offsets,
            ref_offsets,
            ref_ids,
            ref_record_ranges,
        })
    }

    /// Extract `len` ASCII bases starting at joined-text position `joined_pos`.
    pub fn extract_joined(&self, joined_pos: u32, len: u32) -> Vec<u8> {
        let mut out = Vec::with_capacity(len as usize);
        for i in 0..len {
            let p = joined_pos + i;
            let byte = self.packed[(p / 4) as usize];
            // LSB-first: bits [1:0] = char 0, bits [3:2] = char 1, ...
            // (Matches `.1.bt2` BWT packing; the format-spec doc was wrong.)
            let shift = 2 * (p % 4);
            let code = (byte >> shift) & 0b11;
            out.push(b"ACGT"[code as usize]);
        }
        out
    }

    /// Find the (joined offset, record index) covering `(ref_id, ref_off)`,
    /// or `None` if the position falls in an ambiguous gap or out of range.
    /// Uses precomputed `ref_offsets` (absolute stretch starts) so we don't
    /// confuse RefRecord's `off` field (preceding-N count) with an absolute
    /// position.
    ///
    /// O(log K) where K is the number of stretches in the target reference:
    /// per-ref range is an O(1) array lookup, then binary search on the
    /// monotonically-increasing `ref_offsets` slice for that ref.
    pub fn locate(&self, ref_id: u32, ref_off: u32) -> Option<(u32, usize)> {
        let &(lo, hi) = self.ref_record_ranges.get(ref_id as usize)?;
        if lo == hi {
            return None;
        }
        let slice = &self.ref_offsets[lo..hi];
        // partition_point: index of first element > ref_off, so the
        // candidate covering record (if any) is at index-1.
        let idx_in_slice = slice.partition_point(|&v| v <= ref_off);
        if idx_in_slice == 0 {
            return None;
        }
        let i = lo + idx_in_slice - 1;
        let start = self.ref_offsets[i];
        if ref_off < start + self.records[i].len {
            let joined = self.joined_offsets[i] + (ref_off - start);
            Some((joined, i))
        } else {
            None
        }
    }

    /// Extract `len` bases starting at `(ref_id, ref_off)`. Returns `None` if
    /// the requested span runs off the end of the (single) covering record;
    /// multi-record spans (with internal Ns) are not handled here.
    pub fn extract(&self, ref_id: u32, ref_off: u32, len: u32) -> Option<Vec<u8>> {
        let (joined, rec_idx) = self.locate(ref_id, ref_off)?;
        let rec = &self.records[rec_idx];
        let start = self.ref_offsets[rec_idx];
        if ref_off + len > start + rec.len {
            return None;
        }
        Some(self.extract_joined(joined, len))
    }
}

fn path_with_ext(base: &Path, ext: &str) -> std::path::PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_base() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
    }

    fn read_lambda_fasta() -> Vec<u8> {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../validate/fixtures/lambda_virus.fa");
        let text = std::fs::read_to_string(&path).unwrap();
        let mut out = Vec::new();
        for line in text.lines() {
            if !line.starts_with('>') {
                out.extend(line.trim().as_bytes());
            }
        }
        out
    }

    #[test]
    fn extracts_match_fasta() {
        let r = BitPairReference::open(fixture_base()).unwrap();
        let fasta = read_lambda_fasta();
        assert_eq!(r.total_chars as usize, fasta.len());
        assert_eq!(r.records.len(), 1, "lambda has one unambiguous stretch");
        // Spot-check three windows.
        for &(off, len) in &[(0u32, 50u32), (10_000, 100), (48_402, 100)] {
            let got = r.extract(0, off, len).unwrap();
            let want = &fasta[off as usize..(off + len) as usize];
            assert_eq!(
                std::str::from_utf8(&got).unwrap(),
                std::str::from_utf8(want).unwrap(),
                "off={off} len={len}"
            );
        }
    }
}
