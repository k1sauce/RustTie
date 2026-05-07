//! `.bt2` index reader.
//!
//! Spike scope: 32-bit small index, native little-endian only.
//! Loads `.1.bt2` (BWT + occ + ftab + refnames) and `.2.bt2` (SA samples).
//! `.3.bt2` and `.4.bt2` are skipped — not needed for SA-range → ref-pos.

use std::io::{Cursor, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use byteorder::{LittleEndian, ReadBytesExt};

use crate::format::{ENDIAN_NATIVE, ENDIAN_SWAPPED, EbwtParams};

/// Fragment record: contiguous unambiguous stretch in the joined text.
#[derive(Debug, Clone, Copy)]
pub struct RStart {
    /// Offset in joined text where this fragment starts.
    pub joined_off: u32,
    /// Reference sequence id.
    pub ref_id: u32,
    /// Offset within the reference sequence.
    pub ref_off: u32,
}

/// Loaded `.bt2` index, in-memory.
pub struct Bt2Index {
    pub params: EbwtParams,
    pub n_pat: u32,
    pub plen: Vec<u32>,
    pub n_frag: u32,
    pub rstarts: Vec<RStart>,
    /// Raw BWT + per-side checkpoint blocks. Length = `params.ebwt_tot_len`.
    pub ebwt: Vec<u8>,
    /// BWT row of the suffix beginning at text position 0 (the `$`).
    pub z_off: u32,
    /// Cumulative first-character counts: `[$, A, C, G, T]`.
    pub fchr: [u32; 5],
    /// First-level lookup: ftab[i] = first BWT row whose 10-mer prefix == i.
    pub ftab: Vec<u32>,
    /// Extended ftab for short patterns.
    pub eftab: Vec<u32>,
    /// Reference sequence names, in order.
    pub refnames: Vec<String>,
    /// Sampled suffix array (text position for every `off_stride`-th BWT row).
    pub offs: Vec<u32>,
}

impl Bt2Index {
    /// Open the index whose files are `<base>.1.bt2`, `<base>.2.bt2`, etc.
    pub fn open(base: impl AsRef<Path>) -> Result<Self> {
        let base = base.as_ref();
        let p1 = with_ext(base, "1.bt2");
        let p2 = with_ext(base, "2.bt2");

        let mut buf1 = Vec::new();
        std::fs::File::open(&p1)
            .with_context(|| format!("opening {}", p1.display()))?
            .read_to_end(&mut buf1)?;
        let mut buf2 = Vec::new();
        std::fs::File::open(&p2)
            .with_context(|| format!("opening {}", p2.display()))?
            .read_to_end(&mut buf2)?;

        Self::from_bytes(&buf1, &buf2)
    }

    /// Parse the index from in-memory bytes (test-friendly).
    pub fn from_bytes(file1: &[u8], file2: &[u8]) -> Result<Self> {
        let mut c = Cursor::new(file1);

        // ---- Header (.1.bt2) ----
        let endian = c.read_u32::<LittleEndian>()?;
        match endian {
            ENDIAN_NATIVE => {}
            ENDIAN_SWAPPED => bail!("byte-swapped index not supported in spike"),
            other => bail!("bad endianness sentinel: 0x{other:08x}"),
        }
        let len = c.read_u32::<LittleEndian>()?;
        let line_rate = c.read_u32::<LittleEndian>()?;
        let _lines_per_side = c.read_u32::<LittleEndian>()?; // deprecated, always 2
        let off_rate = c.read_u32::<LittleEndian>()?;
        let ftab_chars = c.read_u32::<LittleEndian>()?;
        let flags_raw = c.read_i32::<LittleEndian>()?;
        if flags_raw >= 0 {
            bail!("expected negative flags block, got {flags_raw}");
        }
        let params = EbwtParams::from_header(len, line_rate, off_rate, ftab_chars, flags_raw);
        if params.color {
            bail!("colorspace indexes not supported in spike");
        }

        // ---- nPat + plen ----
        let n_pat = c.read_u32::<LittleEndian>()?;
        let mut plen = Vec::with_capacity(n_pat as usize);
        for _ in 0..n_pat {
            plen.push(c.read_u32::<LittleEndian>()?);
        }

        // ---- nFrag + rstarts ----
        let n_frag = c.read_u32::<LittleEndian>()?;
        let mut rstarts = Vec::with_capacity(n_frag as usize);
        for _ in 0..n_frag {
            let joined_off = c.read_u32::<LittleEndian>()?;
            let ref_id = c.read_u32::<LittleEndian>()?;
            let ref_off = c.read_u32::<LittleEndian>()?;
            rstarts.push(RStart {
                joined_off,
                ref_id,
                ref_off,
            });
        }

        // ---- ebwt ----
        let mut ebwt = vec![0u8; params.ebwt_tot_len as usize];
        c.read_exact(&mut ebwt).context("reading ebwt block")?;

        // ---- zOff ----
        let z_off = c.read_u32::<LittleEndian>()?;
        if z_off >= params.bwt_len {
            bail!("zOff {z_off} >= bwt_len {}", params.bwt_len);
        }

        // ---- fchr (5 entries) ----
        let mut fchr = [0u32; 5];
        for slot in &mut fchr {
            *slot = c.read_u32::<LittleEndian>()?;
        }

        // ---- ftab ----
        let mut ftab = Vec::with_capacity(params.ftab_len as usize);
        for _ in 0..params.ftab_len {
            ftab.push(c.read_u32::<LittleEndian>()?);
        }

        // ---- eftab ----
        let mut eftab = Vec::with_capacity(params.eftab_len as usize);
        for _ in 0..params.eftab_len {
            eftab.push(c.read_u32::<LittleEndian>()?);
        }

        // ---- refnames ----
        let mut refnames = Vec::with_capacity(n_pat as usize);
        let mut current = String::new();
        loop {
            let mut byte = [0u8; 1];
            c.read_exact(&mut byte).context("reading refnames block")?;
            match byte[0] {
                0 => {
                    if !current.is_empty() {
                        refnames.push(std::mem::take(&mut current));
                    }
                    break;
                }
                b'\n' => {
                    refnames.push(std::mem::take(&mut current));
                }
                b => current.push(b as char),
            }
        }

        // ---- .2.bt2: 4-byte endian sentinel, then offs ----
        // (bt2_io.cpp:136-141: only the endian word is rewritten, not full header)
        let mut c2 = Cursor::new(file2);
        let endian2 = c2.read_u32::<LittleEndian>()?;
        if endian2 != ENDIAN_NATIVE {
            bail!("byte-swapped .2.bt2 not supported");
        }
        let mut offs = Vec::with_capacity(params.offs_len as usize);
        for _ in 0..params.offs_len {
            offs.push(c2.read_u32::<LittleEndian>()?);
        }

        Ok(Self {
            params,
            n_pat,
            plen,
            n_frag,
            rstarts,
            ebwt,
            z_off,
            fchr,
            ftab,
            eftab,
            refnames,
            offs,
        })
    }
}

fn with_ext(base: &Path, ext: &str) -> PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_base() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
    }

    #[test]
    fn opens_lambda_index() {
        let idx = Bt2Index::open(fixture_base()).expect("open");
        assert_eq!(idx.params.len, 48502, "lambda phage length");
        assert_eq!(idx.params.bwt_len, 48503);
        assert_eq!(idx.params.line_rate, 6);
        assert_eq!(idx.params.off_rate, 4);
        assert_eq!(idx.params.ftab_chars, 10);
        assert_eq!(idx.n_pat, 1);
        assert_eq!(idx.plen, vec![48502]);
        assert_eq!(idx.refnames.len(), 1);
        // BT2 stores the full FASTA description line, not just the accession.
        assert!(idx.refnames[0].starts_with("gi|9626243|ref|NC_001416.1|"));

        // fchr invariants: fchr[0] == 0, fchr[4] == len.
        // The $ symbol is encoded as 'A' (slot 0), so fchr is over A/C/G/T
        // with fchr[0]=0 (count of chars < A, which is just $) and
        // fchr[4]=len (= total non-$ chars).
        assert_eq!(idx.fchr[0], 0);
        assert_eq!(idx.fchr[4], idx.params.len);

        // offs has one entry per (1 << off_rate) BWT rows.
        assert_eq!(idx.offs.len(), idx.params.offs_len as usize);

        // zOff is a valid BWT row.
        assert!(idx.z_off < idx.params.bwt_len);
    }
}
