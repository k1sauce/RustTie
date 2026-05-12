//! `rusttie-build`: construct a `.bt2` index from a FASTA reference.
//!
//! Uses [`sais_rs`] for SA construction — output is byte-for-byte
//! compatible with libsais (which is what BowTie 2 itself uses), so the
//! per-partition ordering of suffixes matches BT2's convention without any
//! rotation tricks.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use byteorder::{LittleEndian, WriteBytesExt};
use sais_rs::suffix_array;

use crate::format::{ENDIAN_NATIVE, FLAG_ENTIRE_REV, FLAG_VALID, OFF_SIZE, ascii_to_2bit};

/// Default geometry parameters (`bowtie2-build-s` defaults).
pub const DEFAULT_LINE_RATE: u32 = 6; // → side_sz = 64
pub const DEFAULT_OFF_RATE: u32 = 4; // SA sample every 16 rows
pub const DEFAULT_FTAB_CHARS: u32 = 10;

/// Sentinel byte appended before SA construction. BT2 uses `(char)127`
/// (`vendor/bowtie2/blockwise_sa.h:204`) — strictly larger than any
/// ACGT byte, so the sentinel-only suffix sorts last.
const SENTINEL: u8 = 127;

/// Build a forward `.bt2` index for `fasta_path`, writing files prefixed
/// by `out_basename` (e.g., `out_basename.1.bt2`).
pub fn build_index(fasta_path: &Path, out_basename: &Path) -> Result<()> {
    let (refnames, refseqs) = read_fasta(fasta_path)?;
    if refnames.is_empty() {
        bail!("FASTA contains no sequences");
    }
    if refseqs.iter().any(|s| s.is_empty()) {
        bail!("FASTA sequence is empty");
    }

    // Walk each FASTA sequence to identify unambiguous (ACGT) stretches.
    // For each stretch we emit a `.3.bt2` record `(off, len, first)`:
    //   `off` = number of preceding ambiguous bases (Ns) in this same
    //   stretch group (since the previous emit, or since the start of the
    //   reference if `first`),
    //   `len` = number of unambiguous bases,
    //   `first` = 1 iff this is the first record of a new reference.
    // Joined text concatenates only the unambiguous stretches.
    // `rstarts` tracks the (joined_off, ref_id, ref_off) of each stretch.
    let plen: Vec<u32> = refseqs.iter().map(|s| s.len() as u32).collect();
    let mut joined = Vec::new();
    let mut rstarts: Vec<(u32, u32, u32)> = Vec::new();
    let mut records3: Vec<(u32, u32, u8)> = Vec::new();
    for (ref_id, seq) in refseqs.iter().enumerate() {
        let mut first_in_ref = true;
        let mut i = 0usize;
        let n = seq.len();
        while i < n {
            // Count leading Ns since the previous emit / start of ref.
            let mut n_count: u32 = 0;
            while i < n && ascii_to_2bit(seq[i]).is_none() {
                n_count += 1;
                i += 1;
            }
            let stretch_ref_off = i as u32;
            let stretch_joined_off = joined.len() as u32;
            let mut stretch_len: u32 = 0;
            while i < n && ascii_to_2bit(seq[i]).is_some() {
                joined.push(seq[i]);
                stretch_len += 1;
                i += 1;
            }
            // Emit a record:
            //  - whenever there's an ACGT run, OR
            //  - when this is the first record of a new reference (preserves
            //    boundaries even for all-N references), OR
            //  - when this run had leading Ns but no following ACGT (i.e.,
            //    trailing Ns at end of reference) — BT2 emits these as
            //    `RefRecord(off, 0, false)` to record the ambiguous tail.
            if stretch_len > 0 || first_in_ref || n_count > 0 {
                records3.push((n_count, stretch_len, first_in_ref as u8));
                if stretch_len > 0 {
                    rstarts.push((stretch_joined_off, ref_id as u32, stretch_ref_off));
                }
                first_in_ref = false;
            }
        }
    }
    if joined.is_empty() {
        bail!("FASTA contains no unambiguous (ACGT) bases");
    }

    // Forward index: build over the joined text as-is.
    let fwd = build_pass(&joined)?;
    let p1 = with_ext(out_basename, "1.bt2");
    write_primary(&p1, &fwd, &refnames, &plen, &rstarts, FLAG_VALID)?;
    let p2 = with_ext(out_basename, "2.bt2");
    write_secondary(&p2, &fwd)?;
    // .3.bt2 + .4.bt2 are over the original text and don't change for the
    // reverse index — emit once.
    let p3 = with_ext(out_basename, "3.bt2");
    write_records3(&p3, &records3)?;
    let p4 = with_ext(out_basename, "4.bt2");
    write_packed_ref(&p4, &joined)?;

    // Reverse index: build over the *entire* concatenated text reversed
    // (BT2's "EBWT_ENTIRE_REV" mode). Stretches appear in reverse order in
    // the reversed text, so rstarts are recomputed: keep `(ref_id, ref_off)`
    // pairs in reverse order, accumulate `joined_off` from the reversed
    // stretch lengths.
    let mut joined_rev = joined.clone();
    joined_rev.reverse();
    let rev = build_pass(&joined_rev)?;
    let rstarts_rev: Vec<(u32, u32, u32)> = {
        let lens: Vec<u32> = rstarts
            .windows(2)
            .map(|w| w[1].0 - w[0].0)
            .chain(std::iter::once(
                joined.len() as u32 - rstarts.last().unwrap().0,
            ))
            .collect();
        let mut out = Vec::with_capacity(rstarts.len());
        let mut acc: u32 = 0;
        for i in (0..rstarts.len()).rev() {
            out.push((acc, rstarts[i].1, rstarts[i].2));
            acc += lens[i];
        }
        out
    };
    let p_rev_1 = with_ext(out_basename, "rev.1.bt2");
    write_primary(
        &p_rev_1,
        &rev,
        &refnames,
        &plen,
        &rstarts_rev,
        FLAG_VALID | FLAG_ENTIRE_REV,
    )?;
    let p_rev_2 = with_ext(out_basename, "rev.2.bt2");
    write_secondary(&p_rev_2, &rev)?;

    Ok(())
}

/// All in-memory data structures derived from one text (used twice: once
/// for the forward index, once for the reversed text).
struct BuildPass {
    geo: Geometry,
    sa: Vec<usize>,
    ebwt: Vec<u8>,
    z_off: u32,
    fchr: [u32; 5],
    ftab: Vec<u32>,
    eftab: Vec<u32>,
}

fn build_pass(text: &[u8]) -> Result<BuildPass> {
    let len = text.len() as u32;
    let bwt_len = len + 1;

    let mut sa_input = text.to_vec();
    sa_input.push(SENTINEL);
    let sa_i32 = suffix_array(&sa_input).context("sais-rs suffix_array")?;
    debug_assert_eq!(sa_i32.len(), bwt_len as usize);
    debug_assert_eq!(
        sa_i32[(bwt_len - 1) as usize] as u32,
        len,
        "sentinel-only suffix should be the largest row"
    );
    let sa: Vec<usize> = sa_i32.iter().map(|&v| v as usize).collect();

    let z_off = sa
        .iter()
        .position(|&p| p == 0)
        .ok_or_else(|| anyhow!("z_off not found in SA"))? as u32;

    let mut bwt_codes: Vec<u8> = Vec::with_capacity(bwt_len as usize);
    #[allow(clippy::needless_range_loop)] // explicit indexing reads cleaner here
    for i in 0..bwt_len as usize {
        let s = sa[i];
        let c = if s == 0 {
            0u8
        } else if s == len as usize {
            ascii_to_2bit(text[(len - 1) as usize]).expect("ACGT pre-checked")
        } else {
            ascii_to_2bit(text[s - 1]).expect("ACGT pre-checked")
        };
        bwt_codes.push(c);
    }

    let geo = Geometry::compute(len, DEFAULT_LINE_RATE, DEFAULT_OFF_RATE, DEFAULT_FTAB_CHARS);
    let ebwt = pack_ebwt(&bwt_codes, z_off, &geo);

    let mut counts = [0u32; 4];
    for (i, &c) in bwt_codes.iter().enumerate() {
        if i as u32 == z_off {
            continue;
        }
        counts[c as usize] += 1;
    }
    let mut fchr = [0u32; 5];
    for c in 0..4 {
        fchr[c + 1] = fchr[c] + counts[c];
    }
    debug_assert_eq!(fchr[4], len);

    let (ftab, eftab) = build_ftab_eftab(
        text,
        &sa,
        geo.ftab_chars as usize,
        geo.ftab_len as usize,
        geo.eftab_len as usize,
    );

    Ok(BuildPass {
        geo,
        sa,
        ebwt,
        z_off,
        fchr,
        ftab,
        eftab,
    })
}

#[derive(Debug, Clone)]
struct Geometry {
    len: u32,
    bwt_len: u32,
    line_rate: u32,
    #[allow(dead_code)]
    line_sz: u32,
    off_rate: u32,
    ftab_chars: u32,
    ftab_len: u32,
    eftab_len: u32,
    side_sz: u32,
    side_bwt_sz: u32,
    side_bwt_len: u32,
    num_sides: u32,
    ebwt_tot_len: u32,
}

impl Geometry {
    fn compute(len: u32, line_rate: u32, off_rate: u32, ftab_chars: u32) -> Self {
        let bwt_len = len + 1;
        let line_sz = 1u32 << line_rate;
        let ftab_len = (1u32 << (ftab_chars * 2)) + 1;
        let eftab_len = ftab_chars * 2;
        let side_sz = line_sz;
        let side_bwt_sz = side_sz - 4 * OFF_SIZE as u32;
        let side_bwt_len = side_bwt_sz * 4;
        let num_sides = bwt_len.div_ceil(side_bwt_len);
        let ebwt_tot_len = num_sides * side_sz;
        Self {
            len,
            bwt_len,
            line_rate,
            line_sz,
            off_rate,
            ftab_chars,
            ftab_len,
            eftab_len,
            side_sz,
            side_bwt_sz,
            side_bwt_len,
            num_sides,
            ebwt_tot_len,
        }
    }
}

/// Pack `bwt_codes` (length `bwt_len`) into the side-based ebwt blob.
/// Each side is `side_sz` bytes: first `side_bwt_sz` bytes hold up to
/// `side_bwt_len = side_bwt_sz*4` chars (LSB-first within each byte),
/// followed by 4 × OFF_SIZE bytes of cumulative A/C/G/T counts up to that
/// side's start. The encoded `$` at `z_off` is NOT counted in the
/// checkpoints (matching `bt2_idx.h:2955-2963`).
fn pack_ebwt(bwt_codes: &[u8], z_off: u32, geo: &Geometry) -> Vec<u8> {
    let mut out = vec![0u8; geo.ebwt_tot_len as usize];
    let mut cum = [0u32; 4]; // cumulative A/C/G/T from the start of the BWT
    for s in 0..geo.num_sides as usize {
        let side_start_byte = s * geo.side_sz as usize;
        // First write the checkpoint for this side: cumulative counts up to
        // (but not including) the chars in this side.
        let occ_off = side_start_byte + geo.side_bwt_sz as usize;
        for c in 0..4 {
            let bytes = cum[c].to_le_bytes();
            out[occ_off + c * OFF_SIZE..occ_off + (c + 1) * OFF_SIZE].copy_from_slice(&bytes);
        }
        // Pack BWT chars into this side, LSB-first within each byte.
        let chars_in_side =
            ((s + 1) as u32 * geo.side_bwt_len).min(geo.bwt_len) - s as u32 * geo.side_bwt_len;
        for k in 0..chars_in_side as usize {
            let bwt_idx = s * geo.side_bwt_len as usize + k;
            let code = bwt_codes[bwt_idx];
            let byte_off = side_start_byte + (k / 4);
            let bp = (k % 4) as u8;
            out[byte_off] |= code << (bp * 2);
            // Update cumulative count, except for the encoded $.
            if bwt_idx as u32 != z_off {
                cum[code as usize] += 1;
            }
        }
    }
    out
}

/// Build ftab + eftab matching BT2's exact encoding (`bt2_idx.h:2970-3160`).
///
/// `ftab[i]` is the first BWT row whose first `ftab_chars` chars equal the
/// integer `i` (interpreted base-4, MSB-first). For 10-mer patterns where
/// "short" suffixes (length < `ftab_chars`) sort immediately after the
/// pattern's last long-suffix row, BT2 stores the row-range in `eftab` and
/// puts an "encoded pointer" into `ftab[i]` (the value `eftab_idx ^ OFF_MASK`,
/// i.e., `!eftab_idx`). Read-side helpers `ftabLo` / `ftabHi` decode either
/// case based on whether `ftab[i] <= len`.
///
/// Algorithm (two passes over SA in row order):
///
/// Pass 1: walk SA, accumulate per-10-mer counts in `ftab[i+1]` and per-
/// transition absorbed-short counts in `absorb_ftab[i]`.
///
/// Pass 2: running prefix sum that resolves to either a row index or an
/// eftab-pointer, with `eftab` populated for absorbed transitions.
fn build_ftab_eftab(
    text: &[u8],
    sa: &[usize],
    ftab_chars: usize,
    ftab_len: usize,
    eftab_len: usize,
) -> (Vec<u32>, Vec<u32>) {
    let n = text.len();
    let off_mask: u32 = u32::MAX;

    // Pass 1: counts + absorption tracking.
    let mut ftab = vec![0u32; ftab_len];
    let mut absorb_ftab = vec![0u32; ftab_len];
    let mut absorb_cnt: u32 = 0;
    for &s in sa.iter() {
        if s + ftab_chars <= n {
            // Long suffix: decode 10-mer key.
            let mut suf_int: u32 = 0;
            for k in 0..ftab_chars {
                suf_int =
                    (suf_int << 2) | ascii_to_2bit(text[s + k]).expect("ACGT pre-checked") as u32;
            }
            ftab[suf_int as usize + 1] += 1;
            if absorb_cnt > 0 {
                absorb_ftab[suf_int as usize] = absorb_cnt;
                absorb_cnt = 0;
            }
        } else {
            // Short suffix (length < ftab_chars). Will be absorbed at the
            // next long-suffix transition.
            absorb_cnt += 1;
        }
    }
    if absorb_cnt > 0 {
        absorb_ftab[ftab_len - 1] = absorb_cnt;
    }

    // Pass 2: running prefix sum + eftab. Reads `ftab[i-1]` AFTER it's been
    // finalized in the previous iteration (so it's a row index or eftab
    // pointer, not the raw count). `ftab_hi` decodes either case.
    let mut eftab = vec![0u32; eftab_len];
    let mut eftab_cur: u32 = 0;

    let ftab_hi = |ftab: &[u32], eftab: &[u32], i: u32| -> u32 {
        let v = ftab[i as usize];
        if v <= n as u32 {
            v
        } else {
            let ef_idx = (v ^ off_mask) as usize;
            eftab[ef_idx * 2 + 1]
        }
    };

    for i in 1..ftab_len {
        let lo = ftab[i] + ftab_hi(&ftab, &eftab, (i - 1) as u32);
        if absorb_ftab[i] > 0 {
            let hi = lo + absorb_ftab[i];
            eftab[(eftab_cur * 2) as usize] = lo;
            eftab[(eftab_cur * 2 + 1) as usize] = hi;
            ftab[i] = eftab_cur ^ off_mask;
            eftab_cur += 1;
        } else {
            ftab[i] = lo;
        }
    }

    (ftab, eftab)
}

fn read_fasta(path: &Path) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    // We parse FASTA manually rather than going through noodles-fasta's
    // `Definition` because BT2's `.1.bt2` refnames region stores the full
    // header line BYTE-FOR-BYTE (including any double-spaces between name
    // and description). noodles' `Definition::Display` normalizes whitespace,
    // breaking byte-equivalence vs `bowtie2-build` on FASTAs like UCSC's
    // hg38.chrM.fa whose header is `>chrM  AC:J01415.2  ...` (double space).
    use std::io::{BufRead, BufReader};
    let file =
        std::fs::File::open(path).with_context(|| format!("opening fasta {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.context("reading fasta line")?;
        if let Some(stripped) = line.strip_prefix('>') {
            if !names.is_empty() {
                seqs.push(std::mem::take(&mut cur_seq));
            }
            names.push(stripped.to_string());
        } else {
            cur_seq.extend(line.trim_end().bytes());
        }
    }
    if !names.is_empty() {
        seqs.push(cur_seq);
    }
    Ok((names, seqs))
}

fn write_primary(
    path: &Path,
    pass: &BuildPass,
    refnames: &[String],
    plen: &[u32],
    rstarts: &[(u32, u32, u32)],
    flags: i32,
) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);
    let geo = &pass.geo;

    // Header.
    w.write_u32::<LittleEndian>(ENDIAN_NATIVE)?;
    w.write_u32::<LittleEndian>(geo.len)?;
    w.write_u32::<LittleEndian>(geo.line_rate)?;
    w.write_u32::<LittleEndian>(2)?; // linesPerSide (deprecated, always 2)
    w.write_u32::<LittleEndian>(geo.off_rate)?;
    w.write_u32::<LittleEndian>(geo.ftab_chars)?;
    // BT2 writes `-flags` (negative-tagged validity block).
    w.write_i32::<LittleEndian>(-flags)?;

    // nPat + plen[nPat].
    w.write_u32::<LittleEndian>(refnames.len() as u32)?;
    for &l in plen {
        w.write_u32::<LittleEndian>(l)?;
    }
    // nFrag + rstarts.
    w.write_u32::<LittleEndian>(rstarts.len() as u32)?;
    for r in rstarts {
        w.write_u32::<LittleEndian>(r.0)?; // joined_off
        w.write_u32::<LittleEndian>(r.1)?; // ref_id
        w.write_u32::<LittleEndian>(r.2)?; // ref_off
    }
    w.write_all(&pass.ebwt)?;
    w.write_u32::<LittleEndian>(pass.z_off)?;
    for v in &pass.fchr {
        w.write_u32::<LittleEndian>(*v)?;
    }
    for v in &pass.ftab {
        w.write_u32::<LittleEndian>(*v)?;
    }
    for v in &pass.eftab {
        w.write_u32::<LittleEndian>(*v)?;
    }
    // refnames: BT2 emits `name1\nname2\n...\nlast\n\0` — newline after
    // every name plus a final null terminator.
    for name in refnames {
        w.write_all(name.as_bytes())?;
        w.write_all(b"\n")?;
    }
    w.write_all(&[0u8])?;
    w.flush()?;
    Ok(())
}

fn write_secondary(path: &Path, pass: &BuildPass) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);
    w.write_u32::<LittleEndian>(ENDIAN_NATIVE)?;
    let stride = 1u32 << pass.geo.off_rate;
    // Sample SA at every `stride`-th BWT row.
    for i in (0..pass.geo.bwt_len).step_by(stride as usize) {
        w.write_u32::<LittleEndian>(pass.sa[i as usize] as u32)?;
    }
    w.flush()?;
    Ok(())
}

fn write_records3(path: &Path, records: &[(u32, u32, u8)]) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);
    w.write_u32::<LittleEndian>(ENDIAN_NATIVE)?;
    w.write_u32::<LittleEndian>(records.len() as u32)?;
    for r in records {
        w.write_u32::<LittleEndian>(r.0)?;
        w.write_u32::<LittleEndian>(r.1)?;
        w.write_all(&[r.2])?;
    }
    w.flush()?;
    Ok(())
}

fn write_packed_ref(path: &Path, joined: &[u8]) -> Result<()> {
    // 2-bit pack, LSB-first within each byte (matches `.4.bt2` empirical
    // packing — see `docs/bt2-format.md` correction note).
    let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);
    let mut buf = vec![0u8; joined.len().div_ceil(4)];
    for (i, &b) in joined.iter().enumerate() {
        let code = ascii_to_2bit(b).expect("ACGT");
        let byte = i / 4;
        let bp = (i % 4) as u8;
        buf[byte] |= code << (bp * 2);
    }
    w.write_all(&buf)?;
    w.flush()?;
    Ok(())
}

fn with_ext(base: &Path, ext: &str) -> std::path::PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    #![allow(dead_code)]
    use super::*;

    // Diagnostic tests removed: see module doc-comment for the SA-ordering
    // mismatch with BT2. Phase 3d is on hold until we swap the SA backend.

    /// Verify rust-bio's `suffix_array` matches the canonical answer on a
    /// 5-character text.
    #[test]
    fn rust_bio_sa_canonical() {
        let s = b"ACGT$";
        let sa = bio::data_structures::suffix_array::suffix_array(s);
        // Sorted suffixes: $, ACGT$, CGT$, GT$, T$
        // Positions:        4, 0,     1,    2,   3
        assert_eq!(sa, vec![4, 0, 1, 2, 3]);
    }

    /// Direct BWT computation for `ACGT$`. SA = [4, 0, 1, 2, 3].
    /// BWT[i] = T[SA[i] - 1] cyclic.
    /// SA[0] = 4 → T[3] = T
    /// SA[1] = 0 → T[-1 mod 5] = T[4] = $ (encoded as A, this is z_off)
    /// SA[2] = 1 → T[0] = A
    /// SA[3] = 2 → T[1] = C
    /// SA[4] = 3 → T[2] = G
    #[test]
    fn bwt_layout_for_acgt() {
        let seq = b"ACGT";
        let mut sa_input = seq.to_vec();
        sa_input.push(b'$');
        let sa = bio::data_structures::suffix_array::suffix_array(&sa_input);
        assert_eq!(sa, vec![4, 0, 1, 2, 3]);

        let mut bwt = Vec::new();
        let mut z = 0u32;
        for (i, &p) in sa.iter().enumerate() {
            if p == 0 {
                bwt.push(0u8);
                z = i as u32;
            } else {
                bwt.push(ascii_to_2bit(seq[p - 1]).unwrap());
            }
        }
        // Expected encoded: [T=3, A=0(z), A=0, C=1, G=2]
        assert_eq!(bwt, vec![3, 0, 0, 1, 2]);
        assert_eq!(z, 1);
    }
}
