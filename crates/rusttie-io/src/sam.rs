//! Minimal SAM output. We hand-format records rather than going through
//! `noodles-sam` records to keep tight control over field/tag order, which
//! matters for byte-for-byte BT2 SAM matching later.
//!
//! For now we emit only the bare-minimum tags needed to be valid SAM.
//! Full BT2 tag fidelity (AS, XN, XM, XO, XG, NM, MD, YT, ...) lands in 2b.

use std::io::{BufWriter, Write};

use anyhow::Result;
use rusttie_index::Bt2Index;

/// SAM FLAG bits we use.
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_REVERSE: u16 = 0x10;

/// Optional `@RG` read-group declaration. `id` is required (becomes
/// `ID:<id>`); `extra_fields` are pre-formatted `KEY:VALUE` entries and
/// emitted in the order given. The same `id` is then stamped as `RG:Z:<id>`
/// on every alignment record, matching BT2's behavior.
#[derive(Debug, Clone, Default)]
pub struct ReadGroup {
    pub id: String,
    pub extra_fields: Vec<String>,
}

pub struct SamWriter<W: Write> {
    inner: BufWriter<W>,
    read_group: Option<ReadGroup>,
}

impl<W: Write> SamWriter<W> {
    pub fn new(w: W) -> Self {
        Self {
            inner: BufWriter::new(w),
            read_group: None,
        }
    }

    /// Attach a read group. Emits `@RG` in the header and appends
    /// `RG:Z:<id>` to every record's tags.
    pub fn set_read_group(&mut self, rg: ReadGroup) {
        self.read_group = Some(rg);
    }

    /// Emit the `@HD` + `@SQ` header. Programs (`@PG`) deferred to 2b.
    pub fn write_header(&mut self, idx: &Bt2Index) -> Result<()> {
        writeln!(self.inner, "@HD\tVN:1.0\tSO:unsorted")?;
        for (i, name) in idx.refnames.iter().enumerate() {
            // BT2 trims at first whitespace for SQ:SN to match standard tools.
            let sn = name.split_whitespace().next().unwrap_or(name);
            let ln = idx.plen[i];
            writeln!(self.inner, "@SQ\tSN:{sn}\tLN:{ln}")?;
        }
        if let Some(rg) = &self.read_group {
            write!(self.inner, "@RG\tID:{}", rg.id)?;
            for f in &rg.extra_fields {
                write!(self.inner, "\t{f}")?;
            }
            writeln!(self.inner)?;
        }
        Ok(())
    }

    fn write_rg_tag(&mut self) -> Result<()> {
        if let Some(rg) = &self.read_group {
            write!(self.inner, "\tRG:Z:{}", rg.id)?;
        }
        Ok(())
    }

    /// One alignment record. `tags` are pre-formatted (`"AS:i:0"`, `"MD:Z:50"`)
    /// and emitted in the supplied order — order matters for byte-for-byte
    /// SAM diffs against BT2.
    #[allow(clippy::too_many_arguments)]
    pub fn write_record(
        &mut self,
        qname: &[u8],
        flag: u16,
        rname: &str,
        pos_1based: u32,
        mapq: u8,
        cigar: &str,
        seq: &[u8],
        qual: &[u8],
        tags: &[String],
    ) -> Result<()> {
        let qname = std::str::from_utf8(qname)?;
        let seq = std::str::from_utf8(seq)?;
        let qual = if qual.is_empty() {
            "*".to_string()
        } else {
            std::str::from_utf8(qual)?.to_string()
        };
        write!(
            self.inner,
            "{qname}\t{flag}\t{rname}\t{pos_1based}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}"
        )?;
        for t in tags {
            write!(self.inner, "\t{t}")?;
        }
        self.write_rg_tag()?;
        writeln!(self.inner)?;
        Ok(())
    }

    /// Paired-end alignment record. Same shape as `write_record` but with
    /// caller-provided RNEXT, PNEXT, TLEN.
    #[allow(clippy::too_many_arguments)]
    pub fn write_paired_record(
        &mut self,
        qname: &[u8],
        flag: u16,
        rname: &str,
        pos_1based: u32,
        mapq: u8,
        cigar: &str,
        rnext: &str,
        pnext: u32,
        tlen: i64,
        seq: &[u8],
        qual: &[u8],
        tags: &[String],
    ) -> Result<()> {
        let qname = std::str::from_utf8(qname)?;
        let seq = std::str::from_utf8(seq)?;
        let qual = if qual.is_empty() {
            "*".to_string()
        } else {
            std::str::from_utf8(qual)?.to_string()
        };
        write!(
            self.inner,
            "{qname}\t{flag}\t{rname}\t{pos_1based}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
        )?;
        for t in tags {
            write!(self.inner, "\t{t}")?;
        }
        self.write_rg_tag()?;
        writeln!(self.inner)?;
        Ok(())
    }

    /// Unmapped record (FLAG=4, RNAME=*, POS=0, MAPQ=0, CIGAR=*).
    /// Emits BT2's standard YT:Z:UU tag for unpaired unmapped reads.
    pub fn write_unmapped(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        let qname = std::str::from_utf8(qname)?;
        let seq = std::str::from_utf8(seq)?;
        let qual = if qual.is_empty() {
            "*".to_string()
        } else {
            std::str::from_utf8(qual)?.to_string()
        };
        write!(
            self.inner,
            "{qname}\t{}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tYT:Z:UU",
            FLAG_UNMAPPED
        )?;
        self.write_rg_tag()?;
        writeln!(self.inner)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}
