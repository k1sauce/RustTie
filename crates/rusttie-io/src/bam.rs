//! BAM output via `noodles-bam`. We get there by converting the SAM text
//! we already produce — a small extra pass that reuses all the SAM-emit
//! logic (FLAG bits, tags, paired fields, MAPQ, etc.) so SAM and BAM stay
//! in lockstep automatically.

use std::io::{BufRead, Write};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

/// Read SAM text from `sam_text` and write equivalent BAM to `out`.
/// `out` is wrapped in BGZF compression by `bam::io::Writer`.
pub fn convert_sam_text_to_bam<W: Write>(sam_text: &[u8], out: W) -> Result<()> {
    let mut sam_reader = sam::io::Reader::new(sam_text);
    let header: sam::Header = sam_reader
        .read_header()
        .context("parsing SAM header for BAM conversion")?;

    let mut bam_writer = bam::io::Writer::new(out);
    bam_writer
        .write_header(&header)
        .context("writing BAM header")?;

    // Read SAM records line-by-line and re-emit as BAM. We use the
    // line-based approach because the records are owned and easy to inspect
    // if there's a parse error.
    let mut line = String::new();
    let inner = sam_reader.get_mut();
    loop {
        line.clear();
        let n = inner.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with('@') || line.trim().is_empty() {
            continue;
        }
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        let raw = sam::Record::try_from(trimmed.as_bytes())
            .with_context(|| format!("parsing SAM line: {trimmed}"))?;
        let record_buf = sam::alignment::RecordBuf::try_from_alignment_record(&header, &raw)
            .context("converting SAM record to RecordBuf")?;
        bam_writer
            .write_alignment_record(&header, &record_buf)
            .context("writing BAM record")?;
    }
    bam_writer.try_finish()?;
    Ok(())
}
