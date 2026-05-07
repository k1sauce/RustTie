//! FASTQ input via `noodles-fastq`. Auto-detects gzipped input by `.gz`
//! extension and wraps the file in a `MultiGzDecoder` (handles concatenated
//! gzip members, which is common in real sequencing data).

use std::fs::File;
use std::io::{BufReader, Read as IoRead};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles_fastq as fastq;

/// One FASTQ record, owned. We make our own simple struct rather than
/// exposing `noodles_fastq::Record` to keep the public API minimal.
#[derive(Debug, Clone)]
pub struct Read {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// Streaming FASTQ reader. Boxed inner reader so the same type works for
/// plain and gzipped inputs.
pub struct FastqReader<R: IoRead> {
    inner: fastq::io::Reader<BufReader<R>>,
    buf: fastq::Record,
    /// Human-readable source identifier (file path or `"<reader>"`) used to
    /// add context to parse errors so users see which file is malformed.
    source: String,
}

impl FastqReader<Box<dyn IoRead + Send>> {
    /// Open a FASTQ file, auto-detecting gzip via `.gz` extension.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let is_gz = path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
        let r: Box<dyn IoRead + Send> = if is_gz {
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        Ok(Self::from_reader_with_source(r, path.display().to_string()))
    }
}

impl<R: IoRead> FastqReader<R> {
    pub fn from_reader(r: R) -> Self {
        Self::from_reader_with_source(r, "<reader>".to_string())
    }

    fn from_reader_with_source(r: R, source: String) -> Self {
        Self {
            inner: fastq::io::Reader::new(BufReader::new(r)),
            buf: fastq::Record::default(),
            source,
        }
    }

    /// Returns the next read, or `Ok(None)` at EOF.
    pub fn next_read(&mut self) -> Result<Option<Read>> {
        let n = self
            .inner
            .read_record(&mut self.buf)
            .with_context(|| format!("parsing FASTQ record from {}", self.source))?;
        if n == 0 {
            return Ok(None);
        }
        Ok(Some(Read {
            name: self.buf.name().to_vec(),
            seq: self.buf.sequence().to_vec(),
            qual: self.buf.quality_scores().to_vec(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    #[test]
    fn parses_minimal_fastq() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nGGCC\n+\n!!!!\n";
        let mut r = FastqReader::from_reader(Cursor::new(&data[..]));
        let r1 = r.next_read().unwrap().unwrap();
        assert_eq!(r1.name, b"read1");
        assert_eq!(r1.seq, b"ACGT");
        assert_eq!(r1.qual, b"IIII");
        let r2 = r.next_read().unwrap().unwrap();
        assert_eq!(r2.name, b"read2");
        assert_eq!(r2.seq, b"GGCC");
        assert!(r.next_read().unwrap().is_none());
    }

    #[test]
    fn parses_gzipped_fastq() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nGGCC\n+\n!!!!\n";
        let mut gz = Vec::new();
        let mut enc = flate2::write::GzEncoder::new(&mut gz, flate2::Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap();

        // Write to a temp file with .gz extension and reopen via FastqReader::open.
        let path = std::env::temp_dir().join("rusttie_io_gz_test.fq.gz");
        std::fs::write(&path, &gz).unwrap();
        let mut r = FastqReader::open(&path).unwrap();
        let r1 = r.next_read().unwrap().unwrap();
        assert_eq!(r1.name, b"read1");
        assert_eq!(r1.seq, b"ACGT");
        let r2 = r.next_read().unwrap().unwrap();
        assert_eq!(r2.name, b"read2");
        assert!(r.next_read().unwrap().is_none());
    }
}
