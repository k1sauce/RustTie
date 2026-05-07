//! Input/output: FASTQ in, SAM/BAM out. Thin wrappers around `noodles`.

pub mod bam;
pub mod fastq;
pub mod sam;

pub use bam::convert_sam_text_to_bam;
pub use fastq::{FastqReader, Read};
pub use sam::{ReadGroup, SamWriter};
