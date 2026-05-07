use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use rusttie_index::build::build_index;

#[derive(Parser, Debug)]
#[command(
    name = "rusttie-build",
    version,
    about = "Build a .bt2 index from a FASTA"
)]
struct Cli {
    /// Path to FASTA reference.
    fasta: PathBuf,
    /// Index basename. Files written: `<base>.{1,2,3,4}.bt2`.
    base: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    build_index(&cli.fasta, &cli.base)?;
    eprintln!(
        "Built index files: {b}.{{1,2,3,4,rev.1,rev.2}}.bt2",
        b = cli.base.display()
    );
    Ok(())
}
