//! Cross-validation against `bowtie2` itself: for each query, run BT2 and
//! compare the alignment position to RustTie's `exact_hits` output.
//!
//! BT2 may report a match with mismatches; RustTie reports exact matches.
//! For queries that occur exactly in the reference, both should agree on
//! at least one position.
//!
//! Skipped automatically if `bowtie2` is not on PATH (e.g. running outside
//! the devbox shell).

use std::path::PathBuf;
use std::process::Command;

use rusttie_index::{Bt2Index, exact_hits};

fn fixture_base() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus")
}

fn read_lambda_seq() -> Vec<u8> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures/lambda_virus.fa");
    let text = std::fs::read_to_string(&path).expect("read lambda fasta");
    let mut seq = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            continue;
        }
        seq.extend(line.trim().as_bytes());
    }
    seq
}

fn bowtie2_available() -> bool {
    Command::new("bowtie2")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Run bowtie2 with a single read; return the 0-based POS reported in the SAM.
/// `None` if BT2 reports unmapped.
fn bowtie2_align_one(query: &[u8]) -> Option<u32> {
    let idx = fixture_base();
    let qstr = std::str::from_utf8(query).expect("ascii query");
    let out = Command::new("bowtie2")
        .args(["-x", idx.to_str().unwrap(), "--no-head", "-c", qstr])
        .output()
        .expect("run bowtie2");
    assert!(
        out.status.success(),
        "bowtie2 failed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let line = stdout.lines().next().expect("bowtie2 SAM line");
    let fields: Vec<&str> = line.split('\t').collect();
    let flag: u32 = fields[1].parse().unwrap();
    if flag & 0x4 != 0 {
        return None; // unmapped
    }
    let pos_1based: u32 = fields[3].parse().unwrap();
    Some(pos_1based - 1)
}

#[test]
fn rusttie_hits_agree_with_bowtie2() {
    if !bowtie2_available() {
        eprintln!("skipping: bowtie2 not on PATH");
        return;
    }
    let idx = Bt2Index::open(fixture_base()).unwrap();
    let seq = read_lambda_seq();

    // 50-bp queries from across the genome.
    let positions: &[u32] = &[100, 5_000, 24_251, 40_000, 48_452];

    for &start in positions {
        let q = &seq[start as usize..start as usize + 50];
        let bt2_pos = bowtie2_align_one(q).expect("bt2 should align an exact match");
        let rusttie = exact_hits(&idx, q);
        let positions: Vec<u32> = rusttie.iter().map(|h| h.ref_off).collect();
        assert!(
            positions.contains(&bt2_pos),
            "BT2 position {bt2_pos} not in RustTie hits {positions:?} for query at {start}"
        );
    }
}
