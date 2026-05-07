//! BAM output round-trip: emit BAM, read it back via `samtools view`, and
//! verify the records match what we'd get from the SAM-text path directly.

use std::path::PathBuf;
use std::process::Command;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures")
}

fn rusttie_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rusttie"))
}

fn tool_available(name: &str) -> bool {
    Command::new("which")
        .arg(name)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[test]
fn bam_matches_sam_via_samtools() {
    if !tool_available("samtools") || !tool_available("wgsim") {
        eprintln!("skipping: samtools or wgsim not on PATH");
        return;
    }
    let tmp = std::env::temp_dir();
    let fq = tmp.join("rusttie_bamrt.fq");
    let r2 = tmp.join("rusttie_bamrt.r2.fq");
    if !fq.exists() {
        let ref_fa = fixture_dir().join("lambda_virus.fa");
        Command::new("wgsim")
            .args([
                "-e",
                "0.005",
                "-r",
                "0",
                "-R",
                "0",
                "-N",
                "200",
                "-1",
                "50",
                "-2",
                "50",
                "-S",
                "42",
                ref_fa.to_str().unwrap(),
                fq.to_str().unwrap(),
                r2.to_str().unwrap(),
            ])
            .output()
            .expect("wgsim");
    }

    let idx = fixture_dir().join("lambda_virus");
    let sam_path = tmp.join("rusttie_bamrt.sam");
    let bam_path = tmp.join("rusttie_bamrt.bam");

    let sam_run = Command::new(rusttie_bin())
        .args([
            "-x",
            idx.to_str().unwrap(),
            "-U",
            fq.to_str().unwrap(),
            "-S",
            sam_path.to_str().unwrap(),
        ])
        .output()
        .expect("rusttie sam");
    assert!(
        sam_run.status.success(),
        "rusttie sam: {}",
        String::from_utf8_lossy(&sam_run.stderr)
    );
    let bam_run = Command::new(rusttie_bin())
        .args([
            "-x",
            idx.to_str().unwrap(),
            "-U",
            fq.to_str().unwrap(),
            "-S",
            bam_path.to_str().unwrap(),
        ])
        .output()
        .expect("rusttie bam");
    assert!(
        bam_run.status.success(),
        "rusttie bam: {}",
        String::from_utf8_lossy(&bam_run.stderr)
    );

    // Convert BAM back to SAM with samtools.
    let st = Command::new("samtools")
        .args(["view", "-h", bam_path.to_str().unwrap()])
        .output()
        .expect("samtools view");
    assert!(
        st.status.success(),
        "samtools view: {}",
        String::from_utf8_lossy(&st.stderr)
    );
    let from_bam = String::from_utf8_lossy(&st.stdout).into_owned();
    let direct_sam = std::fs::read_to_string(&sam_path).unwrap();

    // Strip @HD and @SQ headers from both and compare alignment records.
    // (samtools may reorder header lines or add a @PG, which is fine.)
    fn records(text: &str) -> Vec<&str> {
        text.lines().filter(|l| !l.starts_with('@')).collect()
    }
    let a = records(&direct_sam);
    let b = records(&from_bam);
    assert_eq!(
        a.len(),
        b.len(),
        "record count differs: SAM-direct={} BAM->SAM={}",
        a.len(),
        b.len()
    );
    let mut diffs = 0;
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        if x != y && diffs < 3 {
            eprintln!("record {i} differs:\n  sam: {x}\n  bam: {y}");
            diffs += 1;
        }
    }
    if diffs > 0 {
        panic!("BAM round-trip diverged on {diffs} records");
    }
}
