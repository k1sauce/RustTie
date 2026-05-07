//! Phase 3d validation. Build lambda with `rusttie-build`, run the
//! `rusttie` aligner against our self-built index, and diff SAM against
//! the same aligner on a `bowtie2-build`-produced index. Identical SAM →
//! the index is correct from RustTie's reader's perspective.

use std::path::PathBuf;
use std::process::Command;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../validate/fixtures")
}

fn rusttie_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rusttie"))
}

fn rusttie_build_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rusttie-build"))
}

fn tool_available(name: &str) -> bool {
    Command::new("which")
        .arg(name)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[test]
fn rusttie_build_produces_aligner_compatible_index() {
    if !tool_available("wgsim") {
        eprintln!("skipping: wgsim not on PATH");
        return;
    }
    let tmp = std::env::temp_dir();

    let our_base = tmp.join("rusttie_self_lambda");
    let build_out = Command::new(rusttie_build_bin())
        .arg(fixture_dir().join("lambda_virus.fa"))
        .arg(&our_base)
        .output()
        .expect("rusttie-build");
    assert!(
        build_out.status.success(),
        "rusttie-build failed: {}",
        String::from_utf8_lossy(&build_out.stderr)
    );
    for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2"] {
        let p = format!("{}.{ext}", our_base.display());
        assert!(std::path::Path::new(&p).exists(), "missing {p}");
    }

    let fq = tmp.join("rusttie_self_reads.fq");
    let r2 = tmp.join("rusttie_self_reads.r2.fq");
    if !fq.exists() {
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
                fixture_dir().join("lambda_virus.fa").to_str().unwrap(),
                fq.to_str().unwrap(),
                r2.to_str().unwrap(),
            ])
            .output()
            .expect("wgsim");
    }

    let our_sam = tmp.join("rusttie_self.our.sam");
    let our_run = Command::new(rusttie_bin())
        .args([
            "-x",
            our_base.to_str().unwrap(),
            "-U",
            fq.to_str().unwrap(),
            "-S",
            our_sam.to_str().unwrap(),
            "--no-head",
        ])
        .output()
        .expect("rusttie on self-built");
    assert!(
        our_run.status.success(),
        "aligner on self-built index failed: {}",
        String::from_utf8_lossy(&our_run.stderr)
    );

    let bt_base = fixture_dir().join("lambda_virus");
    let bt_sam = tmp.join("rusttie_self.bt.sam");
    let bt_run = Command::new(rusttie_bin())
        .args([
            "-x",
            bt_base.to_str().unwrap(),
            "-U",
            fq.to_str().unwrap(),
            "-S",
            bt_sam.to_str().unwrap(),
            "--no-head",
        ])
        .output()
        .expect("rusttie on bt2-built");
    assert!(bt_run.status.success());

    let our = std::fs::read_to_string(&our_sam).unwrap();
    let bt = std::fs::read_to_string(&bt_sam).unwrap();
    if our != bt {
        for (i, (a, b)) in our.lines().zip(bt.lines()).enumerate() {
            if a != b {
                panic!(
                    "self-built index produces different SAM at line {i}:\n  ours: {a}\n  bt2 : {b}"
                );
            }
        }
        panic!(
            "SAM line count differs: {} vs {}",
            our.lines().count(),
            bt.lines().count()
        );
    }
}

/// Multi-contig + N-stretch alignment regression: extract a known
/// substring from EACH unambiguous stretch in `multi_n_long.fa` and verify
/// the aligner reports it back at the correct chromosome and position.
/// Stretches in the fixture are ≥30bp so a 22-bp seed always fits. Catches
/// bugs where the aligner only handles the first stretch (the chr22 failure
/// mode).
#[test]
fn rusttie_aligner_handles_all_stretches_in_multi_contig() {
    let tmp = std::env::temp_dir();
    let base = tmp.join("multi_stretch_probe");
    let build_out = Command::new(rusttie_build_bin())
        .arg(fixture_dir().join("multi_n_long.fa"))
        .arg(&base)
        .output()
        .expect("rusttie-build");
    assert!(build_out.status.success());

    // Hand-derived from validate/fixtures/multi_n_long.fa:
    //   chr1: 40 ACGT, 3 N, 35 chars, 4 N, 30 chars  (total 112)
    //   chr2: 2 N, 37 chars, 6 N, 36 chars            (total 81)
    //   chr3: 36 chars                                (total 36)
    let _probes: &[(&str, u32, &str)] = &[
        ("chr1", 0, "ACGTACGTACGTACGTACGTACGTACGT"), // stretch 0 (40bp)
        ("chr1", 43, "CGGAATAGCATGCATGCATGCATGCATGCATGCATGCATG"), // stretch 1 — wait
                                                     // Stretches need careful re-derivation. Instead of pre-computing,
                                                     // read the FASTA directly and slice each stretch.
    ];
    // Re-derive probes programmatically to avoid manual offset errors.
    let fa_text = std::fs::read_to_string(fixture_dir().join("multi_n_long.fa")).unwrap();
    let mut chrs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut cur_name = String::new();
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in fa_text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !cur_name.is_empty() {
                chrs.push((cur_name.clone(), std::mem::take(&mut cur_seq)));
            }
            cur_name = rest.split_whitespace().next().unwrap().to_string();
        } else {
            cur_seq.extend_from_slice(line.as_bytes());
        }
    }
    if !cur_name.is_empty() {
        chrs.push((cur_name, cur_seq));
    }
    // For each stretch (ACGT run >= 22 bp), emit a probe (chrom, start, seq).
    let mut probes: Vec<(String, u32, Vec<u8>)> = Vec::new();
    for (name, seq) in &chrs {
        let mut i = 0usize;
        while i < seq.len() {
            while i < seq.len() && !matches!(seq[i], b'A' | b'C' | b'G' | b'T') {
                i += 1;
            }
            let start = i;
            while i < seq.len() && matches!(seq[i], b'A' | b'C' | b'G' | b'T') {
                i += 1;
            }
            let len = i - start;
            if len >= 22 {
                probes.push((name.clone(), start as u32, seq[start..start + len].to_vec()));
            }
        }
    }
    assert!(
        probes.len() >= 4,
        "need at least 4 probes, got {}",
        probes.len()
    );

    let fq_path = tmp.join("multi_stretch_probe.fq");
    let mut fq = String::new();
    for (i, (chrom, pos, seq)) in probes.iter().enumerate() {
        fq.push_str(&format!(
            "@probe{i}_{chrom}_{pos}\n{}\n+\n{}\n",
            std::str::from_utf8(seq).unwrap(),
            "I".repeat(seq.len())
        ));
    }
    std::fs::write(&fq_path, &fq).unwrap();

    let sam_path = tmp.join("multi_stretch_probe.sam");
    let out = Command::new(rusttie_bin())
        .args([
            "-x",
            base.to_str().unwrap(),
            "-U",
            fq_path.to_str().unwrap(),
            "-S",
            sam_path.to_str().unwrap(),
            "--no-head",
        ])
        .output()
        .expect("rusttie");
    assert!(
        out.status.success(),
        "rusttie failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let sam = std::fs::read_to_string(&sam_path).unwrap();
    let mut sam_rows = std::collections::HashMap::new();
    for line in sam.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        sam_rows.insert(
            f[0].to_string(),
            (f[1].to_string(), f[2].to_string(), f[3].to_string()),
        );
    }

    let mut failed = Vec::new();
    for (i, (chrom, pos, _seq)) in probes.iter().enumerate() {
        let qname = format!("probe{i}_{chrom}_{pos}");
        let (flag, rname, sam_pos) = sam_rows.get(&qname).expect("missing record");
        let expected_pos_1based = (*pos + 1).to_string();
        if flag != "0" || rname != chrom || sam_pos != &expected_pos_1based {
            failed.push(format!(
                "{qname}: expected (flag=0, rname={chrom}, pos={expected_pos_1based}), got (flag={flag}, rname={rname}, pos={sam_pos})"
            ));
        }
    }
    if !failed.is_empty() {
        panic!("multi-stretch alignment failures:\n{}", failed.join("\n"));
    }
}

/// Byte-equivalence regression: every `.bt2` file produced by `rusttie-build`
/// must match `bowtie2-build`'s output exactly. Covers two fixtures: the
/// single-contig ACGT-only lambda phage and a multi-contig fixture with
/// embedded N runs at various positions.
#[test]
fn rusttie_build_byte_identical_to_bowtie2_build() {
    if !tool_available("bowtie2-build") {
        eprintln!("skipping: bowtie2-build not on PATH");
        return;
    }
    let tmp = std::env::temp_dir();

    let cases: &[(&str, &str)] = &[("lambda_virus.fa", "lambda"), ("multi_n.fa", "multi_n")];

    for (fa_name, label) in cases {
        let fa_path = fixture_dir().join(fa_name);
        let bt_base = tmp.join(format!("byte_eq_bt_{label}"));
        let our_base = tmp.join(format!("byte_eq_rt_{label}"));

        let bt = Command::new("bowtie2-build")
            .args([fa_path.to_str().unwrap(), bt_base.to_str().unwrap()])
            .output()
            .expect("bowtie2-build");
        assert!(bt.status.success(), "bowtie2-build failed for {label}");
        let our = Command::new(rusttie_build_bin())
            .args([fa_path.to_str().unwrap(), our_base.to_str().unwrap()])
            .output()
            .expect("rusttie-build");
        assert!(
            our.status.success(),
            "rusttie-build failed for {label}: {}",
            String::from_utf8_lossy(&our.stderr)
        );

        for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"] {
            let our_p = format!("{}.{ext}", our_base.display());
            let bt_p = format!("{}.{ext}", bt_base.display());
            let our_bytes = std::fs::read(&our_p).unwrap();
            let bt_bytes = std::fs::read(&bt_p).unwrap();
            assert_eq!(
                our_bytes,
                bt_bytes,
                "{label} .{ext} differs from bowtie2-build's output ({} vs {} bytes)",
                our_bytes.len(),
                bt_bytes.len(),
            );
        }
    }
}
