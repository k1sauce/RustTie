//! Phase 3c paired-end validation. Generate ~500 paired reads via `wgsim`
//! at fragment length 200 ± 20, run both aligners, diff SAM.

use std::collections::HashMap;
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

#[derive(Debug, PartialEq, Eq, Clone)]
struct Row {
    flag: u16,
    rname: String,
    pos: u32,
    mapq: u8,
    cigar: String,
    rnext: String,
    pnext: u32,
    tlen: i64,
    tags: Vec<(String, String)>,
}

/// Key into the per-record map: (qname, R1 vs R2 from FLAG 0x40 bit).
type Key = (String, bool);

fn parse_sam(text: &str) -> HashMap<Key, Row> {
    let mut out = HashMap::new();
    for line in text.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let qname = f[0].to_string();
        let flag: u16 = f[1].parse().unwrap();
        let mut tags: Vec<(String, String)> = f[11..]
            .iter()
            .map(|t| {
                let (k, v) = t.split_once(':').expect("k:type:v");
                (k.to_string(), v.to_string())
            })
            .collect();
        tags.sort();
        let row = Row {
            flag,
            rname: f[2].to_string(),
            pos: f[3].parse().unwrap(),
            mapq: f[4].parse().unwrap(),
            cigar: f[5].to_string(),
            rnext: f[6].to_string(),
            pnext: f[7].parse().unwrap(),
            tlen: f[8].parse().unwrap(),
            tags,
        };
        let is_r1 = flag & 0x40 != 0;
        out.insert((qname, is_r1), row);
    }
    out
}

#[test]
fn paired_corpus_diff_summary() {
    if !tool_available("bowtie2") || !tool_available("wgsim") {
        eprintln!("skipping: bowtie2 or wgsim not on PATH");
        return;
    }
    let tmp = std::env::temp_dir();
    let r1 = tmp.join("rusttie_pair_r1.fq");
    let r2 = tmp.join("rusttie_pair_r2.fq");
    if !r1.exists() || !r2.exists() {
        let ref_fa = fixture_dir().join("lambda_virus.fa");
        let out = Command::new("wgsim")
            .args([
                "-e",
                "0.005",
                "-r",
                "0",
                "-R",
                "0",
                "-N",
                "500",
                "-1",
                "50",
                "-2",
                "50",
                "-d",
                "200", // mean fragment length
                "-s",
                "20", // stddev
                "-S",
                "42",
                ref_fa.to_str().unwrap(),
                r1.to_str().unwrap(),
                r2.to_str().unwrap(),
            ])
            .output()
            .expect("wgsim");
        assert!(out.status.success());
    }

    let idx = fixture_dir().join("lambda_virus");
    let rt_sam = tmp.join("rusttie_paired.rt.sam");
    let bt_sam = tmp.join("rusttie_paired.bt.sam");
    let rt = Command::new(rusttie_bin())
        .args([
            "-x",
            idx.to_str().unwrap(),
            "-1",
            r1.to_str().unwrap(),
            "-2",
            r2.to_str().unwrap(),
            "-S",
            rt_sam.to_str().unwrap(),
        ])
        .output()
        .expect("rusttie");
    assert!(
        rt.status.success(),
        "rusttie: {}",
        String::from_utf8_lossy(&rt.stderr)
    );
    let bt = Command::new("bowtie2")
        .args([
            "-x",
            idx.to_str().unwrap(),
            "-1",
            r1.to_str().unwrap(),
            "-2",
            r2.to_str().unwrap(),
            "-S",
            bt_sam.to_str().unwrap(),
        ])
        .output()
        .expect("bowtie2");
    assert!(
        bt.status.success(),
        "bowtie2: {}",
        String::from_utf8_lossy(&bt.stderr)
    );

    let rt_rows = parse_sam(&std::fs::read_to_string(&rt_sam).unwrap());
    let bt_rows = parse_sam(&std::fs::read_to_string(&bt_sam).unwrap());

    let total = rt_rows.len();
    let mut pos_match = 0;
    let mut flag_match = 0;
    let mut tlen_match = 0;
    let mut yt_match = 0;
    let mut tag_match = 0;
    let mut tag_diffs: HashMap<&'static str, u32> = HashMap::new();
    let mut sample_diffs: Vec<String> = Vec::new();

    let key_tags = ["AS", "NM", "MD", "XM", "XO", "XG", "YS"];

    for (k, rt_r) in &rt_rows {
        let Some(bt_r) = bt_rows.get(k) else { continue };
        if rt_r.pos == bt_r.pos && rt_r.rname == bt_r.rname {
            pos_match += 1;
        }
        if rt_r.flag == bt_r.flag {
            flag_match += 1;
        }
        if rt_r.tlen == bt_r.tlen {
            tlen_match += 1;
        }
        let yt_rt = rt_r
            .tags
            .iter()
            .find(|(k, _)| k == "YT")
            .map(|(_, v)| v.clone());
        let yt_bt = bt_r
            .tags
            .iter()
            .find(|(k, _)| k == "YT")
            .map(|(_, v)| v.clone());
        if yt_rt == yt_bt {
            yt_match += 1;
        }
        let core_tags_rt: Vec<_> = rt_r
            .tags
            .iter()
            .filter(|(k, _)| key_tags.contains(&k.as_str()))
            .cloned()
            .collect();
        let core_tags_bt: Vec<_> = bt_r
            .tags
            .iter()
            .filter(|(k, _)| key_tags.contains(&k.as_str()))
            .cloned()
            .collect();
        if core_tags_rt == core_tags_bt {
            tag_match += 1;
        } else {
            for k in key_tags {
                let r = rt_r.tags.iter().find(|(kk, _)| kk == k);
                let b = bt_r.tags.iter().find(|(kk, _)| kk == k);
                if r != b {
                    *tag_diffs.entry(k).or_insert(0) += 1;
                }
            }
        }
        if (rt_r != bt_r) && sample_diffs.len() < 3 {
            sample_diffs.push(format!(
                "{} (r1={}):\n  rusttie: flag={} pos={} cigar={} rnext={} pnext={} tlen={} yt={:?}\n  bowtie2: flag={} pos={} cigar={} rnext={} pnext={} tlen={} yt={:?}",
                k.0, k.1,
                rt_r.flag, rt_r.pos, rt_r.cigar, rt_r.rnext, rt_r.pnext, rt_r.tlen, yt_rt,
                bt_r.flag, bt_r.pos, bt_r.cigar, bt_r.rnext, bt_r.pnext, bt_r.tlen, yt_bt,
            ));
        }
    }

    eprintln!("=== paired-end diff ({total} records) ===");
    eprintln!("  position+rname agree: {pos_match}");
    eprintln!("  FLAG agree:           {flag_match}");
    eprintln!("  TLEN agree:           {tlen_match}");
    eprintln!("  YT agree:             {yt_match}");
    eprintln!("  core tags agree:      {tag_match}");
    eprintln!("  tag-diff counts:");
    for k in key_tags {
        let n = tag_diffs.get(k).copied().unwrap_or(0);
        eprintln!("    {k:>3}: {n}");
    }
    if !sample_diffs.is_empty() {
        eprintln!("  sample divergent records:");
        for d in &sample_diffs {
            eprintln!("{d}");
        }
    }

    // Phase 3c MVP gate: ≥98% on every metric. Some reads where RustTie's
    // seed-and-extend misses an alignment BT2 finds are expected — the same
    // edge case shows up as ~4 unmapped reads in the single-end corpus.
    assert!(
        pos_match * 100 / total >= 98,
        "position agreement: {pos_match}/{total}"
    );
    assert!(
        flag_match * 100 / total >= 98,
        "FLAG agreement: {flag_match}/{total}"
    );
    assert!(
        tlen_match * 100 / total >= 98,
        "TLEN agreement: {tlen_match}/{total}"
    );
    assert!(
        yt_match * 100 / total >= 98,
        "YT agreement: {yt_match}/{total}"
    );
    assert!(
        tag_match * 100 / total >= 98,
        "core tags agreement: {tag_match}/{total}"
    );
}
