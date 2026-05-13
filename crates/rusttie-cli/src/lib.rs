//! CLI for `rusttie`. Phase 2a: subset of BT2 flags for single-end exact-match.

use std::fs::File;
use std::io::{BufWriter, Write, stdout};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;
use rayon::prelude::*;
use rusttie_align::{
    AlignResult, Alignment, DESCENT_D_DEFAULT, DESCENT_R_DEFAULT, FRAG_LEN_MAX, FRAG_LEN_MIN,
    PER_SEED_HIT_CAP, PairOutcome, PairType, Scoring, Strand, align_read_with_descent,
    PairCandidate, classify_pair_set, mapq_v2, mate_rescue, reverse_complement,
};
use rusttie_index::{BitPairReference, Bt2Index};
use rusttie_io::{FastqReader, Read, ReadGroup, SamWriter, convert_sam_text_to_bam, sam};

/// Subset of BT2 flags. We mirror BT2's short-flag spelling for drop-in
/// compatibility on the supported subset.
#[derive(Parser, Debug)]
#[command(name = "rusttie", version, about = "Rust port of BowTie 2 (early)")]
pub struct Cli {
    /// Index basename (the prefix of `<base>.1.bt2` etc.).
    #[arg(short = 'x', long = "index", value_name = "BASE")]
    pub index: PathBuf,

    /// Single-end FASTQ input. Mutually exclusive with `-1`/`-2`.
    #[arg(short = 'U', long = "unpaired", value_name = "FILE", conflicts_with_all = ["m1", "m2"])]
    pub unpaired: Option<PathBuf>,

    /// Paired-end mate 1 FASTQ. Used together with `-2`.
    #[arg(short = '1', long = "m1", value_name = "FILE", requires = "m2")]
    pub m1: Option<PathBuf>,

    /// Paired-end mate 2 FASTQ.
    #[arg(short = '2', long = "m2", value_name = "FILE", requires = "m1")]
    pub m2: Option<PathBuf>,

    /// SAM output. Stdout if omitted.
    #[arg(short = 'S', long = "sam", value_name = "FILE")]
    pub sam: Option<PathBuf>,

    /// Suppress SAM header lines (useful for diffing).
    #[arg(long = "no-head")]
    pub no_head: bool,

    /// Number of parallel alignment threads. 0 = use all cores (rayon default).
    /// Mirrors BT2's `-p`.
    #[arg(short = 'p', long = "threads", default_value_t = 0)]
    pub threads: usize,

    /// Suppress SAM records for unaligned reads. Mirrors BT2's `--no-unal`.
    #[arg(long = "no-unal")]
    pub no_unal: bool,

    /// Mismatch penalty bounds, "MX,MN" — Q≥40 → MX, Q=0 → MN. BT2 default
    /// `--mp 6,2`.
    #[arg(long = "mp", value_name = "MX,MN")]
    pub mp: Option<String>,

    /// Read gap open/extend penalty, "OPEN,EXT". BT2 default `--rdg 5,3`.
    #[arg(long = "rdg", value_name = "OPEN,EXT")]
    pub rdg: Option<String>,

    /// Reference gap open/extend penalty, "OPEN,EXT". BT2 default `--rfg 5,3`.
    #[arg(long = "rfg", value_name = "OPEN,EXT")]
    pub rfg: Option<String>,

    /// Score-min function, "L,A,B" → minimum score = A + B * read_len.
    /// BT2 default `--score-min L,-0.6,-0.6`. Only `L` (linear) supported.
    #[arg(long = "score-min", value_name = "FUNC,A,B")]
    pub score_min: Option<String>,

    /// Sensitivity preset (BT2 has these as no-ops at the moment because
    /// our defaults already match BT2's `--sensitive`; flagged here so users
    /// can pass them on the command line without errors). Future iterations
    /// will plumb these through to `-D`/`-R`/`-L`/`-i`.
    #[arg(long = "very-fast", conflicts_with_all = ["fast", "sensitive", "very_sensitive"])]
    pub very_fast: bool,
    #[arg(long = "fast", conflicts_with_all = ["sensitive", "very_sensitive"])]
    pub fast: bool,
    #[arg(long = "sensitive", conflicts_with = "very_sensitive")]
    pub sensitive: bool,
    #[arg(long = "very-sensitive")]
    pub very_sensitive: bool,

    /// Per-seed hit cap. Seeds with more than this many SA-range hits are
    /// skipped as too repetitive; this is the main perf knob. See
    /// `align::PER_SEED_HIT_CAP` for the perf/recall curve.
    #[arg(long = "seed-hit-cap", default_value_t = PER_SEED_HIT_CAP)]
    pub seed_hit_cap: u32,

    /// BT2 `-D`: stop after this many consecutive seed-extension failures.
    /// An extension "fails" if it doesn't yield a new best or new secbest
    /// alignment. Default 15.
    #[arg(short = 'D', long = "descent-budget", default_value_t = DESCENT_D_DEFAULT)]
    pub descent_budget: u32,

    /// BT2 `-R`: maximum re-seedings when the seed set is repetitive
    /// (avg hits per aligned seed > 300). Each re-seeding shifts seed
    /// offsets to expose new candidates. Default 2.
    #[arg(short = 'R', long = "descent-reseed", default_value_t = DESCENT_R_DEFAULT)]
    pub descent_reseed: u32,

    /// Mate-rescue: from each side's top-K anchor alignments, run SW for
    /// the *other* mate within the FR concordance window. Mirrors BT2's
    /// `extendSeedsPaired` mate-find step. Set to 0 to disable mate-rescue.
    /// Default 5.
    #[arg(long = "mate-rescue", default_value_t = MATE_RESCUE_TOP_K_DEFAULT)]
    pub mate_rescue_top_k: u32,

    /// Read-group ID (`ID:` field of `@RG` and `RG:Z:` on every record).
    /// Required when using `--rg`. BT2-compatible.
    #[arg(long = "rg-id", value_name = "ID")]
    pub rg_id: Option<String>,

    /// Additional `@RG` field, in `KEY:VALUE` form (e.g. `--rg SM:NA12878
    /// --rg LB:lib1 --rg PL:ILLUMINA`). Repeatable. Requires `--rg-id`.
    #[arg(long = "rg", value_name = "KEY:VALUE", requires = "rg_id")]
    pub rg: Vec<String>,

    /// BT2 `-k <int>`: report up to N alignments per read (or per pair, for
    /// paired-end). Primary keeps standard FLAGs; secondaries set FLAG
    /// 0x100. Default 1. Conflicts with `-a`.
    #[arg(short = 'k', long = "k", value_name = "N", conflicts_with = "all")]
    pub k: Option<u32>,

    /// BT2 `-a` / `--all`: report all valid alignments. Equivalent to a very
    /// large `-k`.
    #[arg(short = 'a', long = "all")]
    pub all: bool,
}

fn parse_pair_i32(s: &str, name: &str) -> Result<(i32, i32)> {
    let (a, b) = s
        .split_once(',')
        .ok_or_else(|| anyhow::anyhow!("--{name} expects \"A,B\""))?;
    Ok((a.parse()?, b.parse()?))
}

/// Build a [`ReadGroup`] from `--rg-id` + repeated `--rg KEY:VALUE`. Each
/// `--rg` entry is validated to look like `XX:value` (BT2 / SAM-spec form).
/// Returns `Ok(None)` if `--rg-id` was not given.
fn build_read_group(cli: &Cli) -> Result<Option<ReadGroup>> {
    let Some(id) = cli.rg_id.clone() else {
        if !cli.rg.is_empty() {
            bail!("--rg requires --rg-id");
        }
        return Ok(None);
    };
    let mut extra = Vec::with_capacity(cli.rg.len());
    for entry in &cli.rg {
        let (k, v) = entry
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("--rg expects KEY:VALUE form, got {entry:?}"))?;
        if k.is_empty() || v.is_empty() {
            bail!("--rg KEY and VALUE must both be non-empty: {entry:?}");
        }
        if k.eq_ignore_ascii_case("ID") {
            bail!("--rg ID:... conflicts with --rg-id; specify ID via --rg-id only");
        }
        extra.push(entry.clone());
    }
    Ok(Some(ReadGroup {
        id,
        extra_fields: extra,
    }))
}

fn parse_score_min(s: &str) -> Result<(f64, f64)> {
    // "L,A,B" → (A, B). Only the linear form is implemented; reject others.
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        bail!("--score-min expects \"FUNC,A,B\"");
    }
    if !parts[0].eq_ignore_ascii_case("L") {
        bail!("--score-min: only the L (linear) function is implemented");
    }
    Ok((parts[1].parse()?, parts[2].parse()?))
}

fn build_scoring(cli: &Cli) -> Result<Scoring> {
    let mut sc = Scoring::default();
    if let Some(s) = &cli.mp {
        let (mx, mn) = parse_pair_i32(s, "mp")?;
        sc.mp_max = mx;
        sc.mp_min = mn;
    }
    if let Some(s) = &cli.rdg {
        let (o, e) = parse_pair_i32(s, "rdg")?;
        sc.rdg_open = o;
        sc.rdg_extend = e;
    }
    if let Some(s) = &cli.rfg {
        let (o, e) = parse_pair_i32(s, "rfg")?;
        sc.rfg_open = o;
        sc.rfg_extend = e;
    }
    if let Some(s) = &cli.score_min {
        let (a, b) = parse_score_min(s)?;
        sc.score_min_const = a;
        sc.score_min_coeff = b;
    }
    Ok(sc)
}

/// Reads processed per parallel batch. Bigger = better throughput, more
/// peak memory; 4096 reads × ~200B per read ≈ <1MB so this is fine.
const BATCH_SIZE: usize = 4096;

/// Run the aligner with parsed CLI args.
pub fn run(cli: Cli) -> Result<()> {
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok(); // ignore error if already initialized (e.g., in tests)
    }
    let idx = Bt2Index::open(&cli.index)
        .with_context(|| format!("opening index {}", cli.index.display()))?;
    let refs = BitPairReference::open(&cli.index)
        .with_context(|| format!("opening reference (.3/.4.bt2) for {}", cli.index.display()))?;
    let scoring = build_scoring(&cli)?;
    let no_unal = cli.no_unal;
    let seed_hit_cap = cli.seed_hit_cap;
    let descent_budget = cli.descent_budget;
    let descent_reseed = cli.descent_reseed;
    let mate_rescue_top_k = cli.mate_rescue_top_k;
    let read_group = build_read_group(&cli)?;
    // -k <int> / -a: how many alignments to report per read (or per pair).
    // -a wins over -k (clap rejects both) and means "all valid".
    let report_k: u32 = if cli.all {
        u32::MAX
    } else {
        cli.k.unwrap_or(1)
    };
    if report_k == 0 {
        bail!("-k must be at least 1");
    }

    // BAM output is auto-detected by the `.bam` extension on `-S`.
    let bam_mode = cli
        .sam
        .as_ref()
        .is_some_and(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("bam")));
    if bam_mode {
        // Two-pass: emit SAM to memory, then convert to BAM. Header is
        // mandatory in BAM, so --no-head is silently ignored.
        let mut sam_buf = Vec::<u8>::new();
        match (&cli.unpaired, &cli.m1, &cli.m2) {
            (Some(u), None, None) => {
                let mut reader = FastqReader::open(u)
                    .with_context(|| format!("opening fastq {}", u.display()))?;
                run_unpaired(
                    &idx,
                    &refs,
                    &scoring,
                    no_unal,
                    seed_hit_cap,
                    descent_budget,
                    descent_reseed,
                    report_k,
                    read_group.clone(),
                    &mut reader,
                    &mut sam_buf,
                    false,
                )?;
            }
            (None, Some(p1), Some(p2)) => {
                let mut r1 =
                    FastqReader::open(p1).with_context(|| format!("opening {}", p1.display()))?;
                let mut r2 =
                    FastqReader::open(p2).with_context(|| format!("opening {}", p2.display()))?;
                run_paired(
                    &idx,
                    &refs,
                    &scoring,
                    no_unal,
                    seed_hit_cap,
                    descent_budget,
                    descent_reseed,
                    mate_rescue_top_k,
                    report_k,
                    read_group.clone(),
                    &mut r1,
                    &mut r2,
                    &mut sam_buf,
                    false,
                )?;
            }
            _ => bail!("provide -U <reads.fq>  OR  -1 <mate1.fq> -2 <mate2.fq>"),
        }
        // bam_mode = true implies cli.sam.is_some() (we matched on its .bam extension).
        let path = cli.sam.as_ref().expect("bam_mode requires -S <file.bam>");
        let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        convert_sam_text_to_bam(&sam_buf, BufWriter::new(f))
            .with_context(|| format!("writing BAM to {}", path.display()))?;
        return Ok(());
    }

    let writer: Box<dyn Write> = match &cli.sam {
        Some(path) => {
            let f = File::create(path).with_context(|| format!("creating {}", path.display()))?;
            Box::new(BufWriter::new(f))
        }
        None => Box::new(BufWriter::new(stdout().lock())),
    };

    let result = match (&cli.unpaired, &cli.m1, &cli.m2) {
        (Some(u), None, None) => {
            let mut reader =
                FastqReader::open(u).with_context(|| format!("opening fastq {}", u.display()))?;
            run_unpaired(
                &idx,
                &refs,
                &scoring,
                no_unal,
                seed_hit_cap,
                descent_budget,
                descent_reseed,
                report_k,
                read_group,
                &mut reader,
                writer,
                cli.no_head,
            )
        }
        (None, Some(p1), Some(p2)) => {
            let mut r1 =
                FastqReader::open(p1).with_context(|| format!("opening {}", p1.display()))?;
            let mut r2 =
                FastqReader::open(p2).with_context(|| format!("opening {}", p2.display()))?;
            run_paired(
                &idx,
                &refs,
                &scoring,
                no_unal,
                seed_hit_cap,
                descent_budget,
                descent_reseed,
                mate_rescue_top_k,
                report_k,
                read_group,
                &mut r1,
                &mut r2,
                writer,
                cli.no_head,
            )
        }
        _ => bail!("provide -U <reads.fq>  OR  -1 <mate1.fq> -2 <mate2.fq>"),
    };
    if std::env::var_os("RUSTTIE_PROFILE").is_some() {
        rusttie_align::profile::print();
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn run_unpaired<R: std::io::Read, W: Write>(
    idx: &Bt2Index,
    refs: &BitPairReference,
    scoring: &Scoring,
    no_unal: bool,
    seed_hit_cap: u32,
    descent_budget: u32,
    descent_reseed: u32,
    report_k: u32,
    read_group: Option<ReadGroup>,
    reader: &mut FastqReader<R>,
    writer: W,
    no_head: bool,
) -> Result<()> {
    let mut sam_w = SamWriter::new(writer);
    if let Some(rg) = read_group {
        sam_w.set_read_group(rg);
    }
    if !no_head {
        sam_w.write_header(idx)?;
    }

    // Producer / parallel-map / serial-writer pipeline. Reads are pulled
    // into batches of `BATCH_SIZE`, aligned in parallel via rayon, then
    // emitted in original order so the SAM output is deterministic and
    // byte-equivalent to the single-threaded path.
    let mut batch: Vec<Read> = Vec::with_capacity(BATCH_SIZE);
    loop {
        batch.clear();
        for _ in 0..BATCH_SIZE {
            match reader.next_read()? {
                Some(r) => batch.push(r),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        let results: Vec<Option<AlignResult>> = batch
            .par_iter()
            .map(|r| {
                align_read_with_descent(
                    idx,
                    refs,
                    &r.seq,
                    &r.qual,
                    scoring,
                    seed_hit_cap,
                    descent_budget,
                    descent_reseed,
                )
            })
            .collect();
        for (read, result) in batch.iter().zip(results.iter()) {
            emit_sam(
                &mut sam_w,
                idx,
                read,
                result.as_ref(),
                scoring,
                no_unal,
                report_k,
            )?;
        }
    }
    sam_w.flush()?;
    Ok(())
}

/// Default cap on how many top-scoring anchor alignments per side to
/// mate-rescue from. BT2's descent driver mate-rescues from every
/// seed-extended anchor; we cap because most reads have very few
/// alignments and the long tail (highly repetitive reads with many tied
/// anchors) inflates work without changing MAPQ.
///
/// Empirical chr22 sweep (10k pairs, `-p 8`, warm cache):
///
/// | K | Wall | MAPQ agree | Notes |
/// |---|---|---|---|
/// | 0 | 0.59s | 90.0% | mate-rescue disabled |
/// | 1 | 0.47s | 90.5% | best-anchor only |
/// | 3 | 0.81s | 93.8% | knee — captures tied alternates |
/// | 5 | 1.11s | 93.8% | no further MAPQ gain |
/// | 50 | 2.58s | 93.8% | wasted work past K=3 |
pub const MATE_RESCUE_TOP_K_DEFAULT: u32 = 3;

/// Run mate-rescue from each side's top-K anchor alignments and merge
/// any newly-found alignments into the per-mate sets, deduped by
/// `(ref_id, ref_off, strand)`. Mirrors BT2's `extendSeedsPaired` mate-find
/// step (`vendor/bowtie2/aligner_sw_driver.cpp:2226-2347`): finds mate
/// alignments that wouldn't be reachable via pure FM-index seed search
/// (e.g., when the other mate's seeds all hit the cap).
///
/// Returns the explicit list of `(anchor, rescued)` pair candidates produced
/// — these are BT2's `rs1_`/`rs2_` parallel entries for the cases where the
/// paired-mode aligner emitted a pair via mate-find (see
/// `vendor/bowtie2/aln_sink.cpp:1413`). Each rescue is also pushed into the
/// per-mate `r1_alns`/`r2_alns` lists so the unpaired/discordant fallbacks
/// downstream still see them, but the pair-candidate identity is preserved
/// here so `classify_pair_set` can compute pair-secbest without
/// Cartesian-ing.
#[allow(clippy::too_many_arguments)]
fn augment_via_mate_rescue(
    refs: &BitPairReference,
    r1_seq: &[u8],
    r1_qual: &[u8],
    r2_seq: &[u8],
    r2_qual: &[u8],
    scoring: &Scoring,
    top_k: u32,
    r1_alns: &mut Vec<Alignment>,
    r2_alns: &mut Vec<Alignment>,
) -> Vec<PairCandidate> {
    let mut pair_candidates: Vec<PairCandidate> = Vec::new();
    if top_k == 0 || (r1_alns.is_empty() && r2_alns.is_empty()) {
        return pair_candidates;
    }
    let top_k_u = top_k as usize;

    let r1_rc = reverse_complement(r1_seq);
    let r2_rc = reverse_complement(r2_seq);
    let r1_qual_rev: Vec<u8> = r1_qual.iter().rev().copied().collect();
    let r2_qual_rev: Vec<u8> = r2_qual.iter().rev().copied().collect();

    // Pre-build dedup keys for what we already have so rescued alignments
    // that duplicate seed-found ones are dropped.
    use std::collections::HashSet;
    let mut r1_keys: HashSet<(u32, u32, Strand)> = r1_alns
        .iter()
        .map(|a| (a.ref_id, a.ref_off, a.strand))
        .collect();
    let mut r2_keys: HashSet<(u32, u32, Strand)> = r2_alns
        .iter()
        .map(|a| (a.ref_id, a.ref_off, a.strand))
        .collect();

    // Anchor in r1 → rescue r2.
    let r1_anchors: Vec<Alignment> = r1_alns.iter().take(top_k_u).cloned().collect();
    for anchor in &r1_anchors {
        if let Some(rescued) = mate_rescue(
            refs,
            anchor,
            r2_seq,
            r2_qual,
            &r2_rc,
            &r2_qual_rev,
            scoring,
            FRAG_LEN_MIN,
            FRAG_LEN_MAX,
        ) {
            if let Some(cand) = PairCandidate::try_new(anchor.clone(), rescued.clone()) {
                pair_candidates.push(cand);
            }
            let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
            if r2_keys.insert(key) {
                r2_alns.push(rescued);
            }
        }
    }

    // Anchor in r2 → rescue r1.
    let r2_anchors: Vec<Alignment> = r2_alns.iter().take(top_k_u).cloned().collect();
    for anchor in &r2_anchors {
        if let Some(rescued) = mate_rescue(
            refs,
            anchor,
            r1_seq,
            r1_qual,
            &r1_rc,
            &r1_qual_rev,
            scoring,
            FRAG_LEN_MIN,
            FRAG_LEN_MAX,
        ) {
            if let Some(cand) = PairCandidate::try_new(rescued.clone(), anchor.clone()) {
                pair_candidates.push(cand);
            }
            let key = (rescued.ref_id, rescued.ref_off, rescued.strand);
            if r1_keys.insert(key) {
                r1_alns.push(rescued);
            }
        }
    }

    // Re-sort each side score-descending so classify_pair_set's fallback
    // Cartesian (which assumes sorted input) stays consistent.
    r1_alns.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.ref_id.cmp(&b.ref_id))
            .then(a.ref_off.cmp(&b.ref_off))
    });
    r2_alns.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.ref_id.cmp(&b.ref_id))
            .then(a.ref_off.cmp(&b.ref_off))
    });

    pair_candidates
}

/// Strip the conventional pair suffix (`/1`, `/2`) so both mates share QNAME
/// per the SAM spec — matches BT2's behavior.
fn strip_pair_suffix(name: &[u8]) -> &[u8] {
    let n = name.len();
    if n >= 2 && name[n - 2] == b'/' && (name[n - 1] == b'1' || name[n - 1] == b'2') {
        &name[..n - 2]
    } else {
        name
    }
}

/// Paired-end alignment pipeline. Mate FQs are read in lockstep; each pair
/// is aligned independently then classified for concordance. Two SAM
/// records are emitted per pair.
#[allow(clippy::too_many_arguments)]
fn run_paired<R: std::io::Read, W: Write>(
    idx: &Bt2Index,
    refs: &BitPairReference,
    scoring: &Scoring,
    no_unal: bool,
    seed_hit_cap: u32,
    descent_budget: u32,
    descent_reseed: u32,
    mate_rescue_top_k: u32,
    report_k: u32,
    read_group: Option<ReadGroup>,
    r1_reader: &mut FastqReader<R>,
    r2_reader: &mut FastqReader<R>,
    writer: W,
    no_head: bool,
) -> Result<()> {
    let mut sam_w = SamWriter::new(writer);
    if let Some(rg) = read_group {
        sam_w.set_read_group(rg);
    }
    if !no_head {
        sam_w.write_header(idx)?;
    }

    let mut batch1: Vec<Read> = Vec::with_capacity(BATCH_SIZE);
    let mut batch2: Vec<Read> = Vec::with_capacity(BATCH_SIZE);
    loop {
        batch1.clear();
        batch2.clear();
        for _ in 0..BATCH_SIZE {
            match (r1_reader.next_read()?, r2_reader.next_read()?) {
                (Some(mut a), Some(mut b)) => {
                    a.name = strip_pair_suffix(&a.name).to_vec();
                    b.name = strip_pair_suffix(&b.name).to_vec();
                    batch1.push(a);
                    batch2.push(b);
                }
                (None, None) => break,
                _ => bail!("mate FASTQs are different lengths"),
            }
        }
        if batch1.is_empty() {
            break;
        }
        let outcomes: Vec<PairOutcome> = batch1
            .par_iter()
            .zip(batch2.par_iter())
            .map(|(r1, r2)| {
                let res1 = align_read_with_descent(
                    idx,
                    refs,
                    &r1.seq,
                    &r1.qual,
                    scoring,
                    seed_hit_cap,
                    descent_budget,
                    descent_reseed,
                );
                let res2 = align_read_with_descent(
                    idx,
                    refs,
                    &r2.seq,
                    &r2.qual,
                    scoring,
                    seed_hit_cap,
                    descent_budget,
                    descent_reseed,
                );
                let r1_sec = res1.as_ref().and_then(|r| r.secbest_score);
                let r2_sec = res2.as_ref().and_then(|r| r.secbest_score);
                let r1_best = res1.as_ref().map(|r| r.best.clone());
                let r2_best = res2.as_ref().map(|r| r.best.clone());
                let mut r1_alns: Vec<Alignment> = res1.map(|r| r.all).unwrap_or_default();
                let mut r2_alns: Vec<Alignment> = res2.map(|r| r.all).unwrap_or_default();
                // Pair-candidate pool: each entry is one (r1, r2) tuple
                // produced TOGETHER by the paired-mode aligner (matches
                // BT2's `rs1_`/`rs2_` parallel lists). Two sources:
                // (1) primary candidate from independent seed-and-extend
                //     when both mates' best alignments form a concordant
                //     pair; (2) each successful mate-rescue.
                let mut pair_pool: Vec<PairCandidate> = Vec::new();
                if let (Some(b1), Some(b2)) = (r1_best, r2_best)
                    && let Some(cand) = PairCandidate::try_new(b1, b2)
                {
                    pair_pool.push(cand);
                }
                let rescue_cands = augment_via_mate_rescue(
                    refs,
                    &r1.seq,
                    &r1.qual,
                    &r2.seq,
                    &r2.qual,
                    scoring,
                    mate_rescue_top_k,
                    &mut r1_alns,
                    &mut r2_alns,
                );
                pair_pool.extend(rescue_cands);
                // Dedup by (r1 pos+strand, r2 pos+strand): different
                // anchor paths can rediscover the same pair.
                let mut seen = std::collections::HashSet::new();
                pair_pool.retain(|c| {
                    seen.insert((
                        c.r1.ref_id,
                        c.r1.ref_off,
                        c.r1.strand,
                        c.r2.ref_off,
                        c.r2.strand,
                    ))
                });
                classify_pair_set(&pair_pool, &r1_alns, &r2_alns, r1_sec, r2_sec, report_k)
            })
            .collect();
        for ((r1, r2), outcome) in batch1.iter().zip(batch2.iter()).zip(outcomes) {
            emit_pair(&mut sam_w, idx, r1, r2, &outcome, scoring, no_unal)?;
        }
    }
    sam_w.flush()?;
    Ok(())
}

fn emit_sam<W: Write>(
    sam_w: &mut SamWriter<W>,
    idx: &Bt2Index,
    read: &Read,
    result: Option<&AlignResult>,
    scoring: &Scoring,
    no_unal: bool,
    report_k: u32,
) -> Result<()> {
    let Some(result) = result else {
        if !no_unal {
            sam_w.write_unmapped(&read.name, &read.seq, &read.qual)?;
        }
        return Ok(());
    };
    // Primary MAPQ uses the descent driver's per-mate secbest. All reported
    // alignments share this MAPQ since BT2 does the same; secondaries
    // distinguish themselves via FLAG 0x100, not via MAPQ.
    let primary = &result.best;
    let mapq = mapq_v2(
        primary.score,
        result.secbest_score,
        scoring.score_min(primary.read_len),
    );
    let want = (report_k as usize).min(result.all.len()).max(1);
    for (i, a) in result.all.iter().take(want).enumerate() {
        let is_secondary = i > 0;
        let rname = idx.refnames[a.ref_id as usize]
            .split_whitespace()
            .next()
            .unwrap_or(&idx.refnames[a.ref_id as usize])
            .to_string();
        let pos_1based = a.ref_off + 1;
        let cigar = a.cigar.clone();
        let mut flag = match a.strand {
            Strand::Forward => 0u16,
            Strand::Reverse => sam::FLAG_REVERSE,
        };
        if is_secondary {
            flag |= 0x100;
        }
        // Secondary alignments per SAM spec may omit SEQ/QUAL (use *), but
        // BT2 emits them; we mirror BT2 so downstream tools that index by
        // qname find consistent records.
        let (seq, qual) = match a.strand {
            Strand::Forward => (read.seq.clone(), read.qual.clone()),
            Strand::Reverse => (
                reverse_complement(&read.seq),
                read.qual.iter().rev().copied().collect::<Vec<_>>(),
            ),
        };
        let nm = a.mismatches + a.gap_extends;
        let tags = [
            format!("AS:i:{}", a.score),
            "XN:i:0".to_string(),
            format!("XM:i:{}", a.mismatches),
            format!("XO:i:{}", a.gap_opens),
            format!("XG:i:{}", a.gap_extends),
            format!("NM:i:{}", nm),
            format!("MD:Z:{}", a.md),
            "YT:Z:UU".to_string(),
        ];
        sam_w.write_record(
            &read.name, flag, &rname, pos_1based, mapq, &cigar, &seq, &qual, &tags,
        )?;
    }
    Ok(())
}

/// Emit one read of a pair. `is_r1` selects the FLAG 0x40/0x80 bit. `mate`
/// is the other read's alignment, or None if unmapped. `pair_type` drives
/// FLAG 0x2 (concordant) and `YT` tag.
#[allow(clippy::too_many_arguments)]
fn emit_pair<W: Write>(
    sam_w: &mut SamWriter<W>,
    idx: &Bt2Index,
    r1: &Read,
    r2: &Read,
    outcome: &PairOutcome,
    scoring: &Scoring,
    no_unal: bool,
) -> Result<()> {
    let p_type = outcome.pair_type;
    // Concordant pair: compute pair MAPQ once from the sum of mate scores
    // and the second-best concordant *pair* score (BowtieMapq2's paired
    // formula), then apply to both mates. Otherwise fall back to per-mate
    // single-end MAPQ. The pair MAPQ is shared by all reported pairs
    // (primary + secondaries) since they all rank against the same
    // alignment set; secondaries are distinguished by FLAG 0x100.
    let pair_mapq = if p_type == PairType::Concordant
        && let (Some(a1), Some(a2)) = (outcome.r1.as_ref(), outcome.r2.as_ref())
    {
        let pair_score = a1.score + a2.score;
        let pair_smin = scoring.score_min(a1.read_len) + scoring.score_min(a2.read_len);
        Some(mapq_v2(
            pair_score,
            outcome.concordant_pair_secbest,
            pair_smin,
        ))
    } else {
        None
    };

    // Primary first.
    emit_one_of_pair(
        sam_w,
        idx,
        r1,
        true,
        outcome.r1.as_ref(),
        outcome.r2.as_ref(),
        outcome.r1_secbest,
        pair_mapq,
        p_type,
        outcome.frag_len,
        scoring,
        no_unal,
        false,
    )?;
    emit_one_of_pair(
        sam_w,
        idx,
        r2,
        false,
        outcome.r2.as_ref(),
        outcome.r1.as_ref(),
        outcome.r2_secbest,
        pair_mapq,
        p_type,
        outcome.frag_len,
        scoring,
        no_unal,
        false,
    )?;

    // Then each additional concordant pair as secondary alignments
    // (FLAG 0x100). `--no-unal` doesn't apply because both mates are mapped
    // by definition for an additional concordant pair.
    for (a1, a2, frag_len) in &outcome.additional_concordant {
        emit_one_of_pair(
            sam_w,
            idx,
            r1,
            true,
            Some(a1),
            Some(a2),
            outcome.r1_secbest,
            pair_mapq,
            PairType::Concordant,
            *frag_len,
            scoring,
            false,
            true,
        )?;
        emit_one_of_pair(
            sam_w,
            idx,
            r2,
            false,
            Some(a2),
            Some(a1),
            outcome.r2_secbest,
            pair_mapq,
            PairType::Concordant,
            *frag_len,
            scoring,
            false,
            true,
        )?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_one_of_pair<W: Write>(
    sam_w: &mut SamWriter<W>,
    idx: &Bt2Index,
    read: &Read,
    is_r1: bool,
    self_aln: Option<&Alignment>,
    mate_aln: Option<&Alignment>,
    self_secbest: Option<i32>,
    // Pre-computed pair-level MAPQ. Some only when this is a concordant
    // pair; both mates of the pair share the value.
    pair_mapq: Option<u8>,
    pair_type: PairType,
    frag_len: u32,
    scoring: &Scoring,
    no_unal: bool,
    is_secondary: bool,
) -> Result<()> {
    // --no-unal: drop the record entirely if this read is unmapped. Note
    // BT2 still emits the *mapped* mate when the other half is unmapped;
    // we mirror that here.
    if self_aln.is_none() && no_unal {
        return Ok(());
    }
    // Construct FLAG.
    let mut flag: u16 = 0x1; // paired
    flag |= if is_r1 { 0x40 } else { 0x80 };
    if pair_type == PairType::Concordant {
        flag |= 0x2;
    }
    if is_secondary {
        flag |= 0x100;
    }
    match self_aln {
        None => flag |= 0x4, // unmapped
        Some(a) if a.strand == Strand::Reverse => flag |= 0x10,
        Some(_) => {}
    }
    match mate_aln {
        None => flag |= 0x8, // mate unmapped
        Some(a) if a.strand == Strand::Reverse => flag |= 0x20,
        Some(_) => {}
    }

    let yt = match (pair_type, self_aln.is_some(), mate_aln.is_some()) {
        (PairType::Concordant, _, _) => "CP",
        (PairType::Discordant, _, _) => "DP",
        (PairType::Unpaired, _, _) => "UP",
    };

    if let Some(a) = self_aln {
        // Mapped half.
        let rname = idx.refnames[a.ref_id as usize]
            .split_whitespace()
            .next()
            .unwrap_or(&idx.refnames[a.ref_id as usize])
            .to_string();
        let pos_1based = a.ref_off + 1;
        let cigar = a.cigar.clone();
        let (seq, qual) = if a.strand == Strand::Reverse {
            (
                reverse_complement(&read.seq),
                read.qual.iter().rev().copied().collect::<Vec<_>>(),
            )
        } else {
            (read.seq.clone(), read.qual.clone())
        };
        // RNEXT / PNEXT / TLEN.
        let (rnext, pnext, tlen): (String, u32, i64) = match (mate_aln, pair_type) {
            (Some(m), PairType::Concordant) => {
                let mate_rname = if m.ref_id == a.ref_id {
                    "=".to_string()
                } else {
                    idx.refnames[m.ref_id as usize]
                        .split_whitespace()
                        .next()
                        .unwrap_or(&idx.refnames[m.ref_id as usize])
                        .to_string()
                };
                // TLEN: positive for the leftmost mate, negative for rightmost.
                let signed: i64 = if a.ref_off <= m.ref_off {
                    frag_len as i64
                } else {
                    -(frag_len as i64)
                };
                (mate_rname, m.ref_off + 1, signed)
            }
            (Some(m), _) => {
                let mate_rname = if m.ref_id == a.ref_id {
                    "=".to_string()
                } else {
                    idx.refnames[m.ref_id as usize]
                        .split_whitespace()
                        .next()
                        .unwrap_or(&idx.refnames[m.ref_id as usize])
                        .to_string()
                };
                (mate_rname, m.ref_off + 1, 0)
            }
            (None, _) => ("=".to_string(), pos_1based, 0),
        };

        // Use the pair MAPQ for concordant pairs (BowtieMapq2 paired formula);
        // fall back to per-mate single-end MAPQ otherwise.
        let mapq = pair_mapq
            .unwrap_or_else(|| mapq_v2(a.score, self_secbest, scoring.score_min(a.read_len)));
        let nm = a.mismatches + a.gap_extends;
        let mut tags = vec![
            format!("AS:i:{}", a.score),
            "XN:i:0".to_string(),
            format!("XM:i:{}", a.mismatches),
            format!("XO:i:{}", a.gap_opens),
            format!("XG:i:{}", a.gap_extends),
            format!("NM:i:{}", nm),
            format!("MD:Z:{}", a.md),
        ];
        if let Some(m) = mate_aln {
            tags.push(format!("YS:i:{}", m.score));
        }
        tags.push(format!("YT:Z:{yt}"));

        sam_w.write_paired_record(
            &read.name, flag, &rname, pos_1based, mapq, &cigar, &rnext, pnext, tlen, &seq, &qual,
            &tags,
        )?;
    } else {
        // Unmapped half. RNAME/POS conventionally inherit mate's position
        // (so SAM-aware tools group the pair) when the mate is mapped.
        let (rname, pos_1based, rnext, pnext) = match mate_aln {
            Some(m) => {
                let r = idx.refnames[m.ref_id as usize]
                    .split_whitespace()
                    .next()
                    .unwrap_or(&idx.refnames[m.ref_id as usize])
                    .to_string();
                (r, m.ref_off + 1, "=".to_string(), m.ref_off + 1)
            }
            None => ("*".to_string(), 0, "*".to_string(), 0),
        };
        sam_w.write_paired_record(
            &read.name,
            flag,
            &rname,
            pos_1based,
            0,
            "*",
            &rnext,
            pnext,
            0,
            &read.seq,
            &read.qual,
            &[format!("YT:Z:{yt}")],
        )?;
    }
    Ok(())
}
