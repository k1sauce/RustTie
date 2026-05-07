# RustTie

A Rust port of [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/), aiming
for drop-in compatibility on the supported subset: same `.bt2` index format,
same SAM output, same default scoring and end-to-end alignment semantics.

**Status: alpha.** Not for production. The default-config paired-end
end-to-end pipeline matches Bowtie 2 closely enough for diffing experiments
(see [Validation](#validation)), but there are known feature gaps (see
[Limitations](#limitations)).

## What works

- **`.bt2` index format** — `rusttie-build` produces files byte-identical
  to `bowtie2-build` on the test corpora (lambda phage, multi-contig with N
  gaps, human chr22). Existing `bowtie2-build` indexes load and align
  correctly with `rusttie`.
- **Paired-end end-to-end alignment** with BT2-faithful semantics:
  - BT2's default scoring (`--mp 6,2 --rdg 5,3 --rfg 5,3 --score-min L,-0.6,-0.6`).
  - Quality-scaled mismatch penalties.
  - Multi-seed search with seed length 22 and BT2's `S,1,1.15` interval.
  - Descent driver with `-D` failure budget and `-R` re-seeding.
  - **Mate-rescue** for missing-partner alignments (BT2's
    `extendSeedsPaired` mate-find step).
  - **BT2-faithful paired MAPQ** (`BowtieMapq2` with summed pair scores +
    second-best concordant pair).
- **SAM and BAM output** — BAM via `noodles-bam`.
- **Gzipped FASTQ input.**
- **Multi-threading** via `rayon` — one parallel-map per batch, deterministic
  output order.

### Validation

**Synthetic** (10k paired 100bp wgsim reads from human chr22, GRCh38):

| Metric | RustTie vs Bowtie 2 |
|---|---|
| `.bt2` files | byte-identical |
| Reads mapped | 19,983 / 20,000 (BT2: 20,000) |
| Position agreement | 97.8% |
| TLEN agreement | 98.8% |
| CIGAR agreement | 100.0% |
| AS / NM agreement | 99.6% |
| MD agreement | 99.3% |
| MAPQ agreement | 93.8% |
| Wall time (`-p 8`) | 1.0s (BT2: 0.48s) |

**Real-data** (NA12878 mitochondrial paired-end Illumina reads from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets/tree/Sarek),
322,856 reads, hg38 chrM reference):

| Metric | RustTie vs Bowtie 2 |
|---|---|
| `.bt2` files | byte-identical |
| Reads mapped | 1,040 (BT2: 1,047) — 99.3% of BT2's recall |
| Position agreement | 97.5% |
| MAPQ agreement | 97.1% |
| CIGAR / AS / NM / MD | 96.5 – 96.9% |

Reproducible via `scripts/real_data_validate.sh`. Remaining MAPQ
disagreement is mostly `both_low` cases — both tools report low MAPQ on
multi-mappers but land on slightly different exact bin values from
disagreeing alternate sets. See [`rusttie.md`](./rusttie.md) for the full
per-phase development log.

## Quick start

### Build

```bash
# Requires Rust 1.85 (edition 2024). The sais-rs dependency wants nightly;
# pin via rust-toolchain.toml in this repo.
cargo build --release
```

### Index a reference

```bash
target/release/rusttie-build my_reference.fa my_index
# Produces my_index.{1,2,3,4,rev.1,rev.2}.bt2
```

You can also use an existing `bowtie2-build` index — they're byte-compatible.

### Align reads

```bash
# Paired-end
target/release/rusttie -p 8 -x my_index \
    -1 reads_R1.fq -2 reads_R2.fq \
    -S out.sam

# Single-end
target/release/rusttie -p 8 -x my_index \
    -U reads.fq \
    -S out.sam

# Output to BAM (auto-detected by extension)
target/release/rusttie -p 8 -x my_index \
    -1 reads_R1.fq -2 reads_R2.fq \
    -S out.bam
```

## Supported flags

The flag set is a strict subset of Bowtie 2's. Where flags exist they take
the same meaning and defaults as upstream.

| Flag | Meaning |
|---|---|
| `-x <BASE>` | Index basename (`<base>.1.bt2` etc.) |
| `-1` / `-2` | Paired-end mate 1 / mate 2 FASTQ |
| `-U` | Single-end FASTQ |
| `-S` | Output SAM/BAM (extension auto-detected) |
| `-p` / `--threads` | Threads (0 = all cores) |
| `-D` / `--descent-budget` | Consecutive seed-extension failures (default 15) |
| `-R` / `--descent-reseed` | Max re-seedings on repetitive seeds (default 2) |
| `--mate-rescue <K>` | Mate-rescue from top-K anchors per side (default 3, 0 disables) |
| `--seed-hit-cap` | Per-seed hit cap (default 50; tuning knob) |
| `--mp MX,MN` | Mismatch penalty bounds |
| `--rdg O,E` | Read-gap open/extend |
| `--rfg O,E` | Reference-gap open/extend |
| `--score-min L,A,B` | Score-min function (only `L` supported) |
| `--no-unal` | Suppress unmapped records |
| `--no-head` | Suppress SAM header |
| `--very-fast` / `--fast` / `--sensitive` / `--very-sensitive` | Accepted as no-ops (defaults already match `--sensitive`) |

## Limitations

These are *not* implemented and are left as future work:

- **64-bit "large" `.bt2l` indexes** — references larger than ~4 Gbp
  (e.g., wheat, some metagenomes). Human GRCh38 (3.1 Gbp) fits in the
  small index and works.
- **Local alignment** (`--local`). End-to-end only.
- **Mate-pair orientations other than FR.** `--rf` / `--ff` are ignored.
- **Custom fragment range** (`-I` / `-X`). Hardcoded to `[0, 500]`.
- **Mismatches in seeds** (`-N 1`). We always use `-N 0`.
- **Read group tagging** (`--rg` / `--rg-id`). No `@RG` header injection.
- **Multi-alignment reporting** (`-k <int>` / `-a`). Always reports best.
- **Read trimming flags** (`-3` / `-5` / `--trim-to`).
- **Quality encodings other than Phred+33.**
- **Compressed input other than gzip.** No bzip2 / zstd.
- **Stdin input** (`-U -`).
- **Multi-file input** (`-U f1.fq,f2.fq`).
- **Color-space.** BT2 deprecated this; we don't implement it either.

## How this is organized

The workspace has four crates:

- [`rusttie-io`](./crates/rusttie-io) — FASTQ/SAM/BAM I/O.
- [`rusttie-index`](./crates/rusttie-index) — `.bt2` reader, BWT/FM-index,
  SA construction (via [`sais-rs`](https://crates.io/crates/sais-rs)).
- [`rusttie-align`](./crates/rusttie-align) — seed-and-extend, descent
  driver, MAPQ, paired-end logic, mate-rescue.
- [`rusttie-cli`](./crates/rusttie-cli) — the `rusttie` and `rusttie-build`
  binaries.

The full development log is in [`rusttie.md`](./rusttie.md). Bowtie 2's
source is vendored under `vendor/bowtie2/` for reference (excluded from
published crates).

## License

Dual-licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](./LICENSE-APACHE))
- MIT license ([LICENSE-MIT](./LICENSE-MIT))

at your option.

Bowtie 2 (vendored under `vendor/bowtie2/`) is licensed separately under
GPL-3.0; see [`vendor/bowtie2/LICENSE`](./vendor/bowtie2/LICENSE). It is
included here only as reference material for porting and testing — no
GPL-3 code is incorporated into RustTie's published binary or library.
