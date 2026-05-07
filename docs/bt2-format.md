# BowTie 2 Index Format (.bt2) — Complete Byte-Level Specification

This document provides a precise, byte-level specification of the BowTie 2 index file format as implemented in the official BowTie 2 source code. All line numbers refer to `/vendor/bowtie2/` source files.

---

## File Organization

A BowTie 2 index consists of **6 files** (or 3, for unpacked indexes):

1. **`.1.bt2`** (primary file): Header, BWT, suffix array samples, lookup tables
2. **`.2.bt2`** (secondary file): Suffix array samples (offsets)
3. **`.3.bt2`** (reference metadata): RefRecords describing unambiguous stretches
4. **`.4.bt2`** (reference sequence): Packed 2-bit nucleotide reference
5. **`.rev.1.bt2`**, **`.rev.2.bt2`**: Reverse-complement index files (optional, same format as .1 and .2)

---

## Header (Common to .1 and .2 files)

The header is written identically to both `.1.bt2` and `.2.bt2` files at read time (bt2_io.cpp:813–823).

### Header Layout (in disk order)

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | uint32 | endian_hint | Endianness sentinel: `0x00000001` (little-endian) or `0x01000000` (big-endian). Read as `readU<uint32_t>()` at bt2_io.cpp:134. If reads as 1, native endian; if reads as 0x01000000, swap all subsequent reads. |
| 4 | OFF_SIZE* | TIndexOffU | len | Length of the text (number of characters indexed, excluding $). Off_SIZE is 4 bytes (32-bit) in standard mode or 8 bytes (64-bit) in BOWTIE_64BIT_INDEX mode. Read at bt2_io.cpp:159. |
| 4+OFF_SIZE | 4 | int32 | lineRate | log₂(line size). Line size = 2^lineRate bytes. Read at bt2_io.cpp:161. Default is 6 (for 32-bit) or 7 (for 64-bit). |
| 8+OFF_SIZE | 4 | int32 | linesPerSide | Not used; deprecated. Always written as 2 (bt2_idx.h:817). Read and discarded at bt2_io.cpp:163. |
| 12+OFF_SIZE | 4 | int32 | offRate | log₂(suffix array sample stride). One offset every 2^offRate rows. Read at bt2_io.cpp:165. Default is 5 (every 32 rows). |
| 16+OFF_SIZE | 4 | int32 | ftabChars | Number of characters used to index the first-level lookup table (ftab). Typically 10. Read at bt2_io.cpp:170. |
| 20+OFF_SIZE | 4 | int32 | flags | Bit-encoded flags. Negative integer indicates valid flag block. Bit 1 (EBWT_COLOR=2) set if colorspace index. Bit 2 (EBWT_ENTIRE_REV=4) set if reverse complement is entire concatenated reference reversed (as opposed to per-sequence reversal). Read at bt2_io.cpp:174–202. Value written as `-flags` where flags ∈ {1, 1|EBWT_COLOR, 1|EBWT_ENTIRE_REV, 1|EBWT_COLOR|EBWT_ENTIRE_REV} (bt2_idx.h:820–823). |

**Total header size: 24 + OFF_SIZE bytes** (28 bytes in 32-bit mode, 32 bytes in 64-bit mode)

---

## `.1.bt2` — Primary Index File

Written sequentially from header, then data sections. All multi-byte integers are endianness-swapped if `_switchEndian` is true.

### Complete `.1.bt2` Layout (in disk order)

#### Phase 1: Header (offsets 0–27)
As above.

#### Phase 2: Pattern Metadata
Written by `joinToDisk()` (bt2_idx.h:2721–2745) and read by `readIntoMemory()` (bt2_io.cpp:241–274).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 28 | OFF_SIZE | TIndexOffU | nPat | Number of reference sequences. Written at bt2_idx.h:2721. Read at bt2_io.cpp:241. |
| 28+OFF_SIZE | nPat × OFF_SIZE | TIndexOffU[] | plen[] | Pattern lengths. One entry per reference sequence. Written at bt2_idx.h:2735, 2743. Read at bt2_io.cpp:257–268. |

#### Phase 3: Fragment Metadata
Written by `joinToDisk()` and read by `readIntoMemory()`.

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 28+OFF_SIZE+(nPat × OFF_SIZE) | OFF_SIZE | TIndexOffU | nFrag | Number of fragments (unambiguous stretches) in joined reference. Written at bt2_idx.h:2745. Read at bt2_io.cpp:284. |
| previous+OFF_SIZE | nFrag × 3 × OFF_SIZE | TIndexOffU[nFrag×3] | rstarts[] | Fragment metadata: for each fragment, store (offset_in_joined_text, text_id, offset_in_text). Written during join phase and read at bt2_io.cpp:300–315. |

#### Phase 4: BWT and Occurrence Tables
Written by `buildToDisk()` (bt2_idx.h:2829) and read by `readIntoMemory()` (bt2_io.cpp:324–385).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous+nFrag×3×OFF_SIZE | eh._ebwtTotLen | uint8[] | ebwt | Extended BWT data, packed 2 bits per nucleotide (or 4 bits in SIXTY4_FORMAT). Stored in "sides" of size `eh._sideSz` bytes each. Within each side, BWT data occupies `eh._sideBwtSz` bytes, followed by 4 × OFF_SIZE bytes of occurrence counts for A, C, G, T. Each side's occurrence counts store the cumulative count of that nucleotide from the start of the BWT up to that side boundary. Read at bt2_io.cpp:332–385. |

**BWT Packing Details** (bt2_idx.h:2939–3048):
- Each byte holds 4 DNA symbols (2 bits each): bits [7:6] = char 3, bits [5:4] = char 2, bits [3:2] = char 1, bits [1:0] = char 0.
- Forward sides fill left-to-right (LSB first). Reverse sides fill right-to-left (MSB first).
- The $ symbol (end-of-text marker, position 0 in SA) is encoded as 'A' (0) in the BWT but its position (`_zOff`) is recorded separately.

#### Phase 5: Dollar Sign Position
Written by `buildToDisk()` (bt2_idx.h:3099) and read by `readIntoMemory()` (bt2_io.cpp:388–390).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous | OFF_SIZE | TIndexOffU | zOff | BWT index (row number) corresponding to the suffix starting at position 0 in the text (the $ position). This is where the LF mapping resolves position 0. Read at bt2_io.cpp:388. |

#### Phase 6: First-Character Boundary Table (fchr)
Written by `buildToDisk()` (bt2_idx.h:3119–3120) and read by `readIntoMemory()` (bt2_io.cpp:402–410).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous | 5 × OFF_SIZE | TIndexOffU[5] | fchr[5] | Cumulative first-character counts: fchr[c] = number of BWT rows starting with a character < c. Ordered: fchr[$], fchr[A], fchr[C], fchr[G], fchr[T]. fchr[0] = 0 (always). fchr[4] = len (total BWT length, excluding $). Built via exclusive prefix sum of observed character counts (bt2_idx.h:3105–3113). Read at bt2_io.cpp:404–408. |

#### Phase 7: First Lookup Table (ftab)
Written by `buildToDisk()` (bt2_idx.h:3162–3163) and read conditionally by `readIntoMemory()` (bt2_io.cpp:421–440).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous | (2^(2×ftabChars) + 1) × OFF_SIZE | TIndexOffU[ftabLen] | ftab[ftabLen] | First-level lookup table. ftabLen = (1 << (ftabChars × 2)) + 1. Each entry ftab[i] holds the first BWT row where the first ftabChars characters match pattern i (or, if i encodes a short pattern, the first row of a longer pattern). Entries are typically 32-bit row offsets but may encode pointers into eftab (Extended FTab). Built via prefix sum during BWT building (bt2_idx.h:3129). Only loaded if loadFtab=true in readIntoMemory (bt2_io.cpp:421). |

#### Phase 8: Extended Lookup Table (eftab)
Written by `buildToDisk()` (bt2_idx.h:3166–3167) and read conditionally by `readIntoMemory()` (bt2_io.cpp:451–470).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous | (ftabChars × 2) × OFF_SIZE | TIndexOffU[eftabLen] | eftab[eftabLen] | Extended FTab. eftabLen = ftabChars × 2. Stores low/high row ranges for patterns shorter than ftabChars characters. Populated during BWT construction to handle short suffixes that "absorb" into longer ones. Accessed indirectly via ftab entries that encode pointers (negated, XOR'd with OFF_MASK at bt2_idx.h:3153). Read conditionally at bt2_io.cpp:459–470. |

#### Phase 9: Reference Names
Written by `writeFromMemory()` after ftab/eftab and read by `readIntoMemory()` (bt2_io.cpp:496–511).

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| previous | variable | char[] | refnames | Zero-terminated reference sequence names. Format: null byte = end of all names; newline (`\n`) = separator between names. One name per reference sequence. Read character-by-character until a null terminator is encountered (bt2_io.cpp:499–510). Embedded null terminators are possible but treated as name delimiters. Populated by `_refnames` vector during index building. |

---

## `.2.bt2` — Secondary Index File

Written by `writeFromMemory()` (bt2_idx.h:801–864) and read by `readIntoMemory()` (bt2_io.cpp:514–614).

### Complete `.2.bt2` Layout

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | uint32 | endian_hint | Same as `.1.bt2` (bt2_idx.h:814). |
| 4 | 4 | skip | (padding) | Unused bytes (reserved). |
| 8 | offsLenSampled × OFF_SIZE | TIndexOffU[] | offs[] | Suffix array sample (offset) array. One entry for every 2^offRate rows in the BWT. offsLenSampled = (bwtLen + (1 << offRate) - 1) >> offRate. Each entry is the reference position of the sampled BWT row (read at bt2_io.cpp:581–602). May be sampled sparsely if offRate is overridden at load time. Written directly to file at bt2_idx.h:3011. |

**Note:** `.2.bt2` is not loaded if loadSASamp=false in readIntoMemory (bt2_io.cpp:65–70).

---

## `.3.bt2` — Reference Metadata File

Contains RefRecord entries describing unambiguous stretches of the reference. Read by `BitPairReference` constructor (reference.cpp:30–171) and written by `BitPairReference::szsFromFasta()` (reference.cpp:587–646).

### `.3.bt2` Layout

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | 4 | uint32 | endian_hint | Endianness sentinel (1 or 0x01000000). Written at reference.cpp:613. Read at reference.cpp:103. |
| 4 | OFF_SIZE | TIndexOffU | num_records | Number of RefRecord entries. Written at reference.cpp:623, 637. Read at reference.cpp:115. Must be > 0 (reference.cpp:116). |
| 4+OFF_SIZE | num_records × (2×OFF_SIZE + 1) | RefRecord[] | records | Array of RefRecord structures. Each record describes an unambiguous stretch: off (offset in text), len (length), first (boolean flag indicating start of new reference sequence). Written at reference.cpp:625, 638. Read at reference.cpp:131. See RefRecord structure below. |

**RefRecord Structure** (ref_read.h:73–103):
Each RefRecord is packed as: `[off: OFF_SIZE bytes] [len: OFF_SIZE bytes] [first: 1 byte]`

| Field | Type | Size | Notes |
|-------|------|------|-------|
| off | TIndexOffU | OFF_SIZE | Offset (in characters) of this unambiguous stretch within its containing reference sequence |
| len | TIndexOffU | OFF_SIZE | Length (in characters) of this unambiguous stretch |
| first | bool (uint8) | 1 | 1 if this record starts a new reference sequence, 0 otherwise |

Written using `RefRecord::write()` (ref_read.h:94–97) and read via `RefRecord` constructor from FILE (ref_read.h:79–92).

---

## `.4.bt2` — Packed Reference Sequence File

Contains the reference sequences compressed to 2 bits per nucleotide. Written during index building and read by `BitPairReference` constructor (reference.cpp:198–225).

### `.4.bt2` Layout

| Offset | Size | Type | Name | Description |
|--------|------|------|------|-------------|
| 0 | (bufSz + 3) / 4 | uint8[] | packed_ref | Packed 2-bit nucleotide reference. Total unambiguous characters (sum of all RefRecord.len values) is round up to the nearest 4-character boundary. Each byte encodes 4 nucleotides: bits [7:6] = nt[0], bits [5:4] = nt[1], bits [3:2] = nt[2], bits [1:0] = nt[3]. Written by `BitpairOutFileBuf` during `fastaRefReadSizes()` (reference.cpp:621, 636). Read directly into `buf_` at reference.cpp:183–190 or reference.cpp:226–237. |

**Nucleotide Encoding (2 bits):**
- `00` = A
- `01` = C
- `10` = G
- `11` = T

**Packing is LSB-first** (correction to original draft): the first character of the stretch goes in bits 1–0, the second in bits 3–2, etc. Verified empirically against `.4.bt2` for lambda phage (first byte `0x6a` decodes to `GGGC`, matching the FASTA). This matches the `.1.bt2` BWT packing convention.

---

## Reverse-Complement Index Files

Files `.rev.1.bt2` and `.rev.2.bt2` follow the identical format to `.1.bt2` and `.2.bt2`. The entire concatenated reference is reversed (character-reversed, then reverse-complemented) before building the BWT. The `entireReverse` flag (bit 2 of the flags field) indicates this. Colorspace indexes are built from DNA followed by color-space, and `.rev.*` files contain the color-space reverse complement.

---

## Key Parameters Embedded in the Index

The EbwtParams structure (bt2_idx.h:112–254) fully describes the index geometry:

| Parameter | Source | Calculation |
|-----------|--------|-------------|
| `len` | Header | Text length (read from offset 4) |
| `bwtLen` | Computed | = len + 1 (includes $) |
| `lineRate` | Header | = log₂(line_size) |
| `lineSz` | Computed | = 2^lineRate |
| `offRate` | Header | = log₂(sample_stride) |
| `offMask` | Computed | = OFF_MASK << offRate |
| `ftabChars` | Header | Characters used for first lookup |
| `ftabLen` | Computed | = (1 << (ftabChars × 2)) + 1 |
| `eftabLen` | Computed | = ftabChars × 2 |
| `offsLen` | Computed | = (bwtLen + (1 << offRate) - 1) >> offRate |
| `offsSz` | Computed | = offsLen × OFF_SIZE |
| `sideSz` | Computed | = lineSz (one side per line) |
| `sideBwtSz` | Computed | = sideSz - 4×OFF_SIZE (BWT portion of side) |
| `sideBwtLen` | Computed | = sideBwtSz × 4 (characters, since 2-bit packed) |
| `numSides` | Computed | = (bwtSz + (sideBwtSz - 1)) / sideBwtSz |
| `numLines` | Computed | = numSides |
| `ebwtTotLen` | Computed | = numSides × sideSz (total BWT+occurrence storage) |
| `color` | Header (flags) | Colorspace index if set |
| `entireReverse` | Header (flags) | Reverse-entire if set |

The init() method (bt2_idx.h:133–167) populates all derived values from len, lineRate, offRate, ftabChars, color, and entireReverse.

---

## Endianness and Byte-Swapping

All multi-byte integers are subject to endianness. The first 4 bytes of `.1.bt2` and `.2.bt2` are a sentinel:
- `0x00000001` = native endian (little-endian on x86)
- `0x01000000` = opposite endian (will be swapped)

After reading the sentinel, if it does not equal 1 (native) and equals 0x01000000 (opposite), set `_switchEndian = true` and call `endianSwapU()` on all subsequent reads (bt2_io.cpp:134–147, word_io.h).

The `.3.bt2` file has its own endianness sentinel (reference.cpp:103, 613).

---

## Colorspace Indexing

In colorspace mode (color=true in header flags):
- The BWT is built over color-space nucleotides (00=A→C, 01=C→G, etc.), not base-space.
- The `.3.bt2` and `.4.bt2` files store the base-space reference (not colors), for reconstruction.
- The reverse `.rev.1.bt2` / `.rev.2.bt2` store the color-space reverse complement (reference.cpp:614–633).
- Four additional zero-bytes are written at the beginning of the color-space reference (to encode the first color).

---

## Building to Disk vs. In-Memory

BowTie 2 builds the index in phases:
1. **Header**: Written first (bt2_idx.h:813–823, bt2_io.cpp:813–823).
2. **Join** (reference metadata): nPat, plen, nFrag, rstarts (bt2_idx.h:2721–2745). These are small; written early.
3. **Build** (BWT + occurrence data): Streamed to disk side-by-side with suffix-array computation (bt2_idx.h:2829–3174). The ebwt array is written chunk-by-chunk as "sides" are completed; offsets are written to the secondary file in parallel.
4. **Finalize**: zOff, fchr, ftab, eftab (bt2_idx.h:3099–3168). These are computed only after the entire BWT is built.
5. **Reference metadata**: The `.3.bt2` and `.4.bt2` files are written during FASTA reading, in parallel with suffix-array input, before the BWT build (reference.cpp:587–646).

---

## Conditional Loading and Sparse Sampling

The `readIntoMemory()` function accepts flags to load only parts of the index:
- `justHeader=true`: Stop after reading nPat and plen (bt2_io.cpp:282).
- `loadFtab=false`: Skip ftab and eftab, leaving them NULL (bt2_io.cpp:421, 478–487).
- `loadSASamp=false`: Skip offs array and `.2.bt2` file (bt2_io.cpp:65–70, 514).
- `loadRstarts=false`: Skip rstarts, leaving it NULL (bt2_io.cpp:292, 317–322).
- `_overrideOffRate`: Read every 2^(_overrideOffRate) rows instead of 2^offRate. This downsamples the loaded offs array (bt2_io.cpp:222–237, 542–614).

---

## Sanity Checking and Verification

The `writeFromMemory()` function can optionally perform a sanity check by reading the written files back and comparing all fields (bt2_io.cpp:891–911 has a disabled example). When disabled (the default), no re-read verification occurs.

---

## References to Source Code

All claims in this spec reference specific line numbers in the BowTie 2 source:

- Header definition and reading: **bt2_io.cpp:134–174** (read) and **bt2_idx.h:813–823** (write).
- BWT building and packing: **bt2_idx.h:2829–3174** (buildToDisk).
- BWT I/O: **bt2_io.cpp:324–385** (read) and **bt2_idx.h:847–851** (write).
- FTab/EFTab: **bt2_idx.h:3145–3168** (build and write) and **bt2_io.cpp:412–487** (read).
- Offsets (suffix array sample): **bt2_idx.h:3007–3011** (write) and **bt2_io.cpp:581–602** (read).
- Reference metadata: **reference.cpp:613–638** (write) and **reference.cpp:103–131** (read).
- Packed reference: **reference.cpp:609–646** (write) and **reference.cpp:178–225** (read).

---

## Alignment and Padding

- No explicit padding between fields; all structures are byte-aligned.
- When writing an array of fixed-size records, offsets are cumulative and depend on the exact size of TIndexOffU (4 or 8 bytes).
- The `.4.bt2` file may have up to 3 trailing padding bytes if (total_unambiguous_chars % 4) != 0, to round up to the nearest 32-bit word boundary.

---

## Notes for Implementation

1. **Endianness is global per file**: Once the sentinel is read, all subsequent values in that file are endian-swapped uniformly.
2. **The $ symbol is implicit**: It is not stored in the BWT but recorded as `zOff`. The BWT has len+1 rows (including the row for $), but $ itself does not appear in the data stream; it is encoded as 'A' during packing and marked separately.
3. **Fragment to reference mapping**: The rstarts array is crucial for mapping BWT positions back to reference coordinates. Each fragment has 3 values: joined-text offset, reference ID, and offset within reference.
4. **The BWT is side-organized**: Each "side" is a fixed-size block (lineSz bytes) containing BWT data (sideBwtSz bytes) and a 4-entry occurrence checkpoint (4 × OFF_SIZE bytes). This layout enables efficient counting within a side via lookup and addition.
5. **First-level lookup (ftab) is mandatory for search**: While `loadFtab=false` is allowed during reading, the search operations expect ftab to be present. The eftab is used for short patterns that don't have a dedicated ftab entry.
6. **Suffix array sampling is lossy**: offs[] is not a complete suffix array; it is a sparse sample at stride 2^offRate. Full positions are recomputed on demand using the LF mapping during search.

