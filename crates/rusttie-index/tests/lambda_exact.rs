//! Phase 1 acceptance test: exact-match queries against the lambda phage
//! BT2 index produce hits at the known reference positions.

use std::path::PathBuf;

use rusttie_index::{Bt2Index, RefHit, exact_hits};

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

/// For each chosen substring length and start, the substring should be found
/// at exactly that start (and possibly other positions if it isn't unique).
#[test]
fn substrings_resolve_to_their_origin() {
    let idx = Bt2Index::open(fixture_base()).expect("open index");
    let seq = read_lambda_seq();
    assert_eq!(seq.len(), idx.params.len as usize, "fasta vs index length");

    // Sample positions across the genome, varied length.
    let queries = [
        (0, 25),
        (100, 50),
        (1234, 30),
        (10_000, 50),
        (24_251, 50), // middle of genome
        (40_000, 50),
        (48_452, 50), // last 50 bases
    ];

    for (start, len) in queries {
        let q = &seq[start..start + len];
        // Skip queries that contain non-ACGT (lambda is pure ACGT but be safe).
        if q.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let hits = exact_hits(&idx, q);
        assert!(
            !hits.is_empty(),
            "no hits for query at {start} len {len}: {}",
            std::str::from_utf8(q).unwrap()
        );
        let want = RefHit {
            ref_id: 0,
            ref_off: start as u32,
        };
        assert!(
            hits.contains(&want),
            "expected hit at {want:?} not in {hits:?} for query at {start} len {len}"
        );
    }
}

/// Queries that don't exist in the genome should produce no hits.
#[test]
fn nonexistent_query_no_hits() {
    let idx = Bt2Index::open(fixture_base()).expect("open index");
    // 30-mer of all A's is exceedingly unlikely to occur in lambda.
    let q = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let hits = exact_hits(&idx, q);
    assert!(hits.is_empty(), "unexpected hits: {hits:?}");
}
