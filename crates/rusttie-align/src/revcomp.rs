//! Reverse complement for ASCII DNA. N → N.

#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'N' | b'n' => b,
        _ => b'N',
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_simple() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATC"), b"GATC");
        assert_eq!(reverse_complement(b"AANTT"), b"AANTT");
    }
}
