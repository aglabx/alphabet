use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Monomer {
    pub array_id: String,
    pub idx: u32,
    pub length: u32,
    pub period: u16,
    pub source: String,
    pub ed_prev: String,
    pub ed_next: String,
    pub sequence: String,
    pub strand_flipped: bool,
    pub rotation_offset: i32,
    pub anchor_confidence: f64,
}

const COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table
};

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMPLEMENT[b as usize]).collect()
}

pub fn revcomp_str(seq: &str) -> String {
    String::from_utf8(revcomp(seq.as_bytes())).unwrap()
}

/// Homopolymer compression: collapse consecutive runs of the same base.
pub fn hpc(seq: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(seq.len());
    for &b in seq {
        if result.last().copied() != Some(b) {
            result.push(b);
        }
    }
    result
}

/// Encode nucleotide to 0-3 (A=0, C=1, G=2, T=3), N=255
pub fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}

/// Compute k-mer hash for a window. Returns None if any N.
pub fn kmer_hash(seq: &[u8], k: usize) -> Option<u64> {
    let mut hash: u64 = 0;
    for &b in &seq[..k] {
        let enc = encode_base(b);
        if enc == 255 {
            return None;
        }
        hash = (hash << 2) | (enc as u64);
    }
    Some(hash)
}

/// Count occurrences of anchor k-mers in a sequence
pub fn count_anchor_hits(seq: &[u8], anchor_set: &std::collections::HashSet<u64>, k: usize) -> u32 {
    if seq.len() < k {
        return 0;
    }
    let mut count = 0u32;
    for i in 0..=(seq.len() - k) {
        if let Some(h) = kmer_hash(&seq[i..i + k], k) {
            if anchor_set.contains(&h) {
                count += 1;
            }
        }
    }
    count
}

/// Build a set of frequent k-mers from a seed sequence
pub fn build_anchor_kmers(seed: &[u8], k: usize, top_n: usize) -> std::collections::HashSet<u64> {
    use std::collections::HashMap;
    let mut counts: HashMap<u64, u32> = HashMap::new();
    if seed.len() < k {
        return std::collections::HashSet::new();
    }
    for i in 0..=(seed.len() - k) {
        if let Some(h) = kmer_hash(&seed[i..i + k], k) {
            *counts.entry(h).or_insert(0) += 1;
        }
    }
    let mut pairs: Vec<_> = counts.into_iter().collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1));
    pairs.into_iter().take(top_n).map(|(h, _)| h).collect()
}

/// Cyclic rotation of a sequence
pub fn cyclic_rotate(seq: &[u8], offset: usize) -> Vec<u8> {
    let n = seq.len();
    if n == 0 || offset.is_multiple_of(n) {
        return seq.to_vec();
    }
    let off = offset % n;
    let mut result = Vec::with_capacity(n);
    result.extend_from_slice(&seq[off..]);
    result.extend_from_slice(&seq[..off]);
    result
}
