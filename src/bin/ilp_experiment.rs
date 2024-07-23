use std::collections::HashSet;
use std::{fmt::Debug, str::from_utf8};
use std::fmt::{Formatter, Result};
use counter::Counter;
use aho_corasick::AhoCorasick;

#[derive(PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
struct Segment(Vec<u8>);

impl Debug for Segment {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", from_utf8(&self.0).unwrap())?;
        Ok(())
    }
}

fn main() {
    // let seqs = [
    //     b"ACGTGCGTGC".to_vec(),
    //     b"ACACACACAC".to_vec(),
    //     b"ACGTGTGCGC".to_vec()
    // ].into_iter();

    use bio::io::fasta;
    use std::fs::File;
    let filename = "phiX174/phiX174.2line.fna";
    let f = File::open(filename).expect("Unable to read file.");
    let seqs = fasta::Reader::new(f)
        .records()
        .map(|x| x.expect("Incorrect fasta record.").seq().to_vec());

    let counts = get_valid_segments(seqs, 5);
    // println!("{:?}", counts);
    // let mut counts: Vec<_> = counts.iter().collect();
    // counts.sort();
    // println!("{:?}", counts.len());
    for (s, n) in counts {
        println!("{}\t{:?}", n, s);
    }
}

fn get_valid_segments<T>(seqs: T, k: usize) -> Counter<Segment, usize> 
    where T: Iterator<Item = Vec<u8>>
{
    let mut count = Counter::new();
    for s in seqs {
        let x = seq_to_segments(&s, k);
        count.extend(x);
        let x = border_segments(&s, k);
        count.extend(x);
    }
    return count;
}

fn count_kmers(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
    let mut count = HashSet::new();
    for i in 0..seq.len()-k {
        let kmer = seq[i..i+k].to_vec();
        count.insert(kmer);
    }
    return count;
}

fn border_segments(seq: &[u8], k: usize) -> Counter<Segment, usize> {
    let kmers = count_kmers(seq, k);

    let sentinel = vec![b'$'; k];
    let mut seq2 = Vec::new(); // sentinel + seq + sentinel;
    seq2.extend_from_slice(&sentinel);
    seq2.extend_from_slice(seq);
    seq2.extend_from_slice(&sentinel);

    let mut counts = Counter::new();

    for kmer in kmers.iter() {
        let ac = AhoCorasick::new([kmer]).unwrap();
        if let Some(m) = ac.find(&seq2) {
            counts[&Segment(seq2[0..m.end()].to_vec())] += 1;
        };
    }

    let mut seq3 = seq2.clone();
    seq3.reverse();
    let n = seq3.len();
    for kmer in kmers.iter() {
        let mut kmer = kmer.clone();
        kmer.reverse();
        let ac = AhoCorasick::new([kmer]).unwrap();
        if let Some(m) = ac.find(&seq3) {
            let end = n - m.end();
            counts[&Segment(seq2[end..n].to_vec())] += 1;
        };
    }
    counts[&Segment(seq2)] += 1;
    return counts;
}

fn seq_to_segments(seq: &[u8], k: usize) -> Counter<Segment, usize> {
    let mut count = Counter::new();
    let n = seq.len();

    for i in 0..n-k {
        let ac = AhoCorasick::new([&seq[i..i+k]]).unwrap();
        let end = match ac.find(&seq[i+1..]) {
            Some(m) => { i+1+m.end() },
            None => { n },
        };

        for j in i+k..end {
            if seg_is_valid(&seq[i..=j], k) {
                count[&Segment(seq[i..=j].to_vec())] += 1;
            }
        }
    }
    return count;
}

// -- test related stuff
#[cfg(test)]
fn all_segments(seq: &[u8], k: usize) -> Counter<Segment, usize> {
    let mut count = Counter::new();
    let n = seq.len();
    for i in 0..n-k {
        for j in i+k..n {
            count[&Segment(seq[i..=j].to_vec())] += 1;
        }
    }
    return count;
}

fn seg_is_valid(seg: &[u8], k: usize) -> bool {
    let n = seg.len();
    let ac = AhoCorasick::new([&seg[0..k], &seg[n-k..n]]).unwrap();
    match ac.find(&seg[1..n-1]) {
        Some(_) => return false,
        None => return true,
    }
}

#[test]
fn test_seg_validation() {
    assert!(seg_is_valid(b"ACGTAC", 2));
    assert!(!seg_is_valid(b"ACACAC", 2));
}

#[test]
fn test_seq_to_seg() {
    use std::collections::HashMap;
    let seq = b"ACGTGCGTGC";
    let result = all_segments(seq, 2);
    let result: HashMap<_, _> = 
        result.into_iter().filter(|(seg, _)| seg_is_valid(&seg.0, 2)).collect();

    let count = seq_to_segments(seq, 2).into_map();
    assert!(result == count);
}

