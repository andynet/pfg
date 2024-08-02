use std::collections::HashSet;
use std::{fmt::Debug, str::from_utf8};
use std::fmt::{Formatter, Result};
use counter::Counter;
use aho_corasick::AhoCorasick;
use bio::io::fasta;
use itertools::Itertools;
use std::fs::File;
use pfg::pf::load_trigs;
use std::io::{BufWriter, Write};


fn main() {
    let sequences = "phiX174/phiX174.2line.clean.fna";
    let trigger_file = "phiX174/nonshiftable.08mer.txt";
    let output = "phiX174/nonshiftable.08mer.segments.txt";
    // let sequences = "small_seqs/seqs.fna";
    // let trigger_file = "small_seqs/nonshiftable.06mer.txt";

    let f = File::open(sequences).expect("Unable to read file.");
    let seqs = fasta::Reader::new(f).records()
        .map(|x| x.expect("Incorrect fasta record.").seq().to_vec());

    let kmers = load_trigs(trigger_file);
    let k = kmers.first().expect("Needs at least 1 trigger word.").len();

    let mut counter: Counter<Segment, usize> = Counter::new();
    for seq in seqs {
        let seq = [vec![b'$'; k], seq, vec![b'$'; k]].concat();
        counter.extend(get_inner_segments(&seq, &kmers));
        counter.extend(get_border_segments(&seq, &kmers));
        counter[&Segment(seq)] += 1;
    }

    let mut out = BufWriter::new(File::create(output).expect("Cannot open file for writing."));
    for (s, n) in counter {
        writeln!(out, "{}\t{}", n, from_utf8(&s.0).unwrap()).expect("Error writing.");
    }
}

#[derive(PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
struct Segment(Vec<u8>);

impl Debug for Segment {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", from_utf8(&self.0).unwrap())?;
        Ok(())
    }
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

fn get_inner_segments(seq: &[u8], kmers: &[Vec<u8>]) -> Counter<Segment, usize> {
    let k = kmers[0].len();

    let mut count = Counter::new();
    let ac = AhoCorasick::new(kmers).unwrap();
    let matched = ac.find_overlapping_iter(seq).map(|m| m.start()).collect_vec();

    for i in 0..matched.len()-1 {
        let first_kmer = &seq[matched[i]..matched[i]+k];

        let mut seen = HashSet::new();
        for j in i+1..matched.len() {
            let last_kmer = &seq[matched[j]..matched[j]+k];

            // skips segments with structure x..y..y
            if seen.contains(last_kmer) { continue; }
            seen.insert(last_kmer);

            let segment = Segment(seq[matched[i]..matched[j]+k].to_vec());
            count[&segment] += 1;

            // skips segments with structure x..x..y (and terminates next search)
            if last_kmer == first_kmer { break; }
        }
    }
    return count;
}

#[test]
fn test_get_inner_segments() {
    let seq = b"AGCTGATCGTCGCTAGTCAA";
    let res_inner: Counter<Segment, usize> = Counter::from_iter([
        Segment(b"TCGCTAGTC".to_vec()),
        Segment(b"CTGATCGTCGCT".to_vec()),
        Segment(b"TCGTC".to_vec()),
        Segment(b"CTGATC".to_vec()),
        Segment(b"TCGCT".to_vec()),
        Segment(b"CTAGTC".to_vec())
    ]); 

    //  seq = b"  *   *  *  *   *   ";
    //  seq = b"  CTGATC            "  ;
    //  seq = b"      TCGTC         "  ;
    //  seq = b"         TCGCT      "  ;
    //  seq = b"            CTAGTC  "  ;
    //  seq = b"  CTGATCGTCGCT      "  ;
    //  seq = --------TCGTC----------  ;
    //  seq = b"         TCGCTAGTC  "  ;
    //  outer ========================
    //  seq = $$AGCT                "  ;
    //  seq = $$AGCTGATC            "  ;
    //  seq = b"            CTAGTCAA$$";
    //  seq = b"                TCAA$$";
    let kmers = &[b"CT".to_vec(), b"TC".to_vec()];

    let tmp = get_inner_segments(seq, kmers);
    println!("{:?}", tmp);
    let tmp2 = get_border_segments(seq, kmers);
    println!("{:?}", tmp2);
    // Counter { map: {AGCTGATC: 1, TCAA: 1, AGCT: 1, CTAGTCAA: 1}, zero: 0 }

}

fn get_border_segments(seq: &[u8], kmers: &[Vec<u8>]) -> Counter<Segment, usize> {
    let k = kmers[0].len();

    let mut count = Counter::new();
    let ac = AhoCorasick::new(kmers).unwrap();
    let matched = ac.find_overlapping_iter(seq).map(|m| m.start()).collect_vec();

    let mut seen = HashSet::new();
    for i in 0..matched.len() {
        let last_kmer = &seq[matched[i]..matched[i]+k];

        // skips segments with structure x..y..y
        if seen.contains(last_kmer) { continue; }
        seen.insert(last_kmer);

        let segment = Segment(seq[..matched[i]+k].to_vec());
        count[&segment] += 1;

        // we have seen all kmers, no point in continuing
        if seen.len() == kmers.len() { break; } 
    }

    seen.clear();
    for i in (0..matched.len()).rev() {
        let first_kmer = &seq[matched[i]..matched[i]+k];

        // skips segments with structure x..x..y
        if seen.contains(first_kmer) { continue; }
        seen.insert(first_kmer);

        let segment = Segment(seq[matched[i]..].to_vec());
        count[&segment] += 1;

        // we have seen all kmers, no point in continuing
        if seen.len() == kmers.len() { break; } 
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

