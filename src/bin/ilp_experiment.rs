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
    let trigger_file = "phiX174/nonshiftable_v2.06mer.txt";
    let output = "phiX174/nonshiftable_v2.06mer.segments.txt";
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

#[test]
fn test_get_segments() {
    let kmers = &[b"CT".to_vec(), b"TC".to_vec()];
    let seq1 = b"AGCTGATCGTCGCTAGTCAA";
    //         b"  *   *  *  *   *   "
    let res_inner: Counter<Segment, usize> = Counter::from_iter([
        Segment( b"CTGATC".to_vec()),
        Segment( b"CTGATCGTCGCT".to_vec()),
        Segment(     b"TCGTC".to_vec()),
        Segment(        b"TCGCTAGTC".to_vec()),
        Segment(        b"TCGCT".to_vec()),
        Segment(           b"CTAGTC".to_vec())
    ]); 

    let seq2 =  b"$$AGCTGATCGTCGCTAGTCAA$$";
    let res_border: Counter<Segment, usize> = Counter::from_iter([
        Segment(b"$$AGCT".to_vec()),
        Segment(b"$$AGCTGATC".to_vec()),
        Segment(                  b"TCAA$$".to_vec()),
        Segment(              b"CTAGTCAA$$".to_vec()),
    ]);
    let inner = get_inner_segments(seq1, kmers);
    assert_eq!(inner, res_inner);
    let border = get_border_segments(seq2, kmers);
    assert_eq!(border, res_border);
}
