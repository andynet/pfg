#![allow(unused_imports)]

use aho_corasick::AhoCorasick;
use bio::{data_structures::suffix_array::SuffixArray, io::fasta};
use itertools::Itertools;
use std::fs::File;
use std::io::{BufReader, Write};
use pfg::kmap::{KMap, Kmer};
use pfg::pf;
use std::str::from_utf8;

// use rand::{self, seq::IteratorRandom};
// let n = 20;
// let mut rng = rand::thread_rng();
// let v = (&kmap).into_iter().choose_multiple(&mut rng, n);

fn main() {
    let ks = [2, 3, 4, 5, 6, 8, 10, 12, 14, 16];
    let filename = "phiX174/phiX174.2line.clean.fna";
    let output = "phiX174/nonshiftable";

    for k in ks {
        let mut kmap: KMap = KMap::new(k);
        let f = File::open(filename).unwrap();
        let mut records = fasta::Reader::new(f).records();
        while let Some(Ok(record)) = records.next() {
            let seq = record.seq().to_owned();
            kmap.add_kmers(&seq);
        }

        let trigs = (&kmap).into_iter()
            .map(|(k, _)| k.clone()).collect_vec();
        let n_all = trigs.len();

        let nonshiftable = trigs.into_iter()
            .filter(|kmer| ! is_shiftable(&kmap, kmer)).collect_vec();
        let n_nshf = nonshiftable.len();

        println!("k={k}:\t{n_nshf} / {n_all}");
        let out_filename = output.to_string() + &format!(".{k:02}mer.txt");
        let mut f = File::create(out_filename).unwrap();
        for kmer in nonshiftable {
            writeln!(f, "{kmer}").expect("Cannot write to file.");
        }
    }
}

fn is_shiftable(kmap: &KMap, kmer: &Kmer) -> bool {
    let left = kmap.to_left(kmer);
    let right = kmap.to_right(kmer);

    // maybe measure == 1 are shiftable?
    if left.len() != 1 || right.len() != 1 { return false; }
    if kmap.to_right(&left[0]).len() != 1 { return false; }
    if kmap.to_left(&right[0]).len() != 1 { return false; }
    return true
}

#[test]
fn visualize_6mer_contexts() {
    let k = 6;
    let filename = "phiX174/phiX174.2line.clean.fna";

    let mut kmap: KMap = KMap::new(k);
    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq().to_owned();
        kmap.add_kmers(&seq);
    }

    let mut trigs = vec![
        b"TTGCTG".to_vec(), b"TTTCTG".to_vec(), b"TATTGA".to_vec(), b"ACGTTC".to_vec(), 
        b"ATTTTA".to_vec(), b"TTAATG".to_vec(), b"TTACTT".to_vec(), b"ATGGTT".to_vec(), 
        b"AAAATG".to_vec(), b"TGATGT".to_vec(), b"CTAATG".to_vec(), b"AATTGT".to_vec(), 
        b"TGGTGA".to_vec(), b"GATAAA".to_vec(), b"CCGCTT".to_vec(), b"CTTCGG".to_vec(), 
        b"GCGTTG".to_vec(), b"TTCAAC".to_vec(), b"TTTCGT".to_vec(), b"GTTTCC".to_vec(), 
        b"TCCGCA".to_vec(), b"TAACAC".to_vec(), b"AGTTCA".to_vec(), b"TAAGCA".to_vec(), 
        b"GGACTT".to_vec(), b"CTATGT".to_vec(), b"TTGAGA".to_vec(), b"AAACTG".to_vec(), 
        b"GTGTAC".to_vec(), b"TATTAC".to_vec(), b"GTGCGT".to_vec(), b"CCTCTG".to_vec(), 
        b"CACGCC".to_vec(), b"AGGTTA".to_vec(), b"CGAAGT".to_vec(), b"CGTCAG".to_vec(), 
        b"AAAACC".to_vec(), b"GGCTGA".to_vec(), b"CAAGTT".to_vec(), b"TACATC".to_vec(), 
        b"AGGCCG".to_vec(), b"CGCCAA".to_vec(), b"GGCAAC".to_vec(), b"GAGTGC".to_vec(), 
        b"TTGCCT".to_vec(), b"GAGCAT".to_vec(), b"GTCGCC".to_vec(), b"CTATTC".to_vec(), 
        b"TAATTG".to_vec(), b"TGGACA".to_vec() 
    ];

    // without the next stretch 14433
    // TAATCTCTGGGCATCT
    trigs.push(b"TAATCT".to_vec()); // 14463
    // trigs.push(b"AATCTC".to_vec()); // 14463
    // trigs.push(b"ATCTCT".to_vec()); // 14463

    let shiftable = trigs.iter()
        .filter(|kmer| is_shiftable(&kmap, &kmer[..].into())).count();
    println!("There is {} shiftable kmers.", shiftable);

    let f = BufReader::new(File::open(filename).unwrap());
    let graph = pf::PFGraph::from_fasta(f, &trigs);
    println!("Constructed PFG with:");
    println!("Overlaps:     {:>12}", graph.overlaps());
    println!("#paths:       {:>12}", graph.n_paths());
    println!("#segments:    {:>12}", graph.n_segments());
    println!("Segment size: {:>12}", graph.segment_size());
    println!("Path size:    {:>12}", graph.path_size());
    println!("Total size:   {:>12}", graph.segment_size() + graph.path_size());

    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        let (left, right) = kmap.annotate(seq);
        let matched: Vec<_> = {
            let ac = AhoCorasick::new(&trigs).unwrap();
            ac.find_overlapping_iter(seq).map(|m| m.start()).collect()
        };
        //     find_matches(seq, &trigs);
        println!("{}", record.id());
        println!("{}", from_utf8(record.seq()).unwrap());

        print!("-");
        for x in left { print!("{}", x); }
        for _ in 0..k { print!("-"); }
        println!();

        print!("-");
        for x in right { print!("{}", x); }
        for _ in 0..k { print!("-"); }
        println!();

        let mut j = 0;
        for i in 0..seq.len() {
            match i == matched[j] {
                false => { print!(" "); },
                true => { print!("*"); j += 1; j %= matched.len(); },
            }
        }
        println!();
    }
}
