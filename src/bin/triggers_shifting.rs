#![allow(unused_imports)]

use aho_corasick::AhoCorasick;
use bio::{data_structures::suffix_array::SuffixArray, io::fasta};
use itertools::Itertools;
use std::fs::File;
use std::io::{BufReader, Write};
use pfg::kmap::{KMap, Kmer};
use pfg::pf;
use std::str::from_utf8;
use pfg::kmap::is_shiftable;

// use rand::{self, seq::IteratorRandom};
// let n = 20;
// let mut rng = rand::thread_rng();
// let v = (&kmap).into_iter().choose_multiple(&mut rng, n);

fn main() {
    let ks = [2, 3, 4, 5, 6, 8, 10, 12, 14, 16];
    let filename = "phiX174/phiX174.2line.clean.fna";
    let output = "phiX174/nonshiftable";
    // let filename = "small_seqs/seqs.fna";
    // let output = "small_seqs/nonshiftable";

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
