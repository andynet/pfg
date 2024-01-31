use std::collections::{HashMap, HashSet};
use std::io::stdin;
use bio::io::fasta;
use pfg::split_kmers;
use std::env;

fn main() {
    let k: usize = env::args().nth(1).unwrap().parse().unwrap();

    let mut kmers: HashMap<Vec<u8>, (u8, u8)> = HashMap::new();
    let mut branching_kmers: HashSet<Vec<u8>> = HashSet::new();

    let mut records = fasta::Reader::new(stdin()).records();
    while let Some(Ok(record)) = records.next() {
        let mut seq = record.seq().to_owned();
        seq.push(b'$');
        let (nonbranch, branch) = split_kmers(&seq, k);

        for (k, v) in nonbranch.into_iter() {
            if branching_kmers.contains(&k) { continue; }
            if let Some(tuple) = kmers.get(&k) {
                if *tuple != v { branching_kmers.insert(k); }
                else { continue; }
            } else {
                kmers.insert(k, v);
            }
        }
        branching_kmers = &branching_kmers | &branch;
    }

    for kmer in branching_kmers {
        println!("{}", std::str::from_utf8(&kmer).unwrap());
    }
}
