use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use bio::io::fasta;
use std::env;
use counter::Counter;

fn main() {
    println!("k\ti\tpfg_size\tseg_size\tpat_size");
    for k in 4..=32 {
        let filename: String = env::args().nth(1).unwrap();
        let f = File::open(&filename).expect("");
        let kmers = collect_kmers(f, k);

        let mut vc: Vec<_> = kmers.iter().map(|(k, v)| (v.len(), v.total::<usize>(), k)).collect();
        vc.sort();
        vc.reverse();
        let vc = vc.iter().map(|x| x.2.to_owned()).collect::<Vec<_>>();

        let mut i = 1;
        let mut best = 1000000000000;
        let mut fails = 0;
        while fails <= 5 && i <= vc.len() {
            let trigs = &vc[..i];
            let f = BufReader::new(File::open(&filename).expect(""));
            let graph = pfg::pf::PFGraph::from_fasta(f, trigs);
            let seg_size = graph.segment_size();
            let pat_size = graph.path_size();
            let size = seg_size + pat_size;
            println!("{k}\t{i}\t{size}\t{seg_size}\t{pat_size}");
            if best > size { best = size; fails = 0 }
            else { fails += 1 }
            i += 1;
        }
    }

    // for x in &vc {
    //     println!("{}\t{}\t{}", x.0, x.1, std::str::from_utf8(x.2).unwrap())
    // }

    // for (k, v) in kmers {
    //     print!("{}\t", std::str::from_utf8(&k).unwrap());
    //     print!("{}\t", v.total::<usize>());
    //     print!("{}\t", v.len());
    //     for (context, count) in &v {
    //         print!("({}, {}): {}, ", context.0 as char, context.1 as char, count);
    //     }
    //     println!();
    // }


}

fn collect_kmers<T: Read>(f: T, k: usize) -> HashMap<Vec<u8>, Counter<(u8, u8)>> {
    let mut kmers: HashMap<Vec<u8>, Counter<(u8, u8), usize>> = HashMap::new();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let mut seq = record.seq().to_owned();
        seq.push(b'$');

        for i in 1..seq.len()-1-k {
            let kmer = &seq[i..i+k];
            let a = seq[i-1];
            let b = seq[i+k+1];

            if let Some(counter) = kmers.get_mut(kmer) {
                counter[&(a, b)] += 1;
            } else {
                let mut counter = Counter::new();
                counter[&(a, b)] = 1;
                kmers.insert(kmer.to_owned(), counter);
            }
        }
    }
    return kmers;
}
