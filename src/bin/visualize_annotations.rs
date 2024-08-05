use aho_corasick::AhoCorasick;
use bio::io::fasta;
use std::fs::File;
use std::io::BufReader;
use pfg::kmap::KMap;
use pfg::pf::{self, load_trigs};
use std::str::from_utf8;
use pfg::kmap::is_shiftable;

fn main() {
    let sequences = "phiX174/phiX174.2line.clean.fna";
    // let triggers = "phiX174/in_v2_trigs.txt";
    let triggers = "phiX174/nonshiftable_v2.16mer.txt";
    // let triggers = "phiX174/nonshiftable.10mer.manual.txt";
    // let sequences = "middle_seqs/phiX174.2line.clean.fna";
    // let triggers = "middle_seqs/nonshiftable_v2.04mer.txt";
    let trigs = load_trigs(triggers);
    println!("Number of triggers: {}", trigs.len());
    let k = trigs.first().unwrap().len();

    let mut kmap: KMap = KMap::new(k);
    let f = File::open(sequences).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq().to_owned();
        kmap.add_kmers(&seq);
    }

    let shiftable = trigs.iter()
        .filter(|kmer| is_shiftable(&kmap, &kmer[..].into())).count();
    println!("There is {} shiftable kmers.", shiftable);

    let f = BufReader::new(File::open(sequences).unwrap());
    let graph = pf::PFGraph::from_fasta(f, &trigs);
    println!("Constructed PFG with:");
    println!("Overlaps:     {:>12}", graph.overlaps());
    println!("#paths:       {:>12}", graph.n_paths());
    println!("#segments:    {:>12}", graph.n_segments());
    println!("Segment size: {:>12}", graph.segment_size());
    println!("Path size:    {:>12}", graph.path_size());
    println!("Total size:   {:>12}", graph.segment_size() + graph.path_size());

    let f = File::open(sequences).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        let (left, right, contexts) = kmap.annotate(seq);
        let matched: Vec<_> = {
            let ac = AhoCorasick::new(&trigs).unwrap();
            ac.find_overlapping_iter(seq).map(|m| m.start()).collect()
        };
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

        print!("-");
        for x in contexts {
            match x {
                0..=9 => { print!("{}", x); },
                _ => print!("X"),
            }
        }
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

fn measure(m: &[u32; 16]) -> f64 {
    // try all 24 possibilities of independent columns and rows and choose max
    // keywords: assignment problem, Hungarian algorithm
    #[allow(clippy::identity_op)]
    let imax = *[
        m[0+0] + m[4+1] + m[8+2] + m[12+3], m[0+0] + m[4+1] + m[8+3] + m[12+2],
        m[0+0] + m[4+2] + m[8+1] + m[12+3], m[0+0] + m[4+2] + m[8+3] + m[12+1],
        m[0+0] + m[4+3] + m[8+1] + m[12+2], m[0+0] + m[4+3] + m[8+2] + m[12+1],
        m[0+1] + m[4+0] + m[8+2] + m[12+3], m[0+1] + m[4+0] + m[8+3] + m[12+2],
        m[0+1] + m[4+2] + m[8+0] + m[12+3], m[0+1] + m[4+2] + m[8+3] + m[12+0],
        m[0+1] + m[4+3] + m[8+0] + m[12+2], m[0+1] + m[4+3] + m[8+2] + m[12+0],
        m[0+2] + m[4+0] + m[8+1] + m[12+3], m[0+2] + m[4+0] + m[8+3] + m[12+1],
        m[0+2] + m[4+1] + m[8+0] + m[12+3], m[0+2] + m[4+1] + m[8+3] + m[12+0],
        m[0+2] + m[4+3] + m[8+0] + m[12+1], m[0+2] + m[4+3] + m[8+1] + m[12+0],
        m[0+3] + m[4+0] + m[8+1] + m[12+2], m[0+3] + m[4+0] + m[8+2] + m[12+1],
        m[0+3] + m[4+1] + m[8+0] + m[12+2], m[0+3] + m[4+1] + m[8+2] + m[12+0],
        m[0+3] + m[4+2] + m[8+0] + m[12+1], m[0+3] + m[4+2] + m[8+1] + m[12+0],
    ].iter().max().unwrap();
    let msum: u32 = m.iter().sum();
    return f64::from(imax) / f64::from(msum);
}

fn annotate_measure(kmap: &KMap, seq: &[u8]) -> Vec<f64> {
    Vec::new()
}
