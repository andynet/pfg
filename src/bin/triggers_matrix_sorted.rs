use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::stdin;
use std::io::BufWriter;
use std::io::Write;
use bio::io::fasta;
use std::env;

fn main() {
    let k: usize = env::args().nth(1).unwrap().parse().unwrap();
    let outfile: String = env::args().nth(2).unwrap();

    let mut kmap: HashMap<Vec<u8>, Box<[u32; 16]>> = HashMap::new();

    let mut records = fasta::Reader::new(stdin()).records();
    while let Some(Ok(record)) = records.next() {
        let mut seq = record.seq().to_owned();
        seq.push(b'$');
        add_counts(&mut kmap, &seq, k);
    }

    let mut result = Vec::with_capacity(kmap.len());
    for (k, v) in kmap {
        let value = measure(&v);
        result.push((value, k, v));
    }

    result.sort_by(|a, b| a.0.total_cmp(&b.0));
    println!("Number of kmers: {}", result.len());
    println!("Distribution of kmers w.r.t. nonzero contexts: ");
    println!("Distribution of kmers w.r.t. measure: ");

    let mut out = BufWriter::new(File::create(outfile).unwrap());
    for (val, k, v) in result {
        write!(out, "{}\t", std::str::from_utf8(&k).unwrap()).unwrap();
        write!(out, "{:.4}\t", val).unwrap();
        writeln!(out, "{:?}", v).unwrap();
    }
}

fn add_counts(kmap: &mut HashMap<Vec<u8>, Box<[u32; 16]>>, seq: &[u8], k: usize) {
    for i in 1..seq.len()-k-1 {

        let kmer = seq[i..i+k].to_vec();
        let a = seq[i+k]; // after context
        let b = seq[i-1]; // before context
        let idx = if let Some(x) = ptoi((a, b)) { x } else { continue; };

        match kmap.entry(kmer) {
            Entry::Vacant(e) => { 
                let mut v = Box::new([0; 16]);
                v[idx] += 1;
                e.insert(v); 
            }
            Entry::Occupied(mut e) => { 
                let x = e.get_mut(); 
                x[idx] += 1;
            }
        }
    }
}

fn ptoi(pair: (u8, u8)) -> Option<usize> {
    match pair {
        (b'A', b'A') => return Some(0),
        (b'A', b'C') => return Some(1),
        (b'A', b'G') => return Some(2),
        (b'A', b'T') => return Some(3),
        (b'C', b'A') => return Some(4),
        (b'C', b'C') => return Some(5),
        (b'C', b'G') => return Some(6),
        (b'C', b'T') => return Some(7),
        (b'G', b'A') => return Some(8),
        (b'G', b'C') => return Some(9),
        (b'G', b'G') => return Some(10),
        (b'G', b'T') => return Some(11),
        (b'T', b'A') => return Some(12),
        (b'T', b'C') => return Some(13),
        (b'T', b'G') => return Some(14),
        (b'T', b'T') => return Some(15),
        (_, _) => return None,
    }
}

#[test]
fn test_add_kmers() {
    let mut kmap = HashMap::new();
    let seq = b"ACGTGTGCAXGTGCAGATGCTTAGCTTAGCTTAGCCCTAGATATAGCTAGCTAGCTAGCTAG";
    let k = 2;

    add_counts(&mut kmap, seq, k);
    for (k, v) in kmap.iter() {
        print!("{}\t", std::str::from_utf8(k).unwrap());
        print!("{:.4}\t", measure(v));
        println!("{:?}", v);
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
