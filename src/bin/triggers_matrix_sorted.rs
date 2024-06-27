use core::panic;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::stdin;
use std::io::BufWriter;
use std::io::Write;
use bio::io::fasta;
use itertools::Itertools;
use std::env;

// const ALPH_SIZE: usize = 4;

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

    let mut heatmap = [[0; 16];16];
    for v in kmap.values() {
        let nzeros = *v.iter().counts().get(&0).unwrap_or(&0);
        let value = measure(v);
        heatmap[nzeros][ftointerval(value)] += 1;
    }
    print_heatmap(&heatmap);

    let mut result = Vec::with_capacity(kmap.len());
    for (k, v) in kmap {
        let nzeros = *v.iter().counts().get(&0).unwrap_or(&0);
        let value = measure(&v);
        result.push((value, nzeros, k, v));
    }
    result.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out = BufWriter::new(File::create(outfile).unwrap());
    for (val, nz, k, v) in result {
        write!(out, "{}\t", std::str::from_utf8(&k).unwrap()).unwrap();
        write!(out, "{:.4}\t", val).unwrap();
        write!(out, "{}\t", nz).unwrap();
        writeln!(out, "{:?}", v).unwrap();
    }
}

fn print_heatmap(heatmap: &[[i32; 16]; 16]) {
    let sum: i32 = heatmap.iter()
        .map(|x| -> i32 { x.iter().sum() })
        .sum();
    println!("Number of kmers: {}", sum);

    print!("Number of 0:   ");
    for j in 0..16 { print!("{:6} ", j); }
    println!();

    for i in 0..15 {
        print!("[{:.2}, {:.2}) = ", 0.25+0.05*(i as f64), 0.25+0.05*((i+1) as f64));
        for hm in heatmap {
            print!("{:6} ", hm[i]);
        }
        println!();
    }

    print!("[{:.2}, {:.2}] = ", 1.0, 1.0);
    for hm in heatmap {
        print!("{:6} ", hm[15]);
    }
    println!();
}

// fn expand_chain(kmap: &mut HashMap<>, kmer: &[u8]) -> (left, right, ch_len, n_occ) {
//     let right = expand_right(kmap: &mut HashMap<>, kmer: &[u8]);
//     let left = expand_left(kmap: &mut HashMap<>, kmer: &[u8]);
//     return (left, right, ch_len, n_occ);
// }

fn add_counts(kmap: &mut HashMap<Vec<u8>, Box<[u32; 16]>>, seq: &[u8], k: usize) {
    for i in 1..seq.len()-k-1 {

        let kmer = seq[i..i+k].to_vec();
        let a = seq[i+k]; // after context
        let b = seq[i-1]; // before context
        let idx = if let Some(x) = ptoi((b, a)) { x } else { continue; };

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

fn ftointerval(v: f64) -> usize { 
    match v {
        x if x < 0.30 => { return  0; },
        x if x < 0.35 => { return  1; },
        x if x < 0.40 => { return  2; },
        x if x < 0.45 => { return  3; },
        x if x < 0.50 => { return  4; },
        x if x < 0.55 => { return  5; },
        x if x < 0.60 => { return  6; },
        x if x < 0.65 => { return  7; },
        x if x < 0.70 => { return  8; },
        x if x < 0.75 => { return  9; },
        x if x < 0.80 => { return 10; },
        x if x < 0.85 => { return 11; },
        x if x < 0.90 => { return 12; },
        x if x < 0.95 => { return 13; },
        x if x < 1.00 => { return 14; },
        x if x == 1.0 => { return 15; },
        _ => { panic!("value is outside interval [0.25, 1.00]"); },
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
        let x = v.iter().counts()[&0];
        println!("{:?}", x);

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
