use bio::io::fasta;
use pfg::kmap::KMap;
use std::{collections::HashSet, fs::File, str::from_utf8};
use std::io::Write;

fn main() {
    // let filename = "middle_seqs/phiX174.2line.clean.fna";
    // let output = "middle_seqs/nonshiftable_v2";
    let filename = "phiX174/phiX174.2line.clean.fna";
    let output = "phiX174/nonshiftable_v2";
    let ks = [2, 3, 4, 5, 6, 8, 10, 12, 14, 16];
    // let ks = [16];

    for k in ks {

    let mut kmap: KMap = KMap::new(k);
    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq().to_owned();
        kmap.add_kmers(&seq);
    }

    let mut result = HashSet::new();
    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        let (mut l, mut r, mut c) = kmap.annotate(seq);
        l.insert(0, 1); l.push(1);
        r.insert(0, 1); r.push(1);
        c.insert(0, 1); c.push(1);

        for i in 1..seq.len()-k {
            let contexts: u8 = c[i].try_into().unwrap();
            if l[i] != contexts || r[i] != contexts {
                result.insert(seq[i..i+k].to_vec());
                // println!("{}", from_utf8(&seq[i..i+k]).unwrap());

                if l[i] != 1 {
                    let mut j = i-1;
                    let mut cont = c[j].try_into().unwrap();
                    result.insert(seq[j..j+k].to_vec());
                    // println!("{}", from_utf8(&seq[j..j+k]).unwrap());
                    while l[j] != 1 && (l[j] != cont || r[j] != cont) {
                        j -= 1;
                        cont = c[j].try_into().unwrap();
                        result.insert(seq[j..j+k].to_vec());
                        // println!("{}", from_utf8(&seq[j..j+k]).unwrap());
                    }
                }

                if r[i] != 1 {
                    let mut j = i+1;
                    let mut cont = c[j].try_into().unwrap();
                    result.insert(seq[j..j+k].to_vec());
                    // println!("{}", from_utf8(&seq[j..j+k]).unwrap());
                    while r[j] != 1 && (l[j] != cont || r[j] != cont) {
                        j += 1;
                        cont = c[j].try_into().unwrap();
                        result.insert(seq[j..j+k].to_vec());
                        // println!("{}", from_utf8(&seq[j..j+k]).unwrap());
                    }
                }
            }
        }
    }

    println!("k = {k}: {}", result.len());
    let out_filename = output.to_string() + &format!(".{k:02}mer.txt");
    let mut f = File::create(out_filename).unwrap();
    for kmer in result {
        writeln!(f, "{}", from_utf8(&kmer).unwrap()).expect("Cannot write to file.");
    }

    }
}
