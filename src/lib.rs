use std::collections::{HashMap, HashSet};

pub mod pf;

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let mut result = vec![0; n];
    for (i, c) in s.iter().enumerate() {
        match *c {
            b'A' => { result[n-1-i] = b'T'; },
            b'C' => { result[n-1-i] = b'G'; },
            b'G' => { result[n-1-i] = b'C'; },
            b'T' => { result[n-1-i] = b'A'; },
            b'N' => { result[n-1-i] = b'N'; },
            _ => { panic!("Unexpected letter in s."); }
        }
    }
    return result;
}

type Kmer = Vec<u8>;
pub fn split_kmers(seq: &[u8], k: usize)
    -> (HashMap<Kmer, (u8, u8)>, HashSet<Kmer>)
{
    let mut nonbranching = HashMap::new();
    let mut branching = HashSet::new();

    for i in 1..seq.len()-k {
        let kmer = seq[i..i+k].to_vec();

        let x = nonbranching.get(&kmer);
        let y = branching.get(&kmer);

        let b = seq[i-1];   // letter before kmer
        let a = seq[i+k];   // letter after kmer

        match (x, y) {
            (None, None)        => { nonbranching.insert(kmer, (b, a)); },
            (Some(tuple), None) => { if *tuple != (b, a) {
                nonbranching.remove(&kmer);
                branching.insert(kmer);
            } },
            (None, Some(_))     => { continue; },
            (Some(_), Some(_))  => { unreachable!(); },
        }
    }
    (nonbranching, branching)
}

#[test]
fn test_split_to_branching_and_nonbranching() {
    let seq: &[u8] = b"ACTCTACCTCA$";
    let (nonbranching, branching) = split_kmers(seq, 2);
    println!("{:?}", nonbranching);
    println!("{:?}", branching);
}

#[cfg(test)]
mod tests;
