use std::collections::{hash_map, hash_map::Entry, HashMap};
use std::fmt::{Debug, Display, Formatter, Result};
use std::str::from_utf8;

#[allow(non_camel_case_types)]
type alph = u8;

#[derive(PartialEq, Eq, Hash, Clone)]
pub struct Kmer(Vec<alph>);

impl AsRef<[alph]> for Kmer {
    fn as_ref(&self) -> &[alph] { &self.0 }
}

impl From<&[u8]> for Kmer {
    fn from(v: &[u8]) -> Self { Kmer(v.to_vec()) }
}

impl Debug for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "Kmer({})", from_utf8(&self.0).unwrap())?;
        Ok(())
    }
}

impl Display for Kmer {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}", from_utf8(&self.0).unwrap())?;
        Ok(())
    }
}

#[derive(Clone)]
pub struct KMap{
    map: HashMap<Kmer, Box<[u32; 16]>>,
    k: usize,
}

impl KMap {
    pub fn new(k: usize) -> Self {
        Self { map: HashMap::new(), k }
    }

    pub fn len(&self) -> usize { self.map.len() }
    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn entry(&mut self, key: Kmer) -> Entry<Kmer, Box<[u32; 16]>> {
        self.map.entry(key)
    }

    pub fn get(&self, k: &Kmer) -> Option<Box<[u32; 16]>>{
        self.map.get(k).cloned()
    }

    /// seq is a slice from alphabet 0=$,1=A,2=C,3=G,4=T
    pub fn add_kmers(&mut self, seq: &[alph]) {
        let k = self.k;
        for i in 1..seq.len()-k {

            let kmer = Kmer(seq[i..i+k].to_vec());
            let context = (seq[i-1], seq[i+k]);
            let idx = Self::try_ptoi(context).expect("Unexpected context");

            let map = &mut self.map;
            match map.entry(kmer) {
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

    /// returns all kmers to the left in the de Bruijn graph
    pub fn to_right(&self, kmer: &Kmer) -> Vec<Kmer> {
        let mut result = Vec::new();
        let m = match self.get(kmer) {
            None => return result,
            Some(m) => m,
        };

        let mut new_kmer: Vec<alph> = kmer.as_ref()[1..].to_vec();
        let v = [(0, b'A'), (1, b'C'), (2, b'G'), (3, b'T')];
        for (col, char) in v {
            if m[col] + m[col + 4] + m[col + 8] + m[col + 12] > 0 {
                new_kmer.push(char);
                result.push(Kmer(new_kmer.clone()));
                new_kmer.pop();
            }
        }
        return result;
    }

    /// returns all kmers to the left in the de Bruijn graph
    pub fn to_left(&self, kmer: &Kmer) -> Vec<Kmer> {
        let mut result = Vec::new();
        let m = match self.get(kmer) {
            None => return result,
            Some(m) => m,
        };

        let k = self.k;
        let v = [(0, b'A'), (1, b'C'), (2, b'G'), (3, b'T')];
        for (row, char) in v {
            if m[row*4] + m[row*4 + 1] + m[row*4 + 2] + m[row*4 + 3] > 0 {
                let mut new_kmer = vec![char];
                new_kmer.extend_from_slice(&kmer.as_ref()[..k-1]);
                result.push(Kmer(new_kmer));
            }
        }
        return result;
    }

    /// returns 2 vectors representing for each position of the sequence
    /// how many contexts is to the left and to the right of the kmer
    /// at that position
    pub fn annotate(&self, seq: &[alph]) -> (Vec<u8>, Vec<u8>){
        let mut left = Vec::new();
        let mut right = Vec::new();
        let k = self.k;
        for i in 1..seq.len()-k {
            let kmer = &seq[i..i+k].into();
            left.push(self.to_left(kmer).len().try_into().unwrap());
            right.push(self.to_right(kmer).len().try_into().unwrap());
        }
        return (left, right);
    }

    fn try_ptoi(context: (alph, alph)) -> Option<usize> {
        match context {
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
}

impl<'a> IntoIterator for &'a KMap {
    type Item = (&'a Kmer, &'a Box<[u32; 16]>);
    type IntoIter = hash_map::Iter<'a, Kmer, Box<[u32; 16]>>;

    fn into_iter(self) -> Self::IntoIter {
        self.map.iter()
    }
}

impl IntoIterator for KMap {
    type Item = (Kmer, Box<[u32; 16]>);
    type IntoIter = hash_map::IntoIter<Kmer, Box<[u32; 16]>>;

    fn into_iter(self) -> Self::IntoIter {
        self.map.into_iter()
    }
}

impl FromIterator<(Kmer, Box<[u32; 16]>)> for KMap {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = (Kmer, Box<[u32; 16]>)>
    {
        let mut map = HashMap::new();
        map.extend(iter);
        let (kmer, _) = map.iter().next()
            .expect("Cannot collect KMap from empty iterator.");
        let k = kmer.0.len();
        KMap{map, k}
    }
}

#[test]
fn test_basic_functionality() {
    let k = 2;
    let mut kmap = KMap::new(k);
    let seq = b"AGTCTCGATC";
 
    let records = [b"AGTCTCGATC"].iter();
    for record in records {
        let seq = record.to_vec();
        // seq.push(b'$');
        kmap.add_kmers(&seq);
    }

    for (k, v) in &kmap {
        println!("{:?} {:?}", k, v);
    }

    let new_kmer = kmap.to_right(&b"TC"[..].into());
    println!("{:?}", new_kmer);
    let new_kmer = kmap.to_left(&b"TC"[..].into());
    println!("{:?}", new_kmer);

    let (left, right) = kmap.annotate(seq);
    println!("{}", from_utf8(seq).unwrap());
    for x in left { print!("{}", x); }
    println!();
    for x in right { print!("{}", x); }
    println!();

    let x: KMap = kmap.clone().into_iter().collect();

    for (k, v) in x {
        println!("{} {:?}", k, v);
    }
}

#[test]
fn exp_see_contexts() {
    use bio::io::fasta;
    use std::fs::File;

    let k = 12;
    let filename = "phiX174/phiX174.2line.clean.fna";

    let mut kmap: KMap = KMap::new(k);
    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq().to_owned();
        kmap.add_kmers(&seq);
    }

    let f = File::open(filename).unwrap();
    let mut records = fasta::Reader::new(f).records();
    while let Some(Ok(record)) = records.next() {
        let seq = record.seq();
        let (left, right) = kmap.annotate(seq);
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
    }
}
