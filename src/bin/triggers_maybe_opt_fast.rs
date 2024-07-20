use counter::Counter;
use std::{collections::HashSet, fmt::{Debug, Display}, fs::File, str::from_utf8};
use bio::io::fasta;


fn main() {
    let filename = "phiX174/phiX174.2line.fna";
    let f = File::open(filename).expect("Unable to read file.");
    let seqs = fasta::Reader::new(f)
        .records()
        .map(|x| x.expect("Incorrect fasta record.").seq().to_vec());

    // let seqs = [
    //     b"ACGTGCGTGCA".to_vec(),
    //     b"ACACACACACA".to_vec(),
    //     b"ACGTGTGCGCGATATGAG".to_vec()
    // ].into_iter();

    let k = 6;
    let mut count = count_kmers(seqs, k);
    println!("{}", count.len());
    // let v: Vec<_> = count.iter()
    //     .map(|(k, v)| format!("({}, {})", from_utf8(k).unwrap(), v))
    //     .collect();
    // println!("{}", v.join(", "));
    let f = File::open(filename).expect("Unable to read file.");
    let seqs = fasta::Reader::new(f)
        .records()
        .map(|x| x.expect("Incorrect fasta record.").seq().to_vec());


    let mut state = PFGState::from_stringset(seqs, k);
    // println!("{}", state);

    // let mut kmers = Vec::new();
    let n = count.len();
    for i in 1..=n {
        let kmer = find_locally_best(&state, &count);
        println!("{}: {}", i, from_utf8(&kmer).unwrap());

        state = state.break_at_kmer(&kmer, count[&kmer]);
        count.remove(&kmer);
        // println!("{}", state);
        // kmers.push(kmer);
    }
    // println!("{}", state);
    // let v: Vec<_> = kmers.iter().map(|x| from_utf8(x).unwrap()).collect();
    // println!("{}", v.join(", "));
}

fn find_locally_best(state: &PFGState, counts: &Counter<Vec<u8>, usize>) -> Vec<u8> {
    let mut best_pfg_size = usize::MAX;
    let mut best_kmer = Vec::new();
    for (kmer, &n_occ) in counts {
        let new_state = state.break_at_kmer(kmer, n_occ);
        let new_pfg_size = new_state.pfg_size();
        if new_pfg_size < best_pfg_size {
            best_pfg_size = new_pfg_size;
            best_kmer = kmer.to_owned();
        }
    }
    return best_kmer;
}

fn count_kmers<T>(seqs: T, k: usize) -> Counter<Vec<u8>, usize> 
    where T: Iterator<Item = Vec<u8>>
{
    let mut count = Counter::new();
    for s in seqs {
        for i in 1..s.len()-k {
            let kmer = s[i..i+k].to_vec();
            count[&kmer] += 1;
        }
    }
    return count;
}

#[derive(Debug, Default, Clone)]
struct PFGState {
    segments: HashSet<Vec<u8>>,
    path_size: usize,
    k : usize
}

impl PFGState {
    fn from_stringset<T>(ss: T, k: usize) -> Self 
        where T: Iterator<Item = Vec<u8>>
    {
        let mut segments = HashSet::new();
        let mut path_size = 0;
        for mut s in ss {
            s.extend_from_slice(&vec![b'.'; k]);
            segments.insert(s);
            path_size += 1;
        }
        return Self { segments, path_size, k };
    }

    fn break_at_kmer(&self, kmer: &[u8], n_occ: usize) -> Self {
        assert!(kmer.len() == self.k, "Kmer has different length");
        let k = self.k;
        use aho_corasick::AhoCorasick;
        let ac = AhoCorasick::new([kmer]).unwrap();

        let mut segments = HashSet::new();
        for s in &self.segments {
            let mut i = 0;
            let mut j;

            for m in ac.find_overlapping_iter(&s[1..]) {
                j = m.end()+1;
                segments.insert(s[i..j].to_vec());
                i = m.start()+1;
            }
            j = s.len();
            segments.insert(s[i..j].to_vec());
        }
        let path_size = self.path_size + n_occ;

        Self { segments, path_size, k }
    }

    fn segment_size(&self) -> usize {
        self.segments.iter().map(|x| x.len()).sum()
    }

    fn pfg_size(&self) -> usize {
        self.path_size + self.segment_size()
    }
}

impl Display for PFGState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{} {} {} {}", 
            self.k, self.path_size, self.segment_size(), self.pfg_size()
        )?;

        let mut v: Vec<_> = self.segments.iter().collect();
        v.sort();
        for x in v { write!(f, "{} ", from_utf8(x).unwrap())?; }

        Ok(())
    }
}

