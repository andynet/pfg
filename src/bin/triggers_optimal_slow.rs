use itertools::Itertools;
use pfg::pf;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{fs::File, io::{BufReader, BufWriter}, str};
use std::io::Write;

fn main() {
    let k = 2;
    let alph = [b'A', b'C', b'G', b'T'];
    // let filename = "phiX174/phiX174.2line.fna";
    let filename = "small_seqs/seqs.fna";

    let x = vec![alph; k];
    let multi_prod = x.into_iter().multi_cartesian_product();

    (1..=4).collect_vec().par_iter().for_each(|m| {
        let outfile = format!("{}_k{}.{:02}.txt", filename, k, m);
        let mut out = BufWriter::new(File::create(outfile).unwrap());

        for x in multi_prod.clone().combinations(*m) {
            let f = BufReader::new(File::open(filename).unwrap());
            let graph = pf::PFGraph::from_fasta(f, &x);

            let n_seg = graph.n_segments();
            let seg_s = graph.segment_size();
            let path_s = graph.path_size();
            let kmers = x.iter().map(|z| str::from_utf8(z).unwrap()).join(", ");
            writeln!(out, "{:>6}\t{:>6}\t{:>6}\t{:>6}\t{}", 
                n_seg, seg_s, path_s, seg_s + path_s, kmers).unwrap();
        }
    });
}
