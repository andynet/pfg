use itertools::Itertools;
use pfg::pf;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::io::Write;
use std::sync::{Arc, Mutex};
use std::{fs::File, io::{BufReader, BufWriter}, str};

fn main() {
    let sequences = "phiX174/phiX174.2line.clean.fna";
    let trigger_file = "phiX174/nonshiftable.02mer.txt";
    let k = 2;
    let output_base = "phiX174/exhaustive_stats";

    let triggers = pf::load_trigs(trigger_file).into_iter();

    (1..=16).collect_vec().iter().for_each(|m| {
        println!("Selecting {m} kmers.");
        let output_name = format!("{}_k{}.{:02}.txt", output_base, k, m);
        let output = 
            Arc::new(Mutex::new(BufWriter::new(File::create(output_name).unwrap())));

        // how to avoid collecting?
        let comb = triggers.clone().combinations(*m).collect_vec();
        comb.par_iter().for_each(|x| {
            let f = BufReader::new(File::open(sequences).unwrap());
            let graph = pf::PFGraph::from_fasta(f, x);

            let n_seg = graph.n_segments();
            let seg_s = graph.segment_size();
            let path_s = graph.path_size();
            let kmers = x.iter().map(|z| str::from_utf8(z).unwrap()).join(", ");
            let mut out = output.lock().expect("Failed to lock mutex");
            writeln!(out, "{:>6}\t{:>6}\t{:>6}\t{:>6}\t{}", 
                n_seg, seg_s, path_s, seg_s + path_s, kmers).unwrap();
        })
    });
}
