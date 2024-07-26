use itertools::Itertools;
use pfg::pf::{load_trigs, PFGraph};
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::io::Write;
use std::sync::{Arc, Mutex};
use std::{fs::File, io::{BufReader, BufWriter}, str};

fn main() {
    let sequences = "phiX174/phiX174.2line.clean.fna";
    let trigger_file = "phiX174/nonshiftable.02mer.txt";
    let min = 1;
    let max = 16;
    let output_base = "phiX174/exhaustive_stats_02mer";

    let triggers = load_trigs(trigger_file);
    (min..=max).collect_vec().iter().for_each(|m| {
        print!("Selecting {m} kmers... ");
        let results = Arc::new(Mutex::new(Vec::new()));
        triggers.iter().combinations(*m).par_bridge().for_each(|x| {
            let f = BufReader::new(File::open(sequences).unwrap());
            let graph = PFGraph::from_fasta(f, &x);

            let n_seg = graph.n_segments();
            let seg_s = graph.segment_size();
            let path_s = graph.path_size();
            let kmers = x.iter().map(|z| str::from_utf8(z).unwrap()).join(", ");
            let mut res = results.lock().expect("Poisoned lock.");
            res.push((n_seg, seg_s, path_s, seg_s + path_s, kmers));
        });

        let results = Arc::into_inner(results).expect("Still more than one reference.");
        let mut results = Mutex::into_inner(results).expect("Poisoned lock");
        results.sort_by(|x, y| x.3.cmp(&y.3));

        let output_name = format!("{}.{:02}.txt", output_base, m);
        let mut output = BufWriter::new(File::create(output_name).unwrap());
        for (n_seg, seg_s, path_s, total_s, kmers) in results {
            writeln!(output, "{:>6}\t{:>6}\t{:>6}\t{:>6}\t{}", 
                n_seg, seg_s, path_s, total_s, kmers).unwrap();
        }
        println!("done!");
    });
}
