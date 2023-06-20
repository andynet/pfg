use bio::io::fasta;
use clap::Parser;
use std::collections::HashMap;
use pfg::pf;
use std::io::{stdin, stdout};

/// Build prefix-free graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// trigger file, containing one trigger per line
    #[arg(short)]
    trigger_file: String,
}

fn main() { 
    let args = Args::parse();

    let (trigs, trigs_size) = pf::load_trigs(&args.trigger_file);
    let triggers = pf::get_triggers(&trigs, trigs_size);

    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    let mut records = fasta::Reader::new(stdin()).records();
    while let Some(Ok(record)) = records.next() {
        let mut seq = record.seq().to_owned();
        let v = vec![b'.'; trigs_size];
        seq.extend_from_slice(&v);
        pf::split_prefix_free(&seq, &triggers, &mut segments, &mut paths);
    }

    let (segments, paths) = pf::normalize(segments, paths);

    pf::print_gfa(&segments, &paths, trigs_size, stdout())
        .expect("Error writting GFA");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::stdout;

    #[test]
    fn test_only_prefix_free() {
        let triggers: Vec<&[u8]> = vec![b"T"];
        let k = 1;
        let seq1 = b"ATCTGTTAATG$";
        let seq2 = b"AACGTGTACGTACGAAC$";

        let mut segments = HashMap::new();
        let mut paths = Vec::new();

        pf::split_prefix_free(seq1, &triggers, &mut segments, &mut paths);
        pf::split_prefix_free(seq2, &triggers, &mut segments, &mut paths);

        let (segments, paths) = pf::normalize(segments, paths);

        pf::print_gfa(&segments, &paths, k, stdout()).expect("Error writting GFA");
    }
}
