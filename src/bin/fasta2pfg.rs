use clap::Parser;
use pfg::pf;
use std::{fs::File, io::stdin};

/// Build prefix-free graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// trigger file, containing one trigger per line
    #[arg(short)]
    trigger_file: String,

    /// output file
    #[arg(short)]
    output_file: String,
}

fn main() {
    let args = Args::parse();

    let trigs = pf::load_trigs(&args.trigger_file);
    println!("{:?}", trigs);
    let graph = pf::PFGraph::from_fasta(stdin().lock(), &trigs);
    let out = File::create(args.output_file).expect("Cannot open file.");
    graph.write_gfa(out).expect("Error writting GFA");
    println!("Segment size: {}", graph.segment_size());
    println!("Path size: {}", graph.path_size());
}
