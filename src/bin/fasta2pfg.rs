use clap::Parser;
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

    let trigs = pf::load_trigs(&args.trigger_file);
    let graph = pf::PFGraph::from_fasta(stdin().lock(), &trigs);
    graph.write_gfa(stdout()).expect("Error writting GFA");
}
