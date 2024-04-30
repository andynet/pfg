use clap::Parser;
use gfa::parser::GFAParser;
use std::str;

/// Count length of segments and paths
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// gfa file
    #[arg(short)]
    gfa: String,
}

fn main() {
    let args = Args::parse();
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(args.gfa).expect("Error parsing GFA file.");

    let mut segment_size = 0;
    for s in &graph.segments {
        segment_size += s.sequence.len();
    }

    let mut path_size = 0;
    for p in &graph.paths {
        let p = str::from_utf8(&p.segment_names).unwrap();
        let p: Vec<_> = p.split(',').collect();
        path_size += p.len();
    }

    println!("segment size = {segment_size}");
    println!("path size = {path_size}");
    println!("total = {}", segment_size + path_size);
}

