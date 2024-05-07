use clap::Parser;
use gfa::parser::GFAParser;
use gfa::gfa::GFA;
use pfg::pf::{reconstruct_path, split_prefix_free2, PFGraph};
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::str;
use std::collections::HashMap;

/// Find triggers related to arbitrary graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// arbitrary graph
    #[arg(short)]
    arb_gfa: String,

    /// output base
    #[arg(short)]
    output: String,
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(args.arb_gfa).expect("Error parsing GFA file.");

    for k in (8..=32).step_by(2) {
        let triggers = get_triggers(&graph, k);
        let (old_breaks, new_breaks) = calculate_breaks(&graph, k, &triggers);

        let mut segments = HashMap::new();
        let mut paths = Vec::new();

        let triggers: Vec<_> = triggers.into_iter().collect();
        for p in &graph.paths {
            let (mut seq, _) = reconstruct_path(p, &graph);
            let v = vec![b'.'; k];
            seq.extend_from_slice(&v);
            split_prefix_free2(&seq, &triggers, &mut segments, &mut paths);
        }
        let segments: Vec<Vec<u8>> = segments.into_iter().map(|x| x.0).collect();
        let result = PFGraph::new(k, segments, paths);

        let ssize = result.segment_size();
        let psize = result.path_size();
        println!("k={k}\n\tall_breaks={old_breaks}\n\tnew_breaks={new_breaks}");
        println!("\tseg_size={ssize}\n\tpath_size={psize}");

        let output_file = args.output.clone() + &format!("_{:02}.txt", k);
        let mut out = File::create(output_file).expect("Cannot open file.");
        for trig in triggers {
            writeln!(out, "{}", &str::from_utf8(&trig).unwrap())
                .expect("Failed to write trigger words.");
        }
    }
}

fn get_triggers(graph: &GFA<usize, ()>, k: usize) -> HashSet<Vec<u8>> {
    let mut triggers = HashSet::new();
    for p in &graph.paths {
        let (mut seq, breaks) = reconstruct_path(p, graph);
        let v = vec![b'.'; k];
        seq.extend_from_slice(&v);
        for i in breaks {
            triggers.insert(seq[i..i+k].to_vec());
        }
    }
    return triggers;
}

fn calculate_breaks(graph: &GFA<usize, ()>, k: usize, triggers: &HashSet<Vec<u8>>)
    -> (usize, usize)
{
    let mut new_breaks = 0;
    let mut old_breaks = 0;
    for p in &graph.paths {
        let (mut seq, breaks) = reconstruct_path(p, graph);
        let v = vec![b'.'; k];
        seq.extend_from_slice(&v);
        let mut j = 0;
        for i in 0..seq.len()-k {
            let is_brk = {
                if i == breaks[j] { j+=1; true }
                else { false }
            };

            let is_trig = {
                let kmer = seq[i..i+k].to_vec();
                triggers.get(&kmer).is_some()
            };

            match (is_brk, is_trig) {
                (true, true) => { old_breaks += 1;},
                (true, false) => { unreachable!() },
                (false, true) => { new_breaks += 1; },
                (false, false) => { /* good */ },
            }
        }
    }
    return (old_breaks, new_breaks);
}

#[test]
fn test_reconstruct_seq() {
    use pfg::pf::reconstruct_path;
    use std::str;

    let filename = "./example/pangenome.gfa";
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename).expect("Error parsing GFA file.");

    for p in &graph.paths {
        let (seq, breaks) = reconstruct_path(p, &graph);
        println!("{}", str::from_utf8(&seq).unwrap());
        println!("{:?}", breaks);
    }
}
