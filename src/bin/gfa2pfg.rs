use clap::Parser;
use gfa::parser::GFAParser;
use gfa::gfa::GFA;
use pfg::pf::reconstruct_path;
use std::collections::HashSet;
use std::str;

/// Build prefix-free graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    ///
    #[arg(short)]
    arb_gfa: String,
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(args.arb_gfa).expect("Error parsing GFA file.");

    for k in 1..65 {
        let (set, new_breaks) = is_pfg(&graph, k);
        println!("For {k} we only need to add {new_breaks} new breaks to make it PFG.");
        println!("Triggers: {} {:?}", set.len(), set);

    }
}

fn is_pfg(graph: &GFA<usize, ()>, k: usize) -> (HashSet<Vec<u8>>, usize) {
    let mut triggers = HashSet::new();
    let mut new_breaks = 0;

    for p in &graph.paths {
        let (mut seq, breaks) = reconstruct_path(p, graph);
        let v = vec![b'.'; k];
        seq.extend_from_slice(&v);
        for i in breaks {
            triggers.insert(seq[i..i+k].to_vec());
        }
    }

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
                (true, true) => { /* good */ },
                (true, false) => { unreachable!() },
                (false, true) => { new_breaks += 1; },
                (false, false) => { /* good */ },
            }
        }
    }
    return (triggers, new_breaks);
}

#[test]
fn test_is_pfg() {
    let filename = "./example/pangenome.gfa";
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename).expect("Error parsing GFA file.");

    // assert!(!is_pfg(&graph, 1));
    // assert!(is_pfg(&graph, 2));
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
