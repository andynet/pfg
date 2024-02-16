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
    let graph = parser.parse_file(&args.arb_gfa).expect("Error parsing GFA file.");

    let mut k = 1;
    while ! is_pfg(&graph, k) {
        println!("For {} {} is not PFG.", k, args.arb_gfa);
        k *= 2;
    }
    println!("For {} {} is a PFG", k, args.arb_gfa);
    // let k = binary_search(range, is_pfg, args.arb_gfa);
}

fn is_pfg(graph: &GFA<usize, ()>, k: usize) -> bool {
    let mut triggers = HashSet::new();
    let mut noopers = HashSet::new();
    let mut num_breaks = 0;

    for p in &graph.paths {
        let (mut seq, breaks) = reconstruct_path(p, graph);
        let v = vec![b'.'; k];
        seq.extend_from_slice(&v);
        let mut j = 0;
        for i in 0..seq.len()-k {
            let brk = {
                if i == breaks[j] { j+=1; true }
                else { false }
            };
            let kmer = seq[i..i+k].to_vec();

            let trig = triggers.get(&kmer);
            let noop = noopers.get(&kmer);

            match (trig, noop) {
                (None, None) => {
                    if brk { triggers.insert(kmer); }
                    else { noopers.insert(kmer); }
                },
                (Some(x), None) => { if ! brk {
                    // println!("Should be noop: {:?}", x);
                    // println!("Trig: {:?}", triggers);
                    // println!("Noop: {:?}", noopers);
                    // return false; 
                    println!("Introducing a break.");
                    num_breaks += 1;
                } },
                (None, Some(x)) => { if   brk {
                    println!("Should be in trig: {:?}", x);
                    println!("Trig: {:?}", triggers);
                    println!("Noop: {:?}", noopers);
                    return false;
                } },
                (Some(_), Some(_)) => { unreachable!() }
            }
        }
    }
    println!("Trig: {} {:?}", triggers.len(), triggers);
    println!("Introduced {} breaks.", num_breaks);
    return true;
}

#[test]
fn test_is_pfg() {
    let filename = "./example/pangenome.gfa";
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename).expect("Error parsing GFA file.");

    assert!(!is_pfg(&graph, 1));
    assert!(is_pfg(&graph, 2));
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
