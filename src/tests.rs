use std::fs::File;
use std::io::BufReader;

use crate::pf::PFData;
use crate::pf::load_trigs;
use crate::pf::PFGraph;

#[test]
fn iter_works() {
    let pfdata = PFData::from_pfgraph("example/pfg.gfa");

    for (sa, id, pos) in pfdata.iter() {
        println!("{sa}\t{id}\t{pos}");
    }
}

#[test]
fn iter_works_for_pfg_with_overlaps_1() {
    let pfdata = PFData::from_pfgraph("example/pfg1.gfa");

    for (sa, id, pos) in pfdata.iter() {
        println!("{sa}\t{id}\t{pos}");
    }
}

#[test]
fn show_pfdata() {
    let pfdata = PFData::from_pfgraph("example/pfg.gfa");
    pfdata.print();
}

#[test]
fn can_load_arbitrary_graph() {
    let trigs = load_trigs("example/triggers.txt"); 
    let pfdata = PFData::from_graph("example/pangenome.gfa", &trigs);
    pfdata.print();
}

#[test]
fn trig_types_work() {
    let trigs1: Vec<Vec<u8>> = vec![b"AC".to_vec(), b"GT".to_vec()];
    let trigs2: Vec<&[u8]> = vec![b"AC", b"GT"];
    let trigs3: &[Vec<u8>] = &trigs1[..];
    let trigs4: &[&[u8]] = &trigs2[..];

    let f = BufReader::new(File::open("example/pangenome.fna").unwrap());
    let _pfg = PFGraph::from_fasta(f, &trigs1);

    let f = BufReader::new(File::open("example/pangenome.fna").unwrap());
    let _pfg = PFGraph::from_fasta(f, &trigs2);

    let f = BufReader::new(File::open("example/pangenome.fna").unwrap());
    let _pfg = PFGraph::from_fasta(f, trigs3);

    let f = BufReader::new(File::open("example/pangenome.fna").unwrap());
    let _pfg = PFGraph::from_fasta(f, trigs4);
}
