use crate::pf::PFData;

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
    let pfdata = PFData::from_graph(
        "example/pangenome.gfa", "example/triggers.txt"
    );
    pfdata.print();
}
