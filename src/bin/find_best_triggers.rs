use itertools::Itertools;

fn main() {
    let k = 2;
    let alph = [b'A', b'C', b'G', b'T'];

    let vc = alph.iter()
        .combinations_with_replacement(k).collect::<Vec<_>>();
    println!("{:?}", vc);
    let it = vc.iter().powerset().collect::<Vec<_>>();
    // println!("{:?}", it);
}
