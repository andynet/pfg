# Prefix-free graphs

## Install

## Usage
There are two binary crates and one library crate included in this package. Binary crates serve the purpose of creating prefix-free graphs from FASTA and from GFA formats. Library crate provides an interface for working with the prefix-free graphs, mainly an iterator of suffix array in sublinear space and linear time.

The binary crates are `fasta2pfg` and `gfa2pfg` with usage:

```
fasta2pfg < pangenome.fna > pfg.gfa
```

```
gfa2pfg < pangenome.gfa > pfg.gfa
```

The library can be used from within rust programming language as follows:

```
let pfg = PFG::load("pfg.gfa");

for (i, sa_i) in pfg.iter().enumerate() {
    println!("{}\t{}", i, sa_i);
}
```
