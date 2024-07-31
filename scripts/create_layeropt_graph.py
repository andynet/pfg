#!/bin/python

from collections import defaultdict


def parse_line(line: str) -> set[str]:
    return set(line.strip().split(", "))


def read_kmer_sets(filename: str) -> dict[int, set[str]]:
    with open(filename) as f:
        lines = f.readlines()
        lines2 = map(parse_line, lines)
        return dict(enumerate(lines2))


# file = "small_seqs/layer_optimal.02mer.only.txt"
file = "small_seqs/layer_optimal.manual02mer.only.txt"
kmer_sets = read_kmer_sets(file)

print("\\begin{tikzpicture}[thick]")
n_: defaultdict[int, int] = defaultdict(lambda: 0)
for i, kmers in kmer_sets.items():
    layer = len(kmers)
    neigh = n_[layer]
    n_[layer] += 1
    y_pos = (neigh+1)//2 * ((-1)**neigh)
    print(f"\t\\node[circle,draw] ({i}) at ({layer*2}, {y_pos*2}) {{{i}}};")

for i, kmers in kmer_sets.items():
    j = 0
    other = kmer_sets[i+j]
    while len(other) < len(kmers) + 2:
        symdiff = kmers.symmetric_difference(other)
        if len(symdiff) == 1:
            print(f"\t\\draw [->] ({i}) to node[sloped,above,scale=0.5] {{AC}} ({i+j});")
        j += 1
        other2 = kmer_sets.get(i+j)
        if not other2:
            break
        other = other2

print("\\end{tikzpicture}")
