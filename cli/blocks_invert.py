#!/usr/bin/env python3
import argparse


def parse_args():
    """Parse command-line arguments for inverting A→B blocks TSV to B→A."""
    p = argparse.ArgumentParser(
        description=(
            "Invert A→B blocks TSV into B→A by swapping sides and preserving strand."
        ),
    )
    p.add_argument("blocks_ab", help="Input A→B blocks TSV")
    p.add_argument("-o", "--output", default="B_to_A.blocks.tsv", help="Output B→A blocks TSV")
    return p.parse_args()


def invert_blocks(path_in, path_out):
    """Read A→B blocks and write B→A blocks, sorted per contig by start."""
    by_contig = {}
    with open(path_in) as fin:
        header = fin.readline()
        for line in fin:
            contigA, startA, endA, contigB, startB, endB, strand, mapq = line.rstrip().split("\t")
            b = {
                "contigA": contigB,
                "startA": int(startB),
                "endA": int(endB),
                "contigB": contigA,
                "startB": int(startA),
                "endB": int(endA),
                "strand": strand,
                "mapq": int(mapq),
            }
            by_contig.setdefault(b["contigA"], []).append(b)
    for c in by_contig:
        by_contig[c].sort(key=lambda b: b["startA"])  # ensure sorted
    with open(path_out, "w") as fout:
        fout.write("contigA\tstartA\tendA\tcontigB\tstartB\tendB\tstrand\tmapq\n")
        for c in sorted(by_contig.keys()):
            for b in by_contig[c]:
                fout.write(
                    f"{b['contigA']}\t{b['startA']}\t{b['endA']}\t{b['contigB']}\t{b['startB']}\t{b['endB']}\t{b['strand']}\t{b['mapq']}\n"
                )


def main():
    """Entry point: invert A→B blocks TSV to B→A."""
    args = parse_args()
    invert_blocks(args.blocks_ab, args.output)


if __name__ == "__main__":
    main()