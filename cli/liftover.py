#!/usr/bin/env python3
import argparse
import sys


def parse_args():
    """Parse command-line arguments for liftover using a blocks TSV."""
    p = argparse.ArgumentParser(
        description=(
            "Lift coordinates from genome A to B using a precomputed blocks TSV.\n"
            "Input formats: CHR:POS lines (1-based) or BED (0-based half-open, tab-delimited)."
        ),
    )
    p.add_argument("map", help="Blocks TSV file (contigA startA endA contigB startB endB strand mapq)")
    p.add_argument("input", help="Input coordinates file (CHR:POS or BED)")
    p.add_argument("-f", "--format", choices=["chrpos", "bed"], default="chrpos", help="Input format")
    p.add_argument("-o", "--output", default="liftover.out.tsv", help="Output file path")
    p.add_argument("--allow-split", action="store_true", help="Split intervals crossing blocks (BED only)")
    p.add_argument("--strict", action="store_true", help="Reject intervals that cross blocks (BED only)")
    p.add_argument("--stats", action="store_true", help="Print mapping statistics to stderr")
    return p.parse_args()


def load_blocks(path):
    """Load blocks TSV into per-contig arrays sorted by startA."""
    by_contig = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip().split("\t")
            contigA, startA, endA, contigB, startB, endB, strand, mapq = parts
            b = {
                "contigA": contigA,
                "startA": int(startA),
                "endA": int(endA),
                "contigB": contigB,
                "startB": int(startB),
                "endB": int(endB),
                "strand": strand,
                "mapq": int(mapq),
            }
            by_contig.setdefault(contigA, []).append(b)
    for c in by_contig:
        by_contig[c].sort(key=lambda b: b["startA"])  # ensure sorted
    return by_contig


def find_block(blocks, pos):
    """Binary search to find the block containing pos in genome A."""
    lo, hi = 0, len(blocks)
    while lo < hi:
        mid = (lo + hi) // 2
        b = blocks[mid]
        if pos < b["startA"]:
            hi = mid
        elif pos >= b["endA"]:
            lo = mid + 1
        else:
            return mid
    return None


def map_point(block, posA):
    """Map a single position on genome A to genome B using a block."""
    offset = posA - block["startA"]
    if block["strand"] == "+":
        posB = block["startB"] + offset
    else:
        length = block["endA"] - block["startA"]
        posB = block["startB"] + (length - 1 - offset)
    return block["contigB"], posB, block["strand"]


def liftover_chrpos(blocks_by_contig, infile, outfile):
    """Lift CHR:POS lines and write TSV with mapped coordinates and status."""
    total, mapped = 0, 0
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write("contigA\tposA\tcontigB\tposB\tstrand\tstatus\n")
        for line in fin:
            s = line.strip()
            if not s:
                continue
            total += 1
            try:
                contigA, posA = s.split(":")
                posA_one_based = int(posA)
                if posA_one_based < 1:
                    raise ValueError("posA must be >= 1 for CHR:POS format")
                posA_zero_based = posA_one_based - 1
            except Exception:
                fout.write(f"\t\t\t\t\tBAD_INPUT\n")
                continue
            blocks = blocks_by_contig.get(contigA)
            if not blocks:
                fout.write(f"{contigA}\t{posA_one_based}\t\t\t\tNO_CONTIG\n")
                continue
            idx = find_block(blocks, posA_zero_based)
            if idx is None:
                fout.write(f"{contigA}\t{posA_one_based}\t\t\t\tUNMAPPED\n")
                continue
            b = blocks[idx]
            contigB, posB_zero_based, strand = map_point(b, posA_zero_based)
            posB_one_based = posB_zero_based + 1
            mapped += 1
            fout.write(f"{contigA}\t{posA_one_based}\t{contigB}\t{posB_one_based}\t{strand}\tOK\n")
    return total, mapped


def liftover_bed(blocks_by_contig, infile, outfile, allow_split=False, strict=False):
    """Lift BED intervals and write TSV with mapped intervals and status."""
    total, mapped, split = 0, 0, 0
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write("contigA\tstartA\tendA\tcontigB\tstartB\tendB\tstrand\tstatus\n")
        for line in fin:
            if not line.strip():
                continue
            total += 1
            parts = line.rstrip().split("\t")
            contigA, startA, endA = parts[:3]
            startA, endA = int(startA), int(endA)
            blocks = blocks_by_contig.get(contigA)
            if not blocks:
                fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\tNO_CONTIG\n")
                continue
            idx = find_block(blocks, startA)
            if idx is None:
                fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\tUNMAPPED_START\n")
                continue
            b = blocks[idx]
            if endA <= b["endA"]:
                contigB, startB, strand = map_point(b, startA)
                _, endB, _ = map_point(b, endA - 1)
                fout.write(f"{contigA}\t{startA}\t{endA}\t{contigB}\t{startB}\t{endB+1}\t{strand}\tOK\n")
                mapped += 1
            else:
                if strict and not allow_split:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\tCROSSES_BLOCK\n")
                    continue
                if not allow_split:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\tCROSSES_BLOCK\n")
                    continue
                pieces = []
                a = startA
                while a < endA:
                    idx = find_block(blocks, a)
                    if idx is None:
                        pieces.append((None, None, None, a, min(endA, a + 1)))
                        a += 1
                        continue
                    b = blocks[idx]
                    a_end = min(endA, b["endA"])
                    contigB, b_startB, strand = map_point(b, a)
                    _, b_endB_minus1, _ = map_point(b, a_end - 1)
                    pieces.append((contigB, b_startB, b_endB_minus1 + 1, a, a_end))
                    a = a_end
                for seg in pieces:
                    if seg[0] is None:
                        fout.write(f"{contigA}\t{seg[3]}\t{seg[4]}\t\t\t\t\tUNMAPPED_SEG\n")
                    else:
                        fout.write(f"{contigA}\t{seg[3]}\t{seg[4]}\t{seg[0]}\t{seg[1]}\t{seg[2]}\t{strand}\tSPLIT\n")
                        split += 1
                mapped += 1
    return total, mapped, split


def main():
    """Entry point: perform liftover for CHR:POS or BED using blocks TSV."""
    args = parse_args()
    blocks = load_blocks(args.map)
    if args.format == "chrpos":
        total, mapped = liftover_chrpos(blocks, args.input, args.output)
        if args.stats:
            print(f"total={total} mapped={mapped}", file=sys.stderr)
    else:
        total, mapped, split = liftover_bed(blocks, args.input, args.output, args.allow_split, args.strict)
        if args.stats:
            print(f"total={total} mapped={mapped} split={split}", file=sys.stderr)


if __name__ == "__main__":
    main()
