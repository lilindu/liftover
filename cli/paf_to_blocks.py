#!/usr/bin/env python3
import argparse


def parse_args():
    """Parse command-line arguments for converting a PAF file to a block mapping TSV."""
    p = argparse.ArgumentParser(
        description="Convert minimap2 PAF (with cg:Z CIGAR) to liftover blocks TSV",
    )
    p.add_argument(
        "paf",
        help="Input PAF file produced by minimap2 with -c to include cg:Z CIGAR",
    )
    p.add_argument(
        "-o",
        "--output",
        default="A_to_B.blocks.tsv",
        help="Output blocks TSV file (default: A_to_B.blocks.tsv)",
    )
    return p.parse_args()


def parse_cigar(cg):
    """Parse a compact CIGAR (cg:Z) string into a list of (op, length) tuples."""
    ops = []
    num = []
    for ch in cg:
        if ch.isdigit():
            num.append(ch)
        else:
            if not num:
                raise ValueError(f"Malformed CIGAR: {cg}")
            ops.append((ch, int("".join(num))))
            num = []
    return ops


def paf_to_blocks(paf_line):
    """Convert one PAF alignment line to a list of mapping blocks.

    Returns a list of dicts: {
        contigA, startA, endA, contigB, startB, endB, strand, mapq
    }
    Only 'M' match segments emit blocks; gaps are absorbed as coordinate advances.
    """
    parts = paf_line.rstrip().split("\t")
    q_name = parts[0]
    q_len = int(parts[1])
    q_start = int(parts[2])
    q_end = int(parts[3])
    strand = parts[4]
    t_name = parts[5]
    t_len = int(parts[6])
    t_start = int(parts[7])
    t_end = int(parts[8])
    mapq = int(parts[11]) if len(parts) > 11 else 0
    cg = None
    for field in parts[12:]:
        if field.startswith("cg:Z:"):
            cg = field[5:]
            break
    if cg is None:
        raise ValueError("PAF line lacks cg:Z (use minimap2 -c)")
    ops = parse_cigar(cg)

    blocks = []
    qa = q_start
    if strand == "+":
        ta = t_start
        for op, ln in ops:
            if op == "M":
                blocks.append({
                    "contigA": q_name,
                    "startA": qa,
                    "endA": qa + ln,
                    "contigB": t_name,
                    "startB": ta,
                    "endB": ta + ln,
                    "strand": "+",
                    "mapq": mapq,
                })
                qa += ln
                ta += ln
            elif op == "I":
                qa += ln
            elif op == "D":
                ta += ln
            else:
                # Treat other ops (e.g., =/X if present) as match-like by splitting
                # Here we conservatively skip unknown ops
                pass
    else:
        ta = t_end
        for op, ln in ops:
            if op == "M":
                blocks.append({
                    "contigA": q_name,
                    "startA": qa,
                    "endA": qa + ln,
                    "contigB": t_name,
                    "startB": ta - ln,
                    "endB": ta,
                    "strand": "-",
                    "mapq": mapq,
                })
                qa += ln
                ta -= ln
            elif op == "I":
                qa += ln
            elif op == "D":
                ta -= ln
            else:
                pass
    return blocks


def write_blocks(paf_path, out_path):
    """Read a PAF file and write a block mapping TSV to out_path."""
    with open(paf_path) as fin, open(out_path, "w") as fout:
        fout.write("contigA\tstartA\tendA\tcontigB\tstartB\tendB\tstrand\tmapq\n")
        for line in fin:
            if not line.strip():
                continue
            bs = paf_to_blocks(line)
            for b in bs:
                fout.write(
                    f"{b['contigA']}\t{b['startA']}\t{b['endA']}\t{b['contigB']}\t{b['startB']}\t{b['endB']}\t{b['strand']}\t{b['mapq']}\n"
                )


def main():
    """Entry point: convert PAF to blocks TSV using cg:Z CIGAR segments."""
    args = parse_args()
    write_blocks(args.paf, args.output)


if __name__ == "__main__":
    main()