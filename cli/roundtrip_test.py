#!/usr/bin/env python3
import argparse
import sys


def parse_args():
    """Parse command-line arguments for A→B→A round-trip validation."""
    p = argparse.ArgumentParser(
        description=(
            "Validate round-trip liftover (A→B→A) for CHR:POS (1-based) or BED (0-based half-open) inputs using two blocks TSV files."
        ),
    )
    p.add_argument("--blocks-ab", required=True, help="Blocks TSV for A→B")
    p.add_argument("--blocks-ba", required=True, help="Blocks TSV for B→A")
    p.add_argument("--input", required=True, help="Input coordinates file (CHR:POS or BED)")
    p.add_argument("--format", choices=["chrpos", "bed"], default="chrpos", help="Input format")
    p.add_argument("--allow-split", action="store_true", help="Split intervals crossing blocks (BED only)")
    p.add_argument("--strict", action="store_true", help="Reject intervals that cross blocks (BED only)")
    p.add_argument("--out", default="roundtrip.out.tsv", help="Output TSV path")
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


def liftover_point_once(blocks_by_contig, contig, pos):
    """Perform one liftover mapping for a CHR:POS; return (status, contig, pos, strand)."""
    bs = blocks_by_contig.get(contig)
    if not bs:
        return ("NO_CONTIG", None, None, None)
    idx = find_block(bs, pos)
    if idx is None:
        return ("UNMAPPED", None, None, None)
    cB, pB, s = map_point(bs[idx], pos)
    return ("OK", cB, pB, s)


def roundtrip_chrpos(blocks_ab, blocks_ba, infile, outfile):
    """Run A→B→A round-trip on CHR:POS inputs and write results."""
    total, pass_n, fail_n = 0, 0, 0
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write("contigA\tposA\tcontigB\tposB\tcontigA2\tposA2\tstatus\n")
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
                fout.write("\t\t\t\t\t\tBAD_INPUT\n")
                fail_n += 1
                continue
            st1, cB, pB, _ = liftover_point_once(blocks_ab, contigA, posA_zero_based)
            if st1 != "OK":
                fout.write(f"{contigA}\t{posA_one_based}\t\t\t\t\t{st1}\n")
                fail_n += 1
                continue
            st2, cA2, pA2, _ = liftover_point_once(blocks_ba, cB, pB)
            if st2 != "OK":
                fout.write(f"{contigA}\t{posA_one_based}\t{cB}\t{pB + 1}\t\t\t{st2}\n")
                fail_n += 1
                continue
            status = "PASS" if (cA2 == contigA and pA2 == posA_zero_based) else "FAIL"
            if status == "PASS":
                pass_n += 1
            else:
                fail_n += 1
            fout.write(f"{contigA}\t{posA_one_based}\t{cB}\t{pB + 1}\t{cA2}\t{pA2 + 1}\t{status}\n")
    print(f"total={total} pass={pass_n} fail={fail_n}", file=sys.stderr)


def roundtrip_bed(blocks_ab, blocks_ba, infile, outfile, allow_split=False, strict=False):
    """Run A→B→A round-trip on BED intervals. Strict requires single-block mapping; split allows piecewise."""
    total, pass_n, fail_n = 0, 0, 0
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write("contigA\tstartA\tendA\tcontigB\tstartB\tendB\tcontigA2\tstartA2\tendA2\tstatus\n")
        for line in fin:
            if not line.strip():
                continue
            total += 1
            contigA, startA, endA = line.rstrip().split("\t")[:3]
            startA, endA = int(startA), int(endA)
            blocks = blocks_ab.get(contigA)
            if not blocks:
                fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\t\t\tNO_CONTIG\n")
                fail_n += 1
                continue
            # For strict: require entire interval in one block
            if strict and not allow_split:
                idx = find_block(blocks, startA)
                if idx is None or endA > blocks[idx]["endA"]:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\t\t\tCROSSES_BLOCK\n")
                    fail_n += 1
                    continue
                b = blocks[idx]
                cB, sB, _ = map_point(b, startA)
                _, eBminus1, _ = map_point(b, endA - 1)
                eB = eBminus1 + 1
                # Back to A (strict in BA block set)
                blocks2 = blocks_ba.get(cB)
                if not blocks2:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t{cB}\t{sB}\t{eB}\t\t\t\tNO_CONTIG_BA\n")
                    fail_n += 1
                    continue
                idx2 = find_block(blocks2, sB)
                if idx2 is None or eB > blocks2[idx2]["endA"]:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t{cB}\t{sB}\t{eB}\t\t\t\tCROSSES_BLOCK_BA\n")
                    fail_n += 1
                    continue
                b2 = blocks2[idx2]
                cA2, sA2, _ = map_point(b2, sB)
                _, eA2minus1, _ = map_point(b2, eB - 1)
                eA2 = eA2minus1 + 1
                status = "PASS" if (cA2 == contigA and sA2 == startA and eA2 == endA) else "FAIL"
                if status == "PASS":
                    pass_n += 1
                else:
                    fail_n += 1
                fout.write(f"{contigA}\t{startA}\t{endA}\t{cB}\t{sB}\t{eB}\t{cA2}\t{sA2}\t{eA2}\t{status}\n")
            else:
                # Split mode: piecewise A→B, then each piece back B→A; union must match original
                piecesB = []
                a = startA
                while a < endA:
                    idx = find_block(blocks, a)
                    if idx is None:
                        piecesB.append((None, None, None, a, min(endA, a + 1)))
                        a += 1
                        continue
                    b = blocks[idx]
                    a_end = min(endA, b["endA"])
                    cB, sB, _ = map_point(b, a)
                    _, eBminus1, _ = map_point(b, a_end - 1)
                    piecesB.append((cB, sB, eBminus1 + 1, a, a_end))
                    a = a_end
                # Map each B piece back
                reconstructed = []
                failed = False
                merged = []
                for cB, sB, eB, aS, aE in piecesB:
                    if cB is None:
                        failed = True
                        continue
                    blocks2 = blocks_ba.get(cB)
                    if not blocks2:
                        failed = True
                        continue
                    idx2 = find_block(blocks2, sB)
                    if idx2 is None or eB > blocks2[idx2]["endA"]:
                        failed = True
                        continue
                    b2 = blocks2[idx2]
                    cA2, sA2, _ = map_point(b2, sB)
                    _, eA2minus1, _ = map_point(b2, eB - 1)
                    reconstructed.append((cA2, sA2, eA2minus1 + 1))
                # Simple check: if any piece failed, mark fail; else ensure concatenation equals original interval on same contig
                if failed:
                    status = "FAIL"
                    fail_n += 1
                else:
                    # Merge adjacent pieces
                    reconstructed.sort(key=lambda x: x[1])
                    merged = []
                    for seg in reconstructed:
                        if not merged or seg[0] != merged[-1][0] or seg[1] != merged[-1][2]:
                            merged.append(list(seg))
                        else:
                            merged[-1][2] = seg[2]
                    if len(merged) == 1 and merged[0][0] == contigA and merged[0][1] == startA and merged[0][2] == endA:
                        status = "PASS"
                        pass_n += 1
                    else:
                        status = "FAIL"
                        fail_n += 1
                # For reporting, show the first B piece bounds if available, else blanks
                if piecesB and piecesB[0][0] is not None:
                    cB0, sB0, eB0, _, _ = piecesB[0]
                    if merged:
                        cA2, sA2, eA2 = merged[0]
                    else:
                        cA2 = sA2 = eA2 = ""
                    fout.write(f"{contigA}\t{startA}\t{endA}\t{cB0}\t{sB0}\t{eB0}\t{cA2}\t{sA2}\t{eA2}\t{status}\n")
                else:
                    fout.write(f"{contigA}\t{startA}\t{endA}\t\t\t\t\t\t\t{status}\n")
    print(f"total={total} pass={pass_n} fail={fail_n}", file=sys.stderr)


def main():
    """Entry point: run A→B→A round-trip validation for points or BED."""
    args = parse_args()
    blocks_ab = load_blocks(args.blocks_ab)
    blocks_ba = load_blocks(args.blocks_ba)
    if args.format == "chrpos":
        roundtrip_chrpos(blocks_ab, blocks_ba, args.input, args.out)
    else:
        roundtrip_bed(blocks_ab, blocks_ba, args.input, args.out, args.allow_split, args.strict)


if __name__ == "__main__":
    main()
