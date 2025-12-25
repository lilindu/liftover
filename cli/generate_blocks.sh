#!/bin/bash

# Script to generate liftover block files for a pair of genomes
# Usage: ./generate_blocks.sh <genomeA.fa> <genomeB.fa> <output_dir>

set -e

GENOME_A=$1
GENOME_B=$2
OUT_DIR=$3
MINIMAP2="/Users/lilindu/miniconda3/pkgs/minimap2-2.30-hba9b596_0/bin/minimap2"
SCRIPT_DIR=$(dirname "$0")

if [ -z "$GENOME_A" ] || [ -z "$GENOME_B" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: $0 <genomeA.fa> <genomeB.fa> <output_dir>"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "=== Generating Blocks for:"
echo "  Genome A: $GENOME_A"
echo "  Genome B: $GENOME_B"
echo "  Output:   $OUT_DIR"

# 1. Align A -> B
echo "Step 1: Aligning A -> B..."
$MINIMAP2 -cx asm5 -c --secondary=no "$GENOME_B" "$GENOME_A" > "$OUT_DIR/A_to_B.paf"

# 2. Convert A->B PAF to Blocks
echo "Step 2: Converting A->B PAF to Blocks..."
python3 "$SCRIPT_DIR/paf_to_blocks.py" "$OUT_DIR/A_to_B.paf" -o "$OUT_DIR/A_to_B.blocks.tsv"

# 3. Align B -> A
echo "Step 3: Aligning B -> A..."
$MINIMAP2 -cx asm5 -c --secondary=no "$GENOME_A" "$GENOME_B" > "$OUT_DIR/B_to_A.paf"

# 4. Convert B->A PAF to Blocks
echo "Step 4: Converting B->A PAF to Blocks..."
python3 "$SCRIPT_DIR/paf_to_blocks.py" "$OUT_DIR/B_to_A.paf" -o "$OUT_DIR/B_to_A.blocks.tsv"

# Cleanup
rm "$OUT_DIR/A_to_B.paf" "$OUT_DIR/B_to_A.paf"

echo "=== Done! Block files generated in $OUT_DIR ==="
