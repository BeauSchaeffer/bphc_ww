#!/bin/bash
# check_lineage_mutations.sh
# Usage: bash check_lineage_mutations.sh <lineage> <variants_file>
# Example: bash check_lineage_mutations.sh BA.5 freyja_variants/60975000_variants.tsv

BARCODES="/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja/lib/python3.11/site-packages/freyja/data/usher_barcodes.csv"

# --- input handling ---
if [ $# -ne 2 ]; then
    echo "Usage: bash check_lineage_mutations.sh <lineage> <variants_file>"
    echo "Example: bash check_lineage_mutations.sh BA.5 freyja_variants/60975000_variants.tsv"
    exit 1
fi

LINEAGE=$1
VARIANTS=$2

# --- validate inputs ---
if [ ! -f "$VARIANTS" ]; then
    echo "Error: variants file not found: $VARIANTS"
    exit 1
fi

if ! grep -q "^${LINEAGE}," $BARCODES; then
    echo "Error: lineage $LINEAGE not found in barcode file"
    exit 1
fi

# --- extract defining mutations ---
MUT_FILE=$(mktemp /tmp/lineage_muts.XXXX)
grep "^${LINEAGE}," $BARCODES | tr ',' '\n' | \
    paste - <(head -1 $BARCODES | tr ',' '\n') | \
    awk -F'\t' '$1 == "1.0" {print $2}' > $MUT_FILE

TOTAL_DEFINING=$(wc -l < $MUT_FILE)

# --- reconstruct mutation format from variants file and find matches ---
FOUND_FILE=$(mktemp /tmp/found_muts.XXXX)
awk 'NR>1 {mut=$3$2$4; print mut"\t"$0}' $VARIANTS | \
    grep -Ff $MUT_FILE | \
    cut -f2- > $FOUND_FILE

TOTAL_FOUND=$(wc -l < $FOUND_FILE)

# --- find missing mutations ---
FOUND_MUTS=$(awk 'NR>1 {print $3$2$4}' $VARIANTS)
MISSING_FILE=$(mktemp /tmp/missing_muts.XXXX)
while read mut; do
    if ! echo "$FOUND_MUTS" | grep -q "^${mut}$"; then
        echo "$mut"
    fi
done < $MUT_FILE > $MISSING_FILE

TOTAL_MISSING=$(wc -l < $MISSING_FILE)

# --- output ---
echo "========================================"
echo " Lineage:   $LINEAGE"
echo " Sample:    $VARIANTS"
echo "========================================"
echo ""
echo "--- Summary ---"
echo "Defining mutations in barcode:  $TOTAL_DEFINING"
echo "Detected in sample:             $TOTAL_FOUND / $TOTAL_DEFINING"
echo "Not detected:                   $TOTAL_MISSING / $TOTAL_DEFINING"
echo ""

echo "--- Detected mutations (POS | REF | ALT | ALT_FREQ | ALT_DP | PASS) ---"
if [ -s $FOUND_FILE ]; then
    awk -F'\t' 'BEGIN {printf "%-10s %-6s %-6s %-10s %-10s %-6s\n", "POS", "REF", "ALT", "ALT_FREQ", "ALT_DP", "PASS"}
                {printf "%-10s %-6s %-6s %-10s %-10s %-6s\n", $2, $3, $4, $11, $8, $14}' $FOUND_FILE
else
    echo "None detected"
fi

echo ""
echo "--- Missing mutations ---"
if [ -s $MISSING_FILE ]; then
    cat $MISSING_FILE
else
    echo "All defining mutations detected"
fi

# --- cleanup ---
rm -f $MUT_FILE $FOUND_FILE $MISSING_FILE