#!/bin/bash
# check_lineage_mutations.sh
# Usage: bash check_lineage_mutations.sh <lineage> <variants_file>

BARCODES="/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja/lib/python3.11/site-packages/freyja/data/usher_barcodes.csv"

if [ $# -ne 2 ]; then
    echo "Usage: bash check_lineage_mutations.sh <lineage> <variants_file>"
    exit 1
fi

LINEAGE=$1
VARIANTS=$2

if [ ! -f "$VARIANTS" ]; then
    echo "Error: variants file not found: $VARIANTS"
    exit 1
fi

if ! grep -q "^${LINEAGE}," $BARCODES; then
    echo "Error: lineage $LINEAGE not found in barcode file"
    exit 1
fi

# extract defining mutations
MUT_FILE=$(mktemp /tmp/lineage_muts.XXXX)
grep "^${LINEAGE}," $BARCODES | tr ',' '\n' | \
    paste - <(head -1 $BARCODES | tr ',' '\n') | \
    awk -F'\t' '$1 == "1.0" {print $2}' > $MUT_FILE

TOTAL_DEFINING=$(wc -l < $MUT_FILE)

# parse each defining mutation and do exact position + allele matching
FOUND_FILE=$(mktemp /tmp/found_muts.XXXX)
MISSING_FILE=$(mktemp /tmp/missing_muts.XXXX)

while read mut; do
    # extract position (numeric part) and alt allele (trailing letters/symbols)
    pos=$(echo "$mut" | grep -oP '\d+')
    alt=$(echo "$mut" | grep -oP '[A-Z-]+$')

    # exact match on position (col 2) and alt allele (col 4)
    match=$(awk -F'\t' -v p="$pos" -v a="$alt" 'NR>1 && $2==p && $4==a {print}' $VARIANTS)

    if [ -n "$match" ]; then
        echo "$mut"$'\t'"$match" >> $FOUND_FILE
    else
        echo "$mut" >> $MISSING_FILE
    fi
done < $MUT_FILE

TOTAL_FOUND=$(wc -l < $FOUND_FILE 2>/dev/null || echo 0)
TOTAL_MISSING=$(wc -l < $MISSING_FILE 2>/dev/null || echo 0)

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

echo "--- Detected mutations (MUT | POS | REF | ALT | ALT_FREQ | ALT_DP | PASS) ---"
if [ -s "$FOUND_FILE" ]; then
    awk -F'\t' 'BEGIN {printf "%-12s %-10s %-6s %-6s %-10s %-10s %-6s\n", "MUT", "POS", "REF", "ALT", "ALT_FREQ", "ALT_DP", "PASS"}
                {printf "%-12s %-10s %-6s %-6s %-10s %-10s %-6s\n", $1, $3, $4, $5, $12, $9, $15}' $FOUND_FILE
else
    echo "None detected"
fi

echo ""
echo "--- Missing mutations ---"
if [ -s "$MISSING_FILE" ]; then
    cat $MISSING_FILE
else
    echo "All defining mutations detected"
fi

rm -f $MUT_FILE $FOUND_FILE $MISSING_FILE