#!/bin/bash
# check_lineage_mutations.sh
# Usage: bash scripts/check_lineage_mutations.sh <lineage> <variants_file> <depths_file>
# Example: bash scripts/check_lineage_mutations.sh BA.5 freyja_variants/63135871_variants.tsv freyja_variants/63135871_depths.tsv

BARCODES="/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja/lib/python3.11/site-packages/freyja/data/usher_barcodes.csv"

if [ $# -ne 3 ]; then
    echo "Usage: bash check_lineage_mutations.sh <lineage> <variants_file> <depths_file>"
    echo "Example: bash check_lineage_mutations.sh BA.5 freyja_variants/63135871_variants.tsv freyja_variants/63135871_depths.tsv"
    exit 1
fi

LINEAGE=$1
VARIANTS=$2
DEPTHS=$3

if [ ! -f "$VARIANTS" ]; then
    echo "Error: variants file not found: $VARIANTS"
    exit 1
fi

if [ ! -f "$DEPTHS" ]; then
    echo "Error: depths file not found: $DEPTHS"
    exit 1
fi

if ! grep -q "^${LINEAGE}," "$BARCODES"; then
    echo "Error: lineage $LINEAGE not found in barcode file"
    exit 1
fi

# extract defining mutations using column-index matching
# header row starts with a bare comma (no name for lineage column), mutations start at column 2
# values of 1.0 indicate a defining mutation for the lineage
MUT_FILE=$(mktemp /tmp/lineage_muts.XXXX)
awk -F',' -v lineage="$LINEAGE" '
    NR==1 {
        for (i=2; i<=NF; i++) colname[i] = $i
    }
    NR>1 && $1==lineage {
        for (i=2; i<=NF; i++)
            if ($i == "1.0") print colname[i]
    }
' "$BARCODES" > "$MUT_FILE"

TOTAL_DEFINING=$(wc -l < "$MUT_FILE")

if [ "$TOTAL_DEFINING" -eq 0 ]; then
    echo "Error: no defining mutations found for $LINEAGE in barcode file"
    rm -f "$MUT_FILE"
    exit 1
fi

# parse each defining mutation and do exact position + allele matching
FOUND_FILE=$(mktemp /tmp/found_muts.XXXX)
MISSING_FILE=$(mktemp /tmp/missing_muts.XXXX)

while read mut; do
    pos=$(echo "$mut" | grep -oP '\d+')
    alt=$(echo "$mut" | grep -oP '[A-Z-]+$')

    match=$(awk -F'\t' -v p="$pos" -v a="$alt" 'NR>1 && $2==p && $4==a {print}' "$VARIANTS")

    if [ -n "$match" ]; then
        echo "$mut"$'\t'"$match" >> "$FOUND_FILE"
    else
        # look up depth at missing position
        # depths file format: CHROM, POS, REF, DEPTH ($4)
        depth=$(awk -F'\t' -v p="$pos" '$2==p {print $4}' "$DEPTHS")
        depth=${depth:-0}
        echo "$mut"$'\t'"$depth" >> "$MISSING_FILE"
    fi
done < "$MUT_FILE"

TOTAL_FOUND=$(wc -l < "$FOUND_FILE" 2>/dev/null || echo 0)
TOTAL_MISSING=$(wc -l < "$MISSING_FILE" 2>/dev/null || echo 0)

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
                {printf "%-12s %-10s %-6s %-6s %-10s %-10s %-6s\n", $1, $3, $4, $5, $12, $9, $15}' "$FOUND_FILE"
else
    echo "None detected"
fi

echo ""
echo "--- Missing mutations (MUT | DEPTH_AT_SITE) ---"
if [ -s "$MISSING_FILE" ]; then
    awk -F'\t' 'BEGIN {printf "%-12s %-10s\n", "MUT", "DEPTH"}
                {printf "%-12s %-10s\n", $1, $2}' "$MISSING_FILE"
else
    echo "All defining mutations detected"
fi

rm -f "$MUT_FILE" "$FOUND_FILE" "$MISSING_FILE"