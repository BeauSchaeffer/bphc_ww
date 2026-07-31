#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH -p hsph
#SBATCH --mem=2G
#SBATCH --job-name=coverage_qc
#SBATCH --output=logs/coverage_qc_%A_%a.out
#SBATCH --error=logs/coverage_qc_%A_%a.err
#SBATCH --array=0-684   # <-- Update to match number of samples ###** UPDATE **###

# define samples file and generate sample ID from array index
SAMPLES="/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/bphc_ww/data/baseload_batch2_ids.txt"   ###** UPDATE **### batch1 uses samples.txt
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# define input and output paths
depth_file="freyja_variants/${sample}_depths.tsv"
outfile="coverage/${sample}.qc.tsv"

# spike gene coordinates (Wuhan-Hu-1 reference)
SPIKE_START=21563
SPIKE_END=25384

# per-site depth floors defining "covered". 10x is a consensus-calling
# convention and is too permissive for mixture deconvolution: at 10x the finest
# resolvable allele frequency is 10%. Higher floors match the frequency
# resolution Freyja demix actually needs.
FLOORS="10,50,100,200,500"

# make output and log directories if they don't exist
mkdir -p coverage logs

echo "Running coverage QC for $sample..."

# check that depth file exists
if [ ! -f "$depth_file" ]; then
    echo "Error: depth file not found: $depth_file"
    exit 1
fi

# median depth is more meaningful than the mean for tiled-amplicon libraries:
# amplicon dropout makes the depth distribution bimodal, so the mean is pulled
# up by the amplicons that worked while much of the genome sits near zero.
median_depth=$(awk '{print $4}' "$depth_file" | sort -n | \
    awk '{a[NR]=$1} END {if (NR==0) print "NA"; else if (NR%2) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')

spike_median_depth=$(awk -v s="$SPIKE_START" -v e="$SPIKE_END" \
    '$2 >= s && $2 <= e {print $4}' "$depth_file" | sort -n | \
    awk '{a[NR]=$1} END {if (NR==0) print "NA"; else if (NR%2) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2}')

# single pass for everything else
# depth file format: CHROM, POS, REF, DEPTH ($4)
awk -v s="$SPIKE_START" -v e="$SPIKE_END" -v floors="$FLOORS" -v sample="$sample" \
    -v med="$median_depth" -v spike_med="$spike_median_depth" '
BEGIN { nf = split(floors, F, ",") }
{
    d = $4
    g_sum += d; g_n++
    if (d < 10) g_lt10++
    in_spike = ($2 >= s && $2 <= e)
    if (in_spike) { sp_sum += d; sp_n++; if (d < 10) sp_lt10++ }
    for (i = 1; i <= nf; i++) {
        if (d >= F[i]) {
            gcov[i]++
            if (in_spike) scov[i]++
        }
    }
}
END {
    # naming: depth_* = read depth, breadth_dN = fraction of sites at >= N x
    h = "sample\ttotal_sites\tspike_total_sites\tdepth_mean\tdepth_median"
    h = h "\tspike_depth_mean\tspike_depth_median\tsites_lt_d10\tspike_sites_lt_d10"
    for (i = 1; i <= nf; i++) h = h "\tbreadth_d" F[i]
    for (i = 1; i <= nf; i++) h = h "\tspike_breadth_d" F[i]
    print h

    r = sample "\t" g_n+0 "\t" sp_n+0
    r = r "\t" (g_n  > 0 ? sprintf("%.2f", g_sum/g_n)   : "NA") "\t" med
    r = r "\t" (sp_n > 0 ? sprintf("%.2f", sp_sum/sp_n) : "NA") "\t" spike_med
    r = r "\t" g_lt10+0 "\t" sp_lt10+0
    for (i = 1; i <= nf; i++)
        r = r "\t" (g_n  > 0 ? sprintf("%.4f", gcov[i]/g_n)  : "NA")
    for (i = 1; i <= nf; i++)
        r = r "\t" (sp_n > 0 ? sprintf("%.4f", scov[i]/sp_n) : "NA")
    print r
}' "$depth_file" > "$outfile"

echo "Done with $sample: $outfile"
