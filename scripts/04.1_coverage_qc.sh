#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=coverage_qc
#SBATCH --output=logs/coverage_qc_%A_%a.out
#SBATCH --error=logs/coverage_qc_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/ivar

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output paths
depth_file="freyja_variants/${sample}_depths.tsv"
outfile="coverage/${sample}.qc.tsv"

# spike gene coordinates (Wuhan-Hu-1 reference)
SPIKE_START=21563
SPIKE_END=25384

# make output and log directories if they don't exist
mkdir -p coverage logs

echo "Running coverage QC for $sample..."

# check that depth file exists
if [ ! -f "$depth_file" ]; then
    echo "Error: depth file not found: $depth_file"
    exit 1
fi

# whole genome metrics
# depth file format: CHROM, POS, REF, DEPTH ($4)
avg_depth=$(awk '{sum += $4; n++} END {if (n > 0) printf "%.2f", sum/n; else print "NA"}' "$depth_file")

sites_low_depth=$(awk '$4 <= 10 {n++} END {print n+0}' "$depth_file")

total_sites=$(awk 'END {print NR}' "$depth_file")

frac_covered=$(awk -v total="$total_sites" '$4 >= 10 {n++} END {if (total > 0) printf "%.4f", n/total; else print "NA"}' "$depth_file")

# spike gene metrics
spike_avg_depth=$(awk -v s="$SPIKE_START" -v e="$SPIKE_END" \
    '$2 >= s && $2 <= e {sum += $4; n++} END {if (n > 0) printf "%.2f", sum/n; else print "NA"}' "$depth_file")

spike_sites_low_depth=$(awk -v s="$SPIKE_START" -v e="$SPIKE_END" \
    '$2 >= s && $2 <= e && $4 <= 10 {n++} END {print n+0}' "$depth_file")

spike_frac_covered=$(awk -v s="$SPIKE_START" -v e="$SPIKE_END" \
    '$2 >= s && $2 <= e {total++; if ($4 >= 10) n++} END {if (total > 0) printf "%.4f", n/total; else print "NA"}' "$depth_file")

# write header + results to TSV
echo -e "sample\tavg_depth\tsites_depth_lte10\ttotal_sites\tfrac_genome_cov_10x\tspike_avg_depth\tspike_sites_depth_lte10\tspike_frac_cov_10x" > "$outfile"
echo -e "${sample}\t${avg_depth}\t${sites_low_depth}\t${total_sites}\t${frac_covered}\t${spike_avg_depth}\t${spike_sites_low_depth}\t${spike_frac_covered}" >> "$outfile"

echo "Done with $sample: $outfile"