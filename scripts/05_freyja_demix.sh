#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=freyja_demix
#SBATCH --output=logs/freyja_demix_%A_%a.out
#SBATCH --error=logs/freyja_demix_%A_%a.err
#SBATCH --array=0-9   # update to match number of samples

# Activate Freyja conda environment
source ~/.bashrc
conda activate ./env/freyja_new

# Get sample ID from samples.txt
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# File paths
indir="freyja_variants"
variants="${indir}/${sample}_variants.tsv"
depths="${indir}/${sample}_depths.tsv"
outfile="freyja_demix/${sample}_demix.tsv"

# Create output dir if needed
mkdir -p freyja_demix logs

echo "Running freyja demix for $sample..."

freyja demix \
  "$variants" "$depths" \
  --output "$outfile" \
  --eps 0.000001

echo "Done with $sample"
