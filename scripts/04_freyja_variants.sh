#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=6G
#SBATCH --job-name=freyja_variants
#SBATCH --output=logs/freyja_vars_%A_%a.out
#SBATCH --error=logs/freyja_vars_%A_%a.err
#SBATCH --array=0-9   # update to match number of samples

# Activate Freyja conda environment
source ~/.bashrc
conda activate ./env/freyja_new

# Get sample ID from samples.txt
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# File paths
bam="aligned/${sample}.trimmed.bam"
ref="reference/V4.1/SARS-CoV-2.reference.fasta"

# Output directory and filenames
outdir="freyja_variants"
mkdir -p "$outdir" logs

variants_tsv="${outdir}/${sample}_variants"
depths_tsv="${outdir}/${sample}_depths.tsv"

echo "Running freyja variants for $sample..."

freyja variants \
  "$bam" \
  --ref "$ref" \
  --variants "$variants_tsv" \
  --depths "$depths_tsv"

echo "Done with $sample"
