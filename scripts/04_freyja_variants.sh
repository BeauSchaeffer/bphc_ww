#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=6G
#SBATCH --job-name=freyja_variants
#SBATCH --output=logs/freyja_vars_%A_%a.out
#SBATCH --error=logs/freyja_vars_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate ./env/freyja_new

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output filenames and paths
ref="reference/V4.1/SARS-CoV-2.reference.fasta"
bam="aligned/${sample}.trimmed.bam"
variants_tsv="freyja_variants/${sample}_variants"
depths_tsv="freyja_variants/${sample}_depths.tsv"

# make output and log directories if they don't exist
mkdir -p freyja_variants logs

echo "Running freyja variants for $sample..."

# call variants with Freyja
freyja variants \
  "$bam" \
  --ref "$ref" \
  --variants "$variants_tsv" \
  --depths "$depths_tsv"

echo "Done with $sample"