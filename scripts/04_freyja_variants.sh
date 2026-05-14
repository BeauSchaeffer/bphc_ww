#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=32G
#SBATCH --job-name=freyja_variants
#SBATCH --output=logs/freyja_vars_%A_%a.out
#SBATCH --error=logs/freyja_vars_%A_%a.err
#SBATCH --array=0-684   # <-- Update to match number of samples ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja

# define samples file and generate sample ID from array index
SAMPLES="/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/bphc_ww/data/baseload_batch2_ids.txt"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# define input and output filenames and paths
ref="reference/V5.3.2/SARS-CoV-2.reference.fasta"
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