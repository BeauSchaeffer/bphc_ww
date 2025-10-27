#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=freyja_demix
#SBATCH --output=logs/freyja_demix_%A_%a.out
#SBATCH --error=logs/freyja_demix_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output filenames and paths
variants="freyja_variants/${sample}_variants.tsv"
depths="freyja_variants/${sample}_depths.tsv"
outfile="freyja_demix/${sample}_demix.tsv"

# make output and log directories if they don't exist
mkdir -p freyja_demix logs

echo "Running freyja demix for $sample..."

# deconvolution with Freyja to recover lineage relative abundances
  # --eps 0.000001 = lower threshold for dropping low abundance lineages so that RAs sum to 1
freyja demix \
  "$variants" "$depths" \
  --output "$outfile" \
  --eps 0.000001

echo "Done with $sample"
