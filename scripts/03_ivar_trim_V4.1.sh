#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=ivar_trim
#SBATCH --output=logs/ivar_trim_%A_%a.out
#SBATCH --error=logs/ivar_trim_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/ivar

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output filenames and paths
in_bam="aligned/${sample}.sorted.bam"
out_bam="aligned/${sample}.trimmed.bam"
bed_file="reference/V4.1/SARS-CoV-2.primer.bed"
temp_prefix="aligned/${sample}.ivar_trimmed"

# make output and log directories if they don't exist
mkdir -p aligned logs

echo "Running iVar trim on $sample..."

# post-alignment trim with iVar
# produces: aligned/${sample}.ivar_trimmed.bam - temporary file
  # -e = exclude untrimmed reads
ivar trim \
  -i "$in_bam" \
  -b "$bed_file" \
  -p "$temp_prefix" \
  -e

# sort and index with samtools
samtools sort -@ 4 -o "$out_bam" "${temp_prefix}.bam"
samtools index "$out_bam"

# remove temporary file
rm "${temp_prefix}.bam"

echo "Done with $sample"