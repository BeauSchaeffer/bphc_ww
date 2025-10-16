#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=ivar_trim
#SBATCH --output=logs/ivar_trim_%A_%a.out
#SBATCH --error=logs/ivar_trim_%A_%a.err
#SBATCH --array=0-9   # update this to your sample count

source ~/.bashrc
conda activate ./env/ivar

sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

in_bam="aligned/${sample}.sorted.bam"
out_bam="aligned/${sample}.trimmed.bam"
bed_file="reference/V4.1/SARS-CoV-2.primer.bed"
temp_prefix="aligned/${sample}.ivar_trimmed"

# Make sure output folders exist
mkdir -p aligned logs

echo "Running iVar trim on $sample..."

ivar trim \
  -i "$in_bam" \
  -b "$bed_file" \
  -p "$temp_prefix" \
  -e

# This produces: aligned/${sample}.ivar_trimmed.bam
# Now sort and index
samtools sort -@ 4 -o "$out_bam" "${temp_prefix}.bam"
samtools index "$out_bam"
rm "${temp_prefix}.bam"

echo "Done with $sample"
