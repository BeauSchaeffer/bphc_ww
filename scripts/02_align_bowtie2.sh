#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH -p hsph
#SBATCH --mem=8G
#SBATCH --job-name=bt2_align
#SBATCH --output=logs/bt2_align_%A_%a.out
#SBATCH --error=logs/bt2_align_%A_%a.err
#SBATCH --array=0-9   # <-- Update to match number of samples in samples.txt

# Load Conda environment with bowtie2 + samtools
source ~/.bashrc
conda activate ./env/alignment

# Get sample ID from samples.txt based on SLURM array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# Define paths and filenames
r1="trimmed/${sample}_R1.trimmed.fastq.gz"
r2="trimmed/${sample}_R2.trimmed.fastq.gz"
ref="reference/SARS-CoV-2"

# Output locations
bam_out="aligned/${sample}.sorted.bam"
log_out="logs/${sample}_bowtie2.log"

# Create output directories if they don't exist
mkdir -p aligned logs

echo "Aligning $sample..."

# Align with Bowtie2 and process to BAM
bowtie2 -x "$ref" -1 "$r1" -2 "$r2" --threads 4 2> "$log_out" \
  | samtools view -@ 4 -bS - \
  | samtools sort -@ 4 -o "$bam_out"

# Index the BAM
samtools index "$bam_out"

echo "Done with $sample"
