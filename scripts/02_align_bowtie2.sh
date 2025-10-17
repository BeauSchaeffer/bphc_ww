#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH -p hsph
#SBATCH --mem=8G
#SBATCH --job-name=bt2_align
#SBATCH --output=logs/bt2_align_%A_%a.out
#SBATCH --error=logs/bt2_align_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate ./env/alignment

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output filenames and paths
r1="trimmed/${sample}_R1.trimmed.fastq.gz"
r2="trimmed/${sample}_R2.trimmed.fastq.gz"
ref="reference/SARS-CoV-2"
bam_out="aligned/${sample}.sorted.bam"
log_out="logs/${sample}_bowtie2.log"

# make output and log directories if they don't exist
mkdir -p aligned logs

echo "Aligning $sample..."

# align with Bowtie2 and process to BAM
  # -bS = input SAM, output BAM
  # -@ 4 = 4 threads to match above
bowtie2 -x "$ref" -1 "$r1" -2 "$r2" --threads 4 2> "$log_out" \
  | samtools view -@ 4 -bS - \
  | samtools sort -@ 4 -o "$bam_out"

# Index the BAM with samtools
samtools index "$bam_out"

echo "Done with $sample"