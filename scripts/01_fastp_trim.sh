#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=fastp_trim
#SBATCH --output=logs/fastp_trim_%A_%a.out
#SBATCH --error=logs/fastp_trim_%A_%a.err
#SBATCH --array=0-9   # <-- Update to match number of lines in samples.txt

# Load Conda and activate local environment
source ~/.bashrc
conda activate ./env/fastp

# Get sample ID from samples.txt
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# Define input/output paths
r1="data/${sample}_R1.fastq.gz"
r2="data/${sample}_R2.fastq.gz"
out_r1="trimmed/${sample}_R1.trimmed.fastq.gz"
out_r2="trimmed/${sample}_R2.trimmed.fastq.gz"
html="trimmed/${sample}_fastp.html"
json="trimmed/${sample}_fastp.json"

# Make output and log directories if they don't exist
mkdir -p trimmed logs

echo "Trimming $sample..."

fastp \
  -i "$r1" \
  -I "$r2" \
  -o "$out_r1" \
  -O "$out_r2" \
  -w 4 \
  --detect_adapter_for_pe \
  --html "$html" \
  --json "$json"

echo "Done with $sample"

