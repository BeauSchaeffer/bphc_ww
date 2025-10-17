#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=fastp_trim
#SBATCH --output=logs/fastp_trim_%A_%a.out
#SBATCH --error=logs/fastp_trim_%A_%a.err
#SBATCH --array=0-49   # <-- Update to match number of samples in samples.txt ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate ./env/fastp

# generate sample ID from sample.txt and array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# define input and output filenames and paths
r1="data/${sample}_R1.fastq.gz"
r2="data/${sample}_R2.fastq.gz"
out_r1="trimmed/${sample}_R1.trimmed.fastq.gz"
out_r2="trimmed/${sample}_R2.trimmed.fastq.gz"
html="trimmed/${sample}_fastp.html"
json="trimmed/${sample}_fastp.json"

# make output and log directories if they don't exist
mkdir -p trimmed logs

echo "Trimming $sample..."

# pre-alignment trim with fastp
  # --detect_adapter_for_pe = enable adapter detection for PE data to get ultra-clean data.

fastp \
  -i "$r1" \
  -I "$r2" \
  -o "$out_r1" \
  -O "$out_r2" \
  --detect_adapter_for_pe \
  --html "$html" \
  --json "$json"

echo "Done with $sample"