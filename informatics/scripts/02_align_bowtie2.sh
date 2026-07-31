#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH -p hsph
#SBATCH --mem=8G
#SBATCH --job-name=bt2_align
#SBATCH --output=logs/bt2_align_%A_%a.out
#SBATCH --error=logs/bt2_align_%A_%a.err
#SBATCH --array=0-0   # fallback only; pass --array at submit ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/alignment

# ids.txt lives in the batch working directory, so one script serves both
# batches without editing -- see informatics/REPROCESSING.md
SAMPLES="ids.txt"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# guard: array range longer than the sample list, or --array not passed
if [ -z "$sample" ]; then
    echo "Error: no sample at index ${SLURM_ARRAY_TASK_ID} in $SAMPLES"
    exit 1
fi

# define input and output filenames and paths
r1="trimmed/${sample}_R1.trimmed.fastq.gz"
r2="trimmed/${sample}_R2.trimmed.fastq.gz"
ref="reference/V5.3.2/SARS-CoV-2"
bam_out="aligned/${sample}.sorted.bam"
log_out="logs/${sample}_bowtie2.log"

# make output and log directories if they don't exist
mkdir -p aligned logs

echo "Aligning $sample..."

# clean up any leftover temp files and partial outputs from previous runs
rm -f aligned/${sample}.sorted.bam.tmp.*.bam
rm -f aligned/${sample}.sorted.bam
rm -f aligned/${sample}.sorted.bam.bai

# align with Bowtie2 and process to BAM
  # -bS = input SAM, output BAM
  # -@ 4 = 4 threads to match above
bowtie2 -x "$ref" -1 "$r1" -2 "$r2" --threads 4 2> "$log_out" \
  | samtools view -@ 4 -bS - \
  | samtools sort -@ 4 -o "$bam_out"

# Index the BAM with samtools
samtools index "$bam_out"

echo "Done with $sample"