#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=ivar_trim
#SBATCH --output=logs/ivar_trim_%A_%a.out
#SBATCH --error=logs/ivar_trim_%A_%a.err
#SBATCH --array=0-0   # fallback only; pass --array at submit ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/ivar

# ids.txt lives in the batch working directory, so one script serves both
# batches without editing -- see informatics/REPROCESSING.md.
# NOTE: this is the ARTIC v4.1 script -- run it only from a batch 1 directory.
SAMPLES="ids.txt"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# guard: array range longer than the sample list, or --array not passed
if [ -z "$sample" ]; then
    echo "Error: no sample at index ${SLURM_ARRAY_TASK_ID} in $SAMPLES"
    exit 1
fi

# define input and output filenames and paths
in_bam="aligned/${sample}.sorted.bam"
out_bam="aligned/${sample}.trimmed.bam"
bed_file="reference/V4.1/SARS-CoV-2.primer.bed"
temp_prefix="aligned/${sample}.ivar_trimmed"
# remove temp files if job fails mid-run or at end of script
trap "rm -f ${temp_prefix}.bam" EXIT

# make output and log directories if they don't exist
mkdir -p aligned logs coverage

echo "Running iVar trim on $sample..."

# post-alignment trim with iVar
# produces: aligned/${sample}.ivar_trimmed.bam - temporary file
# -e flag: include reads with no primers

ivar trim \
  -i "$in_bam" \
  -b "$bed_file" \
  -p "$temp_prefix" \
  -q 20 \
  -e

# sort and index with samtools
samtools sort -@ 4 -o "$out_bam" "${temp_prefix}.bam"
samtools index "$out_bam"

echo "Done with $sample"