#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=freyja_demix
#SBATCH --output=logs/freyja_demix_%A_%a.out
#SBATCH --error=logs/freyja_demix_%A_%a.err
#SBATCH --array=0-0   # fallback only; pass --array at submit ###** UPDATE **###

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja

# ids.txt lives in the batch working directory, so one script serves both
# batches without editing -- see informatics/REPROCESSING.md
SAMPLES="${SAMPLES:-ids.txt}"
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

# guard: array range longer than the sample list, or --array not passed
if [ -z "$sample" ]; then
    echo "Error: no sample at index ${SLURM_ARRAY_TASK_ID} in $SAMPLES"
    exit 1
fi

# define input and output filenames and paths
variants="freyja_variants/${sample}_variants.tsv"
depths="freyja_variants/${sample}_depths.tsv"
outfile="freyja_demix/${sample}_demix.tsv"

# make output and log directories if they don't exist
mkdir -p freyja_demix logs

echo "Running freyja demix for $sample..."

# deconvolution with Freyja to recover lineage relative abundances
  # --eps 0.000001 = lower threshold for dropping low abundance lineages so that RAs sum to 1
  # --depthcutoff 10 = default 0, was using this to troubleshoot but have resolved with new version
  # --autoadapt = can help a bit to deal with the high error rates with low quality samples
freyja demix \
  "$variants" "$depths" \
  --output "$outfile" \
  --eps 0.000001 \
  --depthcutoff 10 \
  --autoadapt

echo "Done with $sample"
