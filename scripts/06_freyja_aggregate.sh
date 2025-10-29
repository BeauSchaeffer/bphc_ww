#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH -p hsph
#SBATCH --mem=4G
#SBATCH --job-name=freyja_aggregate
#SBATCH --output=logs/freyja_aggregate_%A.out
#SBATCH --error=logs/freyja_aggregate_%A.err

# load Conda and activate local environment
source ~/.bashrc
conda activate /n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/envs/freyja

# define input and output paths
indir="freyja_demix/"
outfile="freyja_aggregate/aggregated.tsv"

# make output and log directories if they don't exist
mkdir -p freyja_aggregate logs

echo "Running freyja aggregate on all demix outputs..."

# aggregate Freyja demix results into one combined table
freyja aggregate "$indir" --output "$outfile"

echo "Aggregation complete: $outfile"