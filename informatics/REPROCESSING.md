# Reprocessing samples on FASRC

How to (re)run the bioinformatics pipeline, from a single step up to a full
rebuild from raw reads. Everything here runs on the Harvard FASRC cluster;
nothing in `analysis/` runs there.

## Layout

Code and data are kept in separate trees, so a corrupted working directory can't
take the git repo with it.

```
/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/
├── envs/                                   # built conda envs (fastp, alignment, ivar, freyja)
├── bphc_ww_repo/                           # the git clone -- CODE ONLY, never write outputs here
│   └── informatics/{scripts,reference,envs,metadata}
└── bphc_ww/                                # WORK dir
    ├── data/                               # batch1_ids.txt, batch2_ids.txt, metadata
    ├── baseload_batch1_outputs/            # ARTIC v4.1  (sampled through 2023-05-04)
    │   ├── ids.txt                         # this batch's sample list
    │   ├── reference -> ../../bphc_ww_repo/informatics/reference
    │   ├── trimmed/ aligned/ freyja_variants/ freyja_demix/
    │   └── freyja_aggregate/ coverage/ logs/
    └── baseload_batch2_outputs/            # ARTIC v5.3.2 (sampled from 2023-05-05)
        └── (identical structure)
```

Raw reads live outside both trees and are shared by the two batches:
`/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/Project_raw_reads/bphc_baseload/`
(`<sample>_R1.fastq.gz` / `_R2.fastq.gz`).

**The one rule:** every script resolves inputs and outputs relative to the
current working directory -- `ids.txt`, `reference/`, `aligned/`, `logs/`, all of
it. The only absolute paths are the conda envs and `RAW_DIR` in
`01_fastp_trim.sh`. So `cd` into a batch directory and submit the script by its
full path in the repo. **The working directory selects the batch.**

## Current state of the batch directories

As of 2026-07-31 both batch dirs hold outputs from the original run of **990
samples (305 batch 1 / 685 batch 2)**. `ids.txt` now lists **1006 (316 / 690)** --
the 16 extra samples have raw reads and metadata but were never processed (see
`notebook.md` for why). Batch 1's `coverage/` did not survive the move off
scratch; `04.1` regenerates it.

This matters because most reruns consume existing outputs, and the existing
outputs cover fewer samples than `ids.txt`.

## Submitting

```bash
STORE=/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer
REPO=$STORE/bphc_ww_repo

cd $STORE/bphc_ww/baseload_batch2_outputs
sbatch --array=0-$(($(wc -l < ids.txt)-1)) $REPO/informatics/scripts/04_freyja_variants.sh
```

Always pass `--array`, computed from whichever list you're running. It overrides
the `#SBATCH --array=0-0` fallback, so the range can never go stale. If you
forget, you get a single task rather than a wrong-sized run. Each script also
aborts with a clear error if its index falls past the end of the list.

`06_freyja_aggregate.sh` is a single job -- submit with no `--array`.

## Pipeline order

| step | script | reads | writes |
|---|---|---|---|
| 1 | `01_fastp_trim.sh` | `$RAW_DIR/*.fastq.gz` | `trimmed/` |
| 2 | `02_align_bowtie2.sh` | `trimmed/` | `aligned/*.sorted.bam` |
| 3 | `03_ivar_trim_V4.1.sh` **or** `_V5.3.2.sh` | `aligned/*.sorted.bam` | `aligned/*.trimmed.bam` |
| 4 | `04_freyja_variants.sh` | `aligned/*.trimmed.bam` | `freyja_variants/*_{variants,depths}.tsv` |
| 4.1 | `04.1_coverage_qc.sh` | `freyja_variants/*_depths.tsv` | `coverage/*.qc.tsv` |
| 5 | `05_freyja_demix.sh` | `freyja_variants/` | `freyja_demix/*_demix.tsv` |
| 6 | `06_freyja_aggregate.sh` | **all of** `freyja_demix/` | `freyja_aggregate/aggregated.tsv` |

Steps 1-5 are per-sample and driven by `ids.txt`. **Step 6 is not** -- it
aggregates every file in `freyja_demix/`, regardless of what's in `ids.txt`.

### Step 3 is the only primer-scheme-dependent step

Batch 1 used ARTIC v4.1, batch 2 v5.3.2. **Only the primer BED files differ**
between the two reference directories -- the genome FASTA and the bowtie2 index
are byte-identical. Steps 2 and 4 are therefore scheme-independent even though
they name `reference/V5.3.2/` in their paths.

- `baseload_batch1_outputs/` -> `03_ivar_trim_V4.1.sh`
- `baseload_batch2_outputs/` -> `03_ivar_trim_V5.3.2.sh`

Nothing enforces this. Getting it wrong soft-clips the wrong coordinates, which
surfaces downstream as degraded coverage rather than an error.

---

# Scenario A: rerun one step over existing outputs

The common case -- e.g. rerunning `04.1_coverage_qc.sh` after changing `FLOORS`.
Nothing upstream is reprocessed.

Because existing outputs cover fewer samples than `ids.txt`, build a list scoped
to the inputs that actually exist:

```bash
cd $STORE/bphc_ww/baseload_batch2_outputs

# samples that have the input this step needs
ls freyja_variants/*_depths.tsv | xargs -n1 basename | sed 's/_depths\.tsv$//' | sort > /tmp/have_depths.txt
comm -12 <(sort ids.txt) /tmp/have_depths.txt > ids_step.txt
wc -l ids_step.txt        # expect 685 for batch 2, 305 for batch 1
```

Then temporarily point the step at it. The scripts read `ids.txt` by name, so
either swap the file or pass a scoped list:

```bash
cp ids.txt ids_full.txt && cp ids_step.txt ids.txt
sbatch --array=0-$(($(wc -l < ids.txt)-1)) $REPO/informatics/scripts/04.1_coverage_qc.sh
# when the job finishes:
mv ids_full.txt ids.txt
```

Running against the full `ids.txt` instead is *safe but noisy* -- the 16 samples
without inputs fail loudly with "depth file not found" and write nothing.

Do this for **both batches** before syncing.

# Scenario B: process samples that were never processed

For the 16 samples that have reads but no outputs. Run steps 1-6 for just those,
into the existing batch directory. Outputs are per-sample filenames, so they land
alongside the existing ones without disturbing them.

```bash
cd $STORE/bphc_ww/baseload_batch1_outputs

# samples in ids.txt with no demix output yet
ls freyja_demix/*_demix.tsv 2>/dev/null | xargs -n1 basename | sed 's/_demix\.tsv$//' | sort > /tmp/done.txt
comm -23 <(sort ids.txt) /tmp/done.txt > ids_new.txt
wc -l ids_new.txt

cp ids.txt ids_full.txt && cp ids_new.txt ids.txt
for s in 01_fastp_trim 02_align_bowtie2 03_ivar_trim_V4.1 04_freyja_variants 04.1_coverage_qc 05_freyja_demix; do
  echo "submit $s with --array=0-$(($(wc -l < ids.txt)-1))"
done
mv ids_full.txt ids.txt
```

Submit each step only after the previous one completes -- or chain them with
`sbatch --dependency=afterok:<jobid>`.

**Then rerun step 6.** It aggregates the whole `freyja_demix/` directory, so the
existing `aggregated.tsv` is stale the moment you add any demix file.

# Scenario C: full reprocess from raw reads

Rebuilds everything for every sample in `ids.txt`, starting from the FASTQs.

Archive or clear the previous outputs first -- otherwise you end up with a
directory holding a mix of vintages and no way to tell them apart:

```bash
cd $STORE/bphc_ww/baseload_batch2_outputs
mkdir -p ../archive_$(date +%F)/batch2
mv trimmed aligned freyja_variants freyja_demix freyja_aggregate coverage ../archive_$(date +%F)/batch2/
mkdir -p trimmed aligned freyja_variants freyja_demix freyja_aggregate coverage logs
```

Then run steps 1-6 in order against the full `ids.txt`, using the batch's
matching step-3 script.

Budget for it: `trimmed/` and `aligned/` dominate disk, and a full rerun means
holding both the archive and the new outputs. **Check quota before starting.**

---

## Rebuilding `ids.txt`

`ids.txt` = samples assigned to the batch by sampling date, intersected with
samples that have reads.

Batch lists are generated locally from
`analysis/data/bphc_biobot_sequence_metadata.csv` (split at 2023-05-04, keeping
only numeric `FASTQ_ID`s) and pushed to `bphc_ww/data/`. Then:

```bash
RAW=/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/Project_raw_reads/bphc_baseload
ls $RAW/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1\.fastq\.gz$//' | sort > /tmp/have_reads.txt

for b in 1 2; do
  d=$STORE/bphc_ww/baseload_batch${b}_outputs
  comm -12 <(sort $STORE/bphc_ww/data/batch${b}_ids.txt) /tmp/have_reads.txt > $d/ids.txt
  echo "batch${b}: $(wc -l < $d/ids.txt) samples"
done
```

Expected: **batch 1 = 316, batch 2 = 690** (1006 total, matching the raw read
count). 20 batch-1 metadata rows carry kit-style placeholders instead of a
numeric `FASTQ_ID` and are excluded; matching them needs `SHIPPING_KIT_ID`
resolved first.

## Setting up a batch directory from scratch

```bash
for b in 1 2; do
  d=$STORE/bphc_ww/baseload_batch${b}_outputs
  mkdir -p $d/{trimmed,aligned,freyja_variants,freyja_demix,freyja_aggregate,coverage,logs}
  ln -sfn $STORE/bphc_ww_repo/informatics/reference $d/reference
done
```

`mkdir -p` is additive and never touches existing contents. `logs/` must exist
**before** submitting: SLURM opens `#SBATCH --output=logs/...` before your script
runs, so the in-script `mkdir -p logs` is too late.

## Getting results back to the laptop

Run locally, not on the cluster:

- `informatics/scripts/sync_storage.sh` -- `aggregated.tsv` (both batches),
  `bphc_ww/data/`, `coverage/` (both batches)
- `informatics/scripts/sync_variants.sh` -- per-sample `freyja_variants/`
  (`*_variants.tsv` + `*_depths.tsv`) for local diversity work. Large (~1-2 GB),
  one-off.

Metadata and ID lists are pushed *up* by hand with `rsync`; log the command in
`notebook.md`.

After syncing coverage, rerun `analysis/bphc_ww_R/data_clean.R` to rebuild
`compiled_coverage.rds` and `meta_clean.rds`.

## Gotchas

- **Rerun `04.1` for both batches together.** Its output schema changes with the
  `FLOORS` variable. A batch left on the old schema gives all-NA breadth columns
  downstream, and a filter on those NAs would silently drop the whole batch.
  `data_clean.R` has a guard that errors on mismatched columns.
- **Step 6 ignores `ids.txt`.** Any partial reprocess makes `aggregated.tsv`
  stale; rerun step 6 afterwards.
- **`trimmed/` and `aligned/` dominate disk** and are regenerable.
  `freyja_variants/`, `freyja_demix/`, `freyja_aggregate/`, and `coverage/` are
  what feed downstream analysis.
- **`sync_draft_figs.sh` is obsolete.** It pulls from the dead netscratch path
  into `analysis/draft_figures/`, which is now written locally by
  `06_plot_results.R` and `qc_filter_viz.R`. Running it overwrites current
  figures with stale ones.
- **Never write outputs into `bphc_ww_repo/`.** Separate trees are what keep
  `git status` usable.
