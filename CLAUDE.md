# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A SARS-CoV-2 wastewater genomic surveillance project for Boston Public Health Commission (BPHC). It has two
distinct, sequentially-dependent pipelines, not a software package: there is no build system, test suite, or
package manifest (no `renv.lock`/`DESCRIPTION`). The repo is split top-level by where each half runs:
**`informatics/` runs only on the Harvard FASRC cluster; `analysis/` runs only locally.** Nothing in
`analysis/` is meant to execute on the cluster (or vice versa) — see Cluster ↔ local sync below for how output
crosses that boundary.

1. **Bioinformatics pipeline** (`informatics/scripts/`, bash + SLURM, cluster-only) — raw paired-end FASTQ →
   trimmed → aligned → primer-trimmed → variant-called → lineage-deconvolved → aggregated relative abundances
   (via [Freyja](https://github.com/andersen-lab/Freyja)).
2. **Statistical modeling pipeline** (`analysis/bphc_ww_R/`, R, local-only) — Freyja output + BPHC metadata →
   cleaned/joined time series per neighborhood → a hierarchical Bayesian "gravity model" (Stan, via `cmdstanr`)
   estimating cross-neighborhood lineage-transmission coupling.

## Sequencing batches (drives which scripts/params to use)

Samples were sequenced in two batches with **different ARTIC primer schemes** — always check which batch a
sample belongs to before touching alignment/trimming:

- **Batch 1**: samples through 2023-05-04, ARTIC v4.1 primers → `informatics/reference/V4.1/`,
  `informatics/scripts/03_ivar_trim_V4.1.sh`
- **Batch 2**: samples after 2023-05-05, ARTIC v5.3.2 primers → `informatics/reference/V5.3.2/`,
  `informatics/scripts/03_ivar_trim_V5.3.2.sh`

`informatics/notebook.md` (gitignored, local-only) is the running lab notebook of what was actually run per
batch, including one-off sample exclusions, failure counts, and troubleshooting commands — check it before
assuming a script ran cleanly for a given batch. `informatics/testing.md` documents expected output counts for a
50-sample pipeline test run.

## Bioinformatics pipeline (`informatics/scripts/`)

Run as SLURM array jobs from the FASRC cluster working directory (not this local clone — see Cluster ↔ local
sync below). Each script:
- Resolves its sample ID for the job from `SLURM_ARRAY_TASK_ID` by indexing a line in a sample-ID list file
  (`baseload_batch2_ids.txt`, or `samples.txt` for batch 1).
- Has a `#SBATCH --array=0-N` line that **must be updated to match the current sample count** (marked
  `###** UPDATE **###` in each script) — a stale array range silently skips or errors on samples.
- Activates a dedicated conda env from `informatics/envs/*.yml` (`fastp`, `alignment`, `ivar`, `freyja`).

Pipeline order and outputs:

| Script | Output dir | Purpose |
|---|---|---|
| `01_fastp_trim.sh` | `trimmed/` | adapter trim (fastp) |
| `02_align_bowtie2.sh` | `aligned/*.sorted.bam` | align to SARS-CoV-2 reference (bowtie2 → samtools sort/index) |
| `03_ivar_trim_V4.1.sh` / `V5.3.2.sh` | `aligned/*.trimmed.bam` | primer trim (ivar, batch-specific BED file) |
| `04_freyja_variants.sh` | `freyja_variants/*_variants.tsv`, `*_depths.tsv` | variant calling (freyja variants) |
| `04.1_coverage_qc.sh` | `coverage/*.qc.tsv` | genome + spike-gene depth/coverage QC from the depths file |
| `05_freyja_demix.sh` | `freyja_demix/*_demix.tsv` | lineage relative-abundance deconvolution |
| `06_freyja_aggregate.sh` | `freyja_aggregate/aggregated.tsv` | combines all samples' demix output into one table |

`informatics/scripts/check_lineage_mutations.sh <lineage> <variants.tsv> <depths.tsv>` is a standalone diagnostic (not part
of the numbered pipeline): cross-references a lineage's Freyja barcode-defining mutations against a sample's
called variants and reports detected vs. missing (with depth-at-site for missing calls) — useful for manually
sanity-checking a demix call.

## Cluster ↔ local sync

Processing happens on FASRC cluster scratch/storage, not locally. `informatics/scripts/sync_*.sh` are the only
scripts that cross the informatics/analysis boundary — they land cluster output directly into `analysis/`,
which is otherwise a purely local, gitignored data directory (not a "clone" of anything cluster-side anymore):
- `informatics/scripts/sync_storage.sh` — pulls aggregated Freyja outputs and `data/` down to `analysis/data/`,
  `analysis/baseload_batch1_outputs/`, `analysis/baseload_batch2_outputs/`.
- `informatics/scripts/sync_draft_figs.sh` — pulls rendered figures down to `analysis/draft_figures/`.
- Metadata/analysis data are pushed back up via `rsync` (commands logged in `informatics/notebook.md`, not
  scripted).

## Statistical modeling pipeline (`analysis/bphc_ww_R/`)

Two layers, both R, run locally only, against `analysis/` as the data source:

**Data prep / exploration** — `data_clean.R`, `pipeline_prep.R`, `meta_exp.R`, `RAs.R`, `flow_model_data.R`.
These are not standalone/idempotent in the same way as the numbered scripts below: they build up intermediate
RDS objects interactively (metadata cleaning, QC joins, lineage-abundance parsing/collapsing, plotting) and are
the source of most files documented in `analysis/data/README.md`. They read/write against `analysis/data/`
directly (a single local path — there's no cluster variant to toggle, since these never run on the cluster).

**Gravity model pipeline** — `01_abs_abundance.R` through `06_plot_results.R`, plus `gravity_model.stan`. Unlike
the layer above, each numbered script is **standalone and idempotent**: it re-reads its inputs fresh from
`analysis/model_data/` and rewrites its outputs there, and is meant to be rerun top-to-bottom rather than have
cells run out of order. Flow:

1. `01_abs_abundance.R` — subsets `ww_collapsed_complete.rds` to the modeled lineage set, builds an
   absolute-abundance proxy (`relative abundance × mean effective copies/L`) → `abund.rds`
2. `02_build_transitions.R` — for each lineage, builds valid week-over-week log-growth transitions (both
   endpoints must be sampled and detected) plus the neighbor-abundance matrix used in the model's gravity sum →
   `transitions.rds`, `neighbor_abund_mat.rds`, `loc_index.rds`, `week_index.rds`, `lineage_index.rds`
3. `03_distance_matrix.R` — sewershed shapefile + population table → centroid-to-centroid distance matrix and
   population vector, aligned to `loc_index` order → `dist_mat.rds`, `pop_vec.rds`, `catchments.rds`
4. `04_stan_data.R` — assembles/centers covariates and rescales abundance into the final Stan data list →
   `stan_data.rds`, `centering.rds`
5. `05_fit.R` — compiles `gravity_model.stan` and fits via `cmdstanr` (4 chains, explicit tight inits — required
   to avoid divergences) → `gravity_fit.rds`, `gravity_summary.rds`
6. `06_plot_results.R` — reads the saved fit and produces 6 diagnostic figures to `draft_figures/` (connectivity
   map, per-lineage growth/coupling, shared-exponent posterior-vs-prior, distance contribution, informative-null
   summary)

`gravity_model.stan`: hierarchical multi-lineage model of log-growth in each neighborhood as a baseline
(lineage-specific `r`) plus a gravity-coupling term from neighboring catchments' donor abundance. `theta1`
(mass exponent) and `theta3` (distance-decay exponent) are shared/fully pooled across lineages; `r` and `theta0`
(coupling magnitude) are lineage-specific. Student-t likelihood for robustness to near-floor transitions.

Indexing convention used throughout the pipeline: `L` = lineage index, `j`/`k` = neighborhood/location index,
`t` = week index — all built once in `02_build_transitions.R` and reused unchanged through fitting and plotting.

## Gotchas

- `analysis/data/clin_collapsed.R` and `ww_collapsed.R` are **not R source** — they are
  `write_rds()` output that got the `.R` extension by mistake (see `RAs.R` around the `write_rds(ww_collapsed,
  "...ww_collapsed.R")` calls). Treat them as binary RDS data (`readRDS()`), not code to read/edit/source.
- The pipeline output directories (`trimmed/`, `aligned/`, `freyja_variants/`, `freyja_demix/`,
  `freyja_aggregate/`, `coverage/`, `logs/`) only ever exist on the cluster now — informatics is cluster-only,
  so these should not reappear at the repo root locally; if they do, they're stale leftovers, not the source of
  truth. `analysis/` (populated via `informatics/scripts/sync_*.sh`) is.
- `analysis/bphc_ww_R/gravity_model` (no extension) is a compiled `cmdstan` binary produced by `05_fit.R`'s
  `cmdstan_model()` call — it's a build artifact, not something to edit or commit.
- File dictionaries for what each RDS/CSV in the data directories actually contains live in
  `informatics/metadata/README.md` and `analysis/data/README.md` — check these before re-deriving a field's
  meaning from scratch.
