# jaccard_analysis.R

# ---- PATHS (edit these) ------------------------------------------------------
storage_dir <- "../data/"

# concord_v1 = no coverage QC, every sequenced sample contributes
# concord_v2 = coverage QC applied, matching RAs.R and the gravity model
concord_ver <- "concord_v2"

out_dir <- paste0(storage_dir, concord_ver, "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(lubridate)

source("lineage_helpers.R")

# ---- Params -------------------------------------------------------------------

# minimum per-sample Freyja abundance for a call to count at all -- excludes
# near-zero barcode calls that are demix noise rather than real signal. This is
# WW's instrument noise floor, not the presence rule below; clinical has no
# analogue because a pangolin call on a consensus genome is an observation, not
# a mixture estimate.
ww_detect_thresh <- 0.01

# coverage QC, matching RAs.R. Unlike the modeling pipeline -- which nulls
# meaneffCopiesL because its completed grid would otherwise read a failure as a
# confident zero -- the monthly compositions here are built only from observed
# samples, so a failing sample is simply dropped and contributes nothing.
qc_breadth_col <- "spike_breadth_d200"
qc_breadth_min <- 0.4
apply_qc <- identical(concord_ver, "concord_v2")

# per-month prevalence threshold for the min05 panel: a lineage counts as
# present in a month if it reaches this share of that month's citywide
# composition. Applied identically to WW and clinical, and evaluated on each
# month alone -- deliberately unlike RAs.R's retrospective "ever clears 0.2"
# rule, which exists so the bar plots show a lineage's early growth.
monthly_min_prop <- 0.05

# ---- Helpers --------------------------------------------------------------

# monthly composition -> one row per month with a list-col of lineages present
sets_from_comp <- function(comp) {
  comp |>
    filter(proportion > 0) |>
    group_by(month) |>
    summarise(lineages = list(unique(lineage)), n_lineages = n(), .groups = "drop")
}

# monthly Jaccard between two set tables, joined against monthly sample-count
# tables purely for plotting context
jaccard_over_time <- function(ww_sets, clin_sets, ww_n_samples, clin_n_samples) {
  months <- sort(unique(c(ww_sets$month, clin_sets$month)))

  tibble(month = months) |>
    left_join(ww_sets   |> select(month, ww_lineages = lineages),   by = "month") |>
    left_join(clin_sets |> select(month, clin_lineages = lineages), by = "month") |>
    rowwise() |>
    mutate(
      ww_lineages   = list(ww_lineages   %||% character(0)),
      clin_lineages = list(clin_lineages %||% character(0)),
      n_ww     = length(ww_lineages),
      n_clin   = length(clin_lineages),
      n_shared = length(intersect(ww_lineages, clin_lineages)),
      n_union  = length(union(ww_lineages, clin_lineages)),
      # undefined, not zero, when either source is absent that month -- a 0
      # would plot as real disagreement
      jaccard  = if_else(n_ww > 0 & n_clin > 0, n_shared / n_union, NA_real_)
    ) |>
    ungroup() |>
    select(-ww_lineages, -clin_lineages) |>
    left_join(ww_n_samples, by = "month") |>
    left_join(clin_n_samples, by = "month") |>
    mutate(
      n_ww_samples   = replace_na(n_ww_samples, 0),
      n_clin_samples = replace_na(n_clin_samples, 0)
    )
}

# ---- 1. Load + parse WW --------------------------------------------------------

ww_metadata <- read_rds(paste0(storage_dir, "meta_clean.rds"))

sublineages.final <- bind_rows(
  parse_freyja_aggregate("../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv"),
  parse_freyja_aggregate("../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv")
)

sublin_meta <- left_join(sublineages.final, ww_metadata, by = c("Sample" = "FASTQ_ID")) |>
  mutate(abundance = as.numeric(abundance),
         month = floor_date(SAMPLING_DATE, "month"))

if (apply_qc) {
  n_before <- n_distinct(sublin_meta$Sample)
  sublin_meta <- sublin_meta |> filter(.data[[qc_breadth_col]] >= qc_breadth_min)
  cat(sprintf("QC: %s >= %s | %d of %d samples retained\n",
              qc_breadth_col, qc_breadth_min, n_distinct(sublin_meta$Sample), n_before))
}

# renormalize after the QC drop, so each retained sample still sums to 1
sublin_meta <- sublin_meta |>
  apply_ww_detection(ww_detect_thresh) |>
  mutate(group = parent_group(sublineage))

# ---- 2. Load clinical (one row per sequence) ----------------------------------

clin_lin <- read_rds(paste0(storage_dir, "clin_lin.rds")) |>
  mutate(month = floor_date(SAMPLING_DATE, "month"),
         group = parent_group(sublineage))

# sample-volume context (not used in Jaccard itself, but needed to interpret it --
# a low Jaccard in a month with 1 clinical sequence means something different
# than in a month with 100)
ww_n_samples   <- sublin_meta |> distinct(month, Sample) |> count(month, name = "n_ww_samples")
clin_n_samples <- clin_lin |> count(month, name = "n_clin_samples")

# ---- 3. Monthly compositions at both resolutions ------------------------------

ww_fine    <- ww_monthly_composition(sublin_meta, sublineage)
ww_group   <- ww_monthly_composition(sublin_meta, group)
clin_fine  <- clin_monthly_composition(clin_lin, sublineage)
clin_group <- clin_monthly_composition(clin_lin, group)

# per-neighborhood WW compositions, compared against the same citywide clinical
# sets below -- the clinical source is state-level, so there is no
# neighborhood-resolved clinical reference to compare against
ww_fine_nb  <- ww_monthly_composition_by_loc(sublin_meta, sublineage)
ww_group_nb <- ww_monthly_composition_by_loc(sublin_meta, group)

ww_n_samples_nb <- sublin_meta |>
  distinct(LOCATION, month, Sample) |>
  count(LOCATION, month, name = "n_ww_samples")

# ---- 4. Panels ----------------------------------------------------------------

run_panel <- function(ww_comp, clin_comp, tag) {
  ww_sets   <- sets_from_comp(ww_comp)
  clin_sets <- sets_from_comp(clin_comp)

  jm <- jaccard_over_time(ww_sets, clin_sets, ww_n_samples, clin_n_samples)
  write_rds(jm, paste0(out_dir, "jaccard_monthly_citywide_", tag, ".rds"))

  # study-period summary: union across months of each source's monthly sets
  ww_all   <- unique(unlist(ww_sets$lineages))
  clin_all <- unique(unlist(clin_sets$lineages))
  cat(sprintf(
    "[%s/%s] Overall citywide WW vs. clinical Jaccard (full study period): %.3f (%d shared / %d union; %d WW-only, %d clinical-only)\n",
    concord_ver, tag,
    length(intersect(ww_all, clin_all)) / length(union(ww_all, clin_all)),
    length(intersect(ww_all, clin_all)), length(union(ww_all, clin_all)),
    length(setdiff(ww_all, clin_all)), length(setdiff(clin_all, ww_all))
  ))

  jm
}

# per-neighborhood: same metric, same clinical sets, one neighborhood at a time.
# The monthly threshold (min05 panel) is applied within neighborhood-month by
# the caller, so "0.05 of the composition" means that neighborhood's own.
run_panel_nb <- function(ww_comp_nb, clin_comp, tag) {
  clin_sets <- sets_from_comp(clin_comp)

  jm <- ww_comp_nb |>
    group_split(LOCATION) |>
    map_dfr(function(d) {
      loc <- unique(d$LOCATION)
      jaccard_over_time(
        sets_from_comp(d |> select(-LOCATION)),
        clin_sets,
        ww_n_samples_nb |> filter(LOCATION == loc) |> select(month, n_ww_samples),
        clin_n_samples
      ) |>
        mutate(LOCATION = loc, .before = 1)
    })

  write_rds(jm, paste0(out_dir, "jaccard_monthly_nb_", tag, ".rds"))

  cat(sprintf("[%s/%s] per-neighborhood median Jaccard (vs. citywide clinical):\n", concord_ver, tag))
  jm |>
    group_by(LOCATION) |>
    summarise(months = sum(!is.na(jaccard)), med = median(jaccard, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(med)) |>
    pwalk(function(LOCATION, months, med) cat(sprintf("    %-24s %2d months  %.3f\n", LOCATION, months, med)))

  jm
}

# A. full lineage resolution, everything present
jaccard_monthly_fine <- run_panel(ww_fine, clin_fine, "fine")

# B. parent-group resolution, everything present
jaccard_monthly_all <- run_panel(ww_group, clin_group, "all")

# C. parent-group resolution, >= monthly_min_prop that month (symmetric)
jaccard_monthly_min05 <- run_panel(
  apply_monthly_threshold(ww_group,   monthly_min_prop),
  apply_monthly_threshold(clin_group, monthly_min_prop),
  "min05"
)

# ---- 5. Same three panels, per neighborhood ----------------------------------

jaccard_monthly_nb_fine  <- run_panel_nb(ww_fine_nb,  clin_fine,  "fine")
jaccard_monthly_nb_all   <- run_panel_nb(ww_group_nb, clin_group, "all")
jaccard_monthly_nb_min05 <- run_panel_nb(
  ww_group_nb |>
    group_split(LOCATION) |>
    map_dfr(~ apply_monthly_threshold(.x |> select(-LOCATION), monthly_min_prop) |>
              mutate(LOCATION = unique(.x$LOCATION), .before = 1)),
  apply_monthly_threshold(clin_group, monthly_min_prop),
  "min05"
)
