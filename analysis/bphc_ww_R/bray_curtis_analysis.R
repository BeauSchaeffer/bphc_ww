# bray_curtis_analysis.R

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

# TOGGLE: how to combine neighborhoods into a single WW citywide composition
# each month.
#   FALSE -> pool all WW samples that month equally (each sample counted
#            once, regardless of neighborhood) -- matches the unweighted
#            "detected anywhere" citywide definition used for the Jaccard
#            index in jaccard_analysis.R
#   TRUE  -> population-weight neighborhoods (each neighborhood's monthly
#            composition weighted by pop_wts.rds, as in RAs.R's
#            ww_citywide_weighted)
ww_pop_weighted <- FALSE

# WW instrument noise floor, paired with a per-sample renormalization inside
# apply_ww_detection() -- BC charges a mass deficit as dissimilarity even when
# the compositions match in shape, so the removed mass has to be restored.
ww_detect_thresh <- 0.01

# coverage QC, matching RAs.R. Unlike the modeling pipeline -- which nulls
# meaneffCopiesL because its completed grid would otherwise read a failure as a
# confident zero -- the monthly compositions here are built only from observed
# samples, so a failing sample is simply dropped and contributes nothing.
qc_breadth_col <- "spike_breadth_d200"
qc_breadth_min <- 0.4
apply_qc <- identical(concord_ver, "concord_v2")

# per-month prevalence threshold for the min05 panel (mirrors
# jaccard_analysis.R), applied identically to WW and clinical
monthly_min_prop <- 0.05

# ---- Helpers --------------------------------------------------------------

# Bray-Curtis dissimilarity between two named proportion vectors (named by
# lineage, each summing to ~1). BC(A,B) = sum(|pA-pB|) / sum(pA+pB).
# Returns NA if either side has zero total mass (composition undefined).
bray_curtis <- function(p_ww, p_clin) {
  lineages <- union(names(p_ww), names(p_clin))
  a <- setNames(rep(0, length(lineages)), lineages)
  b <- a
  a[names(p_ww)]   <- p_ww
  b[names(p_clin)] <- p_clin
  if (sum(a) == 0 || sum(b) == 0) return(NA_real_)
  sum(abs(a - b)) / sum(a + b)
}

# long monthly composition -> one named proportion vector per month
comp_to_vectors <- function(comp) {
  comp |>
    group_by(month) |>
    summarise(composition = list(setNames(proportion, lineage)),
              n_samples = first(n_samples), .groups = "drop")
}

bc_over_time <- function(ww_comp, clin_comp) {
  months <- sort(unique(c(ww_comp$month, clin_comp$month)))

  tibble(month = months) |>
    left_join(ww_comp   |> select(month, ww_composition = composition, n_ww_samples = n_samples),   by = "month") |>
    left_join(clin_comp |> select(month, clin_composition = composition, n_clin_samples = n_samples), by = "month") |>
    rowwise() |>
    mutate(
      ww_composition   = list(ww_composition   %||% setNames(numeric(0), character(0))),
      clin_composition = list(clin_composition %||% setNames(numeric(0), character(0))),
      # named lineages contributing to the month's composition; the "other"
      # bucket in the min05 panel is mass accounting, not a lineage, so it is
      # excluded from the count
      n_ww   = sum(names(ww_composition)   != "other"),
      n_clin = sum(names(clin_composition) != "other"),
      bray_curtis = bray_curtis(ww_composition, clin_composition),
      bc_similarity = 1 - bray_curtis
    ) |>
    ungroup() |>
    select(-ww_composition, -clin_composition) |>
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

pop_wts <- if (ww_pop_weighted) read_rds(paste0(storage_dir, "pop_wts.rds")) else NULL

# ---- 3. Monthly compositions at both resolutions ------------------------------

ww_fine    <- ww_monthly_composition(sublin_meta, sublineage, ww_pop_weighted, pop_wts)
ww_group   <- ww_monthly_composition(sublin_meta, group,      ww_pop_weighted, pop_wts)
clin_fine  <- clin_monthly_composition(clin_lin, sublineage)
clin_group <- clin_monthly_composition(clin_lin, group)

# ---- 4. Panels ----------------------------------------------------------------

# study-period summaries, pooled across all months. The min05 panel applies the
# same threshold rule to the pooled composition that the monthly series applies
# month by month.
pooled_ww   <- function(df, col) { v <- df |> group_by({{ col }}) |> summarise(a = sum(abundance)) |> deframe(); v / sum(v) }
pooled_clin <- function(df, col) { v <- df |> count({{ col }}) |> deframe(); v / sum(v) }

pool_threshold <- function(v, min_prop) {
  keep <- v >= min_prop
  c(v[keep], other = sum(v[!keep]))
}

run_panel <- function(ww_comp, clin_comp, overall_ww, overall_clin, tag) {
  bm <- bc_over_time(comp_to_vectors(ww_comp), comp_to_vectors(clin_comp))
  write_rds(bm, paste0(out_dir, "bc_monthly_citywide_", tag, ".rds"))

  cat(sprintf(
    "[%s/%s] Overall citywide WW vs. clinical Bray-Curtis dissimilarity (full study period): %.3f (similarity = %.3f)\n",
    concord_ver, tag, bray_curtis(overall_ww, overall_clin), 1 - bray_curtis(overall_ww, overall_clin)
  ))

  bm
}

# A. full lineage resolution, everything present
bc_monthly_fine <- run_panel(
  ww_fine, clin_fine,
  pooled_ww(sublin_meta, sublineage), pooled_clin(clin_lin, sublineage), "fine"
)

# B. parent-group resolution, everything present
bc_monthly_all <- run_panel(
  ww_group, clin_group,
  pooled_ww(sublin_meta, group), pooled_clin(clin_lin, group), "all"
)

# C. parent-group resolution, >= monthly_min_prop that month (symmetric).
# Sub-threshold mass is pooled into "other" rather than dropped: BC is
# proportion-weighted, so removing mass without accounting for it reads as
# dissimilarity. Two sources both carrying a large unnamed remainder therefore
# count as agreement, which is a limitation to keep in mind.
bc_monthly_min05 <- run_panel(
  apply_monthly_threshold(ww_group,   monthly_min_prop, other_label = "other"),
  apply_monthly_threshold(clin_group, monthly_min_prop, other_label = "other"),
  pool_threshold(pooled_ww(sublin_meta, group), monthly_min_prop),
  pool_threshold(pooled_clin(clin_lin, group), monthly_min_prop), "min05"
)
