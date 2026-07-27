# bray_curtis_analysis.R

# ---- PATHS (edit these) ------------------------------------------------------
storage_dir <- "../data/"
fig_dir     <- "../draft_figures/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(lubridate)

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

# threshold for the parent-group collapse below (mirrors jaccard_analysis.R)
collapse_threshold <- 0.05

clinical_decline_date <- ymd("2023-03-01")

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

# monthly WW citywide composition: named-vector-per-month (list-col) of
# lineage -> proportion
ww_monthly_composition <- function(df, lineage_col, pop_weighted) {
  df <- df |> rename(lineage = {{ lineage_col }})

  n_samples_month <- df |> distinct(month, Sample) |> count(month, name = "n_samples")

  if (!pop_weighted) {
    prop <- df |>
      group_by(month, lineage) |>
      summarise(abund_sum = sum(abundance, na.rm = TRUE), .groups = "drop") |>
      left_join(n_samples_month, by = "month") |>
      mutate(proportion = abund_sum / n_samples)
  } else {
    pop_wts <- read_rds(paste0(storage_dir, "pop_wts.rds"))

    loc_month_n <- df |> distinct(LOCATION, month, Sample) |> count(LOCATION, month, name = "n_samples_loc")

    lineages_by_month <- df |> distinct(month, lineage)

    # complete grid: every (LOCATION, month) that has WW samples, crossed with
    # every lineage seen citywide that month -- so lineages absent from a
    # given neighborhood contribute an explicit 0, not a missing row
    grid <- loc_month_n |>
      select(LOCATION, month) |>
      inner_join(lineages_by_month, by = "month")

    loc_month_lineage <- df |>
      group_by(LOCATION, month, lineage) |>
      summarise(abund_sum = sum(abundance, na.rm = TRUE), .groups = "drop")

    p_loc <- grid |>
      left_join(loc_month_lineage, by = c("LOCATION", "month", "lineage")) |>
      mutate(abund_sum = replace_na(abund_sum, 0)) |>
      left_join(loc_month_n, by = c("LOCATION", "month")) |>
      mutate(p_loc = abund_sum / n_samples_loc) |>
      left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION")

    prop <- p_loc |>
      group_by(month, lineage) |>
      summarise(proportion = sum(pop * p_loc) / sum(pop), .groups = "drop") |>
      left_join(n_samples_month, by = "month")
  }

  prop |>
    group_by(month) |>
    summarise(composition = list(setNames(proportion, lineage)),
              n_samples = first(n_samples),
              .groups = "drop")
}

clin_monthly_composition <- function(df, lineage_col) {
  df <- df |> rename(lineage = {{ lineage_col }})

  n_samples_month <- df |> count(month, name = "n_samples")

  df |>
    count(month, lineage, name = "n_lineage") |>
    left_join(n_samples_month, by = "month") |>
    mutate(proportion = n_lineage / n_samples) |>
    group_by(month) |>
    summarise(composition = list(setNames(proportion, lineage)),
              n_samples = first(n_samples),
              .groups = "drop")
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

plot_bc_trend <- function(bc_monthly, title, subtitle) {
  p_bc <- bc_monthly |>
    ggplot(aes(x = month, y = bc_similarity)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey50") +
    annotate("text", x = clinical_decline_date, y = 1.02, label = "clinical seq. volume drops",
             hjust = 0, vjust = 0, size = 3, color = "grey40") +
    geom_line(color = "#D95F02") +
    geom_point(aes(size = n_clin_samples), color = "#D95F02") +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    scale_size_continuous(name = "# clinical\nsequences", range = c(0.5, 4)) +
    labs(x = NULL, y = "Bray-Curtis similarity (1 - BC)", title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_samples <- bc_monthly |>
    pivot_longer(c(n_ww_samples, n_clin_samples), names_to = "source", values_to = "n") |>
    mutate(source = recode(source, n_ww_samples = "WW", n_clin_samples = "Clinical")) |>
    ggplot(aes(x = month, y = n, color = source)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey50") +
    geom_line() +
    geom_point(size = 1) +
    scale_color_manual(values = c(WW = "#2C7BB6", Clinical = "#D7191C"), name = NULL) +
    labs(x = "Month", y = "# sequenced samples") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

  patchwork::wrap_plots(p_bc, p_samples, ncol = 1, heights = c(2, 1))
}

# ---- 1. Load + parse WW (full resolution) --------------------------------------

ww_metadata <- read_rds(paste0(storage_dir, "meta_clean.rds"))

parse_freyja_aggregate <- function(path) {
  results <- read.table(path, fill = TRUE, sep = "\t", h = TRUE)
  results <- as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", "")))
  results <- as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x))))

  for (i in 1:nrow(results)) {
    lineages.temp   <- as.data.frame(t(setDT(tstrsplit(as.character(results[i, 3]), " ", fixed = TRUE))[]))
    abundances.temp <- as.data.frame(t(setDT(tstrsplit(as.character(results[i, 4]), " ", fixed = TRUE))[]))
    sample.temp <- rep(results[i, 1], nrow(lineages.temp))
    if (i == 1) {
      sublineages.final <- cbind(sample.temp, lineages.temp, abundances.temp)
    } else {
      sublineages.final <- rbind(sublineages.final, cbind(sample.temp, lineages.temp, abundances.temp))
    }
  }
  names(sublineages.final) <- c("Sample", "sublineage", "abundance")

  sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)

  sublineages.final |>
    mutate(sublineage = str_remove(as.character(sublineage), "-like[0-9]*$"))
}

sublineages.final <- bind_rows(
  parse_freyja_aggregate("../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv"),
  parse_freyja_aggregate("../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv")
)

sublin_meta <- left_join(sublineages.final, ww_metadata, by = c("Sample" = "FASTQ_ID")) |>
  mutate(abundance = as.numeric(abundance),
         month = floor_date(SAMPLING_DATE, "month"))

# ---- 2. Load clinical (full resolution, one row per sequence) -----------------

clin_lin <- read_rds(paste0(storage_dir, "clin_lin.rds")) |>
  mutate(month = floor_date(SAMPLING_DATE, "month"))

# =============================================================================
# A. Uncollapsed (full-resolution) lineages
# =============================================================================

ww_comp_raw   <- ww_monthly_composition(sublin_meta, sublineage, ww_pop_weighted)
clin_comp_raw <- clin_monthly_composition(clin_lin, sublineage)

bc_monthly_raw <- bc_over_time(ww_comp_raw, clin_comp_raw)

write_rds(bc_monthly_raw, paste0(storage_dir, "bc_monthly_citywide_raw.rds"))

overall_ww_raw   <- sublin_meta |> group_by(sublineage) |> summarise(a = sum(abundance)) |> deframe()
overall_clin_raw <- clin_lin |> count(sublineage) |> deframe()
overall_clin_raw <- overall_clin_raw / sum(overall_clin_raw)
overall_ww_raw   <- overall_ww_raw / sum(overall_ww_raw)

cat(sprintf(
  "[raw] Overall citywide WW vs. clinical Bray-Curtis dissimilarity (full study period): %.3f (similarity = %.3f)\n",
  bray_curtis(overall_ww_raw, overall_clin_raw),
  1 - bray_curtis(overall_ww_raw, overall_clin_raw)
))

p_bc_raw <- plot_bc_trend(
  bc_monthly_raw,
  title = "Lineage proportion similarity: wastewater (citywide) vs. clinical (uncollapsed)",
  subtitle = paste0("full-resolution lineages, WW citywide = ",
                     if (ww_pop_weighted) "population-weighted" else "pooled (unweighted)",
                     ", monthly bins")
)

ggsave(paste0(fig_dir, "bc_citywide_raw.jpg"), p_bc_raw, width = 10, height = 8, dpi = 300)

# =============================================================================
# B. Collapsed lineages (parent-group + time-aware threshold, as in RAs.R)
# =============================================================================

force_keep_lineages <- c("BA.2.86")

parent_group <- function(x, extended_groups = force_keep_lineages) {
  default <- str_extract(x, "^[A-Za-z]+\\.[0-9]+")
  default <- ifelse(is.na(default), x, default)  # keep full token if no match (recombinants)

  group <- default
  for (g in extended_groups) {
    is_g <- !is.na(x) & (x == g | startsWith(x, paste0(g, ".")))
    group[is_g] <- g
  }
  group
}

keep_threshold_time <- function(dat, threshold = collapse_threshold,
                                time_col = year_epiweek,
                                suffix = ".x",
                                force_keep = character(0)) {

  dat |>
    mutate(group = parent_group(sublineage),
           abundance = as.numeric(abundance),
           time = {{ time_col }}
    ) |>
    group_by(Sample, time, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    group_by(time, group) |>
    summarise(max_ra_any_sample = max(ra, na.rm = TRUE), .groups = "drop") |>
    group_by(group) |>
    mutate(keep_by_time = any(max_ra_any_sample >= threshold) | group %in% force_keep) |>
    ungroup() |>
    transmute(time, group, keep_by_time,
              label = paste0(group, suffix))
}

# pseudo-aggregate clinical into weekly counts-as-abundance, needed only to
# feed keep_threshold_time() the same (Sample, time, abundance) shape as WW
# -- identical construction to RAs.R's `clinical` object. Composition itself
# is computed from the raw per-sequence clin_lin rows below, not from this.
clinical <- clin_lin |>
  mutate(parentsub = parent_group(sublineage)) |>
  group_by(year_epiweek, parentsub) |>
  summarise(count = n()) |>
  ungroup() |>
  group_by(year_epiweek) |>
  mutate(abundance = count / sum(count)) |>
  ungroup() |>
  rename(sublineage = parentsub) |>
  mutate(year = year_epiweek %/% 100,
         week = year_epiweek %% 100,
         Sample = "clinical_agg",
         LOCATION = "clinical",
         epiweek = week)

combined_for_keep <- bind_rows(
  sublin_meta |> select(Sample, sublineage, abundance, year_epiweek),
  clinical   |> select(Sample, sublineage, abundance, year_epiweek)
)

keep_tbl <- keep_threshold_time(combined_for_keep,
                                threshold = collapse_threshold,
                                time_col = year_epiweek,
                                suffix = ".x",
                                force_keep = force_keep_lineages)

# assign each raw row its collapsed label directly (rather than going through
# collapse_sublineages_timeaware's per-Sample summarise step): summing raw
# per-lineage abundances straight into (month, label) totals is equivalent to
# summing per-sample collapsed totals into the same (sum is associative), and
# this keeps each source's real per-row SAMPLING_DATE/Sample/LOCATION intact
# for ww_monthly_composition()/clin_monthly_composition() to use directly --
# no epiweek-based month reconstruction needed, no pseudo-sample-count bugs.
assign_collapsed_label <- function(df) {
  df |>
    mutate(group = parent_group(sublineage)) |>
    left_join(keep_tbl |> select(time, group, keep_by_time, label),
              by = c("year_epiweek" = "time", "group")) |>
    mutate(
      keep_by_time = as.logical(coalesce(keep_by_time, FALSE)),
      sublin_collapse = if_else(keep_by_time, label, "other")
    )
}

sublin_meta_collapsed <- assign_collapsed_label(sublin_meta)
clin_lin_collapsed    <- assign_collapsed_label(clin_lin)

# unlike the Jaccard version, "other" is kept here rather than dropped: it
# still carries real abundance mass that both systems detected but couldn't
# name individually, and Bray-Curtis (a proportion-weighted metric) needs
# that mass accounted for or the remaining proportions won't reflect true
# composition. It is compared as its own bucket like any other lineage --
# i.e. two systems both having a lot of unnamed "other" mass counts as
# agreement, which is a limitation to keep in mind.
ww_comp_collapsed   <- ww_monthly_composition(sublin_meta_collapsed, sublin_collapse, ww_pop_weighted)
clin_comp_collapsed <- clin_monthly_composition(clin_lin_collapsed, sublin_collapse)

bc_monthly_collapsed <- bc_over_time(ww_comp_collapsed, clin_comp_collapsed)

write_rds(bc_monthly_collapsed, paste0(storage_dir, "bc_monthly_citywide_collapsed.rds"))

overall_ww_collapsed   <- sublin_meta_collapsed |> group_by(sublin_collapse) |> summarise(a = sum(abundance)) |> deframe()
overall_clin_collapsed <- clin_lin_collapsed |> count(sublin_collapse) |> deframe()
overall_ww_collapsed   <- overall_ww_collapsed / sum(overall_ww_collapsed)
overall_clin_collapsed <- overall_clin_collapsed / sum(overall_clin_collapsed)

cat(sprintf(
  "[collapsed] Overall citywide WW vs. clinical Bray-Curtis dissimilarity (full study period): %.3f (similarity = %.3f)\n",
  bray_curtis(overall_ww_collapsed, overall_clin_collapsed),
  1 - bray_curtis(overall_ww_collapsed, overall_clin_collapsed)
))

p_bc_collapsed <- plot_bc_trend(
  bc_monthly_collapsed,
  title = "Lineage proportion similarity: wastewater (citywide) vs. clinical (collapsed)",
  subtitle = paste0("parent-group collapse (threshold = ", collapse_threshold, "), WW citywide = ",
                     if (ww_pop_weighted) "population-weighted" else "pooled (unweighted)",
                     ", monthly bins")
)

ggsave(paste0(fig_dir, "bc_citywide_collapsed.jpg"), p_bc_collapsed, width = 10, height = 8, dpi = 300)
