# jaccard_analysis.R

# ---- PATHS (edit these) ------------------------------------------------------
storage_dir <- "../data/"
fig_dir     <- "../draft_figures/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(lubridate)

# ---- Params -------------------------------------------------------------------

# minimum per-sample Freyja abundance for a lineage to count as "detected" in
# WW at full (uncollapsed) resolution -- excludes near-zero barcode calls that
# are typically demix noise rather than real signal (clinical calls have no
# analogous threshold: each GISAID sequence is a single definitive lineage
# call)
ww_detect_thresh <- 0.01

# threshold for the parent-group collapse below (RAs.R default is 0.2)
collapse_threshold <- 0.05

clinical_decline_date <- ymd("2023-03-01")

# ---- Helpers --------------------------------------------------------------

# distinct (month, lineage) detections -> one row per month with a list-col
# of the lineages seen that month
monthly_lineage_sets <- function(df, month_col, lineage_col) {
  df |>
    transmute(month = {{ month_col }}, lineage = {{ lineage_col }}) |>
    distinct() |>
    group_by(month) |>
    summarise(lineages = list(unique(lineage)), n_lineages = n(), .groups = "drop")
}

# monthly Jaccard between two monthly_lineage_sets() tables, joined against
# monthly sample-count tables (source_n_samples cols) purely for plotting
# context
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
      jaccard  = if_else(n_union > 0, n_shared / n_union, NA_real_)
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

plot_jaccard_trend <- function(jaccard_monthly, title, subtitle) {
  p_jaccard <- jaccard_monthly |>
    ggplot(aes(x = month, y = jaccard)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey50") +
    annotate("text", x = clinical_decline_date, y = 1.02, label = "clinical seq. volume drops",
             hjust = 0, vjust = 0, size = 3, color = "grey40") +
    geom_line(color = "#2C7BB6") +
    geom_point(aes(size = n_clin_samples), color = "#2C7BB6") +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    scale_size_continuous(name = "# clinical\nsequences", range = c(0.5, 4)) +
    labs(x = NULL, y = "Jaccard index", title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_samples <- jaccard_monthly |>
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

  patchwork::wrap_plots(p_jaccard, p_samples, ncol = 1, heights = c(2, 1))
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
  mutate(abundance = as.numeric(abundance))

# ---- 2. Load clinical (full resolution, one row per sequence) -----------------

clin_lin <- read_rds(paste0(storage_dir, "clin_lin.rds"))

# sample-volume context (not used in Jaccard itself, but needed to interpret it --
# a low Jaccard in a month with 1 clinical sequence means something different
# than in a month with 100)
ww_n_samples <- sublin_meta |>
  mutate(month = floor_date(SAMPLING_DATE, "month")) |>
  distinct(month, Sample) |>
  count(month, name = "n_ww_samples")

clin_n_samples <- clin_lin |>
  mutate(month = floor_date(SAMPLING_DATE, "month")) |>
  count(month, name = "n_clin_samples")

# =============================================================================
# A. Uncollapsed (full-resolution) lineages
# =============================================================================

ww_monthly_raw <- sublin_meta |>
  filter(abundance >= ww_detect_thresh) |>
  mutate(month = floor_date(SAMPLING_DATE, "month")) |>
  monthly_lineage_sets(month, sublineage)

clin_monthly_raw <- clin_lin |>
  mutate(month = floor_date(SAMPLING_DATE, "month")) |>
  monthly_lineage_sets(month, sublineage)

jaccard_monthly_raw <- jaccard_over_time(ww_monthly_raw, clin_monthly_raw, ww_n_samples, clin_n_samples)

write_rds(jaccard_monthly_raw, paste0(storage_dir, "jaccard_monthly_citywide_raw.rds"))

ww_all_lineages_raw   <- sublin_meta |> filter(abundance >= ww_detect_thresh) |> pull(sublineage) |> unique()
clin_all_lineages_raw <- clin_lin |> pull(sublineage) |> unique()

overall_jaccard_raw <- length(intersect(ww_all_lineages_raw, clin_all_lineages_raw)) /
  length(union(ww_all_lineages_raw, clin_all_lineages_raw))

cat(sprintf(
  "[raw] Overall citywide WW vs. clinical Jaccard (full study period): %.3f (%d shared / %d union; %d WW-only, %d clinical-only)\n",
  overall_jaccard_raw,
  length(intersect(ww_all_lineages_raw, clin_all_lineages_raw)),
  length(union(ww_all_lineages_raw, clin_all_lineages_raw)),
  length(setdiff(ww_all_lineages_raw, clin_all_lineages_raw)),
  length(setdiff(clin_all_lineages_raw, ww_all_lineages_raw))
))

p_combined_raw <- plot_jaccard_trend(
  jaccard_monthly_raw,
  title = "Lineage set overlap: wastewater (citywide) vs. clinical (uncollapsed)",
  subtitle = paste0("full-resolution lineages, WW detection threshold = ", ww_detect_thresh,
                     " abundance, monthly bins")
)

ggsave(paste0(fig_dir, "diversity_jaccard_citywide_raw.jpg"), p_combined_raw, width = 10, height = 8, dpi = 300)

# =============================================================================
# B. Collapsed lineages (parent-group + time-aware threshold, as in RAs.R)
# =============================================================================

# lineages always retained, and grouped at their own prefix rather than the
# default 2-segment truncation -- identical to RAs.R
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

collapse_sublineages_timeaware <- function(df, threshold = collapse_threshold,
                                           other = "other",
                                           time_col = year_epiweek,
                                           suffix = ".x",
                                           keep_tbl = NULL,
                                           force_keep = character(0)) {

  if (is.null(keep_tbl)) {
    keep_tbl <- keep_threshold_time(df, threshold = threshold,
                                    time_col = {{ time_col }}, suffix = suffix,
                                    force_keep = force_keep)
  }

  df |>
    as_tibble() |>
    mutate(
      sublineage = coalesce(sublineage, ""),
      group = parent_group(sublineage),
      abundance = as.numeric(abundance),
      time = {{ time_col }}
    ) |>
    left_join(keep_tbl |> select(time, group, keep_by_time, label),
              by = c("time", "group")) |>
    mutate(
      keep_by_time = as.logical(coalesce(keep_by_time, FALSE)),
      sublin_collapse = if_else(keep_by_time, label, other),
      sublin_collapse = factor(sublin_collapse)
    ) |>
    group_by(Sample, time, LOCATION, epiweek, year, sublin_collapse) |>
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
}

# pseudo-aggregate clinical into weekly counts-as-abundance so it can run
# through the same per-sample collapse machinery as WW -- identical
# construction to RAs.R's `clinical` object
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

ww_collapsed <- collapse_sublineages_timeaware(sublin_meta,
                                               threshold = collapse_threshold,
                                               time_col = year_epiweek,
                                               suffix = ".x",
                                               keep_tbl = keep_tbl)

clin_collapsed <- collapse_sublineages_timeaware(clinical,
                                                 threshold = collapse_threshold,
                                                 time_col = year_epiweek,
                                                 suffix = ".x",
                                                 keep_tbl = keep_tbl)

# collapse_sublineages_timeaware() summarises away SAMPLING_DATE, so
# reconstruct a month from epiweek/year (same epiweek -> date convention
# RAs.R uses for the clinical pseudo-dates) -- applied identically to both
# sources so the two months line up
epiweek_to_month <- function(year, epiweek) {
  floor_date(ymd(paste0(year, "-01-01")) + weeks(epiweek - 1), "month")
}

ww_collapsed   <- ww_collapsed   |> mutate(month = epiweek_to_month(year, epiweek))
clin_collapsed <- clin_collapsed |> mutate(month = epiweek_to_month(year, epiweek))

# "other" is a catch-all bucket, not a real shared lineage identity -- a WW
# "other" and a clinical "other" in the same month are not the same lineage,
# so exclude it from set membership rather than counting it as overlap
ww_monthly_collapsed <- ww_collapsed |>
  filter(sublin_collapse != "other") |>
  monthly_lineage_sets(month, sublin_collapse)

clin_monthly_collapsed <- clin_collapsed |>
  filter(sublin_collapse != "other") |>
  monthly_lineage_sets(month, sublin_collapse)

jaccard_monthly_collapsed <- jaccard_over_time(ww_monthly_collapsed, clin_monthly_collapsed,
                                               ww_n_samples, clin_n_samples)

write_rds(jaccard_monthly_collapsed, paste0(storage_dir, "jaccard_monthly_citywide_collapsed.rds"))

ww_all_lineages_collapsed   <- ww_collapsed   |> filter(sublin_collapse != "other") |> pull(sublin_collapse) |> unique() |> as.character()
clin_all_lineages_collapsed <- clin_collapsed |> filter(sublin_collapse != "other") |> pull(sublin_collapse) |> unique() |> as.character()

overall_jaccard_collapsed <- length(intersect(ww_all_lineages_collapsed, clin_all_lineages_collapsed)) /
  length(union(ww_all_lineages_collapsed, clin_all_lineages_collapsed))

cat(sprintf(
  "[collapsed] Overall citywide WW vs. clinical Jaccard (full study period): %.3f (%d shared / %d union; %d WW-only, %d clinical-only)\n",
  overall_jaccard_collapsed,
  length(intersect(ww_all_lineages_collapsed, clin_all_lineages_collapsed)),
  length(union(ww_all_lineages_collapsed, clin_all_lineages_collapsed)),
  length(setdiff(ww_all_lineages_collapsed, clin_all_lineages_collapsed)),
  length(setdiff(clin_all_lineages_collapsed, ww_all_lineages_collapsed))
))

p_combined_collapsed <- plot_jaccard_trend(
  jaccard_monthly_collapsed,
  title = "Lineage set overlap: wastewater (citywide) vs. clinical (collapsed)",
  subtitle = paste0("parent-group collapse (RAs.R logic, no manual merges), threshold = ",
                     collapse_threshold, ", 'other' excluded, monthly bins")
)

ggsave(paste0(fig_dir, "diversity_jaccard_citywide_collapsed.jpg"), p_combined_collapsed, width = 10, height = 8, dpi = 300)
