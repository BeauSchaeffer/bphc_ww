# Shared Freyja parsing, parent-group collapse, and monthly composition helpers.

library(data.table)
library(tidyverse)

# lineages grouped at their own prefix rather than the default 2-segment
# truncation (e.g. "BA.2.86" would otherwise fold into "BA.2"). This controls
# grouping only -- retention is a separate `force_keep` argument below.
force_keep_lineages <- c("BA.2.86")

### parse a freyja aggregate table into one row per (Sample, sublineage)

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

  ### code credit for parsing the demix output: https://github.com/a-roguet

  sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)

  ### freyja collapse lineage option yields e.g. "BA.4.6-like" / "BA.5.2-like2"
  sublineages.final |>
    mutate(sublineage = str_remove(as.character(sublineage), "-like[0-9]*$"))
}

### fx - extract parent group

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

### fx - retrospective "ever clears threshold" retention
# Used by the RA bar plots, where a lineage that later grows past the threshold
# should be visible during its early growth. NOT used by the diversity metrics,
# which threshold each month on its own. `threshold` is intentionally required.

keep_threshold_time <- function(dat, threshold,
                                time_col = year_epiweek,
                                suffix = ".x",
                                force_keep = character(0)) {

  dat |>
    mutate(group = parent_group(sublineage), # create parent lineage
           abundance = as.numeric(abundance),
           time = {{ time_col }}
    ) |>
    # sum RA for collapsed lineages within each sample,week
    group_by(Sample, time, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    # for each week,lineage group, take max RA observed in any sample that week
    group_by(time, group) |>
    summarise(max_ra_any_sample = max(ra, na.rm = TRUE), .groups = "drop") |>
    group_by(group) |>
    # uses future info intentionally (retrospective analysis, not live
    # surveillance): a group that ever clears threshold is named for its
    # entire history
    mutate(keep_by_time = any(max_ra_any_sample >= threshold) | group %in% force_keep) |>
    ungroup() |>
    transmute(time, group, keep_by_time,
              label = paste0(group, suffix))
}

### fx - WW detection threshold + per-sample renormalization
# drops near-zero demix noise, then restores each sample's abundances to sum 1
# so proportion-weighted comparisons aren't charged for the removed mass

apply_ww_detection <- function(df, thresh) {
  df |>
    filter(abundance >= thresh) |>
    group_by(Sample) |>
    mutate(abundance = abundance / sum(abundance, na.rm = TRUE)) |>
    ungroup()
}

### fx - monthly citywide WW composition -> long (month, lineage, proportion)
# pop_weighted = FALSE pools every sample that month equally; TRUE weights each
# neighborhood's own monthly composition by population.

ww_monthly_composition <- function(df, lineage_col, pop_weighted = FALSE, pop_wts = NULL) {
  df <- df |> rename(lineage = {{ lineage_col }})

  n_samples_month <- df |> distinct(month, Sample) |> count(month, name = "n_samples")

  if (!pop_weighted) {
    prop <- df |>
      group_by(month, lineage) |>
      summarise(abund_sum = sum(abundance, na.rm = TRUE), .groups = "drop") |>
      left_join(n_samples_month, by = "month") |>
      mutate(proportion = abund_sum / n_samples)
  } else {
    stopifnot(!is.null(pop_wts))

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

    prop <- grid |>
      left_join(loc_month_lineage, by = c("LOCATION", "month", "lineage")) |>
      mutate(abund_sum = replace_na(abund_sum, 0)) |>
      left_join(loc_month_n, by = c("LOCATION", "month")) |>
      mutate(p_loc = abund_sum / n_samples_loc) |>
      left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION") |>
      group_by(month, lineage) |>
      summarise(proportion = sum(pop * p_loc) / sum(pop), .groups = "drop") |>
      left_join(n_samples_month, by = "month")
  }

  prop |> select(month, lineage, proportion, n_samples)
}

### fx - monthly clinical composition -> long (month, lineage, proportion)
# one row per sequence in, so proportion is simply that lineage's share of the
# month's sequences

clin_monthly_composition <- function(df, lineage_col) {
  df <- df |> rename(lineage = {{ lineage_col }})

  n_samples_month <- df |> count(month, name = "n_samples")

  df |>
    count(month, lineage, name = "n_lineage") |>
    left_join(n_samples_month, by = "month") |>
    mutate(proportion = n_lineage / n_samples) |>
    select(month, lineage, proportion, n_samples)
}

### fx - per-month prevalence threshold, applied identically to either source
# other_label = NULL drops sub-threshold lineages (set-membership metrics);
# supplying a label instead pools them into one bucket so the month's
# proportions still sum to 1 (mass-weighted metrics).

apply_monthly_threshold <- function(comp, min_prop, other_label = NULL) {
  if (is.null(other_label)) {
    comp |> filter(proportion >= min_prop)
  } else {
    comp |>
      mutate(lineage = if_else(proportion >= min_prop, lineage, other_label)) |>
      group_by(month, lineage, n_samples) |>
      summarise(proportion = sum(proportion), .groups = "drop") |>
      select(month, lineage, proportion, n_samples)
  }
}
