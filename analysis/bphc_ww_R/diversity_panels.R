# Combined diversity figures. Reads the monthly tables written by
# jaccard_analysis.R and bray_curtis_analysis.R -- run those first if the
# inputs have changed.

# concord_v1 = no coverage QC, concord_v2 = coverage QC applied. Must match the
# version jaccard_analysis.R / bray_curtis_analysis.R last wrote.
concord_ver <- "concord_v2"

storage_dir <- "../data/"
in_dir      <- paste0(storage_dir, concord_ver, "/")
fig_dir     <- paste0("../draft_figures/", concord_ver, "/")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

qc_breadth_col <- "spike_breadth_d200"
qc_breadth_min <- 0.4
apply_qc <- identical(concord_ver, "concord_v2")

library(tidyverse)
library(lubridate)
library(patchwork)

clinical_decline_date <- ymd("2023-03-01")

col_clin <- "#D7191C"
col_jac  <- "#2C7BB6"
col_bc   <- "#D95F02"

# label row sits below the data; y limits are expanded to make room
label_y     <- -0.05
metric_lims <- c(-0.09, 1.05)

# every month tick-labelled, shared by all rows so gridlines line up
x_scale_monthly <- scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m")

qc_subtitle <- if (apply_qc) {
  paste0("coverage QC applied (", qc_breadth_col, " >= ", qc_breadth_min, ")")
} else {
  "no coverage QC -- every sequenced sample contributes"
}

x_text  <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
x_blank <- theme(axis.text.x = element_blank())

tiers <- list(
  fine  = "full resolution, all lineages present",
  all   = "collapsed to parent lineage, all present",
  min05 = "collapsed to parent lineage, >= 0.05 of each month's composition per source"
)

# ---- panel builders -----------------------------------------------------------

# shared decoration for the two metric rows. `lims` is widened only on the row
# carrying the count labels, so the other row has no dead space beneath it.
metric_decor <- function(p, lims) {
  p +
    annotate("segment", x = clinical_decline_date, xend = clinical_decline_date,
             y = 0, yend = 1.05, linetype = "dashed", color = "grey50") +
    scale_y_continuous(limits = lims, breaks = seq(0, 1, 0.25)) +
    x_scale_monthly +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
}

plot_clin_volume <- function(d) {
  d |>
    ggplot(aes(x = month, y = n_clin_samples)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey50") +
    annotate("text", x = clinical_decline_date, y = max(d$n_clin_samples),
             label = "clinical seq. volume drops", hjust = -0.03, vjust = 1,
             size = 3, color = "grey40") +
    geom_line(color = col_clin) +
    geom_point(color = col_clin, size = 1.2) +
    scale_y_continuous() +
    x_scale_monthly +
    labs(x = NULL, y = "# clinical\nsequences") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    x_blank
}

# the two metrics compare the same monthly sets, so the WW/clinical counts are
# identical on both rows -- printed once, with the Jaccard row
plot_jaccard <- function(d, show_x = FALSE, labels = TRUE) {
  p <- d |>
    ggplot(aes(x = month, y = jaccard)) +
    geom_line(color = col_jac) +
    geom_point(color = col_jac, size = 1.5)
  if (labels) {
    p <- p + geom_text(aes(label = paste0(n_ww, "/", n_clin)), y = label_y,
                       size = 2.2, color = "grey35")
  }
  metric_decor(p + labs(x = if (show_x) "Month" else NULL, y = "Jaccard index"),
               lims = if (labels) metric_lims else c(0, 1.05)) +
    (if (show_x) x_text else x_blank)
}

plot_bc <- function(d, show_x = TRUE, labels = FALSE) {
  p <- d |>
    ggplot(aes(x = month, y = bc_similarity)) +
    geom_line(color = col_bc) +
    geom_point(color = col_bc, size = 1.5)
  if (labels) {
    p <- p + geom_text(aes(label = paste0(n_ww, "/", n_clin)), y = label_y,
                       size = 2.2, color = "grey35")
  }
  metric_decor(p + labs(x = if (show_x) "Month" else NULL,
                        y = "Bray-Curtis\nsimilarity (1 - BC)"),
               lims = if (labels) metric_lims else c(0, 1.05)) +
    (if (show_x) x_text else x_blank)
}

# ---- citywide relative abundance (RAs.R collapse) -----------------------------

# Built from ww_collapsed.rds, i.e. RAs.R's labels: retrospective "ever clears
# 0.2" retention plus the manual KP merge. Deliberately NOT the diversity
# scripts' vocabulary -- the time-aware rule is kept for the RA visual only.
# Mirrors RAs.R's ww_citywide_weighted: renormalize within neighborhood-month,
# then population-weight across neighborhoods.
ra_monthly_citywide <- function() {
  ww_collapsed <- read_rds(paste0(storage_dir, "ww_collapsed.rds"))
  meta         <- read_rds(paste0(storage_dir, "meta_clean.rds"))
  pop_wts      <- read_rds(paste0(storage_dir, "pop_wts.rds"))

  d <- ww_collapsed |>
    left_join(meta |> distinct(FASTQ_ID, SAMPLING_DATE, across(all_of(qc_breadth_col))),
              by = c("Sample" = "FASTQ_ID")) |>
    mutate(month = floor_date(SAMPLING_DATE, "month"))

  # same sample-level QC as the metric rows, so both rows of the figure describe
  # the same set of samples
  if (apply_qc) d <- d |> filter(.data[[qc_breadth_col]] >= qc_breadth_min)

  n_loc <- d |> distinct(month, LOCATION, Sample) |> count(month, LOCATION, name = "n")

  loc <- d |>
    group_by(month, LOCATION, sublin_collapse) |>
    summarise(a = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    left_join(n_loc, by = c("month", "LOCATION")) |>
    mutate(p = a / n) |>
    group_by(month, LOCATION) |>
    mutate(p = p / sum(p, na.rm = TRUE)) |>
    ungroup() |>
    left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION")

  # denominator is the whole sampled population that month, not just the
  # neighborhoods carrying a given lineage -- otherwise a lineage seen in one
  # neighborhood is weighted as if it were citywide and the month sums past 1
  pop_month <- loc |>
    distinct(month, LOCATION, pop) |>
    group_by(month) |>
    summarise(pop_tot = sum(pop), .groups = "drop")

  loc |>
    group_by(month, sublin_collapse) |>
    summarise(num = sum(pop * p), .groups = "drop") |>
    left_join(pop_month, by = "month") |>
    mutate(prop = num / pop_tot) |>
    select(month, sublin_collapse, prop)
}

plot_ra <- function(ra) {
  lin_levels <- sort(unique(as.character(ra$sublin_collapse)))
  lin_levels <- c(setdiff(lin_levels, "other"), "other")  # force "other" last
  fill_cols  <- setNames(viridis::viridis(length(lin_levels), option = "turbo"), lin_levels)

  ra |>
    mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels)) |>
    ggplot(aes(x = month, y = prop, fill = sublin_collapse)) +
    geom_col(width = 24) +
    geom_text(aes(label = ifelse(prop > 0.10, as.character(sublin_collapse), "")),
              position = position_stack(vjust = 0.5),
              size = 1.7, color = "white", angle = 90) +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    x_scale_monthly +
    labs(x = NULL, y = "WW relative\nabundance") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    x_blank +
    guides(fill = "none")
}

# ---- A. one 3-row figure per resolution tier ----------------------------------

build_tier <- function(tag, tier_label) {
  jac <- read_rds(paste0(in_dir, "jaccard_monthly_citywide_", tag, ".rds"))
  bc  <- read_rds(paste0(in_dir, "bc_monthly_citywide_",      tag, ".rds"))
  stopifnot(identical(jac$month, bc$month))

  p <- plot_clin_volume(jac) / plot_jaccard(jac) / plot_bc(bc) +
    plot_layout(heights = c(1, 1.7, 1.7)) +
    plot_annotation(
      title = "Wastewater (citywide) vs. clinical lineage agreement",
      subtitle = paste0(tier_label, "; ", qc_subtitle,
                        "\nlabels below points: # lineages WW/clinical, monthly bins"),
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(paste0(fig_dir, "diversity_panel_", tag, ".jpg"), p,
         width = 10, height = 9, dpi = 300)

  cat(sprintf("[%s/%s] months: %d | median Jaccard %.3f | median BC similarity %.3f\n",
              concord_ver, tag, nrow(jac), median(jac$jaccard, na.rm = TRUE),
              median(bc$bc_similarity, na.rm = TRUE)))
}

iwalk(tiers, function(label, tag) build_tier(tag, label))

# ---- B. "all" tier with the RA composition row --------------------------------

# The RA row and the metric row use different lineage vocabularies on purpose
# (see ra_monthly_citywide above) -- called out in the subtitle so the two rows
# are not read as the same categories.
ra <- ra_monthly_citywide()

build_ra_combo <- function(metric) {
  jac <- read_rds(paste0(in_dir, "jaccard_monthly_citywide_all.rds"))
  bc  <- read_rds(paste0(in_dir, "bc_monthly_citywide_all.rds"))

  bottom <- if (metric == "jaccard") {
    plot_jaccard(jac, show_x = TRUE, labels = TRUE)
  } else {
    plot_bc(bc, show_x = TRUE, labels = TRUE)
  }

  p <- plot_clin_volume(jac) / plot_ra(ra) / bottom +
    plot_layout(heights = c(1, 2, 1.7)) +
    plot_annotation(
      title = "Wastewater (citywide) vs. clinical lineage agreement",
      subtitle = paste0(
        "bottom: ", if (metric == "jaccard") "Jaccard" else "Bray-Curtis",
        ", parent groups, all present (>= 0.05 per-month rule not applied)",
        "\nmiddle: RAs.R collapse -- retrospective 'ever clears 0.2' retention + manual KP merge,",
        " a different lineage vocabulary from the row below",
        "\n", qc_subtitle,
        "\nlabels below points: # lineages WW/clinical, monthly bins"),
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(paste0(fig_dir, "diversity_panel_all_ra_", metric, ".jpg"), p,
         width = 10, height = 10, dpi = 300)
}

walk(c("jaccard", "bc"), build_ra_combo)

cat(sprintf("RA row: %d months, %d lineage labels (RAs.R vocabulary)\n",
            n_distinct(ra$month), n_distinct(ra$sublin_collapse)))
