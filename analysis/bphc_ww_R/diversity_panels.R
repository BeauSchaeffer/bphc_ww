# Combined diversity figures. Reads the monthly tables written by
# jaccard_analysis.R and bray_curtis_analysis.R -- run those first if the
# inputs have changed.

# concord_v1 = no coverage QC, concord_v2 = coverage QC applied. Must match the
# version jaccard_analysis.R / bray_curtis_analysis.R last wrote.
concord_ver <- "concord_v2"

library(tidyverse)
library(lubridate)
library(patchwork)

storage_dir <- "../data/"
in_dir      <- paste0(storage_dir, concord_ver, "/")
fig_dir     <- paste0("../draft_figures/", concord_ver, "/agg/")
nb_fig_dir  <- paste0("../draft_figures/", concord_ver, "/nb/")
walk(c(fig_dir, nb_fig_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

qc_breadth_col <- "spike_breadth_d200"
qc_breadth_min <- 0.4
apply_qc <- identical(concord_ver, "concord_v2")

clinical_decline_date <- ymd("2023-03-01")

col_clin <- "#D7191C"
col_jac  <- "#2C7BB6"
col_bc   <- "#D95F02"
col_ww_n <- "#4D4D4D"
col_ref  <- "grey55"   # citywide series when it is the reference, not the subject

# thin translucent per-neighborhood series drawn beneath the citywide one
nb_alpha <- 0.28
nb_lw    <- 0.3
nb_size  <- 0.6

# the clinical reference is state-level for every panel here -- there is no
# neighborhood-resolved clinical data, so each neighborhood's WW composition is
# compared against the same citywide clinical composition
# pre-wrapped: these figures are 9-10 in wide, so lines past ~110 characters clip
nb_caveat <- paste0(
  "every neighborhood is compared against the same citywide (MA) clinical composition -- there is no",
  "\nneighborhood-resolved clinical data. Neighborhood-months rest on 1-5 samples vs. ~30 citywide, so",
  "\ntheir lineage richness -- and with it Jaccard -- is structurally lower.")

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
# `nb` draws the per-neighborhood series as thin translucent lines beneath the
# subject series; `ref` draws a single series as a thin grey reference above the
# nb layer but below the subject. Both optional, so the citywide figures and the
# per-neighborhood figures share one builder each.
metric_layers <- function(p, y_col, nb, ref, col) {
  if (!is.null(nb)) {
    p <- p +
      geom_line(data = nb, aes(y = .data[[y_col]], group = LOCATION),
                color = col, alpha = nb_alpha, linewidth = nb_lw) +
      geom_point(data = nb, aes(y = .data[[y_col]]),
                 color = col, alpha = nb_alpha, size = nb_size)
  }
  if (!is.null(ref)) {
    p <- p + geom_line(data = ref, aes(y = .data[[y_col]]),
                       color = col_ref, linewidth = 0.4)
  }
  p + geom_line(color = col) + geom_point(color = col, size = 1.5)
}

plot_jaccard <- function(d, show_x = FALSE, labels = TRUE, nb = NULL, ref = NULL) {
  p <- metric_layers(ggplot(d, aes(x = month, y = jaccard)), "jaccard", nb, ref, col_jac)
  if (labels) {
    p <- p + geom_text(data = d |> filter(n_ww > 0),
                       aes(label = paste0(n_ww, "/", n_clin)), y = label_y,
                       size = 2.2, color = "grey35")
  }
  metric_decor(p + labs(x = if (show_x) "Month" else NULL, y = "Jaccard index"),
               lims = if (labels) metric_lims else c(0, 1.05)) +
    (if (show_x) x_text else x_blank)
}

plot_bc <- function(d, show_x = TRUE, labels = FALSE, nb = NULL, ref = NULL) {
  p <- metric_layers(ggplot(d, aes(x = month, y = bc_similarity)), "bc_similarity", nb, ref, col_bc)
  if (labels) {
    p <- p + geom_text(data = d |> filter(n_ww > 0),
                       aes(label = paste0(n_ww, "/", n_clin)), y = label_y,
                       size = 2.2, color = "grey35")
  }
  metric_decor(p + labs(x = if (show_x) "Month" else NULL,
                        y = "Bray-Curtis\nsimilarity (1 - BC)"),
               lims = if (labels) metric_lims else c(0, 1.05)) +
    (if (show_x) x_text else x_blank)
}

# that neighborhood's own monthly WW sample count -- 1 vs. 5 samples moves the
# metrics enough that the rows above can't be read without it
plot_ww_volume <- function(d) {
  d |>
    ggplot(aes(x = month, y = n_ww_samples)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey50") +
    geom_col(fill = col_ww_n, width = 20) +
    scale_y_continuous(breaks = scales::pretty_breaks(3),
                       expand = expansion(mult = c(0, 0.05))) +
    x_scale_monthly +
    labs(x = NULL, y = "# WW\nsamples") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    x_blank
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

# read once per tier; the nb tables carry a LOCATION column, the citywide ones
# don't
read_tier <- function(tag) {
  list(
    jac    = read_rds(paste0(in_dir, "jaccard_monthly_citywide_", tag, ".rds")),
    bc     = read_rds(paste0(in_dir, "bc_monthly_citywide_",      tag, ".rds")),
    jac_nb = read_rds(paste0(in_dir, "jaccard_monthly_nb_",       tag, ".rds")),
    bc_nb  = read_rds(paste0(in_dir, "bc_monthly_nb_",            tag, ".rds"))
  )
}

build_tier <- function(tag, tier_label) {
  d <- read_tier(tag)
  stopifnot(identical(d$jac$month, d$bc$month))

  p <- plot_clin_volume(d$jac) /
    plot_jaccard(d$jac, nb = d$jac_nb) /
    plot_bc(d$bc, nb = d$bc_nb) +
    plot_layout(heights = c(1, 1.7, 1.7)) +
    plot_annotation(
      title = "Wastewater (citywide) vs. clinical lineage agreement",
      subtitle = paste0(tier_label, "; ", qc_subtitle,
                        "\nbold line: citywide WW. thin translucent lines: the 12 individual neighborhoods.",
                        "\n", nb_caveat,
                        "\nlabels below points: # lineages WW/clinical, monthly bins (citywide)"),
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(paste0(fig_dir, "diversity_panel_", tag, ".jpg"), p,
         width = 10, height = 9.5, dpi = 300)

  cat(sprintf("[%s/%s] months: %d | citywide median Jaccard %.3f / BC similarity %.3f | across neighborhoods, median of medians %.3f / %.3f\n",
              concord_ver, tag, nrow(d$jac),
              median(d$jac$jaccard, na.rm = TRUE), median(d$bc$bc_similarity, na.rm = TRUE),
              median(tapply(d$jac_nb$jaccard, d$jac_nb$LOCATION, median, na.rm = TRUE), na.rm = TRUE),
              median(tapply(d$bc_nb$bc_similarity, d$bc_nb$LOCATION, median, na.rm = TRUE), na.rm = TRUE)))
}

iwalk(tiers, function(label, tag) build_tier(tag, label))

# ---- B. "all" tier with the RA composition row --------------------------------

# The RA row and the metric row use different lineage vocabularies on purpose
# (see ra_monthly_citywide above) -- called out in the subtitle so the two rows
# are not read as the same categories.
ra <- ra_monthly_citywide()

build_ra_combo <- function(metric) {
  d <- read_tier("all")

  bottom <- if (metric == "jaccard") {
    plot_jaccard(d$jac, show_x = TRUE, labels = TRUE, nb = d$jac_nb)
  } else {
    plot_bc(d$bc, show_x = TRUE, labels = TRUE, nb = d$bc_nb)
  }

  p <- plot_clin_volume(d$jac) / plot_ra(ra) / bottom +
    plot_layout(heights = c(1, 2, 1.7)) +
    plot_annotation(
      title = "Wastewater (citywide) vs. clinical lineage agreement",
      subtitle = paste0(
        "bottom: ", if (metric == "jaccard") "Jaccard" else "Bray-Curtis",
        ", parent groups, all present (>= 0.05 per-month rule not applied)",
        "\nmiddle: RAs.R collapse -- retrospective 'ever clears 0.2' retention + manual KP merge,",
        " a different lineage vocabulary from the row below",
        "\n", qc_subtitle,
        "\nbottom, thin translucent lines: the 12 individual neighborhoods (the RA row stays citywide)",
        "\n", nb_caveat,
        "\nlabels below points: # lineages WW/clinical, monthly bins (citywide)"),
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(paste0(fig_dir, "diversity_panel_all_ra_", metric, ".jpg"), p,
         width = 10, height = 10.5, dpi = 300)
}

walk(c("jaccard", "bc"), build_ra_combo)

cat(sprintf("RA row: %d months, %d lineage labels (RAs.R vocabulary)\n",
            n_distinct(ra$month), n_distinct(ra$sublin_collapse)))

# ---- C. one figure per neighborhood, per tier --------------------------------

# Same rows as the citywide figures plus that neighborhood's WW sample count,
# with the citywide series carried through as a thin grey reference so a
# neighborhood can be read against the aggregate it contributes to.

slug <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]+", "_")

build_nb <- function(tag, tier_label) {
  d <- read_tier(tag)

  locs <- sort(unique(d$jac_nb$LOCATION))

  walk(locs, function(loc) {
    jac_l <- d$jac_nb |> filter(LOCATION == loc)
    bc_l  <- d$bc_nb  |> filter(LOCATION == loc)

    p <- plot_clin_volume(jac_l) /
      plot_ww_volume(jac_l) /
      plot_jaccard(jac_l, labels = TRUE, ref = d$jac) /
      plot_bc(bc_l, show_x = TRUE, ref = d$bc) +
      plot_layout(heights = c(0.8, 0.7, 1.7, 1.7)) +
      plot_annotation(
        title = paste0(loc, ": wastewater vs. clinical lineage agreement"),
        subtitle = paste0(
          tier_label, "; ", qc_subtitle,
          "\nbold lines: this neighborhood. thin grey line: citywide WW (all neighborhoods pooled)",
          "\n", nb_caveat,
          "\nlabels below Jaccard points: # lineages WW/clinical, months with no retained sample left blank"),
        theme = theme(plot.title = element_text(face = "bold"))
      )

    ggsave(paste0(nb_fig_dir, "diversity_panel_", tag, "_", slug(loc), ".jpg"), p,
           width = 9, height = 8, dpi = 200)
  })

  cat(sprintf("[%s/%s] wrote %d per-neighborhood panels\n", concord_ver, tag, length(locs)))
}

iwalk(tiers, function(label, tag) build_nb(tag, label))

# ---- D. all-neighborhood small multiples, both metrics -----------------------

# Data without the facetting variable is drawn into every panel, which is how
# the citywide reference appears once per neighborhood without duplicating rows.

build_nb_grid <- function(tag, tier_label) {
  d <- read_tier(tag)

  metric_levels <- c("Jaccard index", "Bray-Curtis similarity")

  to_long <- function(jac, bc) {
    bind_rows(
      jac |> transmute(across(any_of("LOCATION")), month, value = jaccard,
                       metric = metric_levels[1]),
      bc  |> transmute(across(any_of("LOCATION")), month, value = bc_similarity,
                       metric = metric_levels[2])
    ) |>
      mutate(metric = factor(metric, levels = metric_levels))
  }

  nb_long   <- to_long(d$jac_nb, d$bc_nb)
  city_long <- to_long(d$jac, d$bc)

  # strip labels carry the months-observed count, since a neighborhood with 5
  # months and one with 22 are not comparable at a glance
  n_months <- d$jac_nb |>
    group_by(LOCATION) |>
    summarise(n = sum(!is.na(jaccard)), .groups = "drop") |>
    mutate(label = paste0(LOCATION, " (", n, " mo)"))

  nb_long <- nb_long |> left_join(n_months |> select(LOCATION, label), by = "LOCATION")

  p <- ggplot(nb_long, aes(x = month, y = value)) +
    geom_vline(xintercept = clinical_decline_date, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    geom_line(data = city_long, aes(group = metric), color = col_ref, linewidth = 0.35) +
    geom_line(aes(color = metric), linewidth = 0.5) +
    geom_point(aes(color = metric), size = 0.8) +
    facet_wrap(~ label, ncol = 4) +
    scale_color_manual(values = setNames(c(col_jac, col_bc), metric_levels), name = NULL) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
    scale_x_date(date_breaks = "4 months", date_labels = "%b %y") +
    labs(x = "Month", y = "Agreement with citywide clinical",
         title = "Per-neighborhood WW vs. clinical lineage agreement",
         subtitle = paste0(
           tier_label, "; ", qc_subtitle,
           "\ncolored: that neighborhood. grey: citywide WW, repeated in every panel for reference",
           "\n", nb_caveat)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          plot.title = element_text(face = "bold"))

  ggsave(paste0(nb_fig_dir, "diversity_panel_", tag, "_nb_grid.jpg"), p,
         width = 12, height = 8, dpi = 300)
}

iwalk(tiers, function(label, tag) build_nb_grid(tag, label))
