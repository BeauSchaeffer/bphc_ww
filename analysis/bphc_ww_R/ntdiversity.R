# Per-neighborhood diversity panels over time, after Hill et al. figs. S4-S6.

storage_dir <- "../data/"
div_subdir  <- "AF03DF200"   # threshold set written by 04.2_diversity.sh
fig_dir     <- "../draft_figures/nt_diversity/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(tidyverse)
library(lubridate)
library(patchwork)

source("lineage_helpers.R")

# coverage QC, matching jaccard_analysis.R / RAs.R / the gravity model
qc_breadth_col <- "spike_breadth_d200"
qc_breadth_min <- 0.4
apply_qc <- TRUE

# same WW noise floor as the concordance metrics; the variant-count row uses the
# "fine" tier vocabulary: >= 1% abundance, renormalized, no parent collapse
ww_detect_thresh <- 0.01

# bottom composition row, mirroring draft_figures/RA_plots_qc/
include_ra_panel   <- TRUE
collapse_threshold <- 0.2
manual_collapse_map <- list("KP.x" = c("KP.1.x", "KP.2.x", "KP.3.x"))
label_min <- 0.05

col_fixed    <- "#2C7BB6"
col_callable <- "#D95F02"
col_count    <- "#4D4D4D"

ew_to_date <- function(ew) ymd(paste0(ew %/% 100, "-01-01")) + weeks((ew %% 100) - 1)

# ARTIC primer swap: batch 1 (V4.1) ends 2023-05-01, batch 2 (V5.3.2) begins
# 2023-05-08 -- no overlap, so the line sits in the gap. Batch is confounded with
# time (§12), so this boundary is where a step could be primer or epidemiology.
primer_swap_date <- as.Date("2023-05-05")
primer_vline <- function() {
  geom_vline(xintercept = primer_swap_date, linetype = "dotted",
             color = "grey25", linewidth = 0.3)
}

# ---- 1. diversity tables ------------------------------------------------------

read_diversity <- function(batch) {
  dir <- paste0("../baseload_batch", batch, "_outputs/diversity/", div_subdir, "/")
  files <- list.files(dir, pattern = "[.]diversity[.]tsv$", full.names = TRUE)
  stopifnot(length(files) > 0)
  map_dfr(files, read_tsv, show_col_types = FALSE, progress = FALSE) |>
    mutate(batch = batch)
}

div <- bind_rows(read_diversity(1), read_diversity(2)) |>
  mutate(Sample = as.character(sample))

meta <- read_rds(paste0(storage_dir, "meta_clean.rds")) |>
  mutate(FASTQ_ID = as.character(FASTQ_ID))

div <- div |> inner_join(meta, by = c("Sample" = "FASTQ_ID"))

if (apply_qc) {
  n_before <- n_distinct(div$Sample)
  div <- div |> filter(.data[[qc_breadth_col]] >= qc_breadth_min)
  cat(sprintf("QC: %s >= %s | %d of %d samples retained\n",
              qc_breadth_col, qc_breadth_min, n_distinct(div$Sample), n_before))
}

# ---- 2. Freyja variant count, collapsed to parent lineage ---------------------

sublin_meta <- bind_rows(
  parse_freyja_aggregate("../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv"),
  parse_freyja_aggregate("../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv")
) |>
  mutate(abundance = as.numeric(abundance)) |>
  left_join(meta, by = c("Sample" = "FASTQ_ID")) |>
  filter(!is.na(LOCATION))

# one row per sample: how many lineages clear the detection floor, at full
# Freyja resolution (no parent collapse)
variant_counts <- sublin_meta |>
  apply_ww_detection(ww_detect_thresh) |>
  group_by(Sample) |>
  summarise(n_groups = n_distinct(sublineage), .groups = "drop")

div <- div |> left_join(variant_counts, by = "Sample")

cat(sprintf("samples plotted: %d | missing variant count: %d\n",
            nrow(div), sum(is.na(div$n_groups))))

# ---- 3. composition row, mirroring qc_filter_viz.R's RA_plots_qc --------------
# Lower Roxbury is NOT dropped here (qc_filter_viz.R drops it) so every
# neighborhood with diversity panels also gets a composition row.

if (include_ra_panel) {
  clinical <- read_rds(paste0(storage_dir, "clin_lin.rds")) |>
    mutate(parentsub = parent_group(sublineage)) |>
    count(year_epiweek, parentsub, name = "count") |>
    group_by(year_epiweek) |>
    mutate(abundance = count / sum(count)) |>
    ungroup() |>
    transmute(Sample = "clinical_agg", sublineage = parentsub, abundance, year_epiweek)

  keep_tbl <- keep_threshold_time(
    bind_rows(
      sublin_meta |>
        filter(.data[[qc_breadth_col]] >= qc_breadth_min) |>
        select(Sample, sublineage, abundance, year_epiweek),
      clinical
    ),
    threshold  = collapse_threshold,
    force_keep = force_keep_lineages
  ) |>
    rename(year_epiweek = time)   # shared helper names its time column `time`

  seq_lw <- meta |>
    select(LOCATION, year_epiweek, all_of(qc_breadth_col)) |>
    distinct()
  weeks_seq <- sort(unique(seq_lw$year_epiweek))

  ww_lin <- sublin_meta |>
    mutate(group = parent_group(sublineage)) |>
    left_join(keep_tbl, by = c("year_epiweek", "group")) |>
    mutate(keep_by_time    = as.logical(coalesce(keep_by_time, FALSE)),
           sublin_collapse = if_else(keep_by_time, label, "other")) |>
    mutate(sublin_collapse = reduce(names(manual_collapse_map), \(acc, dest)
             if_else(acc %in% manual_collapse_map[[dest]], dest, acc),
             .init = sublin_collapse)) |>
    group_by(LOCATION, year_epiweek, lin = sublin_collapse) |>
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    filter(year_epiweek %in% weeks_seq) |>
    mutate(date = ew_to_date(year_epiweek))

  # colour order anchored to RAs.R's set so fills match RA_plots/
  canon_levels <- read_rds(paste0(storage_dir, "ww_collapsed_complete.rds")) |>
    pull(sublin_collapse) |> as.character() |> unique() |> sort()
  canon_levels <- union(canon_levels, sort(unique(as.character(ww_lin$lin))))
  lin_levels <- c(setdiff(canon_levels, "other"), "other")
  ww_lin    <- ww_lin |> mutate(lin = factor(lin, levels = lin_levels))
  fill_cols <- setNames(viridis::viridis(length(lin_levels), option = "turbo"),
                        lin_levels)

  ok_primary <- seq_lw |>
    filter(.data[[qc_breadth_col]] >= qc_breadth_min) |>
    distinct(LOCATION, year_epiweek)
  seq_all  <- seq_lw |> distinct(LOCATION, year_epiweek)
  loc_week <- expand_grid(LOCATION = sort(unique(ww_lin$LOCATION)),
                          year_epiweek = weeks_seq) |>
    mutate(date = ew_to_date(year_epiweek))

  ra_kept  <- ww_lin   |> semi_join(ok_primary, by = c("LOCATION", "year_epiweek"))
  ra_fail  <- ww_lin   |> anti_join(ok_primary, by = c("LOCATION", "year_epiweek"))
  ra_unseq <- loc_week |> anti_join(seq_all,    by = c("LOCATION", "year_epiweek"))
}

# ---- 4. panel builders --------------------------------------------------------

# shared x range and per-metric y ranges so the 12 neighborhoods are comparable
# spans both grids: sample dates and the epiweek-start dates the RA row uses
date_lims <- range(c(div$SAMPLING_DATE,
                     if (include_ra_panel) ra_unseq$date else NULL,
                     if (include_ra_panel) ww_lin$date   else NULL), na.rm = TRUE)
# padded by half a bar width: geom_col(width = 6) puts xmin/xmax 3 days either
# side of the week, and the scale drops bars whose extent falls outside limits --
# unpadded, the first and last weeks' bars vanish silently
x_scale <- scale_x_date(limits = date_lims + c(-4, 4), date_breaks = "1 month",
                        date_labels = "%Y-%m")

lim_of <- function(cols) {
  v <- unlist(div[cols], use.names = FALSE)
  range(v[is.finite(v)])
}

base_theme <- theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold"))

x_blank <- theme(axis.text.x = element_blank())
x_text  <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# both denominators on one row -- dropping either would silently pick a side of
# the fixed-vs-callable question these outputs exist to answer
plot_metric <- function(d, fixed_col, callable_col, ylab, ylim, show_x) {
  d |>
    select(SAMPLING_DATE, fixed = all_of(fixed_col), callable = all_of(callable_col)) |>
    pivot_longer(c(fixed, callable), names_to = "denom", values_to = "value") |>
    filter(is.finite(value)) |>
    arrange(SAMPLING_DATE) |>
    ggplot(aes(SAMPLING_DATE, value, color = denom)) +
    primer_vline() +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    geom_point(size = 1.1, alpha = 0.85) +
    scale_color_manual(values = c(fixed = col_fixed, callable = col_callable)) +
    scale_y_continuous(limits = ylim) +
    x_scale +
    labs(x = NULL, y = ylab) +
    base_theme +
    (if (show_x) x_text else x_blank)
}

# y starts at 0 so a low count sits clear of the axis rather than reading as zero
plot_count <- function(d, ymax, show_x) {
  d |>
    filter(!is.na(n_groups)) |>
    arrange(SAMPLING_DATE) |>
    ggplot(aes(SAMPLING_DATE, n_groups)) +
    primer_vline() +
    geom_line(linewidth = 0.25, alpha = 0.7, color = col_count) +
    geom_point(size = 1.1, alpha = 0.85, color = col_count) +
    scale_y_continuous(limits = c(0, ymax), breaks = scales::breaks_width(5)) +
    x_scale +
    labs(x = NULL, y = "Freyja lineages\n(fine resolution)") +
    base_theme +
    (if (show_x) x_text else x_blank)
}

plot_ra <- function(loc) {
  k <- ra_kept  |> filter(LOCATION == loc)
  f <- ra_fail  |> filter(LOCATION == loc)
  u <- ra_unseq |> filter(LOCATION == loc)

  ggplot() +
    geom_col(data = u, aes(x = date, y = 1), fill = "grey85", width = 6) +
    geom_col(data = f, aes(x = date, y = abundance, fill = lin), width = 6, alpha = 0.25) +
    geom_col(data = k, aes(x = date, y = abundance, fill = lin), width = 6) +
    geom_text(data = k,
              aes(x = date, y = abundance,
                  label = ifelse(abundance > label_min, as.character(lin), "")),
              position = position_stack(vjust = 0.5),
              size = 1.5, color = "white", angle = 90) +
    primer_vline() +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    x_scale +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = "none") +
    labs(x = "Sampling date", y = "relative\nabundance") +
    base_theme +
    x_text
}

# ---- 5. one stacked figure per neighborhood x region --------------------------

regions <- list(
  genome = list(pi = c("pi_ww_fixed", "pi_ww_callable"),
                h  = c("h_ww_fixed",  "h_ww_callable"),
                label = "genome-wide"),
  spike  = list(pi = c("pi_spike_fixed", "pi_spike_callable"),
                h  = c("h_spike_fixed",  "h_spike_callable"),
                label = "spike")
)

count_max <- max(div$n_groups, na.rm = TRUE)
locations <- sort(unique(div$LOCATION))

qc_note <- if (apply_qc) {
  paste0("coverage QC applied (", qc_breadth_col, " >= ", qc_breadth_min, ")")
} else {
  "no coverage QC"
}

build_figure <- function(loc, region_tag) {
  r <- regions[[region_tag]]
  d <- div |> filter(LOCATION == loc)
  if (nrow(d) == 0) return(invisible(NULL))

  rows <- list(
    plot_metric(d, r$pi[1], r$pi[2], expression(pi[ww]), lim_of(r$pi), FALSE),
    plot_metric(d, r$h[1], r$h[2], expression(H[ww]), lim_of(r$h), FALSE),
    plot_count(d, count_max, show_x = !include_ra_panel)
  )
  if (include_ra_panel) rows <- c(rows, list(plot_ra(loc)))

  sub <- paste0("blue = fixed denominator (Hill et al.), orange = callable; ",
                div_subdir, "; ", qc_note,
                "\nlineage row: Freyja calls >= ", ww_detect_thresh,
                " abundance, full lineage resolution")
  if (include_ra_panel) {
    sub <- paste0(sub,
                  "\nbottom: RA composition as in RA_plots_qc/ -- faded = QC fail, ",
                  "light grey = not sequenced",
                  "\ndotted vertical line: ARTIC primer swap V4.1 -> V5.3.2 (2023-05-05)")
  }

  p <- wrap_plots(rows, ncol = 1,
                  heights = if (include_ra_panel) c(1, 1, 1, 1.2) else c(1, 1, 1)) +
    plot_annotation(
      title = paste0(loc, " -- ", r$label, " diversity"),
      subtitle = sub,
      theme = theme(plot.title = element_text(face = "bold"),
                    plot.subtitle = element_text(size = 8, color = "grey30"))
    )

  ggsave(paste0(fig_dir, "ntdiv_", region_tag, "_", gsub(" ", "_", loc), ".png"),
         p, width = 11, height = if (include_ra_panel) 10 else 7.5, dpi = 300)
}

walk(locations, function(loc) walk(names(regions), function(r) build_figure(loc, r)))

cat(sprintf("wrote %d figures to %s\n", 2 * length(locations), fig_dir))

# ---- 6. what each panel actually draws ----------------------------------------
# The top three rows are QC-filtered; the RA row deliberately is not (it fades
# failures instead of dropping them). This table is the audit of that gap.

div_all <- bind_rows(read_diversity(1), read_diversity(2)) |>
  mutate(Sample = as.character(sample)) |>
  inner_join(meta, by = c("Sample" = "FASTQ_ID")) |>
  left_join(variant_counts, by = "Sample") |>
  mutate(qc_pass = .data[[qc_breadth_col]] >= qc_breadth_min)

# full LOCATION x week grid, so weeks with no sequencing at all are counted
week_grid <- expand_grid(LOCATION = sort(unique(meta$LOCATION)),
                         year_epiweek = sort(unique(meta$year_epiweek)))

week_status <- week_grid |>
  left_join(div_all |> group_by(LOCATION, year_epiweek) |>
              summarise(any_pass = any(qc_pass), .groups = "drop"),
            by = c("LOCATION", "year_epiweek")) |>
  mutate(status = factor(case_when(is.na(any_pass) ~ "unsequenced",
                                   any_pass        ~ "qc_pass",
                                   TRUE            ~ "qc_fail"),
                         levels = c("unsequenced", "qc_fail", "qc_pass"))) |>
  count(LOCATION, status, .drop = FALSE) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0)

sample_status <- div_all |>
  group_by(LOCATION) |>
  summarise(
    n_seq          = n(),
    n_qc_pass      = sum(qc_pass),
    pi_ww_fix      = sum(qc_pass & is.finite(pi_ww_fixed)),
    pi_ww_call     = sum(qc_pass & is.finite(pi_ww_callable)),
    pi_spike_fix   = sum(qc_pass & is.finite(pi_spike_fixed)),
    pi_spike_call  = sum(qc_pass & is.finite(pi_spike_callable)),
    freyja_count   = sum(qc_pass & !is.na(n_groups)),
    .groups = "drop")

panel_qc <- week_status |>
  full_join(sample_status, by = "LOCATION") |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  relocate(LOCATION, unsequenced, qc_fail, qc_pass)

panel_qc <- bind_rows(panel_qc,
                      panel_qc |> summarise(across(where(is.numeric), sum)) |>
                        mutate(LOCATION = "TOTAL"))

write_csv(panel_qc, paste0(fig_dir, "panel_qc_summary.csv"))
cat("\n--- rows drawn per panel (weeks | samples) ---\n")
print(as.data.frame(panel_qc), row.names = FALSE)
cat("\nweeks: unsequenced/qc_fail/qc_pass -- RA row draws all three,",
    "\n       the pi/H/lineage rows draw qc_pass only\n")

# ---- 7. detection footing per panel -------------------------------------------
# Floors read from the data where they are recorded (the diversity TSV carries
# its own af_floor/depth_floor); the two Freyja values come from
# informatics/scripts/05_freyja_demix.sh and are hardcoded here.

demix_eps        <- 1e-06   # 05_freyja_demix.sh --eps
demix_depthcut   <- 10      # 05_freyja_demix.sh --depthcutoff
af_floor_used    <- unique(div$af_floor)
depth_floor_used <- unique(div$depth_floor)
stopifnot(length(af_floor_used) == 1, length(depth_floor_used) == 1)

footing <- tibble::tribble(
  ~panel,                  ~signal,              ~abundance_floor,             ~depth_floor,                        ~sample_qc,       ~grouping,
  "1-2 pi/H  fixed",       "allele freq",        paste0("AF >= ", af_floor_used), "none (any depth counts)",        "breadth >= 0.4", "-",
  "1-2 pi/H  callable",    "allele freq",        paste0("AF >= ", af_floor_used), paste0(depth_floor_used, "x per site"), "breadth >= 0.4", "-",
  "3 lineage count",       "Freyja lineage RA",  paste0("RA >= ", ww_detect_thresh, ", renorm"), paste0(demix_depthcut, "x (demix)"), "breadth >= 0.4", "none (fine)",
  "4 RA composition",      "Freyja lineage RA",  paste0("RA >= ", demix_eps, " (demix eps)"),    paste0(demix_depthcut, "x (demix)"), "none (faded)",   paste0("keep>=", collapse_threshold, " else 'other'")
)

write_csv(footing, paste0(fig_dir, "panel_detection_footing.csv"))
cat("\n--- detection footing per panel ---\n")
print(as.data.frame(footing), row.names = FALSE)
cat("\nNOT equal footing. Three depth floors in play (none / ", demix_depthcut,
    "x / ", depth_floor_used, "x) and abundance floors spanning ",
    demix_eps, " to ", af_floor_used, ".\n",
    "freyja variants (04_freyja_variants.sh) is run at defaults -- no AF or depth\n",
    "floor -- so every floor above is applied downstream of an unfiltered call set.\n",
    sep = "")

# ---- 8. citywide population-weighted figure -----------------------------------
# Weighting decisions, all confirmed before writing:
#   pi/H        pop-weighted, renormalized over neighborhoods observed that week
#               (Hill Eq. 4). A full-city denominator would multiply the mean by
#               the coverage fraction and reproduce the fixed-denominator
#               pathology this section criticizes. Every week is drawn; the
#               support row carries coverage.
#   lineages    citywide union AND pop-weighted mean of per-sample counts, both
#               plotted -- their divergence is the sampling-effort effect
#   RA          pop-weighted, renormalized, QC-passing samples only
#   support     share of catchment population with >= 1 QC-passing sample
#   no smoothing; raw weekly, matching the per-neighborhood panels

city_drop_locs <- "Lower Roxbury"   # 20 samples, all after 2024-03-11; including
                                    # it steps the denominator 5.5% mid-series

pop_wts <- read_rds(paste0(storage_dir, "pop_wts.rds")) |>
  select(LOCATION, pop) |>
  filter(!LOCATION %in% city_drop_locs)
pop_city <- sum(pop_wts$pop)

div_city <- div |> filter(!LOCATION %in% city_drop_locs)

# neighborhood-week means first (Hill average replicate samples within a site),
# then population-weight across neighborhoods
nb_week <- div_city |>
  group_by(LOCATION, year_epiweek) |>
  summarise(across(c(starts_with("pi_"), starts_with("h_"), n_groups),
                   ~mean(.x, na.rm = TRUE)), .groups = "drop") |>
  inner_join(pop_wts, by = "LOCATION")

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

city <- nb_week |>
  group_by(year_epiweek) |>
  summarise(across(c(starts_with("pi_"), starts_with("h_")), ~wmean(.x, pop)),
            n_lin_wmean = wmean(n_groups, pop),
            pop_share   = sum(pop) / pop_city,
            n_loc       = n_distinct(LOCATION),
            .groups = "drop")

# citywide union: distinct lineages seen anywhere that week, QC-passing samples
lin_union <- sublin_meta |>
  filter(.data[[qc_breadth_col]] >= qc_breadth_min,
         !LOCATION %in% city_drop_locs) |>
  apply_ww_detection(ww_detect_thresh) |>
  group_by(year_epiweek) |>
  summarise(n_lin_union = n_distinct(sublineage), .groups = "drop")

city <- city |>
  left_join(lin_union, by = "year_epiweek") |>
  mutate(date = ew_to_date(year_epiweek))

# citywide RA: mean composition per neighborhood-week, pop-weighted, renormalized
ra_city <- sublin_meta |>
  filter(.data[[qc_breadth_col]] >= qc_breadth_min,
         !LOCATION %in% city_drop_locs) |>
  mutate(group = parent_group(sublineage)) |>
  left_join(keep_tbl, by = c("year_epiweek", "group")) |>
  mutate(keep_by_time = as.logical(coalesce(keep_by_time, FALSE)),
         lin = if_else(keep_by_time, label, "other")) |>
  mutate(lin = reduce(names(manual_collapse_map), \(acc, dest)
           if_else(acc %in% manual_collapse_map[[dest]], dest, acc), .init = lin)) |>
  group_by(LOCATION, year_epiweek, lin) |>
  summarise(a = sum(abundance, na.rm = TRUE), .groups = "drop") |>
  left_join(div_city |> distinct(LOCATION, year_epiweek, Sample) |>
              count(LOCATION, year_epiweek, name = "n_samp"),
            by = c("LOCATION", "year_epiweek")) |>
  filter(!is.na(n_samp)) |>
  mutate(p = a / n_samp) |>                       # mean RA per neighborhood-week
  inner_join(pop_wts, by = "LOCATION") |>
  group_by(year_epiweek, lin) |>
  summarise(num = sum(pop * p), .groups = "drop") |>
  group_by(year_epiweek) |>
  mutate(prop = num / sum(num)) |>                # renormalize over observed
  ungroup() |>
  mutate(date = ew_to_date(year_epiweek),
         lin  = factor(lin, levels = lin_levels))

cat(sprintf("\ncitywide: %d weeks | pop share median %.2f (min %.2f) | %d neighborhoods\n",
            nrow(city), median(city$pop_share), min(city$pop_share), nrow(pop_wts)))

plot_city_metric <- function(fixed_col, callable_col, ylab) {
  city |>
    select(date, fixed = all_of(fixed_col), callable = all_of(callable_col)) |>
    pivot_longer(c(fixed, callable), names_to = "denom", values_to = "value") |>
    filter(is.finite(value)) |>
    ggplot(aes(date, value, color = denom)) +
    primer_vline() +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    geom_point(size = 1.1, alpha = 0.85) +
    scale_color_manual(values = c(fixed = col_fixed, callable = col_callable)) +
    x_scale + labs(x = NULL, y = ylab) + base_theme + x_blank
}

plot_city_lineages <- function() {
  city |>
    select(date, union = n_lin_union, wmean = n_lin_wmean) |>
    pivot_longer(c(union, wmean), names_to = "agg", values_to = "value") |>
    filter(is.finite(value)) |>
    ggplot(aes(date, value, color = agg)) +
    primer_vline() +
    geom_line(linewidth = 0.25, alpha = 0.7) +
    geom_point(size = 1.1, alpha = 0.85) +
    scale_color_manual(values = c(union = "#4D4D4D", wmean = "#9ECAE1")) +
    scale_y_continuous(limits = c(0, NA)) +
    x_scale +
    labs(x = NULL, y = "Freyja lineages\n(dark: union, light: wtd mean)") +
    base_theme + x_blank
}

plot_city_ra <- function() {
  ggplot(ra_city, aes(date, prop, fill = lin)) +
    geom_col(width = 6) +
    primer_vline() +
    geom_text(aes(label = ifelse(prop > label_min, as.character(lin), "")),
              position = position_stack(vjust = 0.5),
              size = 1.5, color = "white", angle = 90) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    x_scale + coord_cartesian(ylim = c(0, 1)) + guides(fill = "none") +
    labs(x = NULL, y = "relative\nabundance") + base_theme + x_blank
}

plot_city_support <- function() {
  city |>
    ggplot(aes(date, pop_share)) +
    geom_col(width = 6, fill = "grey55") +
    primer_vline() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30", linewidth = 0.3) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    x_scale +
    labs(x = "Sampling week", y = "pop. share\npassing QC") +
    base_theme + x_text
}

build_city_figure <- function(region_tag) {
  r <- regions[[region_tag]]
  p <- plot_city_metric(r$pi[1], r$pi[2], expression(pi[ww])) /
    plot_city_metric(r$h[1], r$h[2], expression(H[ww])) /
    plot_city_lineages() / plot_city_ra() / plot_city_support() +
    plot_layout(heights = c(1, 1, 1, 1.2, 0.5)) +
    plot_annotation(
      title = paste0("Citywide (population-weighted) -- ", r$label, " diversity"),
      subtitle = paste0(
        "blue = fixed denominator (Hill et al.), orange = callable; ", div_subdir, "; ", qc_note,
        "\npi/H and RA are population-weighted over neighborhoods observed each week ",
        "(Hill Eq. 4), renormalized; raw weekly, no smoothing",
        "\nbottom: share of catchment population with >= 1 QC-passing sample -- ",
        "dashed line at 0.5. ", city_drop_locs, " excluded (see script)",
        "\ndotted vertical line: ARTIC primer swap V4.1 -> V5.3.2 (2023-05-05)"),
      theme = theme(plot.title = element_text(face = "bold"),
                    plot.subtitle = element_text(size = 8, color = "grey30")))

  ggsave(paste0(fig_dir, "ntdiv_citywide_", region_tag, ".png"),
         p, width = 11, height = 11, dpi = 300)
}

walk(names(regions), build_city_figure)
cat("wrote 2 citywide figures\n")
