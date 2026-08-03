# Visualize effect of spike-coverage QC filters on lineage composition.

# ---- PATHS ------------------------------------------------------------------
data_dir <- "../data/"
fig_dir  <- "../draft_figures/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(patchwork)

# ---- Params -----------------------------------------------------------------
# breadth = fraction of spike at >= Nx per-site depth; 200x matches demix needs
qc_col              <- "spike_breadth_d200"
qc_thresholds       <- c(0.4, 0.5)
drop_locs           <- "Lower Roxbury"   # dropped in 01_abs_abundance.R (n=20)

# same earn-the-name scheme as RAs.R, recomputed over QC-passing samples only
collapse_threshold  <- 0.2
force_keep_lineages <- c("BA.2.86")
manual_collapse_map <- list("KP.x" = c("KP.1.x", "KP.2.x", "KP.3.x"))
label_min           <- 0.05              # label stacked segments above this RA

# ---- 1. Inputs --------------------------------------------------------------
meta     <- read_rds(paste0(data_dir, "meta_clean.rds"))
pop_wts  <- read_rds(paste0(data_dir, "pop_wts.rds"))
clin_lin <- read_rds(paste0(data_dir, "clin_lin.rds"))
conc_wts <- read_rds(paste0(data_dir, "conc_wts.rds"))

ew_to_date <- function(ew) {
  ymd(paste0(ew %/% 100, "-01-01")) + weeks((ew %% 100) - 1)
}

seq_lw <- meta |>
  select(LOCATION, year_epiweek, all_of(qc_col), breadth_d10, depth_mean) |>
  distinct() |>
  filter(!LOCATION %in% drop_locs)

# ---- 1b. Re-parse Freyja aggregate at full lineage resolution ---------------
parse_freyja_aggregate <- function(path) {
  r <- read.table(path, fill = TRUE, sep = "\t", h = TRUE)
  r <- as.data.frame(sapply(r, function(x) str_replace_all(x, "[',()\\]\\[]", "")))
  r <- as.data.frame(sapply(r, function(x) trimws(gsub("\\s+", " ", x))))
  map_dfr(seq_len(nrow(r)), function(i) {
    lin <- str_split(as.character(r[i, 3]), " ")[[1]]
    ab  <- str_split(as.character(r[i, 4]), " ")[[1]]
    n   <- min(length(lin), length(ab))
    tibble(Sample     = as.character(r[i, 1]),
           sublineage = lin[seq_len(n)],
           abundance  = as.numeric(ab[seq_len(n)]))
  }) |>
    mutate(Sample     = sub("^([0-9]+).*", "\\1", Sample),
           sublineage = str_remove(sublineage, "-like[0-9]*$"))
}

sublin_meta <- bind_rows(
  parse_freyja_aggregate("../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv"),
  parse_freyja_aggregate("../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv")
) |>
  left_join(meta, by = c("Sample" = "FASTQ_ID")) |>
  filter(!is.na(LOCATION), !LOCATION %in% drop_locs)

# ---- 1c. Earn-the-name collapse, computed on QC-passing samples only -------
parent_group <- function(x, extended_groups = force_keep_lineages) {
  default <- str_extract(x, "^[A-Za-z]+\\.[0-9]+")
  default <- ifelse(is.na(default), x, default)   # recombinants keep full token
  group <- default
  for (g in extended_groups) {
    is_g <- !is.na(x) & (x == g | startsWith(x, paste0(g, ".")))
    group[is_g] <- g
  }
  group
}

keep_threshold_time <- function(dat, threshold, force_keep = character(0)) {
  dat |>
    mutate(group = parent_group(sublineage),
           abundance = as.numeric(abundance)) |>
    group_by(Sample, year_epiweek, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    group_by(year_epiweek, group) |>
    summarise(max_ra_any_sample = max(ra, na.rm = TRUE), .groups = "drop") |>
    group_by(group) |>
    # retrospective: named for its whole history if it ever clears threshold
    mutate(keep_by_time = any(max_ra_any_sample >= threshold) | group %in% force_keep) |>
    ungroup() |>
    transmute(year_epiweek, group, keep_by_time, label = paste0(group, ".x"))
}

# clinical unfiltered: QC is a wastewater-sequencing property
clinical <- clin_lin |>
  mutate(parentsub = parent_group(sublineage)) |>
  count(year_epiweek, parentsub, name = "count") |>
  group_by(year_epiweek) |>
  mutate(abundance = count / sum(count)) |>
  ungroup() |>
  transmute(Sample = "clinical_agg", sublineage = parentsub, abundance, year_epiweek)

keep_tbl <- keep_threshold_time(
  bind_rows(
    sublin_meta |>
      filter(.data[[qc_col]] >= qc_thresholds[1]) |>   # QC-passing WW only
      select(Sample, sublineage, abundance, year_epiweek),
    clinical
  ),
  threshold  = collapse_threshold,
  force_keep = force_keep_lineages
)

ww_lin <- sublin_meta |>
  mutate(group = parent_group(sublineage)) |>
  left_join(keep_tbl, by = c("year_epiweek", "group")) |>
  mutate(keep_by_time    = as.logical(coalesce(keep_by_time, FALSE)),
         sublin_collapse = if_else(keep_by_time, label, "other")) |>
  mutate(sublin_collapse = reduce(names(manual_collapse_map), \(acc, dest)
           if_else(acc %in% manual_collapse_map[[dest]], dest, acc),
           .init = sublin_collapse)) |>
  group_by(LOCATION, year_epiweek, lin = sublin_collapse) |>
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

weeks_seq <- sort(unique(seq_lw$year_epiweek))
ww_lin <- ww_lin |> filter(year_epiweek %in% weeks_seq) |>
  mutate(date = ew_to_date(year_epiweek))

# level order anchored to RAs.R's set so colours match RA_plots/ (ordering only)
canon_levels <- read_rds(paste0(data_dir, "ww_collapsed_complete.rds")) |>
  pull(sublin_collapse) |> as.character() |> unique() |> sort()

qc_levels <- sort(unique(as.character(ww_lin$lin)))
extra     <- setdiff(qc_levels, canon_levels)
if (length(extra) > 0) {
  warning("QC-named level(s) absent from the canonical set, appending: ",
          paste(extra, collapse = ", "))
  canon_levels <- c(canon_levels, extra)
}

lin_levels <- c(setdiff(canon_levels, "other"), "other")
ww_lin     <- ww_lin |> mutate(lin = factor(lin, levels = lin_levels))
fill_cols  <- setNames(viridis::viridis(length(lin_levels), option = "turbo"),
                       lin_levels)

cat("named lineage groups (QC-passing, threshold", collapse_threshold, "):",
    length(qc_levels), "of", length(lin_levels), "canonical levels\n")

# ---- 2. QC status per LOCATION-week -----------------------------------------
lab_hi   <- paste0("pass ", qc_thresholds[2])
lab_lo   <- paste0("pass ", qc_thresholds[1], " only")
lab_fail <- paste0("fail ", qc_thresholds[1])

loc_week <- expand_grid(LOCATION = sort(unique(ww_lin$LOCATION)),
                        year_epiweek = weeks_seq) |>
  left_join(seq_lw, by = c("LOCATION", "year_epiweek")) |>
  mutate(
    date = ew_to_date(year_epiweek),
    status = case_when(
      is.na(.data[[qc_col]])          ~ "not sequenced",
      .data[[qc_col]] >= qc_thresholds[2] ~ lab_hi,
      .data[[qc_col]] >= qc_thresholds[1] ~ lab_lo,
      TRUE                                ~ lab_fail
    ),
    status = factor(status, levels = c(lab_hi, lab_lo, lab_fail, "not sequenced"))
  )

status_cols <- setNames(c("#1a7f5a", "#8fce9b", "#D7191C", "grey85"),
                        c(lab_hi, lab_lo, lab_fail, "not sequenced"))

p_status <- loc_week |>
  ggplot(aes(x = date, y = fct_rev(factor(LOCATION)), fill = status)) +
  geom_tile(width = 6, height = 0.85) +
  scale_fill_manual(values = status_cols, name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = NULL,
       title = "Spike-gene coverage QC status by neighborhood-week",
       subtitle = paste0(qc_col, " thresholds; grey = no sequencing attempted")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = "bottom")

p_cov <- seq_lw |>
  mutate(date = ew_to_date(year_epiweek)) |>
  ggplot(aes(x = date, y = .data[[qc_col]])) +
  geom_hline(yintercept = qc_thresholds[1], linetype = "dashed", color = "#8fce9b") +
  geom_hline(yintercept = qc_thresholds[2], linetype = "dotted", color = "#1a7f5a") +
  geom_point(aes(color = .data[[qc_col]] >= qc_thresholds[1]), size = 0.7, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#1a7f5a", `FALSE` = "#D7191C"),
                     guide = "none") +
  facet_wrap(~ LOCATION, ncol = 4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %y") +
  labs(x = NULL, y = qc_col,
       title = "Spike coverage over time by neighborhood",
       subtitle = "dashed = 0.5, dotted = 0.7") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 7))

ggsave(paste0(fig_dir, "qc_status_grid.jpg"), p_status, width = 12, height = 5, dpi = 300)
ggsave(paste0(fig_dir, "qc_coverage_by_nb.jpg"), p_cov, width = 12, height = 8, dpi = 300)

# ---- 3. Neighborhood composition under each filter --------------------------
# dropped weeks drawn as a grey band, never as zeros ("absent" != "not measured")
make_nb_plot <- function(thresh, label) {
  ok <- if (is.null(thresh)) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(.data[[qc_col]] >= thresh) |> distinct(LOCATION, year_epiweek)
  }

  kept    <- ww_lin |> semi_join(ok, by = c("LOCATION", "year_epiweek"))
  dropped <- loc_week |> anti_join(ok, by = c("LOCATION", "year_epiweek"))

  ggplot() +
    geom_col(data = dropped, aes(x = date, y = 1),
             fill = "grey40", width = 6, inherit.aes = FALSE) +
    geom_col(data = kept, aes(x = date, y = abundance, fill = lin), width = 6) +
    geom_text(data = kept,
              aes(x = date, y = abundance,
                  label = ifelse(abundance > label_min, as.character(lin), "")),
              position = position_stack(vjust = 0.5),
              size = 1.4, color = "white", angle = 90) +
    facet_wrap(~ LOCATION, ncol = 4) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b %y") +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = "none") +
    labs(x = NULL, y = "relative abundance",
         title = paste0("Neighborhood lineage composition - ", label),
         subtitle = paste0(nrow(ok), " of ", nrow(seq_lw),
                           " sequenced neighborhood-weeks retained; ",
                           "grey = not measured")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 7))
}

specs <- list(
  list(NULL, "unfiltered", "unfiltered"),
  list(qc_thresholds[1], paste0(qc_col, " >= ", qc_thresholds[1]),
       paste0(qc_col, "_", qc_thresholds[1])),
  list(qc_thresholds[2], paste0(qc_col, " >= ", qc_thresholds[2]),
       paste0(qc_col, "_", qc_thresholds[2])))

# ---- 3b. Per-neighborhood RA plots at the chosen threshold ------------------
# mirrors RA_plots/ but greys out QC-failing weeks; RA_plots/ left untouched
ra_qc_dir <- paste0(fig_dir, "RA_plots_qc/")
dir.create(ra_qc_dir, showWarnings = FALSE, recursive = TRUE)

ok_primary <- seq_lw |>
  filter(.data[[qc_col]] >= qc_thresholds[1]) |>
  distinct(LOCATION, year_epiweek)

# failing weeks keep their composition + labels, washed out -- we still want to
# see what Freyja called even where the call is disregarded
seq_all  <- seq_lw |> distinct(LOCATION, year_epiweek)
ra_kept  <- ww_lin   |> semi_join(ok_primary, by = c("LOCATION", "year_epiweek"))
ra_fail  <- ww_lin   |> anti_join(ok_primary, by = c("LOCATION", "year_epiweek"))
ra_unseq <- loc_week |> anti_join(seq_all,    by = c("LOCATION", "year_epiweek"))

conc_label <- conc_wts |>
  filter(!LOCATION %in% drop_locs, year_epiweek %in% weeks_seq) |>
  mutate(date = ew_to_date(year_epiweek),
         conc_quartile = cut(logeff,
                             breaks = quantile(logeff, 0:4 / 4, na.rm = TRUE),
                             labels = c("Q1", "Q2", "Q3", "Q4"),
                             include.lowest = TRUE))

conc_cols <- c(Q1 = "#2C7BB6", Q2 = "#ABD9E9", Q3 = "#FDAE61", Q4 = "#D7191C")

for (loc in sort(unique(ww_lin$LOCATION))) {
  k      <- ra_kept    |> filter(LOCATION == loc)
  f      <- ra_fail    |> filter(LOCATION == loc)
  u      <- ra_unseq   |> filter(LOCATION == loc)
  cq     <- conc_label |> filter(LOCATION == loc)
  n_pass <- sum(ok_primary$LOCATION == loc)
  n_seq  <- sum(seq_lw$LOCATION == loc)

  p_loc <- ggplot() +
    geom_col(data = u, aes(x = date, y = 1),
             fill = "grey85", width = 6, inherit.aes = FALSE) +
    geom_col(data = f, aes(x = date, y = abundance, fill = lin),
             width = 6, alpha = 0.25) +
    geom_text(data = f,
              aes(x = date, y = abundance,
                  label = ifelse(abundance > label_min, as.character(lin), "")),
              position = position_stack(vjust = 0.5),
              size = 2, color = "grey25", angle = 90) +
    geom_col(data = k, aes(x = date, y = abundance, fill = lin), width = 6) +
    geom_text(data = k,
              aes(x = date, y = abundance,
                  label = ifelse(abundance > label_min, as.character(lin), "")),
              position = position_stack(vjust = 0.5),
              size = 2, color = "white", angle = 90) +
    geom_point(data = cq, aes(x = date, y = 1.05, color = conc_quartile),
               size = 1.4, inherit.aes = FALSE) +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    scale_color_manual(values = conc_cols, name = "PCR conc\nquartile",
                       na.translate = FALSE) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    coord_cartesian(ylim = c(0, 1.08)) +
    guides(fill = "none") +
    labs(x = NULL, y = "relative abundance",
         title = paste0(loc, " - lineage composition, QC-filtered"),
         subtitle = paste0(qc_col, " >= ", qc_thresholds[1], "  |  ",
                           n_pass, " of ", n_seq,
                           " sequenced weeks pass; faded = QC fail, ",
                           "light grey = not sequenced")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "right")

  ggsave(paste0(ra_qc_dir, "RAxNB_", gsub(" ", "_", loc), ".jpg"),
         p_loc, width = 14, height = 6, dpi = 300)
}
cat("Wrote", n_distinct(ww_lin$LOCATION), "per-neighborhood QC plots to", ra_qc_dir, "\n")

for (s in specs) {
  p <- make_nb_plot(s[[1]], s[[2]])
  ggsave(paste0(fig_dir, "qc_nb_composition_", s[[3]], ".jpg"),
         p, width = 13, height = 8.5, dpi = 300)
}

# ---- 4. Population-weighted citywide under each filter ----------------------
# weights renormalised over the neighborhoods passing QC each week
city_pop_total <- sum(pop_wts$pop[!pop_wts$LOCATION %in% drop_locs])

citywide_weighted <- function(thresh, label) {
  ok <- if (is.null(thresh)) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(.data[[qc_col]] >= thresh) |> distinct(LOCATION, year_epiweek)
  }

  dat <- ww_lin |> semi_join(ok, by = c("LOCATION", "year_epiweek"))

  # explicit zeros: else sum(pop) covers only neighborhoods that detected it
  comp <- dat |>
    group_by(year_epiweek, date) |>
    complete(LOCATION, lin, fill = list(abundance = 0)) |>
    ungroup() |>
    left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION") |>
    group_by(year_epiweek, date, lin) |>
    summarise(abundance = sum(abundance * pop) / sum(pop), .groups = "drop")

  support <- dat |>
    distinct(year_epiweek, date, LOCATION) |>
    left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION") |>
    group_by(year_epiweek, date) |>
    summarise(n_loc = n(), pop_frac = sum(pop) / city_pop_total, .groups = "drop")

  comp |>
    left_join(support, by = c("year_epiweek", "date")) |>
    mutate(filter = label)
}

city_labs <- c("unfiltered",
               paste0(">= ", qc_thresholds[1]),
               paste0(">= ", qc_thresholds[2]))
city <- bind_rows(
  citywide_weighted(NULL,             city_labs[1]),
  citywide_weighted(qc_thresholds[1], city_labs[2]),
  citywide_weighted(qc_thresholds[2], city_labs[3])
) |>
  mutate(filter = factor(filter, levels = city_labs))

p_city <- city |>
  ggplot(aes(x = date, y = abundance, fill = lin)) +
  geom_col(width = 6) +
  facet_wrap(~ filter, ncol = 1) +
  geom_text(aes(label = ifelse(abundance > label_min, as.character(lin), "")),
            position = position_stack(vjust = 0.5),
            size = 1.6, color = "white", angle = 90) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  guides(fill = "none") +
  labs(x = NULL, y = "population-weighted relative abundance",
       title = paste0("Citywide population-weighted lineage composition - ", qc_col),
       subtitle = "weights renormalised each week over neighborhoods passing QC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# how much of the city population the estimate rests on
p_support <- city |>
  distinct(date, filter, n_loc, pop_frac) |>
  ggplot(aes(x = date, y = pop_frac, color = filter)) +
  geom_step() +
  scale_color_manual(values = setNames(c("grey40", "#8fce9b", "#1a7f5a"), city_labs),
                     name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "fraction of city pop.\ncontributing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

ggsave(paste0(fig_dir, "qc_citywide_weighted.jpg"),
       p_city / p_support + plot_layout(heights = c(3, 1)),
       width = 12, height = 11, dpi = 300)

# ---- 5. Console summary -----------------------------------------------------
cat("\n--- neighborhood-weeks retained ---\n")
for (s in specs) {
  ok <- if (is.null(s[[1]])) seq_lw else filter(seq_lw, .data[[qc_col]] >= s[[1]])
  cat(sprintf("%-24s %3d / %3d (%2.0f%%) | locations with >=20 weeks: %2d\n",
              s[[2]], nrow(ok), nrow(seq_lw), 100 * nrow(ok) / nrow(seq_lw),
              sum(count(ok, LOCATION)$n >= 20)))
}

cat("\n--- citywide population support (fraction of city pop) ---\n")
print(city |> distinct(date, filter, pop_frac) |>
        group_by(filter) |>
        summarise(median_pop_frac = round(median(pop_frac), 2),
                  min_pop_frac    = round(min(pop_frac), 2),
                  weeks           = n(), .groups = "drop"))

cat("\n--- named lineage groups detected under each filter ---\n")
for (s in specs) {
  ok <- if (is.null(s[[1]])) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(.data[[qc_col]] >= s[[1]]) |> distinct(LOCATION, year_epiweek)
  }
  n_grp <- ww_lin |>
    semi_join(ok, by = c("LOCATION", "year_epiweek")) |>
    filter(abundance > 0) |>
    pull(lin) |> as.character() |> unique() |> length()
  cat(sprintf("%-30s %3d groups detected\n", s[[2]], n_grp))
}

cat("\nFigures written to", fig_dir, "\n")
