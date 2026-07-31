# Visualize effect of spike-coverage QC filters on lineage composition.

# ---- PATHS ------------------------------------------------------------------
data_dir <- "../data/"
fig_dir  <- "../draft_figures/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(patchwork)

# ---- Params -----------------------------------------------------------------
qc_thresholds  <- c(0.5, 0.7)          # spike_breadth_d10 cutoffs to compare
lineages_model <- c("BA.5.x", "BQ.1.x", "XBB.1.x", "EG.5.x",
                    "HV.1.x", "JN.1.x", "KP.x")   # matches 01_abs_abundance.R
drop_locs      <- "Lower Roxbury"       # dropped in 01_abs_abundance.R (n=20)

# ---- 1. Inputs --------------------------------------------------------------
ww      <- read_rds(paste0(data_dir, "ww_collapsed_complete.rds"))
meta    <- read_rds(paste0(data_dir, "meta_clean.rds"))
pop_wts <- read_rds(paste0(data_dir, "pop_wts.rds"))

ew_to_date <- function(ew) {
  ymd(paste0(ew %/% 100, "-01-01")) + weeks((ew %% 100) - 1)
}

# one sample per LOCATION-week, so distinct() is lossless
seq_lw <- meta |>
  distinct(LOCATION, year_epiweek, spike_breadth_d10, breadth_d10, depth_mean) |>
  filter(!LOCATION %in% drop_locs)

# collapse to the modelled lineage set; everything else -> "other"
lin_levels <- c(lineages_model, "other")
ww_lin <- ww |>
  filter(!LOCATION %in% drop_locs) |>
  mutate(lin = if_else(as.character(sublin_collapse) %in% lineages_model,
                       as.character(sublin_collapse), "other")) |>
  group_by(LOCATION, year_epiweek, lin) |>
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") |>
  mutate(lin = factor(lin, levels = lin_levels))

# restrict to weeks that were actually sequenced somewhere (ww is a full grid)
weeks_seq <- sort(unique(seq_lw$year_epiweek))
ww_lin <- ww_lin |> filter(year_epiweek %in% weeks_seq) |>
  mutate(date = ew_to_date(year_epiweek))

# "not measured" is a fill level in its own right so it can never be confused
# with the "other" lineage bucket
plot_levels <- c(lineages_model, "other", "not measured")
fill_cols <- setNames(
  c(viridis::viridis(length(lineages_model), option = "turbo"),
    "grey55", "#f7dcd9"),
  plot_levels
)

# ---- 2. QC status per LOCATION-week -----------------------------------------
loc_week <- expand_grid(LOCATION = sort(unique(ww_lin$LOCATION)),
                        year_epiweek = weeks_seq) |>
  left_join(seq_lw, by = c("LOCATION", "year_epiweek")) |>
  mutate(
    date = ew_to_date(year_epiweek),
    status = case_when(
      is.na(spike_breadth_d10)   ~ "not sequenced",
      spike_breadth_d10 >= 0.7   ~ "pass 0.7",
      spike_breadth_d10 >= 0.5   ~ "pass 0.5 only",
      TRUE                        ~ "fail 0.5"
    ),
    status = factor(status, levels = c("pass 0.7", "pass 0.5 only",
                                       "fail 0.5", "not sequenced"))
  )

status_cols <- c("pass 0.7" = "#1a7f5a", "pass 0.5 only" = "#8fce9b",
                 "fail 0.5" = "#D7191C", "not sequenced" = "grey85")

p_status <- loc_week |>
  ggplot(aes(x = date, y = fct_rev(factor(LOCATION)), fill = status)) +
  geom_tile(width = 6, height = 0.85) +
  scale_fill_manual(values = status_cols, name = NULL) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = NULL,
       title = "Spike-gene coverage QC status by neighborhood-week",
       subtitle = "spike_breadth_d10 thresholds; grey = no sequencing attempted") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = "bottom")

p_cov <- seq_lw |>
  mutate(date = ew_to_date(year_epiweek)) |>
  ggplot(aes(x = date, y = spike_breadth_d10)) +
  geom_hline(yintercept = qc_thresholds[1], linetype = "dashed", color = "#8fce9b") +
  geom_hline(yintercept = qc_thresholds[2], linetype = "dotted", color = "#1a7f5a") +
  geom_point(aes(color = spike_breadth_d10 >= 0.5), size = 0.7, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#1a7f5a", `FALSE` = "#D7191C"),
                     guide = "none") +
  facet_wrap(~ LOCATION, ncol = 4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %y") +
  labs(x = NULL, y = "spike breadth at >=10x",
       title = "Spike coverage over time by neighborhood",
       subtitle = "dashed = 0.5, dotted = 0.7") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 7))

ggsave(paste0(fig_dir, "qc_status_grid.jpg"), p_status, width = 12, height = 5, dpi = 300)
ggsave(paste0(fig_dir, "qc_coverage_by_nb.jpg"), p_cov, width = 12, height = 8, dpi = 300)

# ---- 3. Neighborhood composition under each filter --------------------------
# Dropped weeks are rendered as an explicit grey "no data" band, NOT as zeros --
# collapsing them to 0 would read as "lineage absent" rather than "not measured".
make_nb_plot <- function(thresh, label) {
  ok <- if (is.null(thresh)) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(spike_breadth_d10 >= thresh) |> distinct(LOCATION, year_epiweek)
  }

  kept <- ww_lin |> semi_join(ok, by = c("LOCATION", "year_epiweek"))
  dropped <- loc_week |>
    anti_join(ok, by = c("LOCATION", "year_epiweek")) |>
    transmute(LOCATION, year_epiweek, date, lin = "not measured", abundance = 1)

  bind_rows(kept, dropped) |>
    mutate(lin = factor(as.character(lin), levels = plot_levels)) |>
    ggplot(aes(x = date, y = abundance, fill = lin)) +
    geom_col(width = 6) +
    facet_wrap(~ LOCATION, ncol = 4) +
    scale_fill_manual(values = fill_cols, drop = FALSE, name = NULL) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b %y") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "relative abundance",
         title = paste0("Neighborhood lineage composition - ", label),
         subtitle = paste0(nrow(ok), " of ", nrow(seq_lw),
                           " sequenced neighborhood-weeks retained")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 7),
          legend.position = "bottom")
}

specs <- list(list(NULL, "unfiltered", "unfiltered"),
              list(0.5, "spike coverage >= 0.5", "spike0.5"),
              list(0.7, "spike coverage >= 0.7", "spike0.7"))

for (s in specs) {
  p <- make_nb_plot(s[[1]], s[[2]])
  ggsave(paste0(fig_dir, "qc_nb_composition_", s[[3]], ".jpg"),
         p, width = 13, height = 8.5, dpi = 300)
}

# ---- 4. Population-weighted citywide under each filter ----------------------
# Weights are renormalised over the neighborhoods passing QC that week, so a
# dropped neighborhood does not drag the citywide composition toward zero.
citywide_weighted <- function(thresh, label) {
  ok <- if (is.null(thresh)) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(spike_breadth_d10 >= thresh) |> distinct(LOCATION, year_epiweek)
  }

  ww_lin |>
    semi_join(ok, by = c("LOCATION", "year_epiweek")) |>
    left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION") |>
    group_by(year_epiweek, date, lin) |>
    summarise(abundance = sum(abundance * pop) / sum(pop),
              n_loc = n_distinct(LOCATION),
              pop_frac = sum(pop) / sum(pop_wts$pop[!pop_wts$LOCATION %in% drop_locs]),
              .groups = "drop") |>
    mutate(filter = label)
}

city <- bind_rows(
  citywide_weighted(NULL, "unfiltered"),
  citywide_weighted(0.5,  "spike >= 0.5"),
  citywide_weighted(0.7,  "spike >= 0.7")
) |>
  mutate(filter = factor(filter, levels = c("unfiltered", "spike >= 0.5", "spike >= 0.7")))

p_city <- city |>
  ggplot(aes(x = date, y = abundance, fill = lin)) +
  geom_col(width = 6) +
  facet_wrap(~ filter, ncol = 1) +
  scale_fill_manual(values = fill_cols, drop = FALSE, name = "lineage") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = NULL, y = "population-weighted relative abundance",
       title = "Citywide population-weighted lineage composition under QC filters",
       subtitle = "weights renormalised each week over neighborhoods passing QC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

# how much of the city's population the weighted estimate actually rests on
p_support <- city |>
  distinct(date, filter, n_loc, pop_frac) |>
  ggplot(aes(x = date, y = pop_frac, color = filter)) +
  geom_step() +
  scale_color_manual(values = c("unfiltered" = "grey40",
                                "spike >= 0.5" = "#8fce9b",
                                "spike >= 0.7" = "#1a7f5a"), name = NULL) +
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
  ok <- if (is.null(s[[1]])) seq_lw else filter(seq_lw, spike_breadth_d10 >= s[[1]])
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

cat("\n--- named lineage groups surviving QC (all groups, pre-'other' collapse) ---\n")
for (s in specs) {
  ok <- if (is.null(s[[1]])) {
    seq_lw |> distinct(LOCATION, year_epiweek)
  } else {
    seq_lw |> filter(spike_breadth_d10 >= s[[1]]) |> distinct(LOCATION, year_epiweek)
  }
  n_grp <- ww |>
    semi_join(ok, by = c("LOCATION", "year_epiweek")) |>
    filter(abundance > 0) |>
    pull(sublin_collapse) |> unique() |> length()
  cat(sprintf("%-24s %3d groups detected\n", s[[2]], n_grp))
}

cat("\nFigures written to", fig_dir, "\n")
