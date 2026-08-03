# =============================================================================
# lineage_timing.R
# Explore when each modelled lineage first arrives in each neighborhood, and
# how much arrival timing varies across neighborhoods (which catchments lead,
# which lag) -- both per lineage (fig_lineage_arrival_lag.png) and collapsed
# to a per-neighborhood leader/laggard summary (fig_neighborhood_detection_
# timing.png). Standalone: reads abund.rds + shared indices from the gravity-
# model pipeline (01/02), does not modify their outputs. Rerun top-to-bottom.
# =============================================================================

# ---- PATHS (edit these) ------------------------------------------------------
model_ver <- "gravity_v2"
arrival_ver <- "arrival_v2"
data_dir <- paste0("../model_data/", model_ver, "/")
fig_dir  <- paste0("../draft_figures/", arrival_ver, "/")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries ----------------------------------------------------------------
library(tidyverse)

# ---- 1. Read input --------------------------------------------------------
abund         <- read_rds(paste0(data_dir, "abund.rds"))
loc_index     <- read_rds(paste0(data_dir, "loc_index.rds"))
week_index    <- read_rds(paste0(data_dir, "week_index.rds"))
lineage_index <- read_rds(paste0(data_dir, "lineage_index.rds"))

lin_levels <- lineage_index$lineage[order(lineage_index$L)]  # chronological

# ---- 2. First-detection week per lineage x neighborhood -----------------------
# "Arrival" = first sampled week (abs_abundance finite, i.e. the site was
# actually sampled that week) where the lineage's relative abundance clears
# detect_thresh. A hard 0 cutoff is noisy at low sequencing depth, so this
# requires a minimum relative-abundance call rather than any nonzero call.
detect_thresh <- 0.05

arrival <- abund |>
  filter(is.finite(abs_abundance)) |>
  left_join(week_index, by = "time") |>
  left_join(loc_index, by = "LOCATION") |>
  group_by(sublin_collapse, LOCATION, j) |>
  summarise(
    first_t = if (any(abundance > detect_thresh)) min(t[abundance > detect_thresh]) else NA_integer_,
    .groups = "drop"
  ) |>
  rename(lineage = sublin_collapse) |>
  filter(lineage %in% lin_levels) |>
  mutate(lineage = factor(lineage, levels = lin_levels)) |>
  left_join(week_index, by = c("first_t" = "t"))

# citywide arrival = earliest first_t across neighborhoods, per lineage
citywide <- arrival |>
  filter(!is.na(first_t)) |>
  group_by(lineage) |>
  summarise(citywide_t = min(first_t), .groups = "drop")

arrival <- arrival |>
  left_join(citywide, by = "lineage") |>
  mutate(lag_weeks = first_t - citywide_t)

# ---- 3. Diagnostics ------------------------------------------------------------
cat("\n--- citywide first-detection week, by lineage ---\n")
citywide |>
  left_join(week_index, by = c("citywide_t" = "t")) |>
  select(lineage, citywide_t, week_label) |>
  print()

cat("\n--- neighborhoods that never detected a given lineage ---\n")
arrival |> filter(is.na(first_t)) |> count(lineage, name = "n_locations_missed") |> print()

cat("\n--- arrival lag (weeks after citywide first detection), by lineage ---\n")
arrival |>
  filter(!is.na(lag_weeks)) |>
  group_by(lineage) |>
  summarise(mean_lag = mean(lag_weeks), max_lag = max(lag_weeks), .groups = "drop") |>
  print()

cat("\n--- neighborhoods ranked by mean arrival lag (leaders first) ---\n")
arrival |>
  filter(!is.na(lag_weeks)) |>
  group_by(LOCATION) |>
  summarise(mean_lag = mean(lag_weeks), .groups = "drop") |>
  arrange(mean_lag) |>
  print()

# ---- 4. Plot: arrival-lag heatmap ----------------------------------------------
loc_order <- arrival |>
  filter(!is.na(lag_weeks)) |>
  group_by(LOCATION) |>
  summarise(mean_lag = mean(lag_weeks), .groups = "drop") |>
  arrange(mean_lag) |>
  pull(LOCATION)

fig1 <- arrival |>
  mutate(LOCATION = factor(LOCATION, levels = rev(loc_order))) |>
  ggplot(aes(x = lineage, y = LOCATION, fill = lag_weeks)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(lag_weeks), "-", lag_weeks)), size = 3) +
  scale_fill_gradient(low = "grey95", high = "#d7191c", na.value = "grey40",
                      name = "weeks after\ncitywide first\ndetection") +
  labs(title = "Lineage arrival lag by neighborhood",
       subtitle = "0 = first neighborhood to detect that lineage; grey = never detected",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(paste0(fig_dir, "fig_lineage_arrival_lag.png"), fig1,
       width = 8, height = 6, dpi = 200, bg = "white")

# ---- 5. Save --------------------------------------------------------------------
write_rds(arrival, paste0(data_dir, "lineage_arrival.rds"))
cat("\nSaved lineage_arrival.rds and fig_lineage_arrival_lag.png to", fig_dir, "\n")

# ---- 6. Per-neighborhood detection timing (leader/laggard summary) -------------
# Collapses lag_weeks to one row per neighborhood: is this catchment
# consistently early across lineages, or does its low mean lag just reflect
# missing the slow-to-arrive lineages entirely? n_missed distinguishes the two.
loc_summary <- arrival |>
  group_by(LOCATION) |>
  summarise(
    n_detected = sum(!is.na(lag_weeks)),
    n_missed   = sum(is.na(lag_weeks)),
    mean_lag   = mean(lag_weeks, na.rm = TRUE),
    se_lag     = if (n_detected > 1) sd(lag_weeks, na.rm = TRUE) / sqrt(n_detected) else NA_real_,
    .groups = "drop"
  ) |>
  arrange(mean_lag)

cat("\n--- per-neighborhood detection summary (leaders first) ---\n")
print(loc_summary, n = Inf)

fig2 <- ggplot() +
  geom_jitter(
    data = arrival |> filter(!is.na(lag_weeks)) |> mutate(LOCATION = factor(LOCATION, levels = rev(loc_order))),
    aes(x = lag_weeks, y = LOCATION),
    width = 0, height = 0.15, alpha = 0.35, size = 1.8, color = "#2c7fb8"
  ) +
  geom_errorbarh(
    data = loc_summary |> mutate(LOCATION = factor(LOCATION, levels = rev(loc_order))),
    aes(xmin = mean_lag - se_lag, xmax = mean_lag + se_lag, y = LOCATION),
    height = 0, color = "#d7191c"
  ) +
  geom_point(
    data = loc_summary |> mutate(LOCATION = factor(LOCATION, levels = rev(loc_order))),
    aes(x = mean_lag, y = LOCATION, shape = factor(n_missed)),
    color = "#d7191c", size = 3, stroke = 1.2
  ) +
  scale_shape_manual(name = "# lineages\nnever detected", values = c("0" = 16, "1" = 4)) +
  labs(title = "Neighborhood detection timing across lineages",
       subtitle = "Blue dots = per-lineage lag (weeks after citywide first detection); red = mean ± SE",
       x = "arrival lag (weeks)", y = NULL) +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_neighborhood_detection_timing.png"), fig2,
       width = 8, height = 6, dpi = 200, bg = "white")

write_rds(loc_summary, paste0(data_dir, "neighborhood_detection_timing.rds"))
cat("\nSaved neighborhood_detection_timing.rds and fig_neighborhood_detection_timing.png to", fig_dir, "\n")
