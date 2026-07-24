# =============================================================================
# 06_plot_results.R
# Plot the hierarchical multi-lineage fit. Standalone: reads the saved fit +
# intermediates, writes PNGs to draft_figures/. Rerun top to bottom.
#
# Figures:
#   1. Connectivity map (null-aware edges)        fig_connectivity_map.png
#   2. Baseline growth r by lineage               fig_r_by_lineage.png
#   3. Coupling magnitude theta0 by lineage       fig_theta0_by_lineage.png
#   4. Shared exponents: posterior vs prior       fig_shared_exponents.png
#   5. Per-pair distance contribution             fig_distance_contribution.png
# =============================================================================

# ---- PATHS (edit these) -----------------------------------------------------
out_dir <- "../model_data/"
fig_dir <- "../draft_figures/"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries --------------------------------------------------------------
library(cmdstanr)     # to use methods on the saved fit object
library(tidyverse)
library(sf)
library(posterior)

# ---- 1. Read fit + intermediates --------------------------------------------
fit           <- read_rds(paste0(out_dir, "gravity_fit.rds"))
loc_index     <- read_rds(paste0(out_dir, "loc_index.rds"))
lineage_index <- read_rds(paste0(out_dir, "lineage_index.rds"))
dist_mat      <- read_rds(paste0(out_dir, "dist_mat.rds"))
pop_vec       <- read_rds(paste0(out_dir, "pop_vec.rds"))
catchments    <- read_rds(paste0(out_dir, "catchments.rds"))
cent_const    <- read_rds(paste0(out_dir, "centering.rds"))
# from lineage_timing.R; mean_lag = mean arrival lag in weeks (across lineages)
# after that lineage's citywide first detection, per catchment
detection_timing <- read_rds(paste0(out_dir, "neighborhood_detection_timing.rds"))
catchments <- catchments |>
  left_join(detection_timing |> select(LOCATION, mean_lag), by = "LOCATION")

n_loc     <- nrow(loc_index)
n_lineage <- nrow(lineage_index)
lin_levels <- lineage_index$lineage[order(lineage_index$L)]   # chronological

dr <- as_draws_df(fit$draws())

# Centered covariates, reconstructed from raw inputs + saved constants
log_pop_c  <- log(pop_vec) - cent_const$log_pop_mean              # named, loc order
log_dist_c <- log(dist_mat) - cent_const$log_dist_mean            # n_loc x n_loc

# Helper: median + 90% CI for an indexed parameter vector (e.g. r[1..n])
param_ci <- function(base, labels) {
  map_dfr(seq_along(labels), function(i) {
    x <- dr[[sprintf("%s[%d]", base, i)]]
    tibble(L = i, lineage = labels[i],
           median = median(x),
           q5 = quantile(x, 0.05), q95 = quantile(x, 0.95))
  }) |>
    mutate(lineage = factor(lineage, levels = lin_levels))
}

# =============================================================================
# FIGURE 1 — Connectivity map (null-aware edges)
# =============================================================================
# Edge quantity is the SHARED log-kernel for each pair:
#   logK_jk = theta1*(p_j + p_k) - theta3*d_jk
# "Resolved" = 90% CI of logK_jk excludes 0 (the no-modulation / neutral value
# K=1). With theta1, theta3 ~ 0, no edge is resolved -> all dashed/faint, which
# is the honest visual statement of the null.
th1 <- dr$theta1
th3 <- dr$theta3

centroids <- st_centroid(catchments)
nc <- st_coordinates(centroids)
node <- tibble(LOCATION = catchments$LOCATION, x = nc[, 1], y = nc[, 2],
               pop = as.numeric(pop_vec[catchments$LOCATION]))

prs <- which(upper.tri(dist_mat), arr.ind = TRUE)
edges <- map_dfr(seq_len(nrow(prs)), function(p) {
  jj <- prs[p, 1]; kk <- prs[p, 2]
  logK <- th1 * (log_pop_c[jj] + log_pop_c[kk]) - th3 * log_dist_c[jj, kk]
  tibble(j = jj, k = kk,
         med_logK = median(logK),
         q5 = quantile(logK, 0.05), q95 = quantile(logK, 0.95))
}) |>
  mutate(resolved = q5 > 0 | q95 < 0,
         x    = node$x[j],    y    = node$y[j],
         xend = node$x[k],    yend = node$y[k])

n_resolved <- sum(edges$resolved)

fig1 <- ggplot() +
  geom_sf(data = catchments, aes(fill = mean_lag), color = "grey80", linewidth = 0.3) +
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend,
                   color = med_logK, linetype = resolved, alpha = resolved,
                   linewidth = resolved)) +
  geom_point(data = node, aes(x = x, y = y, size = pop),
             color = "grey20", fill = "white", shape = 21, stroke = 0.6) +
  geom_text(data = node, aes(x = x, y = y, label = LOCATION),
            size = 2.6, vjust = -1.1, color = "grey20") +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                      name = "mean detection\nlag (weeks)") +
  scale_color_gradient2(low = "#2c7bb6", mid = "grey75", high = "#d7191c",
                        midpoint = 0, name = "median log-kernel") +
  scale_linetype_manual(values = c(`FALSE` = "dashed", `TRUE` = "solid"),
                        name = "coupling 90% CI\nexcludes neutral") +
  scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 0.95), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 1.1), guide = "none") +
  scale_size_area(max_size = 7, name = "catchment pop.") +
  coord_sf(datum = 26986) +
  labs(title = "Inter-catchment coupling (shared gravity kernel)",
       subtitle = paste0(n_resolved, " of ", nrow(edges),
                         " pairs resolved away from neutral coupling"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), panel.grid = element_blank())

ggsave(paste0(fig_dir, "fig_connectivity_map.png"), fig1,
       width = 8, height = 7, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 2 — Baseline growth rate r by lineage
# =============================================================================
r_df <- param_ci("r", lineage_index$lineage[order(lineage_index$L)])

fig2 <- ggplot(r_df, aes(x = median, y = fct_rev(lineage))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_pointrange(aes(xmin = q5, xmax = q95)) +
  labs(title = "Lineage-specific baseline log-growth (r)",
       subtitle = "median + 90% credible interval",
       x = "r (log-growth per week)", y = NULL) +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_r_by_lineage.png"), fig2,
       width = 7, height = 4, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 3 — Coupling magnitude theta0 by lineage
# =============================================================================
t0_df <- param_ci("theta0", lineage_index$lineage[order(lineage_index$L)])

fig3 <- ggplot(t0_df, aes(x = median, y = fct_rev(lineage))) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_pointrange(aes(xmin = q5, xmax = q95)) +
  labs(title = "Lineage-specific coupling magnitude (theta0)",
       subtitle = "median + 90% CI; intervals crossing 0 = no detectable coupling",
       x = "theta0", y = NULL) +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_theta0_by_lineage.png"), fig3,
       width = 7, height = 4, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 4 — Shared exponents: posterior vs prior
# =============================================================================
# The crux plot: posterior of theta1 / theta3 against their N(0,1) prior.
# Posterior sitting on (and tightening within) the prior is the null.
shared <- tibble(theta1 = dr$theta1, theta3 = dr$theta3) |>
  pivot_longer(everything(), names_to = "param", values_to = "value")

prior_df <- tibble(value = seq(-3.5, 3.5, length.out = 400),
                   dens  = dnorm(value, 0, 1))

fig4 <- ggplot(shared, aes(value)) +
  geom_density(aes(fill = param), color = NA, alpha = 0.6) +
  geom_line(data = prior_df, aes(value, dens),
            linetype = 2, color = "grey30") +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  facet_wrap(~ param, ncol = 1) +
  scale_fill_manual(values = c(theta1 = "#7570b3", theta3 = "#d95f02"),
                    guide = "none") +
  labs(title = "Shared gravity exponents: posterior vs. N(0,1) prior",
       subtitle = "dashed = prior; posteriors track the prior (mass & distance unresolved)",
       x = "exponent value", y = "density") +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_shared_exponents.png"), fig4,
       width = 7, height = 6, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 5 — Per-pair distance contribution
# =============================================================================
# exp(-theta3 * centered-log-dist) per catchment pair vs raw km, median + 90%CI.
dist_contrib <- map_dfr(seq_len(nrow(prs)), function(p) {
  jj <- prs[p, 1]; kk <- prs[p, 2]
  mult <- exp(-th3 * log_dist_c[jj, kk])
  tibble(dist_km = dist_mat[jj, kk],
         median = median(mult),
         lo = quantile(mult, 0.05), hi = quantile(mult, 0.95))
})

fig5 <- ggplot(dist_contrib, aes(dist_km, median)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey60") +
  geom_vline(xintercept = exp(cent_const$log_dist_mean),
             linetype = 3, color = "steelblue") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3, alpha = 0.7) +
  scale_y_log10() +
  labs(title = "Per-pair distance contribution to coupling",
       subtitle = "exp(-theta3 * centered log-dist); dashed = no effect (1), dotted = mean-dist pair",
       x = "catchment-pair distance (km)", y = "distance multiplier (90% CI)") +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_distance_contribution.png"), fig5,
       width = 7, height = 5, dpi = 200, bg = "white")

# =============================================================================
# FIGURE 6 — Informative null: coupling posteriors vs. prior
# =============================================================================
# All coupling params share a Normal(0,1) prior. Posterior intervals (point =
# median, thick = 50%, thin = 90%) against the grey band = prior central 90%.
# Every coupling parameter concentrates at ~0, inside the prior -> the data
# actively constrain coupling toward zero, not merely fail to find it.
null_map <- tibble(
  col   = c("theta1", "theta3", sprintf("theta0[%d]", seq_len(n_lineage))),
  label = c("theta1  (mass, shared)",
            "theta3  (distance, shared)",
            paste0("theta0 - ", lin_levels))
)

null_df <- pmap_dfr(null_map, function(col, label) {
  x <- dr[[col]]
  tibble(label  = label,
         median = median(x),
         lo90 = quantile(x, 0.05), hi90 = quantile(x, 0.95),
         lo50 = quantile(x, 0.25), hi50 = quantile(x, 0.75))
}) |>
  mutate(label = factor(label, levels = rev(null_map$label)))

prior_q90 <- qnorm(c(0.05, 0.95))   # central 90% of N(0,1) = +/- 1.645

fig6 <- ggplot(null_df, aes(y = label)) +
  annotate("rect", xmin = prior_q90[1], xmax = prior_q90[2],
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_linerange(aes(xmin = lo90, xmax = hi90)) +
  geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.5) +
  geom_point(aes(x = median), size = 2) +
  labs(title = "Coupling posteriors vs. Normal(0,1) prior",
       subtitle = "grey = prior central 90%; bars = posterior 90% (thin) / 50% (thick)",
       x = "parameter value", y = NULL) +
  theme_minimal(base_size = 11)

ggsave(paste0(fig_dir, "fig_informative_null.png"), fig6,
       width = 7.5, height = 5, dpi = 200, bg = "white")

# ---- Done -------------------------------------------------------------------
cat("Saved 6 figures to", fig_dir, "\n")
cat("Connectivity map:", n_resolved, "of", nrow(edges), "pairs resolved.\n")
