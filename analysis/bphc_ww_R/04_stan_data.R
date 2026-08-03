# =============================================================================
# 04_stan_data.R
# Assemble the hierarchical Stan data list from 02/03 intermediates.
# Standalone: reads fresh, writes stan_data.rds. Rerun top to bottom.
# =============================================================================

# ---- PATHS (edit these) -----------------------------------------------------
model_ver <- "gravity_v2"
out_dir <- paste0("../model_data/", model_ver, "/")

# ---- Libraries --------------------------------------------------------------
library(tidyverse)

# ---- 1. Read intermediates --------------------------------------------------
transitions        <- read_rds(paste0(out_dir, "transitions.rds"))
neighbor_abund_mat <- read_rds(paste0(out_dir, "neighbor_abund_mat.rds"))
loc_index          <- read_rds(paste0(out_dir, "loc_index.rds"))
lineage_index      <- read_rds(paste0(out_dir, "lineage_index.rds"))
pop_vec            <- read_rds(paste0(out_dir, "pop_vec.rds"))
dist_mat           <- read_rds(paste0(out_dir, "dist_mat.rds"))

n_loc     <- nrow(loc_index)
n_lineage <- nrow(lineage_index)
n_obs     <- nrow(transitions)

# ---- 2. Alignment checks ----------------------------------------------------
stopifnot(
  ncol(neighbor_abund_mat) == n_loc,
  nrow(neighbor_abund_mat) == n_obs,
  length(pop_vec) == n_loc,
  nrow(dist_mat) == n_loc, ncol(dist_mat) == n_loc,
  identical(names(pop_vec), loc_index$LOCATION),
  identical(rownames(dist_mat), loc_index$LOCATION),
  identical(colnames(dist_mat), loc_index$LOCATION),
  all(transitions$j >= 1 & transitions$j <= n_loc),
  all(transitions$L >= 1 & transitions$L <= n_lineage),
  all(is.finite(transitions$log_growth)),
  all(is.finite(neighbor_abund_mat))
)

# ---- 3. Centered covariates + abundance rescale -----------------------------
# Centering log_pop / log_dist makes unit-scale priors on theta1/theta3
# weakly informative (raw exponent args would be O(8) / O(2)).
log_pop          <- log(pop_vec)
log_pop_mean     <- mean(log_pop)
log_pop_c        <- log_pop - log_pop_mean

offdiag          <- dist_mat[upper.tri(dist_mat)]
log_dist_mean    <- mean(log(offdiag))
log_dist_c       <- log(dist_mat) - log_dist_mean
diag(log_dist_c) <- 0          # placeholder; never read in the gravity sum

# Rescale donor abundance to O(1) so theta0 ~ normal(0,1) is sensible (raw
# neighbor_abund ~ 1e4-1e6 forces theta0 ~ 1e-7 and freezes HMC). One global
# scale across all lineages keeps the lineage-specific theta0_l comparable.
neighbor_scale   <- mean(neighbor_abund_mat[neighbor_abund_mat > 0])
neighbor_abund_s <- neighbor_abund_mat / neighbor_scale

# ---- 4. Assemble stan_data --------------------------------------------------
stan_data <- list(
  n_obs          = n_obs,
  n_loc          = n_loc,
  n_lineage      = n_lineage,
  L              = transitions$L,              # lineage index per obs
  j              = transitions$j,              # focal location per obs
  log_growth     = transitions$log_growth,
  neighbor_abund = neighbor_abund_s,           # rescaled; self/missing = 0
  log_pop_c      = as.vector(log_pop_c),
  log_dist_c     = log_dist_c
)

# ---- 5. Report --------------------------------------------------------------
cat("--- stan_data ---\n")
cat("n_obs:", n_obs, "| n_loc:", n_loc, "| n_lineage:", n_lineage, "\n")
cat("obs per lineage:\n")
print(transitions |> count(L, lineage) |> arrange(L))
cat("log_growth range: [",
    round(min(stan_data$log_growth), 2), ",",
    round(max(stan_data$log_growth), 2), "]\n")
cat("log_pop_c range:  [",
    round(min(stan_data$log_pop_c), 2), ",",
    round(max(stan_data$log_pop_c), 2), "]  (mean removed:",
    round(log_pop_mean, 3), ")\n")
cat("log_dist_c offdiag range: [",
    round(min(log_dist_c[upper.tri(log_dist_c)]), 2), ",",
    round(max(log_dist_c[upper.tri(log_dist_c)]), 2), "]  (mean removed:",
    round(log_dist_mean, 3), ")\n")
cat("neighbor_scale (nonzero mean):", round(neighbor_scale),
    "| rescaled nonzero max:", round(max(neighbor_abund_s), 2), "\n")

# ---- 6. Save ----------------------------------------------------------------
write_rds(stan_data, paste0(out_dir, "stan_data.rds"))
write_rds(list(log_pop_mean = log_pop_mean, log_dist_mean = log_dist_mean,
               neighbor_scale = neighbor_scale),
          paste0(out_dir, "centering.rds"))
cat("\nSaved stan_data / centering to", out_dir, "\n")
