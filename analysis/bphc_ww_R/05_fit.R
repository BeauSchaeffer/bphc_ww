# =============================================================================
# 05_fit.R
# Compile gravity_model.stan and fit to stan_data.rds (hierarchical, multi-
# lineage). Standalone: reads fresh, writes fit + summary. Rerun top to bottom.
# =============================================================================

# ---- PATHS (edit these) -----------------------------------------------------
out_dir   <- "../model_data/"
stan_file <- "gravity_model.stan"

# ---- Libraries --------------------------------------------------------------
library(cmdstanr)
library(tidyverse)   # read_rds / write_rds

# ---- 1. Read data -----------------------------------------------------------
stan_data <- read_rds(paste0(out_dir, "stan_data.rds"))
n_lineage <- stan_data$n_lineage

# ---- 2. Compile -------------------------------------------------------------
mod <- cmdstan_model(stan_file)

# ---- 3. Inits ---------------------------------------------------------------
# Explicit, tight inits (wide unconstrained defaults cause divergences).
# r / theta0 are per-lineage vectors now.
init_fun <- function() list(
  r      = rep(0,    n_lineage),
  theta0 = rep(0.01, n_lineage),
  theta1 = 0,
  theta3 = 1,
  sigma  = 1,
  nu     = 4
)

# ---- 4. Fit -----------------------------------------------------------------
fit <- mod$sample(
  data            = stan_data,
  init            = init_fun,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 1,
  adapt_delta     = 0.95,
  refresh         = 200
)

# ---- 5. Diagnostics + summary ----------------------------------------------
fit$cmdstan_diagnose()
print(fit$summary(c("r", "theta0", "theta1", "theta3", "sigma", "nu")),
      n = Inf)

# ---- 6. Save ----------------------------------------------------------------
fit$save_object(paste0(out_dir, "gravity_fit.rds"))
write_rds(fit$summary(), paste0(out_dir, "gravity_summary.rds"))
cat("\nSaved gravity_fit / gravity_summary to", out_dir, "\n")
