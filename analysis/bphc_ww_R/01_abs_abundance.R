# =============================================================================
# 01_abs_abundance.R
# Build absolute-abundance proxy for all modelled lineages from
# ww_collapsed_complete. Standalone: reads raw data, writes abund.rds.
# Rerun top-to-bottom.
# =============================================================================

# ---- PATHS ------------------------------------------------------------------
data_dir <- "../data/"        # raw inputs
out_dir  <- "../model_data/"  # intermediate objects
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Lineages to pool, in chronological order (defines the L index in 02).
lineages  <- c("BA.5.x", "BQ.1.x", "XBB.1.x", "EG.5.x", "HV.1.x", "JN.1.x")
drop_locs <- "Lower Roxbury"   # no detections across sampled weeks -> dropped

# ---- Libraries --------------------------------------------------------------
library(tidyverse)

# ---- 1. Read input ----------------------------------------------------------
ww_collapsed_complete <- read_rds(paste0(data_dir, "ww_collapsed_complete.rds"))

# ---- 2. Subset to lineages, drop excluded locations, build proxy ------------
abund <- ww_collapsed_complete |>
  filter(sublin_collapse %in% lineages,
         !LOCATION %in% drop_locs) |>
  mutate(abs_abundance = abundance * meaneffCopiesL)
# abs_abundance is NA where meaneffCopiesL is NA (week not sampled) and 0 where
# sampled but the lineage was not detected. The not-sampled mask is lineage-
# independent (concentration is), so it is identical across lineages; detection
# (abundance > 0) is lineage-specific. 02 keys off both.

# Guard: confirm every requested lineage is actually present.
missing_lin <- setdiff(lineages, unique(abund$sublin_collapse))
if (length(missing_lin) > 0)
  warning("requested lineage(s) not found in data: ",
          paste(missing_lin, collapse = ", "))

# ---- 3. Diagnostics per lineage (printed, not saved) ------------------------
cat("\n--- counts by lineage ---\n")
abund |>
  group_by(sublin_collapse) |>
  summarise(
    n_total       = n(),
    n_not_sampled = sum(is.na(meaneffCopiesL)),
    n_zero_detect = sum(abundance == 0 & !is.na(meaneffCopiesL)),
    n_detected    = sum(abundance > 0 & !is.na(meaneffCopiesL)),
    frac_detected = n_detected / (n_detected + n_zero_detect),
    .groups = "drop"
  ) |>
  print()
# n_not_sampled should be identical across lineages (shared sampling structure).

# ---- 4. Save ----------------------------------------------------------------
write_rds(abund, paste0(out_dir, "abund.rds"))
cat("\nSaved:", paste0(out_dir, "abund.rds"),
    "| rows:", nrow(abund),
    "| lineages:", n_distinct(abund$sublin_collapse),
    "| locations:", n_distinct(abund$LOCATION), "\n")
