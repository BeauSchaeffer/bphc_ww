# =============================================================================
# 02_build_transitions.R
# Build valid log-growth transitions + neighbor-abundance matrix for every
# lineage, stacked with a lineage index L. Standalone: reads abund.rds, writes
# transitions / neighbor_abund_mat / loc_index / week_index / lineage_index.
# Rerun top-to-bottom.
# =============================================================================

# ---- PATHS (edit these) -----------------------------------------------------
out_dir <- "../model_data/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Libraries --------------------------------------------------------------
library(tidyverse)

# ---- 1. Read input ----------------------------------------------------------
abund <- read_rds(paste0(out_dir, "abund.rds"))

# ---- 2. Indices: lineage, location, week (all global) -----------------------
# Lineage order (= L) is fixed chronologically; keep only those present.
lineage_index <- tibble(
  lineage = c("BA.5.x", "BQ.1.x", "XBB.1.x", "EG.5.x", "HV.1.x", "JN.1.x", "KP.2.x")
) |>
  filter(lineage %in% abund$sublin_collapse) |>
  mutate(L = row_number())

# Location and week indices are shared across lineages (same sampling structure).
loc_index <- abund |>
  distinct(LOCATION) |>
  arrange(LOCATION) |>
  mutate(j = row_number())
n_loc <- nrow(loc_index)

week_index <- abund |>
  distinct(time) |>
  arrange(time) |>
  mutate(t = row_number(), week_label = as.character(time))
n_weeks <- nrow(week_index)

# ---- 3. Per-lineage transition builder --------------------------------------
# Same logic as the single-lineage build: a transition is valid only if the
# focal location is finite and > 0 at BOTH t and t+1 (is.finite excludes NA /
# NaN / +/-Inf; > 0 excludes not-detected and the 0<->+ edge transitions).
# Indexed vectors are pre-extracted before the tibble to avoid tibble()'s
# sequential column masking. Neighbor abundances at t_now: non-finite/self -> 0.
build_lineage <- function(lin) {
  d <- abund |>
    filter(sublin_collapse == lin) |>
    select(time, LOCATION, abs_abundance) |>
    left_join(week_index, by = "time") |>
    left_join(loc_index,  by = "LOCATION")

  A_mat <- matrix(NA_real_, nrow = n_weeks, ncol = n_loc)
  A_mat[cbind(d$t, d$j)] <- d$abs_abundance

  a_now  <- A_mat[1:(n_weeks - 1), , drop = FALSE]
  a_next <- A_mat[2:n_weeks,       , drop = FALSE]
  valid  <- is.finite(a_now) & is.finite(a_next) & a_now > 0 & a_next > 0
  idx    <- which(valid, arr.ind = TRUE)        # cols: row (= t), col (= j)

  a_now_v  <- a_now[valid]
  a_next_v <- a_next[valid]
  tr <- tibble(
    lineage    = lin,
    j          = idx[, "col"],
    t_now      = idx[, "row"],
    t_next     = idx[, "row"] + 1L,
    a_now      = a_now_v,
    a_next     = a_next_v,
    log_growth = log(a_next_v / a_now_v)
  )

  Nsrc <- a_now
  Nsrc[!is.finite(Nsrc)] <- 0
  nb <- Nsrc[idx[, "row"], , drop = FALSE]
  nb[cbind(seq_len(nrow(idx)), idx[, "col"])] <- 0

  list(tr = tr, nb = nb)
}

# ---- 4. Build and stack -----------------------------------------------------
res <- map(lineage_index$lineage, build_lineage)
# bind_rows and rbind preserve lineage order, so rows stay aligned.
transitions <- bind_rows(map(res, "tr")) |>
  left_join(lineage_index, by = "lineage") |>
  relocate(lineage, L)
neighbor_abund_mat <- do.call(rbind, map(res, "nb"))

# ---- 5. Sanity checks -------------------------------------------------------
stopifnot(
  nrow(neighbor_abund_mat) == nrow(transitions),
  ncol(neighbor_abund_mat) == n_loc,
  all(is.finite(transitions$log_growth)),
  all(transitions$j >= 1 & transitions$j <= n_loc),
  all(transitions$L >= 1 & transitions$L <= nrow(lineage_index))
)

cat("\n--- transition build ---\n")
cat("n_loc:", n_loc, "| n_weeks:", n_weeks,
    "| n_lineage:", nrow(lineage_index), "\n")
cat("n_transitions (total):", nrow(transitions), "\n")
cat("log_growth range: [",
    round(min(transitions$log_growth), 3), ",",
    round(max(transitions$log_growth), 3), "]\n")
cat("transitions per lineage:\n")
print(transitions |> count(L, lineage) |> arrange(L))

# ---- 6. Save ----------------------------------------------------------------
write_rds(transitions,        paste0(out_dir, "transitions.rds"))
write_rds(neighbor_abund_mat, paste0(out_dir, "neighbor_abund_mat.rds"))
write_rds(loc_index,          paste0(out_dir, "loc_index.rds"))
write_rds(week_index,         paste0(out_dir, "week_index.rds"))
write_rds(lineage_index,      paste0(out_dir, "lineage_index.rds"))
cat("\nSaved transitions / neighbor_abund_mat / loc_index / week_index /",
    "lineage_index to", out_dir, "\n")
