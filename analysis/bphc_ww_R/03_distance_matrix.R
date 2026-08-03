# =============================================================================
# 03_distance_matrix.R
# Sewershed shapefiles (+ population table) -> dist_mat.rds, pop_vec.rds.
# Both are aligned to loc_index$j order from 02. Standalone; rerun top to bottom.
# =============================================================================

# ---- PATHS / PARAMS (edit these) --------------------------------------------
model_ver   <- "gravity_v2"
out_dir     <- paste0("../model_data/", model_ver, "/")
shared_dir  <- "../model_data/"   # shapefiles, shared across versions
shp_path    <- paste0(shared_dir, "bphc_genomic_sites/genomic_sites.shp")   # EDIT
pop_path <- "../data/pop_wts.rds"             # EDIT
proj_crs <- 26986       # NAD83 / Mass Mainland (meters); correct for Boston

# population table is expected to have a name column and a population column;
# set these to the actual column names in pop_path:
pop_name_col <- "LOCATION"
pop_val_col  <- "pop"

# ---- Libraries --------------------------------------------------------------
library(tidyverse)
library(sf)

# ---- 1. Read indices + shapefile --------------------------------------------
loc_index <- read_rds(paste0(out_dir, "loc_index.rds"))
n_loc <- nrow(loc_index)
shp   <- st_read(shp_path, quiet = TRUE)

# ---- 2. Inspect (printed) ---------------------------------------------------
# Run once to discover the right name_col / CRS, then set params above.
cat("--- shapefile structure ---\n")
cat("features:", nrow(shp),
    "| geometry:", as.character(st_geometry_type(shp, by_geometry = FALSE)[1]), "\n")
cat("CRS:", st_crs(shp)$input, "\n")
cat("columns:", paste(names(shp), collapse = ", "), "\n")
char_cols <- names(st_drop_geometry(shp))[
  map_lgl(st_drop_geometry(shp), is.character)]
cat("name-candidate (character) columns:", paste(char_cols, collapse = ", "), "\n")

# ---- 3. Crosswalk shapefile codes -> LOCATION -------------------------------
# `nh` holds abbreviation codes, not full names. Map explicitly (auditable).
crosswalk <- tibble::tribble(
  ~nh,        ~LOCATION,
  "A/B",      "Allston-Brighton",
  "BB",       "Back Bay",
  "CH",       "Charlestown",
  "DOR 2224", "Dorchester",
  "EB",       "East Boston",
  "HP",       "Hyde Park",
  "JP",       "Jamaica Plain",
  "MT",       "Mattapan",
  "RS/WR",    "Roslindale West Roxbury",
  "RX",       "Roxbury",
  "SB",       "South Boston"
)
# "LR" (Lower Roxbury) intentionally absent -- dropped in 01 (no detections).

# Join codes -> LOCATION (drops LR and any unmatched site), then dissolve so
# multi-polygon sites collapse to one geometry per LOCATION. Mattapan has two
# sub-catchments (site 09 A/B); st_union merges them. No-op for singletons.
shp <- shp |>
  mutate(nh = as.character(nh)) |>
  inner_join(crosswalk, by = "nh") |>
  group_by(LOCATION) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

# Coverage: every modelled LOCATION must have exactly one polygon.
missing_loc <- setdiff(loc_index$LOCATION, shp$LOCATION)
if (length(missing_loc) > 0) {
  stop("No polygon for: ", paste(missing_loc, collapse = ", "),
       " -- resolve before computing distances ",
       "(supply a centroid, use a newer shapefile, or drop the location).")
}
stopifnot(nrow(shp) == n_loc)

# ---- 4. Align polygon rows to loc_index order -------------------------------
shp_aligned <- shp[match(loc_index$LOCATION, shp$LOCATION), ]
stopifnot(identical(shp_aligned$LOCATION, loc_index$LOCATION))

# ---- 5. Distance matrix (km) ------------------------------------------------
if (is.na(st_crs(shp_aligned)))
  stop("shapefile CRS is undefined; set it (st_set_crs) before projecting")
cent <- shp_aligned |> st_transform(proj_crs) |> st_centroid()
dist_mat <- units::drop_units(st_distance(cent)) / 1000   # m -> km
dimnames(dist_mat) <- list(loc_index$LOCATION, loc_index$LOCATION)

# ---- 6. Population vector ----------------------------------------------------
# Separate table joined by name. To instead pull from a shapefile column:
#   pop_vec <- as.numeric(shp_aligned[["POP_COLUMN"]]); names(pop_vec) <- loc_index$LOCATION
pop_raw <- readRDS(pop_path) |>
  rename(LOCATION = all_of(pop_name_col), population = all_of(pop_val_col))
pop_vec <- pop_raw$population[match(loc_index$LOCATION, pop_raw$LOCATION)]
names(pop_vec) <- loc_index$LOCATION

# ---- 7. Sanity checks -------------------------------------------------------
stopifnot(
  nrow(dist_mat) == n_loc, ncol(dist_mat) == n_loc,
  all(diag(dist_mat) == 0),
  isSymmetric(unname(dist_mat)),
  all(dist_mat[upper.tri(dist_mat)] > 0),
  length(pop_vec) == n_loc,
  all(is.finite(pop_vec)), all(pop_vec > 0)
)

cat("\n--- dist_mat (km) ---\n")
print(round(dist_mat, 1))
cat("\n--- pop_vec ---\n")
print(pop_vec)

# ---- 8. Save ----------------------------------------------------------------
write_rds(dist_mat, paste0(out_dir, "dist_mat.rds"))
write_rds(pop_vec,  paste0(out_dir, "pop_vec.rds"))
write_rds(shp_aligned, paste0(out_dir, "catchments.rds"))  # dissolved polygons, loc_index order
cat("\nSaved dist_mat / pop_vec / catchments to", out_dir, "\n")

# -----------------------------------------------------------------------------
# NOTE — `nh` codes are crosswalked, not name-matched. Multi-polygon sites are
# dissolved by LOCATION (step 3), so adding/removing a site is a crosswalk edit.
# Distances are centroid-to-centroid in EPSG:26986 (planar, km). Swap candidates
# for later: shared-boundary length, or commute/mobility travel time.
# -----------------------------------------------------------------------------
