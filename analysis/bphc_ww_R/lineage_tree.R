# Illustrative lineage tree, rooted at BA.2. Backbone branch lengths are
# real barcode mutation distances between named lineages. Triangle sizes
# are currently fixed (real observed-diversity counting is disabled below,
# not deleted -- uncomment to restore it).


library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(ggrepel)


# Config ----------------------------------------------------------------

# edit to change which lineages get their own triangle; matches 01_abs_abundance.R,
# plus BA.2.86 (JN.1's direct ancestor, kept small/faded since it's not a modeled lineage)
major_lineages <- c("BA.5.x", "BQ.1.x", "XBB.1.x", "EG.5.x", "HV.1.x", "JN.1.x", "KP.x", "BA.2.86")
small_buckets <- c("BA.2.86")  # smaller triangle + lower alpha

root_lineage <- "BA.2"

manual_collapse_map <- list("KP.x" = c("KP.1.x", "KP.2.x", "KP.3.x"))

# real nesting (checked via barcode distance): child bucket -> parent bucket;
# only controls draw order (parents drawn under children) now that triangles
# don't overlap, not fill color
nested_bucket_parent <- c(
  "BQ.1.x" = "BA.5.x",
  "EG.5.x" = "XBB.1.x",
  "HV.1.x" = "EG.5.x",
  "KP.x"   = "JN.1.x"
)

# snapshot of RAs.R's fill_cols (viridis turbo over its full historical lin_levels)
# for these buckets, so this figure's colors match the RA plots. Re-extract from
# RAs.R (stop before its write_rds() calls) if that lineage set changes.
bucket_colors <- c(
  "BA.5.x"    = "#392B74FF",
  "BQ.1.x"    = "#4772EAFF",
  "XBB.1.x"   = "#E03F08FF",
  "EG.5.x"    = "#35ACF8FF",
  "HV.1.x"    = "#D9E436FF",
  "JN.1.x"    = "#EECF3AFF",
  "KP.x"      = "#FDAE35FF",
  "BA.2.86"   = "#341C53FF"
)

barcode_path <- "../data/reference/usher_barcodes.csv"

triangle_depth <- 6    # fixed x-extent, in mutation-distance units
triangle_width <- 0.3  # fixed y half-width, in tip-row units (spacing = 1)

# manual_exclude_lineages <- c()        # disabled below with the bucketing filter
# storage_dir <- "../data/"             # disabled below with the Freyja parsing
# triangle_size_cap_n <- 60             # disabled: was for count-based sizing
# triangle_depth_range <- c(3, 10)      # disabled: was for count-based sizing
# triangle_width_range <- c(0.15, 0.42) # disabled: was for count-based sizing


# Observed fine lineages (Freyja calls) -----------------------------------
# DISABLED -- real per-bucket diversity counting; uncomment to restore

# parse_freyja_aggregate <- function(path) {
#   results <- read.table(path, fill = TRUE, sep = "\t", h = TRUE)
#   results <- as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", "")))
#   results <- as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x))))
#   for (i in 1:nrow(results)) {
#     lineages.temp <- as.data.frame(t(setDT(tstrsplit(as.character(results[i, 3]), " ", fixed = TRUE))[]))
#     sample.temp <- rep(results[i, 1], nrow(lineages.temp))
#     if (i == 1) sublineages.final <- cbind(sample.temp, lineages.temp)
#     else sublineages.final <- rbind(sublineages.final, cbind(sample.temp, lineages.temp))
#   }
#   names(sublineages.final) <- c("Sample", "sublineage")
#   sublineages.final |> mutate(sublineage = str_remove(as.character(sublineage), "-like[0-9]*$"))
# }
#
# sublineages_observed <- bind_rows(
#   parse_freyja_aggregate(paste0(storage_dir, "../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv")),
#   parse_freyja_aggregate(paste0(storage_dir, "../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv"))
# ) |>
#   distinct(sublineage) |>
#   filter(!is.na(sublineage), sublineage != "")


# Bucket assignment (matches RAs.R convention) ----------------------------
# DISABLED -- uncomment to restore real per-bucket diversity counting

# parent_group <- function(x, extended_groups = backbone_lineages) {
#   default <- str_extract(x, "^[A-Za-z]+\\.[0-9]+")  # e.g. JN.1.4.3 -> JN.1
#   default <- ifelse(is.na(default), x, default)
#   group <- default
#   for (g in extended_groups) {
#     is_g <- !is.na(x) & (x == g | startsWith(x, paste0(g, ".")))
#     group[is_g] <- g
#   }
#   group
# }
#
# apply_manual_collapse <- function(labels, collapse_map) {
#   for (dest in names(collapse_map)) labels[labels %in% collapse_map[[dest]]] <- dest
#   labels
# }
#
# lineage_buckets <- sublineages_observed |>
#   mutate(
#     group = parent_group(sublineage),
#     bucket = paste0(group, ".x"),
#     bucket = if_else(group %in% backbone_lineages, group, bucket),
#     bucket = apply_manual_collapse(bucket, manual_collapse_map)
#   ) |>
#   filter(bucket %in% major_lineages)


# Barcode distances --------------------------------------------------------

barcodes <- fread(barcode_path, header = TRUE)
setnames(barcodes, 1, "lineage")
barcodes <- barcodes[!duplicated(lineage)]  # file carries duplicate rows for some aliases

get_prototypes <- function(bucket) {
  sources <- if (bucket %in% names(manual_collapse_map)) manual_collapse_map[[bucket]] else bucket
  str_remove(sources, "\\.x$")
}
bucket_prototypes <- setNames(lapply(major_lineages, get_prototypes), major_lineages)
all_prototypes <- unique(unlist(bucket_prototypes))
bucket_anchor <- sapply(bucket_prototypes, `[`, 1)  # single representative prototype per bucket, for tree placement

missing_prototypes <- setdiff(all_prototypes, barcodes$lineage)
if (length(missing_prototypes) > 0) {
  stop("Prototype(s) missing from barcode file: ", paste(missing_prototypes, collapse = ", "))
}

candidate_lineages <- intersect(unique(c(root_lineage, all_prototypes)), barcodes$lineage)
cand_sub <- barcodes[lineage %in% candidate_lineages]
cand_mat <- as.matrix(cand_sub[, -"lineage"])
rownames(cand_mat) <- cand_sub$lineage

# cand_dist <- as.matrix(dist(cand_mat, method = "manhattan"))  # disabled: only needed by the nearest-prototype filter below


# Nearest-prototype filter: assign each observed lineage to its closest bucket
# DISABLED -- uncomment to restore real per-bucket diversity counting

# nearest_is_own <- function(fine_lineage, own_bucket) {
#   own_protos <- intersect(bucket_prototypes[[own_bucket]], rownames(cand_dist))
#   other_protos <- intersect(setdiff(all_prototypes, own_protos), rownames(cand_dist))
#   if (length(own_protos) == 0 || !(fine_lineage %in% rownames(cand_dist))) return(NA)
#   min(cand_dist[fine_lineage, own_protos]) <= min(cand_dist[fine_lineage, other_protos])
# }
#
# lineage_buckets <- lineage_buckets |>
#   filter(sublineage %in% rownames(cand_dist)) |>
#   rowwise() |>
#   mutate(nearest_own = nearest_is_own(sublineage, bucket)) |>
#   ungroup() |>
#   filter(nearest_own, !sublineage %in% manual_exclude_lineages) |>
#   select(sublineage, bucket)
#
# bucket_counts <- lineage_buckets |> count(bucket, name = "n_tips")
# cat("Observed lineage counts per bucket:\n")
# print(bucket_counts)


# Backbone tree: root + backbone lineages + one anchor prototype per bucket

tree_lineages <- intersect(unique(c(root_lineage, bucket_anchor)), rownames(cand_mat))
tree_mat <- cand_mat[tree_lineages, , drop = FALSE]
tree <- nj(dist(tree_mat, method = "manhattan"))
tree <- root(tree, outgroup = root_lineage, resolve.root = TRUE)


# Plot ----------------------------------------------------------------------

p <- ggtree(tree, layout = "rectangular")
tree_data <- p$data

nesting_depth <- function(bucket) {
  depth <- 0
  current <- bucket
  while (current %in% names(nested_bucket_parent)) {
    current <- nested_bucket_parent[[current]]
    depth <- depth + 1
  }
  depth
}

# scale_capped <- function(n, cap, out_range) {  # disabled: was for count-based sizing
#   size <- sqrt(pmin(pmax(n, 1), cap))
#   size_norm <- pmin(pmax((size - 1) / (sqrt(cap) - 1), 0), 1)
#   out_range[1] + size_norm * (out_range[2] - out_range[1])
# }

bucket_geom <- tibble(bucket = major_lineages) |>
  # left_join(bucket_counts, by = "bucket") |>  # disabled: real counts
  mutate(
    depth = sapply(bucket, nesting_depth),
    anchor = bucket_anchor[bucket],
    is_small = bucket %in% small_buckets,
    x_depth = if_else(is_small, triangle_depth / 2, triangle_depth),
    y_half  = if_else(is_small, triangle_width / 2, triangle_width),
    alpha   = if_else(is_small, 0.4, 1)
  ) |>
  left_join(tree_data |> filter(isTip) |> select(label, x, y), by = c("anchor" = "label")) |>
  arrange(depth)

bucket_triangles <- list()
bucket_labels <- list()
for (i in seq_len(nrow(bucket_geom))) {
  row <- bucket_geom[i, ]
  x1 <- row$x + row$x_depth
  bucket_triangles[[row$bucket]] <- tibble(
    bucket = row$bucket, depth = row$depth, alpha = row$alpha,
    x = c(row$x, x1, x1), y = c(row$y, row$y - row$y_half, row$y + row$y_half)
  )
  bucket_labels[[row$bucket]] <- tibble(
    bucket = row$bucket, depth = row$depth, x = x1, y = row$y
  )
}

# draw parents (shallow depth) first so nested children stay visible on top;
# `bucket` must be an explicit factor since geom_polygon's `group` splits by
# factor level (alphabetical for a plain character), not row order
triangle_df <- bind_rows(bucket_triangles) |> arrange(depth)
triangle_df <- triangle_df |> mutate(bucket = factor(bucket, levels = unique(bucket)))
label_df <- bind_rows(bucket_labels) |> arrange(depth)

p <- p +
  geom_polygon(
    data = triangle_df, aes(x = x, y = y, group = bucket, fill = bucket, alpha = alpha),
    color = "grey20", linewidth = 0.3, show.legend = FALSE
  ) +
  scale_fill_manual(values = bucket_colors) +
  scale_alpha_identity() +
  geom_text_repel(
    data = label_df, aes(x = x, y = y, label = bucket),
    size = 3.2, fontface = "bold", direction = "y", hjust = 0,
    nudge_x = 2, segment.size = 0.3, min.segment.length = 0, box.padding = 0.3
  ) +
  geom_tiplab(
    data = td_filter(isTip & label == root_lineage),
    aes(label = label), size = 3, fontface = "bold", hjust = -0.1
  ) +
  theme_tree2() +
  labs(
    title = "Illustrative lineage tree (rooted at BA.2)",
    x = "Barcode mutation distance"
  ) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 90, 5.5, 5.5))

ggsave("../draft_figures/lineage_tree.jpg", p, width = 10, height = 7, dpi = 300)
