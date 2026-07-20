# Visualizing relative abundances


# Libraries ---------------------------------------------------------------


library(data.table)
library(tidyverse)
library(lubridate)
library(runner)
library(viridis)


# Load data ---------------------------------------------------------------

storage_dir <- "../data/"

ww_metadata <- read_rds(paste0(storage_dir, "meta_clean.rds"))

clin_lin <- read_rds(paste0(storage_dir, "clin_lin.rds"))

pop_wts <- read_rds(paste0(storage_dir, "pop_wts.rds"))

conc_wts <- read_rds(paste0(storage_dir, "conc_wts.rds"))


# Parse demix -------------------------------------------------------------


### Code credit for parsing the demix output: https://github.com/a-roguet

results<-read.table("/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/bphc_ww/freyja_aggregate/aggregated.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

### Summarized data
summarized<-as.data.frame(setDT(tstrsplit(as.character(results$summarized), " ", fixed=TRUE))[]) # Extract Summarized data
summarized$sample<-results$X

for(i in 1:((ncol(summarized)-1)/2)){
  if(i==1){
    summarized.final<-summarized[,c(ncol(summarized),1:2)]
  } else {
    start=i*2-1; end=i*2
    summarized.final<-rbind(summarized.final, setNames(summarized[,c(ncol(summarized), start:end)], names(summarized.final)))
  }
}
summarized.final<-summarized.final[complete.cases(summarized.final), ]
names(summarized.final)<-c("Sample", "lineage", "abundance")


### Sublineages data
for(i in 1:nrow(results)){
  lineages.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 3]), " ", fixed=TRUE))[]))
  abundances.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 4]), " ", fixed=TRUE))[]))
  sample.temp<-rep(results[i, 1], nrow(lineages.temp))
  if(i==1){
    sublineages.final<-cbind(sample.temp, lineages.temp, abundances.temp)
  } else {
    sublineages.final<-rbind(sublineages.final, cbind(sample.temp, lineages.temp, abundances.temp))
  }
}
names(sublineages.final)<-c("Sample", "sublineage", "abundance")

### Code credit for parsing the demix output: https://github.com/a-roguet

### clean up parsed output

rm(abundances.temp, lineages.temp, results, summarized, summarized.final, end, i, sample.temp, start)

### files created: summarized.final, sublineages.final

### clean up sample names
sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)

### Freyja collapse lineage option results in classifications like "BA.4.6-like" or "BA.5.2-like2"
### remove "-likeN" suffix

sublineages.final <- sublineages.final %>%
  mutate(sublineage = str_remove(as.character(sublineage), "-like[0-9]*$"))


# Parse demix as function for multi-batch ---------------------------------

parse_freyja_aggregate <- function(path){

  results<-read.table(path, fill = TRUE, sep = "\t", h=T)
  results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
  results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

  ### Sublineages data
  for(i in 1:nrow(results)){
    lineages.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 3]), " ", fixed=TRUE))[]))
    abundances.temp<-as.data.frame(t(setDT(tstrsplit(as.character(results[i, 4]), " ", fixed=TRUE))[]))
    sample.temp<-rep(results[i, 1], nrow(lineages.temp))
    if(i==1){
      sublineages.final<-cbind(sample.temp, lineages.temp, abundances.temp)
    } else {
      sublineages.final<-rbind(sublineages.final, cbind(sample.temp, lineages.temp, abundances.temp))
    }
  }
  names(sublineages.final)<-c("Sample", "sublineage", "abundance")

  ### Code credit for parsing the demix output: https://github.com/a-roguet

  ### clean up sample names
  sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)

  ### Freyja collapse lineage option results in classifications like "BA.4.6-like" or "BA.5.2-like2"
  ### remove "-likeN" suffix

  sublineages.final <- sublineages.final %>%
    mutate(sublineage = str_remove(as.character(sublineage), "-like[0-9]*$"))

  sublineages.final
}

sublineages.final <- bind_rows(
  parse_freyja_aggregate("../baseload_batch1_outputs/freyja_aggregate/aggregated.tsv"),
  parse_freyja_aggregate("../baseload_batch2_outputs/freyja_aggregate/aggregated.tsv")
)


# Merge with metadata -----------------------------------------------------


sublin_meta <- left_join(sublineages.final, ww_metadata, by=c("Sample"="FASTQ_ID")) |> 
  mutate(abundance=as.numeric(abundance))

# sublin_meta$spike_frac_cov_10x |> hist()
# sublin_meta$frac_genome_cov_10x |> hist()

### Apply QC filters at this stage

# sublin_meta <- sublin_meta |> 
#   filter(spike_frac_cov_10x >= 0.75)


qc_flags_nb <- ww_metadata |>
  mutate(time = as.character(year_epiweek),
         low_spike_qc = spike_frac_cov_10x <= 0.5 | is.na(spike_frac_cov_10x)) |>
  group_by(LOCATION, time) |>
  summarise(low_spike_qc = any(low_spike_qc), .groups = "drop")

qc_flags_nb <- qc_flags_nb |> filter(time %in% sublin_meta$year_epiweek)


# Demix progress report ---------------------------------------------------


# p_demix_progress <- ww_metadata |>
#   mutate(freyja = case_when(
#     FASTQ_ID %in% sublineages.final$Sample ~ 1,
#     !FASTQ_ID %in% sublineages.final$Sample ~ 0
#   )) |>
#   arrange(year_epiweek, LOCATION) |>
#   mutate(year_epiweek = factor(year_epiweek)) |>
#   ggplot(aes(x=year_epiweek, y=LOCATION, fill=as.factor(freyja))) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values = c("#4C9AED", "#654CED"), name="demix")
# p_demix_progress
#
# ggsave("../figures/p_demix_progress.jpg", p_demix_progress, scale = 1.5) # 2025-11-10


# Collapsing sublineages --------------------------------------------------


### fx - extract parent group
### this fx used within keep_threshold_time()

parent_group <- function(x) {
  out <- str_extract(x, "^[A-Za-z]+\\.[0-9]+")
  ifelse(is.na(out), x, out)  # keep full token if no match (recombinants)
}

### fx - only keep parent groups that exceed a threshold
### used within collapse_sublineages_timeaware()

keep_threshold_time <- function(dat, threshold = 0.1, 
                                time_col = year_epiweek,
                                suffix = ".x") {
  dat |> 
    mutate(group = parent_group(sublineage), # create parent lineage
           # format variables
           abundance = as.numeric(abundance),
           time = {{ time_col }}
    ) |>
    # sum RA for collapsed lineages within each sample,week
    group_by(Sample, time, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |>
    # for each week,lineage group, take max RA observed in any sample that week
    group_by(time, group) |> 
    summarise(max_ra_any_sample = max(ra, na.rm = TRUE), .groups = "drop") |> 
    arrange(group, time) |> 
    group_by(group) |> 
    # once lineages reach threshold RA anywhere, flag with keep_by_time==TRUE moving forward
    mutate(keep_by_time = cummax(max_ra_any_sample >= threshold)) |>
    # once a lineage reaches the threshold it is named moving forward, rather than "other"
    ungroup() |> 
    transmute(time, group, keep_by_time,
              label = paste0(group, suffix))
  
  # output is a "keep table" with a row for every group/time combo
  # keep_by_time==TRUE starting the first year_epiweek the group reaches threshold,
  # and for all later weeks where the group is observed
  
}

# view keepers
# tholds.t <- keep_threshold_time(sublin_meta) |> arrange(group,time)


### fx - collapse using time-aware threshold
# use the keep table keyed by (time, group), and collapse to "other"
# until keep_by_time becomes TRUE.

# keep_tbl from combined data can be supplied 
# or
# keep_tbl can be computed internally 

collapse_sublineages_timeaware <- function(df, threshold = 0.1,
                                           other = "other",
                                           time_col = year_epiweek,
                                           suffix = ".x",
                                           keep_tbl = NULL) {
  
  # if no keep_tbl provided, generate one
  if (is.null(keep_tbl)) {
    keep_tbl <- keep_threshold_time(df, threshold = threshold,
                                    time_col = {{ time_col }}, suffix = suffix)
  }
  
  df |>
    as_tibble() |> 
    mutate(
      sublineage = coalesce(sublineage, ""),
      group = parent_group(sublineage),
      abundance = as.numeric(abundance),
      time = {{ time_col }}
    ) |> 
    
    # assign keep flag by epiweek + parent group
    left_join(keep_tbl |> select(time, group, keep_by_time, label),
              by = c("time", "group")) |> 
    mutate(
      # safety check
      keep_by_time = as.logical(coalesce(keep_by_time, FALSE)),
      # if group has not yet reached threshold, grouped into other
      sublin_collapse = if_else(keep_by_time, label, other),
      sublin_collapse = factor(sublin_collapse)
    ) |> 
    # sum RAs after collapsing (within sample,week)
    group_by(Sample, time, LOCATION, epiweek, year, sublin_collapse) |> 
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
}



# Calculate RAs for clinical ----------------------------------------------

clinical <- clin_lin |> 
  mutate(parentsub = parent_group(sublineage)) |> 
  group_by(year_epiweek, parentsub) |> 
  summarise(count = n())  |> 
  ungroup() |> 
  group_by(year_epiweek) |> 
  mutate(abundance=count/sum(count)) |> 
  ungroup() |> 
  rename(sublineage=parentsub) |> 
  mutate(year = year_epiweek %/% 100,
         week = year_epiweek %% 100,
         year_first = ymd(paste0(year, "-01-01")) + weeks(week - 1),
         SAMPLING_DATE = floor_date(year_first, "week", week_start = 7),
         
         # placeholder columns so collapse_sublineages_timeaware() functions properly
         LOCATION = NA, 
         Sample=NA,
         epiweek=week)

### remove later - filter to relevant time that samples have been processed
clinical <- clinical |> filter(year_epiweek %in% sublin_meta$year_epiweek)

### pseudo variables

clinical <- clinical |>
  mutate(
    Sample = "clinical_agg",
    LOCATION = "clinical"
  )



# Run functions -----------------------------------------------------------

combined_for_keep <- bind_rows(
  sublin_meta |> select(Sample, sublineage, abundance, year_epiweek),
  clinical   |> select(Sample, sublineage, abundance, year_epiweek)
)


keep_tbl <- keep_threshold_time(combined_for_keep,
                                threshold = 0.2,
                                time_col = year_epiweek,
                                suffix = ".x")

ww_collapsed <- collapse_sublineages_timeaware(sublin_meta,
                                               threshold = 0.2,
                                               time_col = year_epiweek,
                                               suffix = ".x",
                                               keep_tbl = keep_tbl)

clin_collapsed <- collapse_sublineages_timeaware(clinical,
                                                 threshold = 0.2,
                                                 time_col = year_epiweek,
                                                 suffix = ".x",
                                                 keep_tbl = keep_tbl)

# write_rds(ww_collapsed, "../data/ww_collapsed.R")
# write_rds(clin_collapsed, "../data/clin_collapsed.R")

# ww_collapsed_complete <- ww_collapsed |>
#   mutate(time = as.character(time)) |>
#   complete(
#     time = year_epiweek_levels,
#     LOCATION,
#     sublin_collapse,
#     fill = list(abundance = 0)
#   ) |>
#   mutate(time = factor(time, levels = year_epiweek_levels),
#          year_epiweek = as.numeric(as.character(time))) |> 
#   left_join(conc_wts, by = c("LOCATION", "year_epiweek"))
# write_rds(ww_collapsed_complete, "../data/ww_collapsed_complete.rds")
# rm(ww_collapsed_complete)

# stable lineage order for all plots
lin_levels <- sort(unique(c(
  as.character(ww_collapsed$sublin_collapse),
  as.character(clin_collapsed$sublin_collapse)
)))
# force "other" to last
lin_levels <- c(setdiff(lin_levels, "other"), "other")

# fixed color mapping for all lineage levels (stable across plots)
fill_cols <- setNames(viridis::viridis(length(lin_levels), option = "turbo"),
                      lin_levels)


# apply ordering
ww_collapsed  <- ww_collapsed  |> mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels))
clin_collapsed <- clin_collapsed |> mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels))


# Plot clinical -----------------------------------------------------------

### identify set of weeks for plotting
year_epiweek_levels <- sort(unique(c(
  ww_collapsed$time,
  clin_collapsed$time
))) |> as.character()


### create a "complete" dataset such that each location has a placeholder for all samples and lineage
clin_collapsed_complete <- clin_collapsed |>
  mutate(time = as.character(time)) |>
  complete(
    time = year_epiweek_levels,
    sublin_collapse,
    fill = list(abundance = 0)
  ) |>
  mutate(time = factor(time, levels = year_epiweek_levels))

clin_collapsed_complete <- clin_collapsed_complete |>
  mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels))


### plot the "complete" data for clinical
p_RAxClin <- clin_collapsed_complete |>
  ggplot(aes(x = time, y = abundance, fill = sublin_collapse)) +
  geom_col() +
  scale_x_discrete(drop = FALSE) +  # keep empty weeks
  labs(
    x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition - Clinical"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  geom_text(
    aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
    position = position_stack(vjust = 0.5),
    size = 2.25, color = "white", angle=90
  ) +
  guides(fill="none")

p_RAxClin
# ggsave("../draft_figures/p_RAxClin.jpg", p_RAxClin, width = 12, height = 6, dpi = 300)
# ggsave("../draft_figures/RAxClin.jpg", p_RAxClin, width = 12, height = 6, dpi = 300)



# Plot WW -----------------------------------------------------------------

ww_collapsed_complete <- ww_collapsed |>
  mutate(time = as.character(time)) |>
  complete(
    time = year_epiweek_levels,
    LOCATION,
    sublin_collapse,
    fill = list(abundance = 0)
  ) |>
  mutate(time = factor(time, levels = year_epiweek_levels))


ww_collapsed_complete <- ww_collapsed_complete |>
  mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels))

### plot the "complete" data for ww

  ### adding pcr concentration data

conc_label <- conc_wts |>
  mutate(time = as.character(year_epiweek)) |>
  filter(LOCATION != "Lower Roxbury",
         time %in% year_epiweek_levels) |>
  select(LOCATION, time, logeff) |>
  mutate(
    conc_quartile = cut(logeff, 
                        breaks = quantile(logeff, probs = 0:4/4, na.rm = TRUE),
                        labels = c("Q1", "Q2", "Q3", "Q4"),
                        include.lowest = TRUE))


p_RAxNB <- ww_collapsed_complete |>
  ggplot(aes(x = time, y = abundance, fill = sublin_collapse)) +
  geom_col() +
  geom_col(
    data = qc_flags_nb |> filter(low_spike_qc),
    aes(x = time, y = 1),
    fill = "grey40", alpha = 1,
    inherit.aes = FALSE
  ) +
  # geom_text(
  #   data = conc_label |> mutate(time = factor(time, levels = year_epiweek_levels)),
  #   aes(x = time, y = 1.00, label = round(logeff, 1)),
  #   inherit.aes = FALSE,
  #   size = 1.5, angle = 90, hjust = 0
  # ) +
  geom_point(
    data = conc_label |> mutate(time = factor(time, levels = year_epiweek_levels)),
    aes(x = time, y = 1.03, color = conc_quartile),
    inherit.aes = FALSE,
    size = 1.2
  ) +
  scale_color_manual(
    values = c("Q1" = "#2C7BB6", "Q2" = "#ABD9E9", "Q3" = "#FDAE61", "Q4" = "#D7191C"),
    name = "PCR conc\nquartile"
  ) +
  facet_wrap(~ LOCATION, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +  # keep empty weeks
  labs(
    x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition - Wastewater",
    subtitle = "2026-06-13 batch 1+2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=4),
    panel.spacing = unit(0.75, "lines"),
    legend.position = "right"
  ) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  geom_text(
    aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
    position = position_stack(vjust = 0.5),
    size = 1.5, color = "white", angle=90
  ) +
  guides(
    fill = "none",
    color = guide_legend(title = "PCR conc\nquartile")
  )

p_RAxNB
# ggsave("../draft_figures/p_RAxNB.jpg", p_RAxNB, width = 22, height = 14, dpi = 300)

### separate figures rather than multipanel

locations <- unique(ww_collapsed_complete$LOCATION)

plots_by_location <- setNames(
  lapply(locations, function(loc) {
    ww_collapsed_complete |>
      filter(LOCATION == loc) |>
      ggplot(aes(x = time, y = abundance, fill = sublin_collapse)) +
      geom_col() +
      geom_col(
        data = qc_flags_nb |> filter(low_spike_qc, LOCATION == loc),
        aes(x = time, y = 1),
        fill = "grey40", alpha = 1,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = conc_label |>
          filter(LOCATION == loc) |>
          mutate(time = factor(time, levels = year_epiweek_levels)),
        aes(x = time, y = 1.03, color = conc_quartile),
        inherit.aes = FALSE,
        size = 1.2
      ) +
      scale_color_manual(
        values = c("Q1" = "#2C7BB6", "Q2" = "#ABD9E9", "Q3" = "#FDAE61", "Q4" = "#D7191C"),
        name = "PCR conc\nquartile"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_fill_manual(values = fill_cols, drop = FALSE) +
      geom_text(
        aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
        position = position_stack(vjust = 0.5),
        size = 2, color = "white", angle = 90
      ) +
      labs(
        x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
        title = paste("SARS-CoV-2 Weekly Lineage Composition -", loc),
        subtitle = "2026-06-13 batch 1+2"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "right"
      ) +
      guides(fill = "none", color = guide_legend(title = "PCR conc\nquartile"))
  }),
  locations
)

plots_by_location[[locations[12]]]


lapply(locations, function(loc) {
  ggsave(
    paste0("../draft_figures/RA_plots/RAxNB_", gsub(" ", "_", loc), ".jpg"),
    plots_by_location[[loc]],
    width = 14, height = 6, dpi = 300
  )
})


### concentration weighted

plots_by_location_concwt <- setNames(
  lapply(locations, function(loc) {
    ww_collapsed_complete |>
      mutate(year_epiweek = paste0(year, sprintf("%02d", epiweek)),
             year_epiweek = ifelse(year_epiweek=="NANA",NA,year_epiweek),
             year_epiweek=as.numeric(year_epiweek)) |>
      left_join(conc_wts, by=c("LOCATION", "year_epiweek")) |>
      mutate(abundance_LEwtd = abundance*logeff,
             abundance_Ewtd = abundance*meaneffCopiesL) |> 
      filter(LOCATION == loc) |>
      ggplot(aes(x = time, y = abundance_LEwtd, fill = sublin_collapse)) +
      geom_col() +
      geom_col(
        data = qc_flags_nb |> filter(low_spike_qc, LOCATION == loc),
        aes(x = time, y = Inf),
        fill = "grey40", alpha = 1,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = conc_label |>
          filter(LOCATION == loc) |>
          mutate(time = factor(time, levels = year_epiweek_levels)),
        aes(x = time, color = conc_quartile),
        y=Inf,
        position = position_nudge(y = -0.5),
        inherit.aes = FALSE,
        size = 1.2
      ) +
      scale_color_manual(
        values = c("Q1" = "#2C7BB6", "Q2" = "#ABD9E9", "Q3" = "#FDAE61", "Q4" = "#D7191C"),
        name = "PCR conc\nquartile"
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_fill_manual(values = fill_cols, drop = FALSE) +
      geom_text(
        aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
        position = position_stack(vjust = 0.5),
        size = 2, color = "white", angle = 90
      ) +
      labs(
        x = "year_epiweek", y = "Absolute-ish Abundance", fill = "sublineage",
        title = paste("SARS-CoV-2 Weekly Lineage Composition -", loc),
        subtitle = "2026-06-13 batch 1+2"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "right"
      ) +
      guides(fill = "none", color = guide_legend(title = "PCR conc\nquartile"))
  }),
  locations
)


plots_by_location_concwt[[1]]

lapply(locations, function(loc) {
  ggsave(
    paste0("../draft_figures/RA_concwt_plots/RAxNB_", gsub(" ", "_", loc), ".jpg"),
    plots_by_location_concwt[[loc]],
    width = 14, height = 6, dpi = 300
  )
})





### concentration weighted

ww_collapsed_complete |>
  mutate(year_epiweek = paste0(year, sprintf("%02d", epiweek)),
         year_epiweek = ifelse(year_epiweek=="NANA",NA,year_epiweek),
         year_epiweek=as.numeric(year_epiweek)) |>
  left_join(conc_wts, by=c("LOCATION", "year_epiweek")) |>
  mutate(abundance_LEwtd = abundance*logeff,
         abundance_Ewtd = abundance*meaneffCopiesL) |>

  ggplot(aes(x = time, y = abundance_LEwtd, fill = sublin_collapse)) +
  geom_col() +
  facet_wrap(~ LOCATION) +
  scale_x_discrete(drop = FALSE) +  # keep empty weeks
  labs(
    x = "year_epiweek", y = "Absolute-ish Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition - Wastewater"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.75, "lines")
  ) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  geom_text(
    aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
    position = position_stack(vjust = 0.5),
    size = 2.25, color = "white"
  ) +
  guides(
    fill = "none",
    color = guide_legend(title = "PCR conc\nquartile")
  )




# Plot WW citywide --------------------------------------------------------


# check for multiple samples per week
ww_collapsed |>
  count(time, LOCATION, sublin_collapse) |>
  summarise(max_dups = max(n))

ww_locweek <- ww_collapsed |>
  # only needed if multiple samples/week somewhere
  # group_by(time, LOCATION, sublin_collapse) |>
  # summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") |>
  mutate(sublin_collapse = factor(sublin_collapse, levels = lin_levels),
         time = as.character(time))

# for every time,location pair
# create rows for every sublin_collapse in lin_levels
# abundance set to 0 for missing lin_levels
ww_locweek_complete <- ww_locweek |>
  tidyr::complete(
    # nesting - only expand lineages within time-location pairs that already exist in the data.
    tidyr::nesting(time, LOCATION),
    sublin_collapse = factor(lin_levels, levels = lin_levels),
    fill = list(abundance = 0)
  )

# # within-neighborhood, within-week renormalization
# # ensures RAs sum exactly to 1 before weighting
# # probably unnecessary but freyja outputs can sometimes fail to sum to 1
ww_locweek_complete <- ww_locweek_complete |>
  group_by(time, LOCATION) |>
  mutate(
    denom = sum(abundance, na.rm = TRUE),
    abundance = dplyr::if_else(denom > 0, abundance / denom, NA_real_)
  ) |>
  select(-denom) |>
  ungroup()

  # make sure all locations in RA data appear in pop_wts
  # anti_join(ww_locweek_complete, pop_wts, by = "LOCATION") |> distinct(LOCATION)

ww_citywide_weighted <- ww_locweek_complete |>
  left_join(pop_wts |> select(LOCATION, pop), by = "LOCATION") |>
  mutate(pop = dplyr::if_else(is.na(abundance), NA_real_, pop)) |>
  group_by(time, sublin_collapse) |>
  summarise(
    abundance = sum(pop * abundance, na.rm = TRUE) / sum(pop, na.rm = TRUE),
    pop_covered = sum(pop, na.rm = TRUE),
    n_locations = n_distinct(LOCATION),
    .groups = "drop"
  ) |>
  mutate(
    time = as.numeric(time),
    LOCATION = "WW_citywide_weighted",
    sublin_collapse = factor(sublin_collapse, levels = lin_levels)
  )

# checks
ww_citywide_weighted |>
  group_by(time) |>
  summarise(sum_abundance = sum(abundance), .groups = "drop") |>
  summarise(min = min(sum_abundance), max = max(sum_abundance))

# complete data for plotting
ww_citywide_weighted_complete <- ww_citywide_weighted |>
  mutate(time = as.character(time)) |>
  tidyr::complete(
    time = year_epiweek_levels,
    sublin_collapse = factor(lin_levels, levels = lin_levels),
    fill = list(abundance = 0, pop_covered = NA_real_, n_locations = NA_integer_, LOCATION = "WW_citywide_weighted")
  ) |>
  mutate(
    time = factor(time, levels = year_epiweek_levels),
    sublin_collapse = factor(sublin_collapse, levels = lin_levels)
  )


p_RAxCity <- ww_citywide_weighted_complete |>
  ggplot(aes(x = time, y = abundance, fill = sublin_collapse)) +
  geom_col() +
  # geom_col(
  #   data = qc_flags_city |> filter(low_qc),
  #   aes(x = time, y = 1),
  #   fill = "grey40", alpha = 0.4,
  #   inherit.aes = FALSE
  # ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition - Wastewater Citywide",
    subtitle = "population weighted"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = fill_cols, drop = FALSE) +
  geom_text(
    aes(label = ifelse(abundance > 0.05, as.character(sublin_collapse), "")),
    position = position_stack(vjust = 0.5),
    size = 2.25, color = "white", angle=90
  ) +
  guides(fill="none")

p_RAxCity

ggsave("../draft_figures/p_RAxCity.jpg", p_RAxCity, width = 12, height = 6, dpi = 300)



###

p_RAxClin
p_RAxNB
p_RAxCity




