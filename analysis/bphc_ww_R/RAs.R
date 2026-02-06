# Visualizing relative abundances

# Libraries ---------------------------------------------------------------

library(data.table)
library(tidyverse)
library(lubridate)
library(runner)
library(viridis)


# Load data ---------------------------------------------------------------

metadata <- read_rds("../data/meta_clean.rds")

clin_lin <- read_rds("../data/clin_lin.rds")


# Parse demix -------------------------------------------------------------


### Code credit for parsing the demix output: https://github.com/a-roguet

results<-read.table("../../freyja_aggregate/aggregated.tsv", fill = TRUE, sep = "\t", h=T)
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

rm(abundances.temp, lineages.temp, results, summarized, end, i, sample.temp, start)

### files created: summarized.final, sublineages.final

### clean up sample names
sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)

### Freyja collapse lineage option results in classifications like "BA.4.6-like" or "BA.5.2-like2"
  ### remove "-likeN" suffix

sublineages.final <- sublineages.final %>%
  mutate(sublineage = str_remove(sublineage, "-like[0-9]*$"))


# Merge with metadata -----------------------------------------------------


sublin_meta <- left_join(sublineages.final, metadata, by=c("Sample"="FASTQ_ID"))


# Demix progress report ---------------------------------------------------

# exclude <- metadata |>
#   filter(LOCATION=="Charlestown" & year_epiweek=="202334" | 
#            LOCATION=="Lower_Roxbury" & year_epiweek=="202410") |> 
#   select(FASTQ_ID)
# 
# p_demix_progress <- metadata |> 
#   filter(!FASTQ_ID %in% exclude$FASTQ_ID) |> 
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

### fx - only keep parent groups that are greater than RA threshold in any sample
  ### this fx used within collapse_sublineages_timeaware()

keep_threshold_time <- function(dat, threshold = 0.1, 
                                date_col = SAMPLING_DATE,
                                suffix = ".x") {
  dat |> 
    mutate(group = parent_group(sublineage), # create parent lineage
           abundance = as.numeric(abundance),
           date = as.Date({{ date_col }})
           ) |>
    group_by(Sample, date, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |> # sum RAs for collapsed lineages
    group_by(date, group) |> 
    summarise(max_ra_any_sample = max(ra, na.rm = TRUE), .groups = "drop") |> 
    arrange(group, date) |> 
    group_by(group) |> 
    mutate(keep_by_date = cummax(max_ra_any_sample >= threshold)) |> 
    ungroup() |> 
    transmute(date, group, keep_by_date,
              label = paste0(group, suffix))
}

  # view keepers
  # tholds.t <- keep_threshold_time(sublin_meta)

### fx - collapse using time-aware threshold
  # if using combined dataset, keep_tbl from combined data can be supplied 
  # if using single dataset, keep_tbl can be computed internally

collapse_sublineages_timeaware <- function(df, threshold = 0.1,
                                           other = "other",
                                           date_col = SAMPLING_DATE,
                                           suffix = ".x",
                                           keep_tbl = NULL) {
  
  # if no keep_tbl provided, generate one
  if (is.null(keep_tbl)) {
    keep_tbl <- keep_threshold_time(df, threshold = threshold,
                                    date_col = {{ date_col }}, suffix = suffix)
  }
  
  df |> 
    as_tibble() |> 
    mutate(
      sublineage = coalesce(sublineage, ""),
      group = parent_group(sublineage),
      abundance = as.numeric(abundance),
      date = as.Date({{ date_col }})
    ) |> 
    left_join(keep_tbl |> select(date, group, keep_by_date, label),
              by = c("date", "group")) |> 
    mutate(
      sublin_collapse = if_else(keep_by_date==TRUE, label, other),
      sublin_collapse = factor(sublin_collapse)
    ) |> 
    group_by(Sample, SAMPLING_DATE, LOCATION, epiweek, year, year_epiweek, sublin_collapse) |> 
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
}


# Weekly Lineage Composition x Neighborhood -------------------------------


collapse_meta <- collapse_sublineages_timeaware(
  sublin_meta,
  threshold = 0.2,
  date_col = SAMPLING_DATE,
  suffix = ".x",
  other = "other"
)

### identify set of weeks for plotting
year_epiweek_levels <- collapse_meta |> 
  mutate(year_epiweek = as.character(year_epiweek)) |> 
  distinct(year_epiweek) |> 
  arrange(as.integer(year_epiweek)) |> 
  pull(year_epiweek)

### create a "complete" dataset such that each location has a placeholder for all samples and lineage
collapse_meta_complete <- collapse_meta |> 
  mutate(year_epiweek = factor(as.character(year_epiweek), levels = year_epiweek_levels)) |> 
  group_by(LOCATION) |> 
  complete(year_epiweek = year_epiweek_levels, sublin_collapse, fill = list(abundance = 0)) |> 
  ungroup()

### plot the "complete" data x neighborhood
p_RAxNB <- collapse_meta_complete |> 
  ggplot(aes(x = year_epiweek, y = abundance, fill = sublin_collapse)) +
  geom_col() +
  facet_wrap(~ LOCATION) +
  scale_x_discrete(drop = FALSE) +  # keep empty weeks
  labs(
    x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition x Neighborhood"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(option="turbo", discrete = T) +
  geom_text(aes(label = ifelse(abundance > 0.05, 
                               as.character(sublin_collapse), "")),
            position = position_stack(vjust = 0.5),
            size = 2.25, color = "white")
p_RAxNB

ggsave("../figures/p_RAxNB.jpg", p_RAxNB, scale = 1.25) # 2025-02-05


# Incorporate clinical data -----------------------------------------------

### calculate weekly RAs in clinical data

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
         LOCATION = NA,
         Sample=NA,
         epiweek=NA)

### filter to relevant time

clinical <- clinical |> filter(year_epiweek %in% sublin_meta$year_epiweek)

### combine datasets

clin_ww_comb <- bind_rows(
  sublin_meta |> mutate(source = "wastewater", abundance=as.numeric(abundance)),
  clinical |> mutate(source = "clinical")
)

### create keep table from combined dataset

keep_tbl_comb <- keep_threshold_time(clin_ww_comb, threshold = 0.1)

### apply combined keep table to ww and clinical data individually 

collapse_ww <- collapse_sublineages_timeaware(sublin_meta, keep_tbl = keep_tbl_comb)
collapse_clin <- collapse_sublineages_timeaware(clinical, keep_tbl = keep_tbl_comb)

### this again

### identify set of weeks for plotting
year_epiweek_levels <- collapse_clin |> 
  mutate(year_epiweek = as.character(year_epiweek)) |> 
  distinct(year_epiweek) |> 
  arrange(as.integer(year_epiweek)) |> 
  pull(year_epiweek)

### create a "complete" dataset such that each location has a placeholder for all samples and lineage
collapse_clin_complete <- collapse_clin |> 
  mutate(year_epiweek = factor(as.character(year_epiweek), levels = year_epiweek_levels)) |> 
  group_by(LOCATION) |> 
  complete(year_epiweek = year_epiweek_levels, sublin_collapse, fill = list(abundance = 0)) |> 
  ungroup()

### plot the "complete" data x neighborhood
collapse_clin_complete |> 
  ggplot(aes(x = year_epiweek, y = abundance, fill = sublin_collapse)) +
  geom_col() +
  facet_wrap(~ LOCATION) +
  scale_x_discrete(drop = FALSE) +  # keep empty weeks
  labs(
    x = "year_epiweek", y = "Relative Abundance", fill = "sublineage",
    title = "SARS-CoV-2 Weekly Lineage Composition - Clinical"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(option="turbo", discrete = T) +
  geom_text(aes(label = ifelse(abundance > 0.05, 
                               as.character(sublin_collapse), "")),
            position = position_stack(vjust = 0.5),
            size = 2.25, color = "white")








