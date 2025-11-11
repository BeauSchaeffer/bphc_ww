# Visualizing relative abundances

# Libraries ---------------------------------------------------------------

library(data.table)
library(tidyverse)
library(lubridate)
library(runner)
library(viridis)


# Load data ---------------------------------------------------------------

metadata <- read_rds("../data/meta_clean.rds")


# Parse demix -------------------------------------------------------------

### Code credit for parsing the demix output: https://github.com/a-roguet

results<-read.table("../../freyja_aggregate/aggregated.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

# Summarized data
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


# Sublineages data
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


# Merge with metadata -----------------------------------------------------

sublin_meta <- left_join(sublineages.final, metadata, by=c("Sample"="FASTQ_ID"))


# Demix progress report ---------------------------------------------------

p_demix_progress <- metadata |> 
  filter(!FASTQ_ID %in% exclude$FASTQ_ID) |> 
  mutate(freyja = case_when(
    FASTQ_ID %in% sublineages.final$Sample ~ 1,
    !FASTQ_ID %in% sublineages.final$Sample ~ 0
  )) |> 
  arrange(year_epiweek, LOCATION) |>
  mutate(year_epiweek = factor(year_epiweek)) |> 
  ggplot(aes(x=year_epiweek, y=LOCATION, fill=as.factor(freyja))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#4C9AED", "#654CED"), name="demix")
p_demix_progress

ggsave("../figures/p_demix_progress.jpg", p_demix_progress, scale = 1.5) # 2025-11-10

# Collapsing sublineages --------------------------------------------------

### fx - extract parent group

parent_group <- function(x) {
  out <- str_extract(x, "^[A-Za-z]+\\.[0-9]+")
  ifelse(is.na(out), x, out)  # keep full token if no match (recombinants)
}

### fx - only keep parent groups that are greater than RA threshold in any sample

keep_threshold <- function(dat, threshold = 0.05) {
  dat |> 
    mutate(group = parent_group(sublineage), # create parent lineage
           abundance = as.numeric(abundance)) |>
    group_by(Sample, group) |>
    summarise(ra = sum(abundance, na.rm = TRUE), .groups = "drop") |> # sum RAs for collapsed lineages
    filter(ra >= threshold) |>
    distinct(group) |>
    pull(group)
}

### fx - format dots in parent group names

escape_dots <- function(x) gsub("\\.", "\\\\.", x)

### fx - create "rules" for collapsing sublineages

rules_from_keepers <- function(keepers, suffix = ".x") {
  tibble(
    pattern = paste0("^", escape_dots(keepers), "(?:\\.|$)"),
    label   = paste0(keepers, suffix)
  )
}

### fx - collapse sublineages based on rules above

collapse_sublineages <- function(df, rules, other = "other", ignore_case = FALSE) {
  out <- as_tibble(df)
  out$sublineage <- coalesce(out$sublineage, "") # replace NA with empty string
  out$sublin_collapse <- NA_character_ # initialize new col with NAs
  
  for (i in seq_len(nrow(rules))) {
    rx <- rules$pattern[i] # identify current regex pattern
    idx <- is.na(out$sublin_collapse) & str_detect(out$sublineage, rx) # identifies which rows to fill for a given rule
    out$sublin_collapse[idx] <- rules$label[i] # assign corresponding label
  }
  
  out$sublin_collapse[is.na(out$sublin_collapse)] <- other # if still NA, assign other
  
  out$sublin_collapse <- factor(out$sublin_collapse,
                                levels = c(unique(rules$label), other)) # make factor
  
  # sum RA values over collapsed lineages
  out <- out |> 
    group_by(Sample, SAMPLING_DATE, LOCATION, epiweek, year, year_epiweek, sublin_collapse) |> 
    summarise(abundance = sum(as.numeric(abundance), na.rm = TRUE), .groups = "drop")
  
  out
}

keepers <- keep_threshold(sublin_meta, threshold = 0.2)
rules   <- rules_from_keepers(keepers)
collapse_meta <- collapse_sublineages(sublin_meta, rules)

# Weekly Lineage Composition x Neighborhood -------------------------------

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
  scale_fill_viridis(option="turbo", discrete = T)
p_RAxNB

ggsave("../figures/p_RAxNB.jpg", p_RAxNB) # 2025-11-11




