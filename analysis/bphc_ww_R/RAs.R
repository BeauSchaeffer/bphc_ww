# Visualizing relative abundances

# Libraries ---------------------------------------------------------------

library(data.table)
library(tidyverse)
library(lubridate)
library(runner)
library(viridis)


# Load data ---------------------------------------------------------------

metadata <- read_csv("../data/meta_clean.csv")


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

sublineages.final
setDF(sublineages.final)

### clean up sample names
sublineages.final$Sample <- sub("^([0-9]+).*", "\\1", sublineages.final$Sample)


# Merge with metadata -----------------------------------------------------

sublin_meta <- left_join(sublineages.final, metadata, by=c("Sample"="FASTQ_ID"))

### manually collapsing lineages
sublin_meta <- sublin_meta |> 
  as_tibble() |> 
  mutate(abundance = as.numeric(abundance)) |> 
  
  mutate(sublineage_new = case_when(
    
    str_detect(sublineage, "B\\.1\\.1\\.529") ~ "B.1.1.529.x", #
    str_detect(sublineage, "BA\\.1") ~ "BA.1.x", #
    
    str_detect(sublineage, "BA\\.2\\.12\\.1") ~ "BA.2.12.1.x", #
    str_detect(sublineage, "BA\\.2\\.75") ~ "BA.2.75.x", #
    str_detect(sublineage, "CH\\.1\\.1") ~ "CH.1.1.x", #
    str_detect(sublineage, "BA\\.2") ~ "BA.2.x", #
    
    str_detect(sublineage, "BA\\.2\\.86") ~ "BA.2.86.x", # 
    str_detect(sublineage, "JN\\.1\\.11\\.1") ~ "JN.1.11.1.x", # 
    str_detect(sublineage, "JN\\.1") ~ "JN.1.x", # 
    
    str_detect(sublineage, "XBB\\.1\\.5\\.70") ~ "XBB.1.5.70.x", #
    str_detect(sublineage, "XBB\\.1\\.5") ~ "XBB.1.5.x", #
    str_detect(sublineage, "XBB\\.1\\.16") ~ "XBB.1.16.x", #
    str_detect(sublineage, "XBB\\.1\\.9") ~ "XBB.1.9.x", #
    str_detect(sublineage, "XBB\\.2\\.3") ~ "XBB.2.3.x", #
    str_detect(sublineage, "EG\\.5\\.1") ~ "EG.5.1.x", #
    str_detect(sublineage, "HK\\.3") ~ "HK.3.x", #
    str_detect(sublineage, "XBB") ~ "XBB.x", # 
    
    str_detect(sublineage, "BA\\.4") ~ "BA.4.x", #
    
    str_detect(sublineage, "BQ\\.1") ~ "BQ.1.x", #
    str_detect(sublineage, "BA\\.5") ~ "BA.5.x", #
    
    # TRUE ~ sublineage
    TRUE ~ "other"
  )) |> 
  mutate(sublineage_new=factor(sublineage_new)) |>
  mutate(sublineage_new = factor(sublineage_new,
                            levels = c("B.1.1.529.x",
                                       "BA.1.x",
                                       "BA.2.x",
                                       "BA.2.86.x", "JN.1.x", "JN.1.11.1.x",
                                       "BA.2.75.x", "CH.1.1.x",
                                       "XBB.x", "XBB.1.5.x", "XBB.1.5.70.x",
                                       "XBB.1.16.x",
                                       "XBB.1.9.x", "EG.5.1.x", "HK.3.x",
                                       "XBB.2.3.x",
                                       "BA.2.12.1.x",
                                       "BA.5.x",
                                       "BQ.1.x",
                                       "BA.4.x",
                                       "other")))

### Demix progress report

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

# Weekly Lineage Composition x Neighborhood -------------------------------

year_epiweek_levels <- sublin_meta |> 
  mutate(year_epiweek = as.character(year_epiweek)) |> 
  distinct(year_epiweek) |> 
  arrange(as.integer(year_epiweek)) |> 
  pull(year_epiweek)

sublin_meta_complete <- sublin_meta |> 
  mutate(year_epiweek = factor(as.character(year_epiweek), levels = year_epiweek_levels)) |> 
  group_by(LOCATION) |> 
  complete(year_epiweek = year_epiweek_levels, sublineage_new, fill = list(abundance = 0)) |> 
  ungroup()

p_RAxNB <- sublin_meta_complete |> 
  ggplot(aes(x = year_epiweek, y = abundance, fill = sublineage_new)) +
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

ggsave("../figures/p_RAxNB.jpg", p_RAxNB) # 2025-11-10


