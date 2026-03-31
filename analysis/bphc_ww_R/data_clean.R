# Metadata cleaning


# Libraries ---------------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

ww_metadata <- read_csv("../../metadata/bphc_biobot_sequence_metadata.csv")

baseload_ids <- read_table("../../metadata/baseload_ids.txt",col_names = F)

clin_lin <- read_tsv("../data/var_surv_less_filt.tsv")

ww_pcr <- read_csv("../data/pcr_final_baseload.csv")

# coverage_dir <- "../../coverage/"
# coverage <- list.files(coverage_dir, pattern = "\\.qc\\.tsv$", full.names = TRUE) |>
#   map(read_tsv, show_col_types = FALSE) |>
#   list_rbind()
# coverage <- coverage |> mutate(sample=as.character(sample))
# 
# write_rds(coverage, "../data/coverage.rds") # 2026-03-17
# rm(coverage_dir, coverage)


coverage_qc <- readRDS("../data/coverage.rds")


# Wastewater metadata -----------------------------------------------------


ww_metadata <- ww_metadata[,1:6] |> 
  # format dates
  rename(SAMPLING_DATE = SAMPLING_DATE_FOR_REPORT,
         LOCATION = LOCATION_NAME) |> 
  mutate(SAMPLING_DATE = as_date(mdy(SAMPLING_DATE)),
         DATE_SENT_SEQ = as_date(mdy(DATE_SENT_SEQ)),
         BATCH_ID=as.character(BATCH_ID),
         LOCATION=as_factor(LOCATION),
         FASTQ_ID=as.character(FASTQ_ID)) |> 
  mutate(LOCATION=case_when(
    LOCATION=="Dorchester_2224" ~ "Dorchester",
    LOCATION=="Roslindale_WestRoxbury" ~ "Roslindale West Roxbury",
    LOCATION=="BackBay" ~ "Back Bay",
    LOCATION=="Lower_Roxbury" ~ "Lower Roxbury",
    T ~ LOCATION
  ))

baseload_ids <- baseload_ids |> 
  rename(FQ_ID=X1) |> 
  mutate(FQ_ID=as.character(FQ_ID))

### cross checks

  # fastq files without metadata?
baseload_ids$FQ_ID %in% ww_metadata$FASTQ_ID |> table()
  # all have metadata

  # metadata without fastq?
ww_metadata$FASTQ_ID[!ww_metadata$FASTQ_ID %in% baseload_ids$FQ_ID]
  # 20 samples formatted like 2022-1021-91H06WW

ww_metadata <- ww_metadata |> 
  # drop metadata with no fastq file
  filter(FASTQ_ID %in% baseload_ids$FQ_ID)

### sampling by neighborhood

  # add year_epiweek and sort
ww_metadata <- ww_metadata |>
  mutate(epiweek = epiweek(SAMPLING_DATE),
         year = year(SAMPLING_DATE),
         year_epiweek = paste0(year, sprintf("%02d", epiweek))) |>
  mutate(year_epiweek=as.numeric(year_epiweek)) |> 
  arrange(year_epiweek, LOCATION)

counts <- ww_metadata |> 
  count(LOCATION, year_epiweek) |>
  arrange(year_epiweek, LOCATION) |>
  mutate(year_epiweek = factor(year_epiweek))

  # notes RE Lower_Roxbury and South Boston
  # Lower_Roxbury - sampling did not start here until 2024
  # South Boston - sampling ended after 2024
  # sewer was under construction, moved sampling from S Bos to L Rox
  # L Rox and Rox in same zip 
  # for concentration data, compared L Rox to Rox to see if similar
  # could think about a weighted average
  
counts |> filter(n>1) 
  # Charlestown 202334
  # Lower_Roxbury 202410

ww_metadata |>
  filter(LOCATION=="Charlestown" & year_epiweek=="202334" | 
           LOCATION=="Lower Roxbury" & year_epiweek=="202410")

  # Charlestown 202334
  # could be mislabling, biobot did not have info
  # EXCLUDE BOTH

  # Lower_Roxbury 202410
  # 

exclude <- ww_metadata |>
  filter(LOCATION=="Charlestown" & year_epiweek=="202334" | 
           LOCATION=="Lower Roxbury" & year_epiweek=="202410") |> 
  select(FASTQ_ID)

ww_metadata <- ww_metadata |> filter(!FASTQ_ID %in% exclude$FASTQ_ID)

p_counts <- ww_metadata |>
  count(LOCATION, year_epiweek) |>
  arrange(year_epiweek, LOCATION) |>
  mutate(year_epiweek = factor(year_epiweek)) |>
  ggplot(aes(x=year_epiweek, y=LOCATION, fill=as.factor(n))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#4C9AED", "#654CED"), name="sample count")
p_counts

# ggsave("../figures/p_counts.jpg", p_counts, scale = 1.5)

  # RE missing weeks
  # vendor may have not sent samples if below quality threshold
  # check concentration metadata to see if they exist there
  # seq data may have also not been sent for non-detect weeks
  # should be able to see in concentration data
  # nondetects =/= below LOD
  # check PCR folder
  # link with KIT ID

### join with QC data

ww_metadata <- ww_metadata |> 
  left_join(coverage_qc, by=c("FASTQ_ID"="sample"))

# write_rds(ww_metadata, "../data/meta_clean.rds") # 2026-03-17

  # check which have failed QC

ww_metadata |>
  arrange(year_epiweek, LOCATION) |>
  mutate(year_epiweek = factor(year_epiweek)) |>
  mutate(QC_ov_depth10 = ifelse(frac_genome_cov_10x < 0.75, 1,0),
         QC_spk_depth10 = ifelse(spike_frac_cov_10x < 0.75, 1,0),
         fail_either75 = case_when(
           QC_ov_depth10 + QC_spk_depth10 > 1 ~ 1,
           QC_ov_depth10 + QC_spk_depth10 <= 1 ~ 0
         )) |> 
  drop_na(QC_ov_depth10) |> 
  ggplot(aes(x=year_epiweek, y=LOCATION, fill=as.factor(fail_either75))) +
  geom_tile() +
  geom_text(aes(label=round(spike_frac_cov_10x, digits = 2))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("lightgreen", "red"), name="fail either\n 75% @ 10x ?") +
  ggtitle("QC report", subtitle = "labels: spike coverage > 10x")

ww_metadata |> 
  select(frac_genome_cov_10x, spike_frac_cov_10x) |> 
  plot() # pretty correlated


# Clinical data from GISAID -----------------------------------------------


### format to match ww metadata

clin_lin <- clin_lin |>
  mutate(Sample = `Accession ID`,
         sublineage = `Pango lineage`,
         SAMPLING_DATE = `Collection date`,
         .keep = "none") |> 
  mutate(epiweek = epiweek(SAMPLING_DATE),
         year = year(SAMPLING_DATE),
         year_epiweek = paste0(year, sprintf("%02d", epiweek))) |>
  mutate(year_epiweek=as.numeric(year_epiweek))

# write_rds(clin_lin, "../data/clin_lin.rds") # 2026-02-05


# WW PCR data -------------------------------------------------------------


ww_pcr <- ww_pcr |> 
  rename(effCopiesL=biobot_effective_sarscov2_concentration_copies_per_liter,
         rawCopiesL=biobot_raw_sarscov2_concentration_copies_per_liter,
         LOCATION=wwtp_name,
         SAMPLING_DATE=sample_collect_date) |> 
  mutate(SAMPLING_DATE = as.Date(trimws(SAMPLING_DATE), format="%m/%d/%Y"),
         epiweek = epiweek(SAMPLING_DATE),
         year = year(SAMPLING_DATE),
         year_epiweek = paste0(year, sprintf("%02d", epiweek))) |>
  mutate(year_epiweek=as.numeric(year_epiweek)) |> 
  mutate(LOCATION=case_when(
    LOCATION=="Dorchester_2224" ~ "Dorchester",
    LOCATION=="Roslindale_WestRoxbury" ~ "Roslindale West Roxbury",
    LOCATION=="BackBay" ~ "Back Bay",
    LOCATION=="Lower_Roxbury" ~ "Lower Roxbury",
    T ~ LOCATION
  )) |> 
  select(-hum_frac_mic_conc, -hum_frac_mic_unit, -rec_eff_percent, -rec_eff_target_name)

# write_rds(ww_pcr, "../data/ww_pcr.rds") # 2026-02-10

### sites might have 1 or 2 samples per week (1 case of 3 samples)
### 

# create population weights

pop_wts <- ww_pcr |> 
  select(LOCATION, population_served) |> 
  group_by(LOCATION) |> 
  summarise(pop=max(population_served)) |> 
  ungroup() |> 
  mutate(popwt = pop/sum(pop))


# write_rds(pop_wts, "../data/pop_wts.rds") # 2026-02-10

rm(list=setdiff(ls(), c("clin_lin","ww_metadata","ww_pcr","pop_wts")))


# create general weights

conc_wts <- ww_pcr |> 
  group_by(LOCATION, year_epiweek) |> 
  summarise(
    meaneffCopiesL = mean(effCopiesL),
    meanrawCopiesL = mean(rawCopiesL),
    .groups = "drop"
  ) |> 
  
  # replace 0 values with 1s before taking log
  mutate(meaneffCopiesL = ifelse(meaneffCopiesL<=0,1,meaneffCopiesL),
         meanrawCopiesL = ifelse(meanrawCopiesL<=0,1,meanrawCopiesL)) |> 
  
  mutate(logeff = log(meaneffCopiesL),
         lograw = log(meanrawCopiesL))

# write_rds(conc_wts, "../data/conc_wts.rds") # 2026-03-09



