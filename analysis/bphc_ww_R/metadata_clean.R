# Metadata cleaning


# Libraries ---------------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

metadata <- read_csv("../../metadata/bphc_biobot_sequence_metadata.csv")

baseload_ids <- read_table("../../metadata/baseload_ids.txt",col_names = F)

# Cleaning ----------------------------------------------------------------

metadata <- metadata[,1:6] |> 
  # format dates
  rename(SAMPLING_DATE = SAMPLING_DATE_FOR_REPORT,
         LOCATION = LOCATION_NAME) |> 
  mutate(SAMPLING_DATE = as_date(mdy(SAMPLING_DATE)),
         DATE_SENT_SEQ = as_date(mdy(DATE_SENT_SEQ))) |> 
  # batch ID as character
  mutate(BATCH_ID=as.character(BATCH_ID)) |> 
  # locations as factors
  mutate(LOCATION=as_factor(LOCATION))

baseload_ids <- baseload_ids |> 
  rename(FQ_ID=X1) |> 
  mutate(FQ_ID=as.character(FQ_ID))

# Cross checks ------------------------------------------------------------

  # fastq files without metadata?
baseload_ids$FQ_ID %in% metadata$FASTQ_ID |> table()
  # all have metadata

  # metadata without fastq?
metadata$FASTQ_ID[!metadata$FASTQ_ID %in% baseload_ids$FQ_ID]
  # 20 samples formatted like 2022-1021-91H06WW

metadata <- metadata |> 
  # drop metadata with no fastq file
  filter(FASTQ_ID %in% baseload_ids$FQ_ID)


# Sampling by neighborhood ------------------------------------------------

counts <- metadata |>
  mutate(epiweek = epiweek(SAMPLING_DATE),
         year = year(SAMPLING_DATE),
         year_epiweek = paste(year,epiweek,sep = "")) |>
  mutate(year_epiweek=as.numeric(year_epiweek)) |>
  count(LOCATION, year_epiweek) |>
  arrange(year_epiweek, LOCATION) |>
  mutate(year_epiweek = factor(year_epiweek)) |>
  filter(LOCATION!="Lower_Roxbury", # sampling did not start here until 2024
         LOCATION!="South Boston") # sampling ended after 2024

counts |> filter(n>1) # Charlestown 202334

metadata |>
  mutate(epiweek = epiweek(SAMPLING_DATE),
         year = year(SAMPLING_DATE),
         year_epiweek = paste(year,epiweek,sep = "")) |>
  mutate(year_epiweek=as.numeric(year_epiweek)) |> 
  filter(LOCATION=="Charlestown", year_epiweek==202334)
  # 2 samples from 2023-08-21
  # drop one?

plot_counts <- counts |>
  ggplot(aes(x=year_epiweek, y=LOCATION, fill=as.factor(n))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#4C9AED", "#654CED"), name="sample count")
plot_counts

