# Pipeline prep


# Libraries ---------------------------------------------------------------


library(tidyverse)
library(lubridate)


# Load data ---------------------------------------------------------------

storage_dir <- "/n/holylfs05/LABS/hanage_lab/Lab/hsphfs1/bschaeffer/bphc_ww/data/"

meta_clean <- readRDS(paste0(storage_dir,"meta_clean.rds"))

## Baseload batch 1 IDs

baseload_batch1_ids <- meta_clean |> 
  filter(SAMPLING_DATE < as.Date("2023-05-05")) |> 
  pull(FASTQ_ID)

## Baseload batch 2 IDs

baseload_batch2_ids <- meta_clean |> 
  filter(SAMPLING_DATE > as.Date("2023-05-05")) |> 
  pull(FASTQ_ID)

baseload_batch2_ids |> 
  write_lines(file = paste0(storage_dir,"baseload_batch2_ids.txt"))


