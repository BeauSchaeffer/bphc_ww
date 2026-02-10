# WW PCR concentration data


# Libraries ---------------------------------------------------------------


library(tidyverse)


# Load data ---------------------------------------------------------------


ww_pcr <- read_csv("../data/pcr_final_baseload.csv")


ww_pcr |> 
  mutate(sample_collect_date = as.Date(trimws(sample_collect_date), format="%m/%d/%Y"),
         epiwk = epiweek(sample_collect_date),
         year = year(sample_collect_date)) |> 
  
  mutate(sample_collect_date)
