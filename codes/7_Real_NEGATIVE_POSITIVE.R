## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project

# 1.12.25 Madelaine also wonders whether this script is needed at all?

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)

## Importing the Biochemical_VITEK_Combined
Biochemical_VITEK_Combined <- readRDS("data/CLEANED_DATA/Biochemical_VITEK_Combined.rds")

# Select relevant columns
Biochemical_VITEK_Combined <- Biochemical_VITEK_Combined %>%
  select(Isolate_ID,VITEK,TVLA_ID,VITEK_MS_Results,ESC)

# Importing the Cleaned_Lab data_Original
Cleaned_Labdata_Original <- readRDS("data/CLEANED_DATA/1_Cleaned_Labdata_Original-2025-10-07.rds")

## Join the Cleaned_Lab data_Original and Biochemical_VITEK_Combined
Real_NEGATIVE_POSITIVE <- left_join(Cleaned_Labdata_Original, Biochemical_VITEK_Combined, by = "Isolate_ID")





