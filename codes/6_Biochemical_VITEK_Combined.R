## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project


#1-12-25_ Madelaine is unsure if this step is needed at all? Will not all datasets be joined in 9. anyway?

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)

## Importing the VITEK RESULTS
AST_ISOLATES <- readRDS("data/CLEANED_DATA/FINAL_AST_RESULTS.rds")

## Select relevant columns
AST_ISOLATES <- AST_ISOLATES %>%
  select(Isolate_ID, TVLA_ID, VITEK_MS_Results, ESC, ESBL_Selection,ESBL_Presumptivefinal,Isolate,Isolatnr,Maldi_Tof_25_2_26_tmf)

# AST_ISOLATES <- AST_ISOLATES %>%
#   select(Isolate_ID, TVLA_ID, VITEK_MS_Results, ESC, ESBL)

## Importing the Biochemical Results ( The FINAL_VITEK_MS_SELECTED-2025-10-09.rds" have been created in 1.2 from joining the Isolates with the biochemical results (VITEK_MS_SELECTION <- 
#   left_join(Cleaned_Labdata_Original, 
#             Selection %>%) and is named VITEK_MS_SELECTED before written out, not sure this had to be done!!!!!!!!!!!!!!!!!!!!! )
VITEK_MS_SELECTED <- readRDS ("data/CLEANED_DATA/FINAL_VITEK_MS_SELECTED-2025-10-09.rds")

# Extract part before first "-" We need to make sure we have a correct INIKA_ID for all, it was an "I" before some of the IDs in this file above
VITEK_MS_SELECTED<-VITEK_MS_SELECTED%>%
  mutate(INIKA_ID =  sub("_.*", "", Isolate_ID))



 

# Note Isolate_ID is not always correct- need to be corrected above! MN
## Join the VITEK_MS_SELECTED with AST_ISOLATES
Biochemical_VITEK_Combined <- left_join(VITEK_MS_SELECTED, AST_ISOLATES, by = "Isolate_ID")

# Check which ones not joined
#Biochemical_VITEK_Combined_anti<- anti_join(VITEK_MS_SELECTED, AST_ISOLATES, by = "Isolate_ID")


# check dublicates 
Biochemical_VITEK_Combined_duplicates <- Biochemical_VITEK_Combined$Isolate_ID[duplicated(Biochemical_VITEK_Combined$Isolate_ID)]
view(Biochemical_VITEK_Combined_duplicates)
# No duplicates!

## Save the file
write_tsv(Biochemical_VITEK_Combined,"data/CLEANED_DATA/Biochemical_VITEK_Combined.tsv")
saveRDS(Biochemical_VITEK_Combined,"data/CLEANED_DATA/Biochemical_VITEK_Combined.rds")

