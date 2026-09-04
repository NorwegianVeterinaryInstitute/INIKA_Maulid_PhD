
#########################################################################
## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project


## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)


## Importing the DRY_RAINY_MWZ_KILIMANJARO_AST2 file ----
# Might have been checked , but Madelaine like to check again later! 17.10.25

DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  read_excel("data/DRY_RAINY_MWZ_KILIMANJARO.xlsx", 
             sheet = "DRY_RAINY_MWZ_KILIMANJARO_AST2",
             na = c("", "NA")
  )


# Selecting relevant columns
DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  select(-`S/No`, -`SAMPLE TYPE`,
         -`Season Included in Number? ( evt. how?)`,-Comment_)



## Replacing the NIKA_ with INIKA_ and INIKA_ to " " in INIKA_ID column

DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  mutate(
    INIKA_ID = str_replace_all(INIKA_ID, "NIKA_", "INIKA_"), # First replacement
    INIKA_ID = str_replace_all(INIKA_ID, "INIKA_", "") , # Second replacement on the already modified column
    INIKA_ID = str_replace_all(INIKA_ID, "I", "")  
  )


DRY_RAINY_MWZ_KILIMANJARO_AST2%>%
  mutate_if(is.character, factor) %>%
  lapply(levels)

DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  mutate(
    # Extract first part of INIKA_ID
    INIKA_prefix = str_extract(INIKA_ID, "^[^_]+"),
    
    # Remove first part of Isolate_ID
    Isolate_suffix = str_replace(Freezing_ID_No, "^[^_]+_", ""),
    
    # Combine to make corrected ID
    Isolate_ID = paste0(INIKA_prefix, "_", Isolate_suffix)
  )


#Check for duplicates
duplicated_test<-DRY_RAINY_MWZ_KILIMANJARO_AST2$Isolate_ID[duplicated(DRY_RAINY_MWZ_KILIMANJARO_AST2$Isolate_ID)]
print(duplicated_test)



# 1. Group the data by Isolate_ID

DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
group_by(Isolate_ID) %>%
  # 2. Count the number of times each Isolate_ID appears
  mutate(Isolate_ID_repeat_count = n()) %>%
  # 3. Assign a sequential count (1, 2, 3...) within each Isolate_ID group
  mutate(Isolate_ID_row_count = row_number()) %>%
  # 4. Remove the grouping
  ungroup() %>%
  # 5. Modify the Isolate_ID based on the counts
  mutate(Isolate_ID = case_when(
    # The condition: If Isolate_ID appeared exactly 2 times AND it's the second row
    Isolate_ID_repeat_count == 2 & Isolate_ID_row_count == 2 ~
      # Action: Append the row count (2) to the Isolate_ID
      paste(Isolate_ID, Isolate_ID_row_count, sep = "_"),
    
    # The default: Otherwise, keep the original Isolate_ID
    TRUE ~ Isolate_ID
  ))




## Re labeling the data in the columns
# List the columns to be modified

cols_to_relabel <- c(
  "COLONY MORPHOLOGY ON C3GR", 
  "COLONY MORPHOLOGY ON CARBA", 
  "COLONY MORPHOLOGY ON XLD", 
  "COLONY MORPHOLOGY ON BGA"
) 

# Apply the replacement to the selected columns
DRY_RAINY_MWZ_KILIMANJARO_AST2 <- DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  mutate(
    across(
      .cols = all_of(cols_to_relabel),
      .fns = ~ replace(., . == "NBG", "No growth")
    )
  )

## Check the homogeneity of the data registering ----
# Trick to select all the columns of type character in DRY_RAINY_MWZ_KILIMANJARO_AST2 dataset
DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  select(where(is.character)) %>%
  # I modify all the columns from character to factor - because I want to see the "levels" : categories
  # to identify if you wrote homogeneously 
  mutate_all(factor) %>%
  # Select one column and look at the levels
  lapply(levels)

## Further cleaning for homogeneity entries
DRY_RAINY_MWZ_KILIMANJARO_AST2 <- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  # `SAMPLE TYPE` is only Stool BUT I do not want to remove because I want to combine with other dataset. 
  # `SAMPLE FROM` : Fix syntax  "Adult outPatient" "Adult outpatient"
  mutate(`SAMPLE FROM` = str_replace(`SAMPLE FROM`, "Adult outPatient", "Adult outpatient")) %>%
  # `NAME OF SCHOOL/HF` Homeogeneize :  "Pasua HC"  "PasuaHC"   From To 
  mutate(`NAME OF SCHOOL/HF` = str_replace(`NAME OF SCHOOL/HF`, "PasuaHC", "Pasua HC")) %>%
  #  `COLONY MORPHOLOGY ON C3GR` "Pinkish/Reddish" <- keep this one "pinkish/Reddish"
  mutate(`COLONY MORPHOLOGY ON C3GR` = str_replace(`COLONY MORPHOLOGY ON C3GR`, "pinkish/Reddish", "Pinkish/Reddish")) %>%
  # `COLONY MORPHOLOGY ON XLD`  "Black colonies"       "Black colony"         "Blackish colonies"  need to be identical
  ## "No black colony"      "No blackish colonies" must be identical - Fixed with several choices of equality
  mutate(`COLONY MORPHOLOGY ON XLD` = 
           case_when(
             `COLONY MORPHOLOGY ON XLD`=="No blackish colonies" ~ "No black colony",
             `COLONY MORPHOLOGY ON XLD`=="Black colonies" ~ "Black colony",
             `COLONY MORPHOLOGY ON XLD`=="Blackish colonies" ~ "Black colony",
             TRUE ~ `COLONY MORPHOLOGY ON XLD`)
  ) 
#NB: $CITRATE  "+" "-" "d" ->  intermediary : its a standard code : we keep it like that 

## Renaming the data entries in Isolate column
DRY_RAINY_MWZ_KILIMANJARO_AST2 <- DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  mutate(
    Isolate = case_when(
      Isolate == "S.typhimariam" ~ "S.typhimurium",
      # Default: Keep all other values in the 'Isolate' column as they are
      TRUE ~ Isolate 
    )
  )


### Inserting new column "VITEK" ----
# Create the new column 'VITEK' and assign the value 1 to all rows, this means these have been selected for VITEK confirmation
Selection<- 
  DRY_RAINY_MWZ_KILIMANJARO_AST2 %>%
  mutate(VITEK = 1) %>%
  # Relocate the new column immediately after 'Isolate_ID'
  relocate(VITEK, .after = Isolate_ID)

# TODO : add comment of what VITEK 1 means  
#26.11.2025 Check not writing out new files
## Saving the file as tsv to be opened in excel format
#write_tsv(Selection, "data/CLEANED_DATA/Selection-2026-19-07.tsv")
write_tsv(Selection, "data/CLEANED_DATA/Selection-2026-03-11.tsv")

# Save as rds to be used further work in R
#saveRDS(Selection, "data/CLEANED_DATA/Selection-2025-10-07.rds")
saveRDS(Selection, "data/CLEANED_DATA/Selection-2025-03-11.rds")
# This one we use for the join in script number 9. Joining..
############################################################################
# ### Now join "Selection" back with " Cleaned_Labdata_Original"
# dim(Cleaned_Labdata_Original)
# dim(Selection)
# names(Cleaned_Labdata_Original)
# names(Selection)
# 
# 
# VITEK_MS_SELECTION <- 
#   left_join(Cleaned_Labdata_Original, 
#             Selection %>%
#               select(-INIKA_ID, -SEASON, -`SAMPLE FROM`, -`NAME OF SCHOOL/HF`, 
#                      -REGION, -DISTRICT, -AR_RESIDUES_SAMPLE, 
#                      -`COLONY MORPHOLOGY ON C3GR`, -`COLONY MORPHOLOGY ON CARBA`,
#                      -`COLONY MORPHOLOGY ON XLD`, -`COLONY MORPHOLOGY ON BGA`, 
#                      -`TSI Media_Slope`, -`TSI Media_Butt`, -`TSI Media_Gas`, 
#                      -'TSIMedia_H2S, -`SIM Media_H2S`, -`SIM Media Indole`, 
#                      -`SIM Media Motility`, -CITRATE, -UREASE, -Isolate),
#             by = "Isolate_ID"
#   ) 
# 
# ## Creating new column as "Not selected" (0), when the value in "VITEK" column is not selected as 1
# 
# VITEK_MS_SELECTED <- 
#   VITEK_MS_SELECTION %>%
#   mutate(VITEK = as.character(VITEK)) %>%
#   mutate(VITEK = case_when(
#     is.na(VITEK) ~ "0",
#     TRUE ~ VITEK
#   ))
# 
# 
# # Select the specified isolates from the VITEK_MS_SELECTED dataset
# VITEK_MS_SELECTED <-
#   VITEK_MS_SELECTED %>%
#   # Filter rows where the value in the 'Isolate' column is one of the listed organisms
#   filter(Isolate %in% c("E.coli", "K.pneumoniae", "S.typhimurium")) %>%
#   filter(VITEK == "1" )
# 
# ### Checking further
# 
# # Filter the data for S.typhimurium 53
# S_typhimurium_subset <- 
#   Selection[Selection$Isolate == "S.typhimurium", ]
# # Find which INIKA_ID values in the subset are duplicates
# duplicate_ids <- 
#   S_typhimurium_subset$INIKA_ID[
#     duplicated(S_typhimurium_subset$INIKA_ID)
#   ]
# 
# #"I23272" "I23346"
# # MN I see you have seen that you have IDs starting with I here.
# # HELPER to clean dataset 
# S_typhimurium_subset %>%
#   filter(INIKA_ID %in% duplicate_ids) %>%
#   View()
# ## 23272 and 23346 ## These two numbers are all the same, one per each number should be selected
# 
# # Selects one unique row for each INIKA_ID, keeping the first occurrence.
# S_typhimurium_subset <- 
#   S_typhimurium_subset %>%
#   distinct(INIKA_ID, .keep_all = TRUE)  ## 51
# # Here its work ok, but I would specify (best practice)
# # Here it removes 23272_2_R  and 23346_2_R that are exacly the same isolates as .._1_R
# 
# 
# ## E.coli Checking for duplication in INIKA_ID
# # Filter the data for e.coli # 196
# E.coli_subset <- 
#   Selection[Selection$Isolate == "E.coli", ]
# # Find which INIKA_ID values in the subset are duplicates
# duplicate_ids <- 
#   E.coli_subset$INIKA_ID[
#     duplicated(E.coli_subset$INIKA_ID)
#   ]
# # Filter the original subset to keep only the rows with duplicate IDs
# duplicates_E.coli_base <- E.coli_subset[
#   E.coli_subset$INIKA_ID %in% duplicate_ids,
# ]
# 
# # Filter the original subset to keep only the rows with duplicate IDs
# duplicates_E.coli_base <- 
#   E.coli_subset[E.coli_subset$INIKA_ID %in% duplicate_ids, ] ## 0
# 
# 
# 
# ## K.pneumoniae Checking for the duplicates in INKA_ID
# # Filter the data for K.pneumoniae # 162
# K.pneumoniae_subset <- 
#   Selection[Selection$Isolate == "K.pneumoniae", ]
# # Find which INIKA_ID values in the subset are duplicates
# duplicate_ids <- 
#   K.pneumoniae_subset$INIKA_ID[ duplicated(K.pneumoniae_subset$INIKA_ID)]
# # Filter the original subset to keep only the rows with duplicate IDs
# duplicates_K.pneumoniae_base <- K.pneumoniae_subset[
#   K.pneumoniae_subset$INIKA_ID %in% duplicate_ids,
# ]
# 
# # Filter the original subset to keep only the rows with duplicate IDs
# duplicates_K.pneumoniae_base <- K.pneumoniae_subset[
#   K.pneumoniae_subset$INIKA_ID %in% duplicate_ids,
# ] ## 0
# 
# ## Now join the three sub data sets for the VITEK_MS confirmation
# # Combine all three data sets
# VITEK_MS_SELECTED <- rbind(E.coli_subset, K.pneumoniae_subset, S_typhimurium_subset)
# 
# ## Saving the file
# write_tsv(VITEK_MS_SELECTED, "data/CLEANED_DATA/FINAL_VITEK_MS_SELECTED-2025-10-09.tsv")
# saveRDS(VITEK_MS_SELECTED, "data/CLEANED_DATA/FINAL_VITEK_MS_SELECTED-2025-10-09.rds")
# #######################################################################################