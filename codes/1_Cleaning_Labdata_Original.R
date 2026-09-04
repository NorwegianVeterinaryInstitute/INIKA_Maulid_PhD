## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project


## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)

## Importing the original lab data
DRY_RAINY_MWZ_KILIMANJARO <- read_excel("data/DRY_RAINY_MWZ_KILIMANJARO.xlsx",
                                        na = c("", "NA")
                                        )



# glimpse(DRY_RAINY_MWZ_KILIMANJARO)
# str(DRY_RAINY_MWZ_KILIMANJARO)
# names(DRY_RAINY_MWZ_KILIMANJARO)

## Part1: Data cleaning
# Selecting relevant columns
Cleaned_Labdata_Original <- 
  DRY_RAINY_MWZ_KILIMANJARO %>%
  select(-`S/No`, -`SAMPLE TYPE`,
         -`Season Included in Number? ( evt. how?)`,-Comment_)

## Replacing the NIKA_ with INIKA_ and INIKA_ to " " in INIKA_ID column

Cleaned_Labdata_Original <-
  Cleaned_Labdata_Original %>%
  mutate(INIKA_ID = str_replace_all(INIKA_ID, "(NIKA_)|(INIKA_)", ""))


# Some Freezing_IDs have been misspelled
 
# Don't think this applies anymore after the cleaning process; We detected when trying to join the labdata with whonetdata that it has been wrong written in the Freexing_ID
#3 isolates  11321_1_R,  22362_2_R 

  

# We detected that there are spelling mistakes in some of the Isolates_Ids as they are not corresponding to the correct INIKA_IDs.
# Solution: We make a new Isolate_ID calling it Isolate_ID_Corrected, but need to keep the wrong one as well in order to manage to join the data from TVLA (VITEK)


Cleaned_Labdata_Original<-Cleaned_Labdata_Original%>%
  mutate(
    # Extract first part of INIKA_ID
    INIKA_prefix = str_extract(INIKA_ID, "^[^_]+"),
    
    # Remove first part of Isolate_ID
    Isolate_suffix = str_replace(Freezing_ID_No, "^[^_]+_", ""),
    
    # Combine to make corrected ID
    Isolate_ID = paste0(INIKA_prefix, "_", Isolate_suffix)
  )

## Renaming the columns
# Cleaned_Labdata_Original <- 
# Cleaned_Labdata_Original %>%
 #rename(Isolate_ID = Isolate_ID_Corrected) 

# Freezing_IDS still not corrected, needs to be kept, as TVLA has used them
Cleaned_Labdata_Original %>%
  mutate_if(is.character, factor) %>%
  lapply(levels)

## Re labeling the data in the columns
# List the columns to be modified

cols_to_relabel <- c(
  "COLONY MORPHOLOGY ON C3GR", 
  "COLONY MORPHOLOGY ON CARBA", 
  "COLONY MORPHOLOGY ON XLD", 
  "COLONY MORPHOLOGY ON BGA"
) 

# Apply the replacement to the selected columns
Cleaned_Labdata_Original <- 
  Cleaned_Labdata_Original %>%
  mutate(
    across(
      .cols = all_of(cols_to_relabel),
      .fns = ~ replace(., . == "NBG", "No growth")
    )
  )

## Check the homogeneity of the data registering ----
# Trick to select all the columns of type character in Cleaned_Labdata_Original dataset
Cleaned_Labdata_Original %>%
  select(where(is.character)) %>%
  # I modify all the columns from character to factor - because I want to see the "levels" : categories
  # to identify if you wrote homogeneously 
  mutate_all(factor) %>%
  # Select one column and look at the levels
  lapply(levels)

## Further cleaning for homogeneity entries
Cleaned_Labdata_Original <- 
  Cleaned_Labdata_Original %>%
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
Cleaned_Labdata_Original <- Cleaned_Labdata_Original %>%
  mutate(
    Isolate = case_when(
      Isolate == "S.paratyphi A (Except Indole)" ~ "S.paratyphi A",
      Isolate == "S.typhi (Except Gas)" ~ "S.typhi",
      Isolate == "S.typhimariam" ~ "S.typhimurium",
      Isolate == "Shigella spp ( Except Gas)" ~ "Shigella spp",
      # Default: Keep all other values in the 'Isolate' column as they are
      TRUE ~ Isolate 
    )
  )

Cleaned_Labdata_Original %>%
str(.)


## Replacing NA in Isolate column with No growth
Cleaned_Labdata_Original$Isolate[is.na(Cleaned_Labdata_Original$Isolate)] <- "No growth"
#Check for duplicates
duplicated_test<-Cleaned_Labdata_Original$Isolate_ID[duplicated(Cleaned_Labdata_Original$Isolate_ID)]
print(duplicated_test)

#6.12.25 Maulid, I don't see what have been changed? or do you mean # infront?
Cleaned_Labdata_Original <- Cleaned_Labdata_Original %>%
  mutate(
    Isolate_ID = case_when(
      Isolate_ID == "1392_1_D" & Isolate == "K.aerogenes" ~ "1392_2_D",
      Isolate_ID == "1395_1_D" & Isolate == "P.aeruginosa" ~ "1395_2_D",
      Isolate_ID == "13150_2_D" & Isolate == "P.aeruginosa" ~ "13150_3_D",
      Isolate_ID == "21345_1_R" & Isolate == "K.aerogenes" ~ "21345_2_R",
      Isolate_ID == "22287_1_R" & Isolate == "K.aerogenes" ~ "22287_2_R",
      Isolate_ID == "2374_1_D" & Isolate == "C.freundii" ~ "2374_2_D",
      Isolate_ID == "23130_1_D" & Isolate == "C.freundii" ~ "23130_2_D",
      TRUE ~ Isolate_ID
    )) %>%
  mutate(
    INIKA_ID = case_when(
      # Rename 'INIKA_22315' to 'INIKA_22316' only when Freezing_ID_No is '22316_1_R'
      INIKA_ID == "22315" & Freezing_ID_No == "22316_1_R" ~ "22316",
      TRUE ~ INIKA_ID
    )) %>%
  mutate(
    Isolate_ID = case_when(
      Isolate_ID == "22315_1_R" &
        INIKA_ID == "22316" &
        Freezing_ID_No == "22316_1_R" ~ "22316_1_R",
      TRUE ~ Isolate_ID
    ))
      
      
# Get duplicated IDs (including first occurrences)
duplicated_test<-Cleaned_Labdata_Original$Isolate_ID[duplicated(Cleaned_Labdata_Original$Isolate_ID, fromLast = TRUE)]

# Remove duplicates from the list itself
unique_duplicated_ids <- unique(duplicated_test)

# Collapse into a single text string
text_string <- paste(unique_duplicated_ids, collapse = ", ")

# View result
print(text_string) # "21303_2_R"
# 11227_1_R, 11360_1_R, 11366_1_R, 126_1_D, 12134_1_D, 1392_1_D, 1395_1_D, 13150_2_D, 21113_1_D, 21303_2_R, 21345_1_R, 22262_2_R, 22287_1_R, 2374_1_D, 23130_1_D"

CheckingDuplicates <- Cleaned_Labdata_Original%>%
  filter(Isolate_ID %in% 
           c("21303_2_R"))%>%
  select(INIKA_ID, Isolate_ID, Freezing_ID_No,Isolate,REGION, SEASON, `SAMPLE FROM` ) # "1392_1_D", "1395_1_D", "13150_2_D", , "21345_1_R", "22287_1_R", "22315_1_R", "2374_1_D", "23130_1_D"

# Ok_ 8.12.25 Madelaine: Now I see these changes, needs to be included in the rmd script as well!, and make sure it is correct for the joining with Malditof results.
# We detected that there are spelling mistakes in some of the Isolates_Ids as they are not corresponding to the correct INIKA_IDs.
# Solution: We make a new Isolate_ID calling it Isolate_ID_Corrected, but need to keep the wrong one as well in order to manage to join the data from TVLA (VITEK)



# Get duplicated IDs (including first occurrences)
duplicated_test<-Cleaned_Labdata_Original$Isolate_ID[duplicated(Cleaned_Labdata_Original$Isolate_ID, fromLast = TRUE)]

# Remove duplicates from the list itself
unique_duplicated_ids <- unique(duplicated_test)
# Collapse into a single text string
text_string <- paste(unique_duplicated_ids, collapse = ", ")

# View result
print(text_string)


CheckingDuplicates <- Cleaned_Labdata_Original%>%
  filter(Isolate_ID %in% 
           c("21303_2_R"))%>%
  select(INIKA_ID, Isolate_ID,Isolate,REGION, SEASON, `SAMPLE FROM`,`COLONY MORPHOLOGY ON BGA`, `COLONY MORPHOLOGY ON CARBA`, `COLONY MORPHOLOGY ON XLD`, `COLONY MORPHOLOGY ON C3GR` ) # "1392_1_D", "1395_1_D", "13150_2_D", "21345_1_R", "22287_1_R", "22315_1_R", "2374_1_D", "23130_1_D" 




Cleaned_Labdata_Original <- Cleaned_Labdata_Original %>%
  # 1. Group the data by Isolate_ID
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


# We have now new Isolate_IDs as some were duplicates, but we will remove one further below.
Test<- Cleaned_Labdata_Original %>% 
  select(INIKA_ID, Isolate_ID,Isolate, Isolate_ID_repeat_count, Isolate_ID_row_count,
         REGION, SEASON, `SAMPLE FROM`,`COLONY MORPHOLOGY ON BGA`, `COLONY MORPHOLOGY ON CARBA`, `COLONY MORPHOLOGY ON XLD`, `COLONY MORPHOLOGY ON C3GR` ) %>% 
  filter(Isolate_ID_repeat_count=="2")
## Found one true duplicate, remove it below
# Not Madelaine changed below as the Freezing_ID_No=="2146-1_D" should be "2156_1_D", due to mispellings, no 2146 in the selected!

Cleaned_Labdata_Original <- Cleaned_Labdata_Original%>% 
  mutate(Isolate_ID =case_when(Freezing_ID_No=="21161-1_D" ~ "21161_1_D",
                               Freezing_ID_No=="2146-1_D" ~ "2156_1_D",
                               TRUE ~ Isolate_ID,
                               ))%>% 
    filter(Isolate_ID!="22315_1_R_2")%>%
  filter(INIKA_ID !="21267") # We remove this here as well as we know this one is a duplicate in the KOBO samplingdata
  

#Cleaned_Labdata_Original <- Cleaned_Labdata_Original %>%
 # mutate(
#    Isolate_ID_Corrected = case_when(
 #     Isolate == "C.freundii" & Isolate_ID == "2374_1_D" ~ "2374_2_D",
     # TRUE ~ Isolate_ID_Corrected
  #  )
 # )


# MN Checking scripts not writing out new files here if not needed
# Note the total number of isolates remaining in the dataset after removing one duplicate is 2865!

## Saving the file as tsv to be opened in excel format
write_tsv(Cleaned_Labdata_Original, "data/CLEANED_DATA/1_Cleaned_Labdata_Original-2026-03-11.tsv")

# Save as rds to be used further work in R
saveRDS(Cleaned_Labdata_Original, "data/CLEANED_DATA/1_Cleaned_Labdata_Original-2025-03-11.rds")
#######################################################################################################

