## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(systemfonts)

# Madelaine checked, seems ok 17.10 

## Importing the original WHONET DATA
WHONET_DATA_9_9_2025 <- read_excel("data/WHONET DATA_9_9_2025.xlsx")


## Looking at the data structure ----
glimpse(WHONET_DATA_9_9_2025)
str(WHONET_DATA_9_9_2025)
names(WHONET_DATA_9_9_2025)

## Selecting the relevant columns
WHONET_Cleaned_data <-  WHONET_DATA_9_9_2025 %>%
  select(-ROW_IDX, -LABORATORY, -COUNTRY_A, -PATIENT_ID, -X_SPECIES, -INSTITUT, -SPEC_TYPE, - X_AST_DATE, -SPEC_CODE, -X_SPEC_HUM, -X_SPEC_ANI, 
         -X_SPEC_ENV, -ISOL_NUM, -ORG_TYPE, -COMMENT, -DATE_DATA, -TYL_ED30, -ORIGIN)

## Changing the numeric into character across the columns
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
  # Use mutate and across to apply the change to multiple columns
  mutate(
    # Select all columns that are currently numeric (is.numeric)
    across(where(is.numeric), ~ as.character(.))
  )

## Renaming the columns names
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
rename(Isolate_ID = `SPEC_NUM`,
       INIKA_ID = `X_INIKA_ID`,
       REGION = `X_REGION`,
       DISTRICT = `X_DISTRICT`,
       SEASON = `X_SEASON`,
       ORIGIN_OF_SAMPLE = `X_ORIGIN_T`,
       PROTOCOL = `X_PROTOCOL`)
 


# We use the one above even if there are duplication on INIKA_ID to be joined, in script 9.Join
## Duplicates checking in INIKA_ID column
WHONET_duplicates <- WHONET_Cleaned_data$INIKA_ID[duplicated(WHONET_Cleaned_data$INIKA_ID)]
print("Actual duplicated Isolate_IDs (second or subsequent appearance):")
print(WHONET_duplicates) ## "11114" "12195" "13113" "22144" "22160" "22247" "23184" "23223" "23248" "23272" "23277" "23298" "23346"
# From the duplicates, "23346" "23272" have the same isolate "S.typhimurium", while the rest duplicated with different Isolates
# "23272" have confirmed as E.coli by VITEK_MS

## Duplicates checking in Isolate_ID column
WHONET_duplicates <- WHONET_Cleaned_data$Isolate_ID[duplicated(WHONET_Cleaned_data$Isolate_ID)]
print("Actual duplicated Isolate_IDs (second or subsequent appearance):")
print(WHONET_duplicates) ## "11348_1_R"

# 5/3/2026 MMJ: Has detected that "11348_1_R" is duplicated on Isolate ID, while
#on the INIKA_ID are different i.e.11348 and 11342.
#Solution:To re write the Isolate_ID "11348_1_R" to correspond with INIKA_ID "11342" 
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
  mutate(Isolate_ID = if_else(INIKA_ID == "11342", 
                              str_glue("{INIKA_ID}_1_R"), 
                              Isolate_ID))

#MMJ 6/12/2025
# We have detected that, "22131_2_D" has mislabeled, instead of "22131_1_D" the one selected for AST
# Solution: We re write this number accordingly
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
  mutate(Isolate_ID = case_when(
    Isolate_ID == "22131_2_D" ~ "22131_1_D",
    TRUE ~ Isolate_ID
  ))

## Under here done on 18/10/2025, by Maulid

# Inserting a column "ESBL_Selection" its contents being the difference between CTC_ED30 and CTX_ED5, If the results will be >= 5 (ESBL) 
# MN Comment, here there is a misunderstanding mm zones should be below 22mm
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
   mutate(
  CTC_ED30 = as.numeric(as.character(CTC_ED30)),
  CTX_ED5 = as.numeric(as.character(CTX_ED5)),
  CRO_ED5 = as.numeric(as.character(CRO_ED30)),
  
   ESBL_Selection = CTC_ED30 - CTX_ED5
   )

# Select the true ESBL_ # MN This is to early to do here!
# You first need to select those that have mm zones for Cefotaxime  below 22mm for E.coli, and 21 mm for Klebsiella pneumoniae
# Madelaine corrected 26.11 below! Note I have new names on the variables so we need to change throughout the scripts to ensure it works

# Note I have renamed this outputfile to remian the same name! Yu don't need a new file here!

WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
  mutate(
    ESCR_ECO_presumptive = case_when(
      ORGANISM == "eco" & CTX_ED5  < 22 ~ 1,
      ORGANISM == "eco" & CRO_ED30 < 23  ~ 1,
      TRUE ~ 0
    ),
    ESCR_KPN_presumptive = case_when(
      ORGANISM == "kpn" & CTX_ED5 < 21 ~ 1,
      ORGANISM == "kpn" & CRO_ED30 < 23  ~ 1,
      TRUE ~ 0
    ),
    ESBL_Presumptivefinal = case_when(
      (ESCR_KPN_presumptive == 1 | ESCR_ECO_presumptive == 1) & ESBL_Selection >= 5 ~ 1,
      TRUE ~ 0
    )
  )


# In total 190 ESCR Eco ( Presumptive), 170 ESBL Eco ( Presumptive), 148 ESCR K.Pneumoniae (Presumptive), 99 ESBL K.Pneumoinae ( Presumoptive),



  # Note this one WHONET_Cleaned_data_ESBL  has not been printed out!- you need to do that!
       
         
  


#26.11.25 Check have not written out the files again # before the write commands
# Saving the file as tsv
#write_tsv(WHONET_Cleaned_data, "data/CLEANED_DATA/4_WHONET_Cleaned_data-2025-10-08.tsv")

# Exporting the file as rds
#saveRDS(WHONET_Cleaned_data,"data/CLEANED_DATA/4_WHONET_Cleaned_data-2025-10-08.rds")

#write.csv(WHONET_Cleaned_data_ESBL,"data/CLEANED_DATA/WHONET_Cleaned_data_ESBL.csv")
#################################################################################

