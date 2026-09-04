## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)

# Madelaine check 17.10 written down
# Madelaine_New check 30.11.25
# Use the new AST_Isolates , Have # the one below to not import the previous dataset

## Importing the original WHONET DATA; !!!!Note this are not all only the one kept from Cleandata original- We should join the total instead?
#AST_ISOLATES <- readRDS("data/CLEANED_DATA/3_AST_ISOLATES-2025-10-08.rds") %>%
  #filter(!Isolate_ID %in% c("13122_1_D", "2395_2_D")) 
  
## Importing the original WHONET DATA- 1.12.25 I don't import now, use the ones I have cleaned in the previous script 4._
#WHONET_Cleaned_data <- readRDS("data/CLEANED_DATA/4_WHONET_Cleaned_data-2025-10-08.rds")

check<-WHONET_Cleaned_data%>%
  filter(Isolate_ID %in% c("13122_1_D", "2395_2_D"))

# These isolate_Ids needs to be removed as they are not in the TVLA box 
# or selected for sequencing at TVLA, Must be some writing mistake at TVLA
# removed below
#AST_ISOLATES <- AST_ISOLATES%>%
#filter(!Isolate_ID %in% c("13122_1_D", "2395_2_D")) 



## Join the two data sets
AST_RESULTS <- left_join(AST_ISOLATES, WHONET_Cleaned_data, by ="Isolate_ID")

check<-AST_RESULTS%>%
  filter(Isolate_ID %in% c("13122_1_D", "2395_2_D"))
# checking for not joining
not_joinedAST_RESULTS<-anti_join(AST_ISOLATES, WHONET_Cleaned_data, 
                                 by ="Isolate_ID") # "22206_1_R", "2395_2_D"

# 5/3/2026 MMJ:Detected that, Isolate_ID "222o6_1_R" has typo error i.e. "o" instead of "0" 
# Solution: To re write the Isolate_ID "222o6_1_R" to numeric instead of letter
WHONET_Cleaned_data <- WHONET_Cleaned_data %>%
  mutate(Isolate_ID = str_replace(Isolate_ID, fixed("222o6_1_R"), "22206_1_R"))

AST_RESULTS <- left_join(AST_ISOLATES, WHONET_Cleaned_data, by ="Isolate_ID")

not_joinedAST_RESULTS<-anti_join(AST_ISOLATES, WHONET_Cleaned_data, 
                                 by ="Isolate_ID") # "2395_2_D"
#The isolate 2395_2_D is not joined above and has to be a spelling misstake by TVLA as this isoalte_ID has not been selected for confirmation

## Duplicates checking in Isolate_ID column


AST_RESULTS_duplicates <- AST_RESULTS$Isolate_ID[duplicated(AST_RESULTS$Isolate_ID)]
print("Actual duplicated Isolate_IDs (second or subsequent appearance):")
print(AST_RESULTS_duplicates)

# 26.11.25 NOTE! This one is not joined above ! 22131_1_D ,
#The TVLA_ID = 22131_2_D need to be changed as this one has the same ID as the TVLA_ID in the WHONET file! Changes needs to be done in the 3.Cleaning....

## "23184" and "22160" have different isolates from the same ID, so keep them all

## "13122_1_D" and "2395_2_D" were not tested for AST, remove them 

## # Create a vector of IDs to remove
# This has now already been performed
#IDs_to_remove <- c("13122_1_D", "2395_2_D")

# Filter the dataset to keep all rows EXCEPT those in the vector
#AST_RESULTS <- AST_RESULTS %>%
#  filter(!Isolate_ID %in% IDs_to_remove)

## Saving the file for further use
write_tsv(AST_RESULTS, "data/CLEANED_DATA/FINAL_AST_RESULTS.tsv")
saveRDS(AST_RESULTS, "data/CLEANED_DATA/FINAL_AST_RESULTS.rds")



## Counting ESBL
ESBL_ECO <- AST_RESULTS%>%
  filter(VITEK_MS_Results == "Escherichia coli", ESBL_Presumptivefinal == 1) %>%
   count() ## 132- Note this changed from 135 to 132 when including more critera!


ESBL_Kleb <- AST_RESULTS%>%
  filter(VITEK_MS_Results == "Klebsiella pneumoniae", ESBL_Presumptivefinal == 1) %>%
  count() ## 54 




## Draw a pivot table

# ----------------------------------------------------------------------
#  Filter the data for the two target organisms and count them
# ----------------------------------------------------------------------
organism_counts <- AST_RESULTS %>%
  # Filter for the VITEK_MS_Results of interest
  filter(VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae")) %>%
  
  # Group the data by the geographic and temporal columns
  group_by(REGION, DISTRICT, SEASON, VITEK_MS_Results, ORIGIN_OF_SAMPLE) %>%
  
  # Count the number of rows (occurrences) in each group
  summarise(Count = n(), .groups = 'drop')




organism_counts_Region <- AST_RESULTS %>%
  # Filter for the VITEK_MS_Results of interest
  filter(VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae")) %>%
  
  # Group the data by the geographic and temporal columns
  group_by(REGION, VITEK_MS_Results) %>%
  
  # Count the number of rows (occurrences) in each group
  summarise(Count = n(), .groups = 'drop')

organism_counts_Origin <- AST_RESULTS %>%
  # Filter for the VITEK_MS_Results of interest
  filter(VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae")) %>%
  
  # Group the data by the geographic and temporal columns
  group_by(ORIGIN_OF_SAMPLE, VITEK_MS_Results) %>%
  
  # Count the number of rows (occurrences) in each group
  summarise(Count = n(), .groups = 'drop')



# ----------------------------------------------------------------------
#  Reshape the data to have organism names as columns
# ----------------------------------------------------------------------
Isolates_Count_table <- organism_counts %>%
  pivot_wider(
    names_from = VITEK_MS_Results,
    values_from = Count,
    values_fill = 0, # Fill any combination with no occurrences with 0
    names_prefix = "Count_"
  )

 #print(Isolates_Count_table)

## Save the file NB! 1.12.25 MN_Not printed out yet!

write_tsv(Isolates_Count_table, "Results/Isolates_Count_table.tsv")
saveRDS(Isolates_Count_table, "Results/Isolates_Count_table.rds")

