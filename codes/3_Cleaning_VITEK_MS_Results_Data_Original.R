## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)

## Importing the original lab MALDITOF_RESULTS
MALDITOF_RESULTS <- read_excel("data/MALDITOF_RESULTS.xlsx")


# Data cleaning ----
## Looking at the data structure ----
glimpse(MALDITOF_RESULTS)
MALDITOF_RESULTS[, "Isolate ID"]
str(MALDITOF_RESULTS)
names(MALDITOF_RESULTS)

## Selecting relevant columns (removing unwanted columns) ----
# Madelaine_ need to ensure we join with the correct and then we need the biochemical results here as well

MALDITOF_RESULTS_cleaned_data <-
  MALDITOF_RESULTS %>%
  select("Isolate ID" , 'Biochemica Id by Maulid',"VITEK MS  Results")
#MMJ # I would suggest that 'Biochemical Id by Maulid' should not be selected, since its data here are not correct, rather we will join the file with correct data

## Renaming the column and name it TVLA_ID
MALDITOF_RESULTS_cleaned_data <-
  MALDITOF_RESULTS_cleaned_data %>%
  rename(VITEK_MS_Results = `VITEK MS  Results`,
         TVLA_ID = `Isolate ID`)




## Relocate the TVLA_ID, not necessary as we only have one TVLA_ID now!
#MALDITOF_RESULTS_cleaned_data <- MALDITOF_RESULTS_cleaned_data %>%
#relocate(TVLA_ID, .after = TVLA_ID)


## Replacing "-" with "_" in TVLA_ID column
MALDITOF_RESULTS_cleaned_data$TVLA_ID <- gsub("-", "_", MALDITOF_RESULTS_cleaned_data$TVLA_ID)



## Re writing the entries in TVLA_ID to match with that in Original Lab data, 
#seems to have been some mispellings in the data from TVLA,
#I suggest to call this corrected_TVLA_ID instead to not loose TVLA IDs, MN found mistakes in the file that have been joined in ). Join,
#so need to correct below 2374_1_D, 11301_1_R and 13122_1_D 
MALDITOF_RESULTS_cleaned_data <-
  MALDITOF_RESULTS_cleaned_data%>%
 # some Isolates removed due to duplications of ID 
filter(!(TVLA_ID == "23372_1_R" & `Biochemica Id by Maulid` == "Klebsiella pneumoniae"))%>%
filter(!(TVLA_ID == "23224_1_R" & VITEK_MS_Results  == "No Identification" ))%>%
# Correcting misspellings and misidentification as far as possible
  mutate(Isolate_ID = case_when(
    TVLA_ID == "1152_1_D G.N" ~ "1152_1_D",
    TVLA_ID== "11321_(i)_R" ~ "11321_1_R",
    TVLA_ID== "11214(i)_R"  ~ "11214_1_R",
    TVLA_ID== "11253_1_R"  ~ "11252_1_R",
    TVLA_ID== "11262(i)_R"  ~ "11262_1_R",
    TVLA_ID== "11348(i)_R"  ~ "11348_1_R",
    TVLA_ID== "122_1_D"  ~ "125_1_D",
    TVLA_ID== "1268_1_R"  ~ "1268_1_D",
    TVLA_ID== "12132(i)_D"  ~ "12132_1_D",
    TVLA_ID== "21337_(i)_R"  ~ "21337_1_R",
    TVLA_ID== "23124_1_R"  ~ "23124_1_D",
    TVLA_ID== "25256_2_R"  ~ "23256_2_R",
    TVLA_ID== "23316_1R"  ~ "23316_1_R",
    TVLA_ID== "23237_R_"  ~ "23237_1_R",
    TVLA_ID== "13235_2_R"  ~ "13235_1_R",
    TVLA_ID== "12310_2_R"  ~ "12310_1_R",
    TVLA_ID== "12285_1_R_"  ~ "12285_1_R",
    TVLA_ID== "13122_1_D"  ~ "13123_1_D", # MN corrected it has to be in this direction
    TVLA_ID== "2337_1_R"  ~ "23370_1_R",
    TVLA_ID== "2216_1_D"  ~ "22161_1_D",
    TVLA_ID== "2369_1_D"  ~ "2366_1_D",
    TVLA_ID== "23347_2_"  ~ "23347_2_R",
    TVLA_ID== "23220_1_C"  ~ "23220_1_R",
    TVLA_ID== "23282_1_"  ~ "23282_1_R",
    TVLA_ID== "23308_1_R"  ~ "23308_2_R",
    TVLA_ID== "23320_2_"  ~ "23320_2_R",
    TVLA_ID== "21220_(i)_R"  ~ "21220_1_R",
    TVLA_ID== "21209_(i)_R"  ~ "21209_1_R",
    TVLA_ID== "2111_(i)_D"  ~ "2111_1_D",
    TVLA_ID== "12195(II)_R"  ~ "12195_2_R",
    TVLA_ID== "12195(i)_R"  ~ "12195_1_R",
    TVLA_ID== "11353_(i)_R"  ~ "11353_1_R",
    TVLA_ID== "11350(ii)"  ~ "11350_2_R",
    TVLA_ID== "11359(i)_R"  ~ "11359_1_R",
    TVLA_ID== "11308(i)_R"  ~ "11308_1_R",
    TVLA_ID== "11191(ii)_R"  ~ "11191_2_R",
    TVLA_ID== "C.1197_1_DMV"  ~ "1197_1_D",
    TVLA_ID== "11154mv"  ~ "11154_1_D",
    TVLA_ID== "1164_1_DMV"  ~ "1164_1_D",
    TVLA_ID== "12125_1_D"  ~ "12125_2_D",
    TVLA_ID== "2291_1_D"  ~ "2291_2_D",
    TVLA_ID== "1166_1_R"  ~ "1166_1_D",
    TVLA_ID== "12125_1_D"  ~ "12125_2_D",
    TVLA_ID== "23142_1_R"  ~ "23142_1_D",
    TVLA_ID== "23131-1-D"  ~ "23131_1_D",
    TVLA_ID== "12311_1_R"  ~ "12311_2_R",
    TVLA_ID== "2234_2_R"  ~ "22347_2_R",
    TVLA_ID== "2236_2_R"  ~ "22362_2_R",
    TVLA_ID== "212348_1_R" ~ "22348_2_R",
    TVLA_ID == "23308_2" ~ "23308_2_R_2",
    TVLA_ID == "12203_2"  ~ "12203_2_R",
    TVLA_ID == "11301-1-R"  ~ "11301_2_R", # MN new correction
    TVLA_ID == "2374-1-D"  ~ "2374_1_D", # MN new correction
    TVLA_ID == "2220_1_R"  ~ "22204_1_R", # MMJ new correction 23/10/25 and below
    TVLA_ID == "232606_2_R"  ~ "23260_2_R",  
    TVLA_ID == "11301_1_R"  ~ "11301_2_R",
    TVLA_ID == "21161_1_D"  ~ "21161_1_D",
    TVLA_ID == "21185_1_R"  ~ "21185_2_R",   
    TVLA_ID == "21294_1_R"  ~ "21294_2_R",
    TVLA_ID == "23286_2_R"  ~ "23286_1_R",
    TVLA_ID == "12300_1_R"  ~ "12300_2_R",
    TVLA_ID == "12324_11_R"  ~ "12324_1_R",
    TVLA_ID == "22285_1_D"  ~ "22285_1_R",
    TVLA_ID == "22131_2_D"  ~ "22131_1_D", # MN Note why was this one corrected it is numbered 22131_2_D in the WHONETdata! # But if this is the case we should keep it as it is!# MMJ; There is typing error, This number two isolates were isolated, Isolate 1 was Kpn and was  selected for AST, Isolate 2 was E.coli which was not selected
    TVLA_ID == "2223_2_R"  ~ "22231_2_R",
    TVLA_ID == "22145_1_D"  ~ "22145_2_D",
    TVLA_ID == "22845_1_R"  ~ "22345_1_R",
    TVLA_ID == "12336_1_R"  ~ "12336_2_R",
    TVLA_ID =="21335_1_R" ~ "21235_1_R", # MN attached some more, correct?
    
    TRUE ~ TVLA_ID  # This is the "else" condition: keep any other values as they are
  ))%>%
  # Remove first part of Isolate_ID to use as the INIKA_ID for possible joining by two columns
  mutate(
    # Matches the first underscore (_) followed by any character (.) 
    # zero or more times (*), and removes the entire match.
    INIKA_ID = str_remove(Isolate_ID, "_.*")
  )

  MALDITOF_RESULTS_cleaned_data<-unique(MALDITOF_RESULTS_cleaned_data)
  ##############################################################################
  ## Importing the MALDITOF results from NVI ## 02.03.2026 MMJ
  # Note here you can not import the first 3 rows!
  NVI_MALDITOF_results <- read_excel("data/MALDITOF_NVI 25feb2026.xlsx", skip=3)
  str(NVI_MALDITOF_results)
  # Selecting relevant columns and rows
  #NVI_MALDITOF_results <- NVI_MALDITOF_results %>%
  #  select(- `Inika, E. coli-stammer mottatt 24/2-26 tmf`,- ...4) %>%
  #  slice(-c(1, 2)) %>%
  #  rename(
   #   Isolatnr = ...2,
    #  Maldi_Tof_25_2_26_tmf = ...3) %>%
    #slice(-c(1))
  ## Inserting the new column as 'Isolate_ID" its contents being as in "Isolatnr"
  #Madelaine:I changed and just renamed the columns
  NVI_MALDITOF_results <- NVI_MALDITOF_results %>%
    rename(Isolate='Maldi-Tof. 25/2-26 tmf',Isolate_ID = Isolatnr) %>%
   mutate(`Isolate_ID` = str_replace_all(`Isolate_ID`, "-", "_"))
# Combine the MALDITOF_RESULTS_cleaned_data and NVI_MALDITOF_results
  combined_MALDITOF_Results <- full_join(MALDITOF_RESULTS_cleaned_data, NVI_MALDITOF_results, by = c("Isolate_ID"))
 
  
  # 05.03.26_ I have not run below as I did not run the previous 1 and 2 scripts
# Joining all MALDITOF results and those selected for AST
  JoinedMALDITOF_Selected<-left_join(combined_MALDITOF_Results,DRY_RAINY_MWZ_KILIMANJARO_AST2,by = "Isolate_ID")
  notjoinedMALDITOF<-anti_join(combined_MALDITOF_Results,DRY_RAINY_MWZ_KILIMANJARO_AST2,by = "Isolate_ID")
  
  

# Joining and checking until we managed to identify almost all that were not joined, used the helpfile ( where results was imported from the files below) in excel for checking

#JoinedTVLA_Selected<-left_join(MALDITOF_RESULTS_cleaned_data,DRY_RAINY_MWZ_KILIMANJARO_AST2,by = "Isolate_ID")
#notjoinedTVLA<-anti_join(MALDITOF_RESULTS_cleaned_data,DRY_RAINY_MWZ_KILIMANJARO_AST2,by = "Isolate_ID")
# 5 not joined # Now are 4. on 6/12/2025. by MMJ

# "1192_1_D" Was not selected and tested for AST
# "12336_1_D" was not selected and not tested for AST
# "2395_2_D" This ID has only one isolate which is "2395_1_D" and was not selected and tested for AST
# "23308_2_R_2" This ID has duplicates, since "23308_2_R" is the only selected and tested for AST
# All of the above four ID's should not be considered for further analyses  
  DuplicationsMALDITOF<-JoinedMALDITOF_Selected$Isolate_ID[duplicated(JoinedMALDITOF_Selected$Isolate_ID)] 
  Duplications<-DuplicationsMALDITOF
  view(Duplications)
#DuplicationsTVLA<-JoinedTVLA_Selected$Isolate_ID[duplicated(JoinedTVLA_Selected$Isolate_ID)] 
#Duplications<-DuplicationsTVLA
#view(Duplications)

Selection<-DRY_RAINY_MWZ_KILIMANJARO_AST2
#JoinedSelected_TVLA<-left_join(DRY_RAINY_MWZ_KILIMANJARO_AST2,MALDITOF_RESULTS_cleaned_data,by = "Isolate_ID")
#FullJoinedSelected_TVLA<-full_join(DRY_RAINY_MWZ_KILIMANJARO_AST2,MALDITOF_RESULTS_cleaned_data,by = "Isolate_ID")%>%
#  select(INIKA_ID.x,TVLA_ID, Isolate_ID, Freezing_ID_No, 'Biochemica Id by Maulid',Isolate,VITEK_MS_Results )
#notjoinedSelected<-anti_join(DRY_RAINY_MWZ_KILIMANJARO_AST2,MALDITOF_RESULTS_cleaned_data,by = "Isolate_ID")
# 5 not joined

MALDITOF <- combined_MALDITOF_Results
names(MALDITOF)

#VITEK<-MALDITOF_RESULTS_cleaned_data
#names(VITEK)

#JoinedSelectedTVLA <- left_join(Selection,VITEK, 
 #                               by = "Isolate_ID")
#noT<-anti_join(Selection,VITEK, 
 #              by= "Isolate_ID")

#noTvs<-anti_join(VITEK,Selection,
  #               by= "Isolate_ID")

#write.csv(noTvs, "data/CLEANED_DATA/notjoined_fromVITEKbacktoSelected23.10.25.csv")

#write.csv(noT, "data/CLEANED_DATA/notjoined_fromSelectedinVITEK23.10.25.csv")
  

## Duplicates checking 
MALDITOF_RESULTS_duplicates <- combined_MALDITOF_Results$Isolate_ID[duplicated(combined_MALDITOF_Results$Isolate_ID)]
print("Actual duplicated Isolate_IDs (second or subsequent appearance):")
print(MALDITOF_RESULTS_duplicates) # "125_1_D"   "21228_1_R" "23248_2_R" "13235_1_R""23370_1_R"

# MALDITOF_RESULTS_duplicates <- MALDITOF_RESULTS_cleaned_data$Isolate_ID[duplicated(MALDITOF_RESULTS_cleaned_data$Isolate_ID)]
# print("Actual duplicated Isolate_IDs (second or subsequent appearance):")
# print(MALDITOF_RESULTS_duplicates) # "125_1_D"   "21228_1_R" "23248_2_R" "23277_2_R" "23308_2_R" "13235_1_R" "23372_1_R" "23370_1_R"
#[9] "23224_1_R"

#Need to remove duplicates if they are true duplicates in combined_MALDITOF_Results
combined_MALDITOF_Results <- unique(combined_MALDITOF_Results)%>%
filter( TVLA_ID!="125_1_D",
        TVLA_ID!="13235_2_R")



# Below we created new IDs for the duplicates so it will be possible to distinguish them, need to be careful about how these can be joined!
# as there will not be another with exactly the same ID in any of the original files, needs to be assigned in those files!

combined_MALDITOF_Results<- combined_MALDITOF_Results %>%
  # 1. Group the data by Isolate_ID
  group_by(Isolate_ID) %>%
  # 2. Count the number of times each Isolate_ID appears
  mutate(Isolate_ID_repeat_count = n()) %>%
  # 3. Assign a sequential count (1, 2, 3...) within each TVLA_ID group
  mutate(Isolate_ID_row_count = row_number()) %>%
  # 4. Remove the grouping
  ungroup() %>%
  # 5. Modify the Isolate_ID based on the counts
  mutate(Isolate_ID = case_when(
    # The condition: If Isolate_ID appeared exactly 2 times AND it's the second row
    Isolate_ID_repeat_count == 2 & VITEK_MS_Results != "Escherichia coli" ~
      # Action: Append the row count (2) to the Isolate_ID
      paste(Isolate_ID, Isolate_ID_row_count, sep = "_"),
    
    # The default: Otherwise, keep the original TVLA_ID
    TRUE ~ Isolate_ID
  ))


# MALDITOF_RESULTS_cleaned_data<- MALDITOF_RESULTS_cleaned_data %>%
#   # 1. Group the data by Isolate_ID
#   group_by(Isolate_ID) %>%
#   # 2. Count the number of times each Isolate_ID appears
#   mutate(Isolate_ID_repeat_count = n()) %>%
#   # 3. Assign a sequential count (1, 2, 3...) within each TVLA_ID group
#   mutate(Isolate_ID_row_count = row_number()) %>%
#   # 4. Remove the grouping
#   ungroup() %>%
#   # 5. Modify the Isolate_ID based on the counts
#   mutate(Isolate_ID = case_when(
#     # The condition: If Isolate_ID appeared exactly 2 times AND it's the second row
#     Isolate_ID_repeat_count == 2 & VITEK_MS_Results != "Escherichia coli" ~
#       # Action: Append the row count (2) to the Isolate_ID
#       paste(Isolate_ID, Isolate_ID_row_count, sep = "_"),
#     
#     # The default: Otherwise, keep the original TVLA_ID
#     TRUE ~ Isolate_ID
#   ))

# Checking for the isolates that remains in the dataset after removing the duplicates, we could have removed the whole observation, but some of these might hav ebeen selected for further confirmation at NVI
Test<-combined_MALDITOF_Results %>%
  filter(Isolate_ID%in% c("21228_1_R", "23248_2_R", "23370_1_R"))

TestINLabdata<-DRY_RAINY_MWZ_KILIMANJARO_AST2%>%
  filter(Isolate_ID%in% c( "21228_1_R", "23248_2_R", "23370_1_R" ))

TESTjoin<-left_join(Test,TestINLabdata, by="Isolate_ID")%>%
  select("INIKA_ID.x","TVLA_ID", "Isolate_ID", "Biochemica Id by Maulid", "Isolate.x", "Isolate.y", "VITEK_MS_Results")

Test_duplicatesinleftjoined<-JoinedMALDITOF_Selected%>%
  filter(Isolate_ID%in% c( "21228_1_R", "23248_2_R", "23370_1_R" ))

# Test<-MALDITOF_RESULTS_cleaned_data %>%
#   filter(Isolate_ID%in% c( "21228_1_R" ,"23248_2_R" ,"23277_2_R", "23308_2_R", "13235_1_R" ,"23372_1_R","23224_1_R","2337,0_1_R"  ))
# 
# TestINLabdata<-DRY_RAINY_MWZ_KILIMANJARO_AST2%>%
#   filter(Isolate_ID%in% c( "21228_1_R" ,"23248_2_R" ,"23277_2_R", "23308_2_R", "13235_1_R" ,"23372_1_R","23224_1_R","2337,0_1_R"  ))
# 
# TESTjoin<-left_join(Test,TestINLabdata, by="Isolate_ID")%>%
#   select(INIKA_ID.x, TVLA_ID, Isolate_ID, 'Biochemica Id by Maulid', Isolate)
# 
# Test_duplicatesinleftjoined<-JoinedTVLA_Selected%>%
#   filter(Isolate_ID%in% c( "21228_1_R" ,"23248_2_R" ,"23277_2_R", "23308_2_R", "13235_1_R" ,"23372_1_R","23224_1_R","2337,0_1_R"  ))
########################22131_2_D 

# NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Note for one Klebsiella pneumoniae the TVLA_ID = 23308_1_R and the created Isolate_ID = 23308_2_R
# Note for one Escherichia coli the TVLA_ID = 13235_2_R and the created Isolate_ID= 13235_1_R
# Need to check that we manage to join these two isolates properly! ( EVt change the Isolate_ID!)


## Creating new column named ESC indicating the isolates for AST
combined_MALDITOF_Results <- combined_MALDITOF_Results %>%
  mutate(
    ESC = if_else(
      grepl("coli|klebsiella", VITEK_MS_Results, ignore.case = TRUE) | 
        grepl("coli|klebsiella", Isolate, ignore.case = TRUE),
      1,
      0,
      missing = 0 # Ensures NAs are marked as 0 instead of staying NA
    )
  )


# MALDITOF_RESULTS_cleaned_data <- MALDITOF_RESULTS_cleaned_data %>%
#   mutate(
#     ESC = if_else(
#       VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae"),
#       1,
#       0
#     )
#   )

# 26.11.25 checked not writing out any new tables if all correct
#write.csv(MALDITOF_RESULTS_cleaned_data, "data/CLEANED_DATA/MALDITOF_RESULTS_cleaned_data.csv")
##############################################################################################
# 3o/10/2025 Corrected one isolate_ID above has not yet run below
# Selecting confirmed Isolate as ESC for AST

AST_ISOLATES <- combined_MALDITOF_Results %>% 
  filter(ESC == 1)

# AST_ISOLATES <- MALDITOF_RESULTS_cleaned_data %>% 
#   filter(ESC == 1)

## Duplicates checking
AST_ISOLATES_duplicates <- AST_ISOLATES$Isolate_ID[duplicated(AST_ISOLATES$Isolate_ID)]
print("Actual duplicated TVLA_IDs (second or subsequent appearance):")
print(AST_ISOLATES_duplicates) ## 0

### Removing the duplicates and maintain the first occurrence, not necessary now, problem solved above!
#AST_ISOLATES <- AST_ISOLATES %>%
 # distinct(Isolate_ID, .keep_all = TRUE)

## Saving the file for further analyses
# Save as tsv file
#write_tsv(AST_ISOLATES, "data/CLEANED_DATA/3_AST_ISOLATES-2025-10-08.tsv")

# Export as rds file
#saveRDS(AST_ISOLATES, "data/CLEANED_DATA/3_AST_ISOLATES-2025-10-08.rds")
###############################################################################