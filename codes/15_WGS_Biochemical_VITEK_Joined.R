# Loading libraries

library(tidyverse) 
library(readxl) 


# Importing the file
UniqueData <- read_csv("data/CLEANED_DATA/UniqueData.csv")

# Selecting relevant columns
UniqueData_Relevant_Columns<-UniqueData %>% 
  select(INIKA_ID.x,TVLA_ID, Age_yrs, GENDER, REGION.x,DISTRICT.x, 
         SEASON.x,ORIGIN_OF_SAMPLE.x,`COLONY MORPHOLOGY ON C3GR.x`, 
         `COLONY MORPHOLOGY ON CARBA.x`,`COLONY MORPHOLOGY ON XLD.x`, 
         `COLONY MORPHOLOGY ON BGA.x`, `TSI Media_Slope.x`, `TSI Media_Butt.x`,
         `TSI Media_Gas.x`, TSIMedia_H2S.x, `SIM Media_H2S.x`, `SIM Media Indole.x`,
         `SIM Media Motility.x`,CITRATE.x, UREASE.x, Isolate.x, VITEK_MS_Results,
         AMX_ED10,AZM_ED15, CIP_ED5, GEN_ED10,CRO_ED30, FLR_ED30,MEM_ED10, POL_ED300,
         DOX_ED30, OXY_ED30,SXT_ED1_2,CTX_ED5, CTC_ED30)
# Renaming the columns
UniqueData_Renaming_Columns<-UniqueData_Relevant_Columns %>% 
  rename(INIKA_ID = INIKA_ID.x, REGION = REGION.x, DISTRICT = DISTRICT.x,
         SEASON = SEASON.x, ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
         C3GR = `COLONY MORPHOLOGY ON C3GR.x`,
         CARBA = `COLONY MORPHOLOGY ON CARBA.x`, 
         XLD = `COLONY MORPHOLOGY ON XLD.x`,
         BGA = `COLONY MORPHOLOGY ON BGA.x`, 
         TSI_slope = `TSI Media_Slope.x`,
         TSI_Butt = `TSI Media_Butt.x`, TSI_Gas = `TSI Media_Gas.x`,
         TSI_H2S = TSIMedia_H2S.x, SIM_H2S = `SIM Media_H2S.x`,
         SIM_Indole = `SIM Media Indole.x`, SIM_Motility = `SIM Media Motility.x`,
         CITRATE = CITRATE.x, UREASE = UREASE.x, Isolate = Isolate.x)

# Importing the E.coli WGS file 
Escherichia_coli_WGS_NVI <- read_excel("data/Escherichia_coli_WGS_NVI.xlsx")

# Selecting relevant columns
Escherichia_coli_WGS_Relevant_columns <-Escherichia_coli_WGS_NVI %>% 
  select(Sample, `ReadsQC-1.0.0/SpeciesName`)

# Renaming the columns
Escherichia_coli_WGS_Renaming_columns <- Escherichia_coli_WGS_Relevant_columns %>% 
  rename(TVLA_ID = Sample, Isolate_name = `ReadsQC-1.0.0/SpeciesName`)

# Checking the duplicates
DuplicateCheck <- Escherichia_coli_WGS_Renaming_columns[
  duplicated(Escherichia_coli_WGS_Renaming_columns[ , c( "TVLA_ID")]), 
] # 0

# String replacing "-MMJ" with ""
Escherichia_coli_WGS_Renaming_columns <- Escherichia_coli_WGS_Renaming_columns %>% 
  mutate(TVLA_ID = str_replace(TVLA_ID, "-MMJ", "")) %>% 
  mutate(TVLA_ID = str_replace_all(TVLA_ID, "-", "_"))

#Importing the Klebsiella WGS file
Klebsiella_pneumoniae_WGS_NVI <- read_excel("data/Klebsiella_pneumoniae_WGS_NVI.xlsx")
# Selecting relevant columns
Klebsiella_pneumoniae_WGS_Relevant_columns <-Klebsiella_pneumoniae_WGS_NVI %>% 
  select(Sample, `ReadsQC-1.0.0/SpeciesName`)

# Renaming the columns
Klebsiella_pneumoniae_WGS_Renaming_columns <- Klebsiella_pneumoniae_WGS_Relevant_columns %>% 
  rename(TVLA_ID = Sample, Isolate_name = `ReadsQC-1.0.0/SpeciesName`)

# Checking the duplicates
DuplicateCheck <- Klebsiella_pneumoniae_WGS_Renaming_columns[
  duplicated(Klebsiella_pneumoniae_WGS_Renaming_columns[ , c( "TVLA_ID")]), 
] # 0

# String replacing "-MMJ" with ""
Klebsiella_pneumoniae_WGS_Renaming_columns <- Klebsiella_pneumoniae_WGS_Renaming_columns %>% 
  mutate(TVLA_ID = str_replace(TVLA_ID, "-MMJ", "")) %>% 
  mutate(TVLA_ID = str_replace_all(TVLA_ID, "-", "_"))

# Joining the three data sets
WGS_Biochemical_VITEK_TVLA <- UniqueData_Renaming_Columns %>% 
  left_join(Escherichia_coli_WGS_Renaming_columns, by = "TVLA_ID") %>%
  left_join(Klebsiella_pneumoniae_WGS_Renaming_columns, by = "TVLA_ID")

# Checking the unmatched
unmatched1<- anti_join(Escherichia_coli_WGS_Renaming_columns,
                      UniqueData_Renaming_Columns, by = "TVLA_ID") # "22161_1_D" "23316_1_R" "11214_1_R" "11321_1_R"
# Re writing the unmatched ID in UniqueData
UniqueData_Renaming_Columns<-UniqueData_Renaming_Columns %>% 
  mutate(TVLA_ID = case_when(TVLA_ID =="2216_1_D" ~"22161_1_D",
                   TVLA_ID =="23316_1R" ~ "23316_1_R",
                   TVLA_ID =="11214(i)_R" ~"11214_1_R",
                   TVLA_ID =="11321_(i)_R" ~ "11321_1_R"))
WGS_Biochemical_VITEK_TVLA2 <- UniqueData_Renaming_Columns %>% 
  left_join(Escherichia_coli_WGS_Renaming_columns, by = "TVLA_ID") %>%
  left_join(Klebsiella_pneumoniae_WGS_Renaming_columns, by = "TVLA_ID")

unmatched3<- anti_join(Escherichia_coli_WGS_Renaming_columns,
                       UniqueData_Renaming_Columns,
                        by = "TVLA_ID")

unmatched2<- anti_join(Klebsiella_pneumoniae_WGS_Renaming_columns, 
                       UniqueData_Renaming_Columns, by = "TVLA_ID") # 0



# Save the file
write.csv (WGS_Biochemical_VITEK_TVLA, "data/CLEANED_DATA/WGS_Biochemical_VITEK_TVLA.csv")



