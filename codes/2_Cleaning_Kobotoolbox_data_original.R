## Making the environment to be reproducible ----
# library(renv) # https://rstudio.github.io/renv/articles/renv.html
# select -> ctrl + Enter to run
# renv::init() # This has to be run once per project

## Loading libraries ----
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)


## Importing the original KOBOTOOLBOX data
KOBOTOOLBOX_INIKA_SAMPLING_ANALYSIS_HUMAN_2025_06_12_18_24_20 <- 
  read_excel("data/KOBOTOOLBOX_INIKA_SAMPLING_ANALYSIS_HUMAN_-_2025-06-12-18-24-20.xlsx",
             na = c("","NA"))
#'Origin of sample' == "Outpatient"
#'Age (yrs)' == "45")
Check_KOBO<-KOBOTOOLBOX_INIKA_SAMPLING_ANALYSIS_HUMAN_2025_06_12_18_24_20 %>%
  filter(Region=="Kilimanjaro", District=="Moshi Rural", Season=="Wet" )

## Selecting relevant columns (removing unwanted columns) ----
Demographic_cleaned_data <-
  KOBOTOOLBOX_INIKA_SAMPLING_ANALYSIS_HUMAN_2025_06_12_18_24_20 %>%
  select(INIKA_OH_TZ_ID,`Age (yrs)`, Season, GPS_coordinates, `_GPS_coordinates_latitude`,  `_GPS_coordinates_longitude`,
         `_GPS_coordinates_precision`,Region, District, `Origin of sample`, `If school-aged children, name of the school*`,
         `Which class/grade are you?*`,`Who is your caretaker?*`,`If others, mention`,`What is your occupation and/or of your caretaker?*`,
         `Have you ever heard about AMR?`, `If yes, how did you get this information?`, `Have you or your children used any antibiotics at any time?`,
         `Have you or your children used any antibiotics in the past six months?`, `If yes, where did you get these drugs from?`,
         `If it was drug sellers or pharmacy; did you have a prescription from the doctor/prescriber?`,
         `If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?`,
         `Do you have farm animals at your home place?`,`If Yes, do you participate in taking care of animals i.e. cleaning, Feeding?`, 
         `Are your animals being treated for any diseases?`, `Do you have a poultry farm at your home?`,
         `If yes, are they for eggs or meat?`,`Do you participate in taking care of Poultry i.e. Cleaning/Feeding?`,
         `Are your Poultry being treated by veterinary doctor?`, `If Yes, does the veterinary doctor tell you about the withdrawal time?`,
         `Do you have a toilet at your home place?`, `Do you wash your hands with soap after every toilet visit?`,Gender,`Origin of sample`,
         - GPS_coordinates) %>%
  rename("Latitude" = `_GPS_coordinates_latitude`,
         "Longitude" = `_GPS_coordinates_longitude`,
         "GPS_Precision" = `_GPS_coordinates_precision`,
         INIKA_ID = INIKA_OH_TZ_ID,
         Age_yrs = `Age (yrs)`,
         GENDER = `Gender`,
         REGION = `Region`,
         DISTRICT = `District`,
         SEASON = `Season`,
         ORIGIN_OF_SAMPLE = `Origin of sample`
         ) %>% 
  # Consider keeping as double
  mutate_at(.vars = "Age_yrs", as.character) %>% 
  #MMJ
  # Ilemela district is mislabeled as Ilemala, We need to rewrite accordingly
  mutate(DISTRICT = case_when(
    DISTRICT == "Ilemala" ~ "Ilemela",
    TRUE ~ DISTRICT
  ))


# Reordering the columns


Demographic_cleaned_data <- 
  Demographic_cleaned_data %>%
  select(INIKA_ID, 
         Age_yrs,
         GENDER, 
         SEASON,,
         REGION, DISTRICT,
         ORIGIN_OF_SAMPLE,
         Latitude, Longitude, GPS_Precision,
         everything())
# MMJ 
# We have detected that Name of schools have written in different cases
# Solution: We rewrite them here under in the same case, to have homogeneity entries
Demographic_cleaned_data <-Demographic_cleaned_data %>%
  mutate(`If school-aged children, name of the school*`=case_when(
    `If school-aged children, name of the school*` =="Pasua Primary School" ~"Pasua primary school",
    `If school-aged children, name of the school*` =="Mirongo Primary school" ~"Mirongo primary school",
    `If school-aged children, name of the school*` =="Mirongo Primary School" ~"Mirongo primary school",
    `If school-aged children, name of the school*` =="Umbwe Primary School" ~"Umbwe primary school",
    TRUE~`If school-aged children, name of the school*`
  ))

Demographic_cleaned_data %>%
  count(GENDER, sort = TRUE)


# Get a logical vector: TRUE for elements that are duplicates (2nd, 3rd, etc. occurrences)
duplicates_logical <- duplicated(Demographic_cleaned_data$INIKA_ID)
#sum(duplicates_logical)

# Extract the actual duplicated ID values
duplicate_INIKA_IDs <- 
  Demographic_cleaned_data$INIKA_ID[duplicates_logical]

duplicate_INIKA_IDs
# Alternative - not better though
# Demographic_cleaned_data$INIKA_ID[
#   which(duplicated(Demographic_cleaned_data$INIKA_ID))
# ]


# View the unique duplicate ID values
unique_duplicate_INIKA_IDs <- unique(duplicate_INIKA_IDs)

# Print the result
print(unique_duplicate_INIKA_IDs)


## First check if the duplicated ID has exactly the same data
Demographic_cleaned_data %>%
  filter(INIKA_ID %in% unique_duplicate_INIKA_IDs) %>%
  arrange(INIKA_ID) %>%
  View()
# Here your duplicates INIKA IDs contain different information : 
# - can we know what is correct or not ? what should we do ? keep in mind
# if one was sent to sequencing ? 

# HERE FIX or make an explicit choice ---- 
#The INIKA_ID 21267 needs to be totally excluded as it is not possibe to know which of the persons that were sampled.

## Remove the duplicates and maintain the fist occurrence
Demographic_cleaned_data <- 
  Demographic_cleaned_data %>%
  filter(INIKA_ID !=21267)%>%
  # Select only unique rows based on the INIKA_ID column
  # .keep_all = TRUE ensures that all other columns from the first 
  # row for each unique INIKA_ID are retained.
  distinct(INIKA_ID, .keep_all = TRUE)

Check_KOBO<-Demographic_cleaned_data %>%
  filter(REGION=="Kilimanjaro", DISTRICT=="Moshi Rural", SEASON=="Wet" )


# Madelaine:_transform INIKA_ID to charachter
Demographic_cleaned_data$INIKA_ID <- as.character(Demographic_cleaned_data$INIKA_ID)



# View a preview of the resulting data frame
head(Demographic_cleaned_data)

# 2204 unique Ids remained after removing one duplicate

### Saving the file
#26-11-25 check so I don't write out new files unless finding a mistake somewhere!
#write_tsv(Demographic_cleaned_data, "data/CLEANED_DATA/2_Demographic_cleaned_data-2025-10-07.tsv")

#saveRDS(Demographic_cleaned_data, "data/CLEANED_DATA/2_Demographic_cleaned_data-2025-10-07.rds")

#write.csv(Demographic_cleaned_data, "data/CLEANED_DATA/2_Demographic_cleaned_data-2025-10-07.csv")

#######################################################################################################
# # Clean up and remove this to a separate code file call (Can have everything in one file)
#  
#  ## Importing the data set for drawing a table
#  Demographic_cleaned_data <- read_csv("data/CLEANED_DATA/2_Demographic_cleaned_data-2025-10-07.csv")
#  
#  # Define Age Bins and Create the Summary Table
#  summary_table <- Demographic_cleaned_data %>%
#    
#    # Create the Age Range column using the 'Age_yrs' column
#    mutate(`Age Range` = cut(Age_yrs,
#                             breaks = c(10, 18, 28, 38, 48, 58, 68, 78, 89), # 89 ensures 78-88 is captured
#                             labels = c('10 - 17', '18 - 27', '28 - 37', '38 - 47', 
#                                        '48 - 57', '58 - 67', '68 - 77', '78 - 88'),
#                             right = FALSE)) %>% # 'right = FALSE' makes the lower bound inclusive
#    
#    # Group by all requested demographic variables and the new age range
#    group_by(REGION, DISTRICT, GENDER, SEASON, ORIGIN_OF_SAMPLE, `Age Range`) %>%
#    
#    # Count the number of observations (rows) for each group,
#    #    as an explicit 'Total' column was not provided.
#    summarise(`Total Count` = n()) %>%
#    
#    # Remove grouping structure
#    ungroup()
#  
#  # Print the resulting table
#  print(summary_table)
#  
#  
#  # Optional: Write the resulting table to a new CSV file
# write_csv(summary_table, 
#           "Results/Demographic_Summary_Table_Final.csv")
# 
# 
# 
# # Compute the Interquartile Range (IQR) 
# # The na.rm = TRUE argument ensures that missing values are ignored in the calculation.
# iqr_age <- IQR(Demographic_cleaned_data$Age_yrs, na.rm = TRUE)
# 
# # Print the result
# print(paste("The Interquartile Range (IQR) of Age_yrs is:", iqr_age))
# 
# 
# 
# ## Install and Load Necessary Libraries 
# 
# install.packages("leaflet")
# 
# library(dplyr)
# library(tidyr)
# library(leaflet)
# 
# 
# 
# 
# #Drawing a map
# 
# map_data <- Demographic_cleaned_data %>%
#   # Rename columns 7 and 8 for clarity in the map code
#   rename(
#     latitude = 7,  # Column number 7 as Latitude
#     longitude = 8  # Column number 8 as Longitude
#   ) %>%
#   # Ensure both columns are treated as numeric
#   mutate(
#     latitude = as.numeric(latitude),
#     longitude = as.numeric(longitude)
#   ) %>%
#   # Remove any rows where conversion failed or coordinates are missing (NA)
#   filter(!is.na(latitude) & !is.na(longitude))
# 
# # 2. Draw the Interactive Map using Leaflet
# map_leaflet <- map_data %>%
#   leaflet() %>%
#   # Add a base map layer (OpenStreetMap is a good default)
#   addTiles() %>%
#   # Add markers for each GPS coordinate in your dataset
#   addMarkers(
#     lat = ~latitude,
#     lng = ~longitude,
#     # Optional: Add a popup showing the coordinates when a marker is clicked
#     popup = ~paste("Lat:", latitude, "<br>Lng:", longitude)
#   )
# 
# # 3. Display the map
# map_leaflet
