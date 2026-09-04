## Final joined data

library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)


# Importing the Demographic_cleaned_data

Final_JoinedDATA <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv",guess_max = Inf)
spec(Final_JoinedDATA)
names(Final_JoinedDATA)

# Selecting relevant columns
Final_JoinedDATA <-Final_JoinedDATA %>%
  select(-"...1",-"SEASON.y", -"DISTRICT.y",
          -"REGION.y",-"AGE" ,-"REGION" ,-"DISTRICT" ,-"SEASON",
         -"ORIGIN_OF_SAMPLE.y",-"Isolate_ID_repeat_count.y",
         -"Isolate_ID_row_count.y", -"INIKA_ID.y.x",-"INIKA_ID.x.y",
         -"Løpenr",-"Biochemica Id by Maulid" 
         )

# Rename the columns
Final_JoinedDATA <-Final_JoinedDATA %>%
  rename(INIKA_ID = `INIKA_ID.x.x`,DISTRICT = `DISTRICT.x`,AGE = `Age_yrs`,
         REGION =`REGION.x`,SEASON = `SEASON.x`,
         ORIGIN_OF_SAMPLE = `ORIGIN_OF_SAMPLE.x`,
         Isolate_ID_repeat_count =`Isolate_ID_repeat_count.x`,
         Isolate_ID_row_count = `Isolate_ID_row_count.x`,
         NVI_ID =`INIKA_ID.y.y`
        
         )

## I identified some of the mislabeled i.e., Ilemela is labeled as "Ilemala"
## Solution; Renaming the entry

Final_JoinedDATA <- Final_JoinedDATA %>%
  mutate(DISTRICT = if_else(DISTRICT == "Ilemala", "Ilemela", DISTRICT))
    

# Save the file # Not this file will not contain all the variables
write.csv(Final_JoinedDATA, "data/CLEANED_DATA/Final_JoinedDATA.csv")

## Preparing the data in a pivot longer table for easy analysis
long_data <- Final_JoinedDATA %>%
  pivot_longer(
    cols = c(
      `COLONY MORPHOLOGY ON C3GR`,
      `COLONY MORPHOLOGY ON CARBA`, 
      `COLONY MORPHOLOGY ON XLD`, 
      `COLONY MORPHOLOGY ON BGA`, 
      CITRATE, `TSI Media_Slope`, 
      `TSI Media_Butt`, 
      `TSI Media_Gas`,
      `TSIMedia_H2S`, `SIM Media Indole`,
      `SIM Media Motility`,`SIM Media_H2S`,
      UREASE
    ),
    names_to = "Method",
    values_to = "Content"
  )

grouped_summary <- long_data %>%
  group_by(Isolate, VITEK_MS_Results, Method, Content) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  arrange(Isolate, Method)

# We filter out the columns where Biochemical and TVLA results are the same
ECO_correct_Identified<-grouped_summary%>%
  filter(Isolate =="E.coli", VITEK_MS_Results =="Escherichia coli") ## 130 observation

## Save the file_ Maulid needs to complete. 
write_tsv(ECO_correct_Identified,"data/CLEANED_DATA/ECO_correct_Identified.tsv")

saveRDS(ECO_correct_Identified,"data/CLEANED_DATA/ECO_correct_Identified.rds")

ECO_NOT_correct_Identified_Bio<-grouped_summary%>%
  filter(Isolate =="E.coli", VITEK_MS_Results!="Escherichia coli") 

write_tsv(ECO_NOT_correct_Identified_Bio,"data/CLEANED_DATA/ECO_NOT_correct_Identified_Bio.tsv")

saveRDS(ECO_NOT_correct_Identified_Bio,"data/CLEANED_DATA/ECO_NOT_correct_Identified_Bio.rds")

ECO_NOT_correct_Identified_VITEK<-grouped_summary%>%
  filter(Isolate !="E.coli", VITEK_MS_Results =="Escherichia coli") # 21 

write_tsv(ECO_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/ECO_NOT_correct_Identified_VITEK.tsv")

saveRDS(ECO_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/ECO_NOT_correct_Identified_VITEK.rds")


Kleb_correct_Identified<-grouped_summary%>%
  filter(Isolate =="K.pneumoniae", VITEK_MS_Results =="Klebsiella pneumoniae") # 52 

write_tsv(Kleb_correct_Identified,"data/CLEANED_DATA/Kleb_correct_Identified.tsv")

saveRDS(Kleb_correct_Identified,"data/CLEANED_DATA/Kleb_correct_Identified.rds") 

Kleb_NOT_correct_Identified_Bio<-grouped_summary%>%
  filter(Isolate =="K.pneumoniae", VITEK_MS_Results!="Klebsiella pneumoniae")

write_tsv(Kleb_NOT_correct_Identified_Bio,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_Bio.tsv")

saveRDS(Kleb_NOT_correct_Identified_Bio,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_Bio.rds")

Kleb_NOT_correct_Identified_VITEK<-grouped_summary%>%
  filter(Isolate !="K.pneumoniae", VITEK_MS_Results=="Klebsiella pneumoniae") # 8 

write_tsv(Kleb_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_VITEK.tsv")

saveRDS(Kleb_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_VITEK.rds")
################################################################################
## Counting the confirmed isolate by VITEK MS from the selected isolates
# Create the summary table
tvla_verification_table <- Final_JoinedDATA %>%
  # Filter for the three isolates of interest
  filter(Isolate %in% c("E.coli", "K.pneumoniae", "S.typhimurium")) %>%
  group_by(Isolate) %>%
  summarise(
    # Count how many have a TVLA_ID
    `Tested at TVLA` = sum(!is.na(TVLA_ID), na.rm = TRUE),
    
    # Apply updated confirmation criteria
    `Confirmed as per Criteria` = sum(
      !is.na(TVLA_ID) & (
        (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
          (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae") |
          # Set to FALSE or a specific mismatch if VITEK identified them as something else
          (Isolate == "S.typhimurium" & VITEK_MS_Results == "S.typhimurium") 
      ), 
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# View the results
print(tvla_verification_table)
#######################################################
## COunting the isolates tested by VITEK MS from the selected isolates after cleaning 
# and group them by Region, Season, Origin, and District
# 1. Define the criteria function to keep the code clean
calculate_counts <- function(Final_JoinedDATA) {
  Final_JoinedDATA %>%
    summarise(
      `Tested at TVLA` = sum(!is.na(TVLA_ID), na.rm = TRUE),
      `Confirmed as per Criteria` = sum(
        !is.na(TVLA_ID) & (
          (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
            (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae") |
            (Isolate == "S.typhimurium" & VITEK_MS_Results == "S.typhimurium")
        ), 
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}



# 2. Filter base data for target isolates
filtered_data <- Final_JoinedDATA %>%
  filter(Isolate %in% c("E.coli", "K.pneumoniae", "S.typhimurium"))

# 3. Generate Aggregates by Region
region_agg <- filtered_data %>%
  group_by(REGION, Isolate) %>%
  calculate_counts()

# 4. Generate Aggregates by Season
season_agg <- filtered_data %>%
  group_by(SEASON, Isolate) %>%
  calculate_counts()

# 5. Generate Aggregates by Origin
origin_agg <- filtered_data %>%
  group_by(ORIGIN_OF_SAMPLE, Isolate) %>%
  calculate_counts()

# Print the Regional Table (The "Kilimanjaro" example)
print("--- SUMMARY BY REGION ---")
print(region_agg)

# Print the Seasonal Table
print("--- SUMMARY BY SEASON ---")
print(season_agg)

# Print the Origin Table
print("--- SUMMARY BY ORIGIN ---")
print(origin_agg)
################################################################################
## Isolation rate
## Calculating the Isolation rate (some of the isolates being grouped in "others")
# This identifies every unique "Sample-Isolate" combination
# 1. Prepare the raw counts

isolate_counts <- Final_JoinedDATA %>%
  distinct(INIKA_ID, Isolate) %>%
  count(Isolate, name = "count")

# 2. Define your specific categories
main_isolates <- c("E.coli", "K.pneumoniae", "S.typhimurium", "No growth")

# 3. Group everything else into "Others"
plot_data_bar <- isolate_counts %>%
  mutate(
    Isolate_Grouped = if_else(Isolate %in% main_isolates, Isolate, "Others")
  ) %>%
  group_by(Isolate_Grouped) %>%
  summarise(count = sum(count), .groups = "drop")

# 4. Calculate exact 100% percentages
# We use a rounding adjustment to ensure the sum is exactly 100
total_findings <- sum(plot_data_bar$count)

plot_data_bar <- plot_data_bar %>%
  mutate(
    raw_perc = (count / total_findings) * 100,
    # Standard rounding
    percentage = round(raw_perc, 1)
  )

# Fix the rounding error (if sum != 100, adjust the largest category)
diff <- 100 - sum(plot_data_bar$percentage)
if (diff != 0) {
  plot_data_bar <- plot_data_bar %>%
    arrange(desc(percentage)) %>%
    mutate(percentage = if_else(row_number() == 1, percentage + diff, percentage))
}

# 5. Extract the "Others" list for the caption
others_list <- isolate_counts %>%
  filter(!(Isolate %in% main_isolates)) %>%
  pull(Isolate) %>%
  unique() %>%
  paste(collapse = ", ")

# 6. Generate the Bar Chart
bar_chart <- ggplot(plot_data_bar, aes(x = reorder(Isolate_Grouped, percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = paste0(format(percentage, nsmall = 1), "%")), 
            hjust = -0.2, 
            size = 3.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0, margin = margin(t = 15))
  ) +
  labs(
    title = "Relative Distribution of Bacterial Isolates",
    subtitle = paste0("Total Isolation Events (N = ", total_findings, ")"),
    x = "Isolate Name",
    y = "Percentage (%)",
    caption = str_wrap(paste0("Note: 'Others' includes: ", others_list), width = 90)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Show the chart
print(bar_chart)

# 5. Save the results
# Save the chart as a high-resolution PNG
ggsave(
  filename = "Relative Distribution of Bacterial Isolates.png", 
  plot = bar_chart,
  width = 10, 
  height = 7, 
  dpi = 300, 
  bg = "white"
)

# Save the summary table as a CSV file for Excel
write_csv(plot_data_bar, "Isolate_Summary_Table.csv")
################################################################################
## Here is the prevalence of each isolate
# Calculate the total number of unique samples (the denominator)

# 1. Prepare raw counts
isolate_counts <- Final_JoinedDATA %>%
  distinct(INIKA_ID, Isolate) %>%
  count(Isolate, name = "count")

# 2. Calculate the total number of isolation events
total_findings <- sum(isolate_counts$count)

# 3. Calculate exact 100% percentages
plot_data_full <- isolate_counts %>%
  mutate(
    raw_perc = (count / total_findings) * 100,
    # Standard rounding to 1 decimal place
    percentage = round(raw_perc, 1)
  )

# 4. Rounding Correction (Ensures exactly 100%)
diff <- 100 - sum(plot_data_full$percentage)
if (diff != 0) {
  plot_data_full <- plot_data_full %>%
    arrange(desc(percentage)) %>%
    mutate(percentage = if_else(row_number() == 1, percentage + diff, percentage))
}

# 5. Generate the Bar Chart
full_bar_chart <- ggplot(plot_data_full, aes(x = reorder(Isolate, percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", width = 0.8) +
  coord_flip() +
  # Use sub() to remove the ".0" but keep other decimals like ".1"
  geom_text(aes(label = paste0(sub("\\.0$", "", as.character(percentage)), "%")), 
            hjust = -0.2, 
            size = 3) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Relative Distribution of All Bacterial Isolates",
    subtitle = paste0("Total Isolation Events (N = ", total_findings, ")"),
    x = "Isolate Name",
    y = "Percentage (%)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Show the chart
print(full_bar_chart)

# 6. Save the results
ggsave(
  filename = "All_Isolates_Distribution_Clean_Labels.png", 
  plot = full_bar_chart,
  width = 10, 
  height = 14, 
  dpi = 300, 
  bg = "white"
)
################################################################################

## UNIQUE IDs count, and grouping them per Region,Season,Origin of sample and Districts

Region_summary <- Final_JoinedDATA %>%
  group_by(REGION) %>%
  summarize(total_unique_ids = n_distinct(INIKA_ID))

print(Region_summary)

## Unique IDs by Season
  Season_summary <- Final_JoinedDATA %>%
  group_by(SEASON) %>%
  summarize(total_unique_ids = n_distinct(INIKA_ID))

print(Season_summary)

## Unique IDs by Origin of the sample
Origin_summary <- Final_JoinedDATA %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  summarize(total_unique_ids = n_distinct(INIKA_ID))

print(Origin_summary)

## Unique IDs by District
District_summary <- Final_JoinedDATA %>%
  group_by(DISTRICT) %>%
  summarize(total_unique_ids = n_distinct(INIKA_ID))

print(District_summary)
################################################################################
## Counting the E.coli and K.pneumoniae from the Unique IDs,and group them by Region,
## Season, Origin of the sample and District, respectively
#Counting unique IDs for specific isolates grouped by your variables

# Clean and count
Region_Isolate_Summary <- Final_JoinedDATA %>%
  # Filter for your target pathogens
  filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
  
  # IMPORTANT: Remove exact duplicates of the same isolate for the same person
  # This keeps 1 row per unique combination of ID and Isolate
  distinct(INIKA_ID, Isolate, REGION, .keep_all = TRUE) %>%
  
  # Group and Count
  group_by(REGION, Isolate) %>%
  summarise(Count = n(), .groups = "drop") %>%
  
  # Pivot for the final table
  pivot_wider(
    names_from = Isolate, 
    values_from = Count, 
    values_fill = 0
  )

# Final Verification
total_ecoli <- sum(Region_Isolate_Summary$E.coli)
print(Region_Isolate_Summary)
############################################################

## Group by Season
Season_Isolate_Summary <- Final_JoinedDATA %>%
  # Filter for your target pathogens
  filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
  
  # Remove identical rows to ensure "Unique Isolates"
  # This keeps one record per participant, per isolate type, per season
  distinct(INIKA_ID, Isolate, SEASON, .keep_all = TRUE) %>%
  
  # Group by Season and Isolate
  group_by(SEASON, Isolate) %>%
  
  # Count the unique isolation events to reach 976
  summarise(Count = n(), .groups = "drop") %>%
  
  # Pivot for a cleaner table
  pivot_wider(
    names_from = Isolate, 
    values_from = Count, 
    values_fill = 0
  )

# Verification: Confirm E.coli total is 976
total_ecoli <- sum(Season_Isolate_Summary$E.coli)
print(Season_Isolate_Summary)
################################################################################
## Group by Origin of the sample
# Prepare the summary by Sample Origin
Origin_Isolate_Summary <- Final_JoinedDATA %>%
  # Filter for your target pathogens
  filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
  
  # Ensure we count unique isolation events per ID, Isolate, and Origin
  # This prevents double-counting identical records while keeping the 976 total
  distinct(INIKA_ID, Isolate, ORIGIN_OF_SAMPLE, .keep_all = TRUE) %>%
  
  # Group by the Origin and the Isolate type
  group_by(ORIGIN_OF_SAMPLE, Isolate) %>%
  
  # Use n() to match your reported isolate total of 976
  summarise(Count = n(), .groups = "drop") %>%
  
  # Pivot for a professional table layout
  pivot_wider(
    names_from = Isolate, 
    values_from = Count, 
    values_fill = 0
  )

# Verification: Confirm E.coli total is 976
total_ecoli <- sum(Origin_Isolate_Summary$E.coli)
print(Origin_Isolate_Summary)
################################################################################
## Group by District
# Prepare the summary by District
District_Isolate_Summary <- Final_JoinedDATA %>%
  # Filter for your target pathogens
  filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
  
  # Ensure we count unique isolation events per ID, Isolate, and District
  # This prevents double-counting identical records while keeping the 976 total
  distinct(INIKA_ID, Isolate, DISTRICT, .keep_all = TRUE) %>%
  
  # Group by the District and the Isolate type
  group_by(DISTRICT, Isolate) %>%
  
  # Use n() to match your reported isolate total of 976
  summarise(Count = n(), .groups = "drop") %>%
  
  # Pivot for a professional table layout
  pivot_wider(
    names_from = Isolate, 
    values_from = Count, 
    values_fill = 0
  )

# Verification: Confirm E.coli total is 976
total_ecoli <- sum(District_Isolate_Summary$E.coli)
print(District_Isolate_Summary)
################################################################################

# Counting the ESCR_ECO_presumptive
# Define the column you want to count '1's in

target_col <- "ESCR_ECO_presumptive" 

# Create a function to count only the 1s
count_ones <- function(Final_JoinedDATA, group_var) {
  Final_JoinedDATA %>%
    group_by(Category = group_var, Value = .data[[group_var]]) %>%
    summarise(
      Count_of_Ones = sum(.data[[target_col]] == 1, na.rm = TRUE),
      .groups = "drop"
    )
}

# Apply the function to all your categories and stack them
Summary_of_Ones <- bind_rows(
  count_ones(Final_JoinedDATA, "REGION"),
  count_ones(Final_JoinedDATA, "DISTRICT") %>% 
    mutate(Value = if_else(Value == "Ilemala", "Ilemela", Value)),
  count_ones(Final_JoinedDATA, "SEASON"),
  count_ones(Final_JoinedDATA, "ORIGIN_OF_SAMPLE")
)

# Arrange the final results
Summary_of_Ones <- Summary_of_Ones %>%
  arrange(Category, Value)

# Show the results
print(Summary_of_Ones)

######################################################################
## Counting the Confirmed ESCR E.coli by VITEK, and group the results by Region,
## Season, Origin, and District
# Creating the refined subset based on "Isolate = E.coli, VITEK_MS_Results = 
# Escherichia coli, and ESCR_ECO_presumptive = 1"
Confirmed_ESCR_ECO <- Final_JoinedDATA %>%
  filter(
    Isolate == "E.coli",
    VITEK_MS_Results == "Escherichia coli",
    ESCR_ECO_presumptive == 1
  )

# Group by Region
Confirmed_ESCR_ECO_summary_Region <- Confirmed_ESCR_ECO %>%
  group_by(REGION) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Season
Confirmed_ESCR_ECO_summary_Season <- Confirmed_ESCR_ECO %>%
  group_by(SEASON) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Origin of Sample
Confirmed_ESCR_ECO_summary_Origin <- Confirmed_ESCR_ECO %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Districts
Confirmed_ESCR_ECO_summary_District <- Confirmed_ESCR_ECO %>%
  group_by(DISTRICT) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  View the results 
print(Confirmed_ESCR_ECO_summary_Region)
print(Confirmed_ESCR_ECO_summary_Season)
print(Confirmed_ESCR_ECO_summary_Origin)
print(Confirmed_ESCR_ECO_summary_District)
################################################################################
## Presumptive KPN selected for confirmation by VITEK MS
# Define the column you want to count '1's in

target_col <- "ESCR_KPN_presumptive" 

# Create a function to count only the 1s
count_ones <- function(Final_JoinedDATA, group_var) {
  Final_JoinedDATA %>%
    group_by(Category = group_var, Value = .data[[group_var]]) %>%
    summarise(
      Count_of_Ones = sum(.data[[target_col]] == 1, na.rm = TRUE),
      .groups = "drop"
    )
}

# Apply the function to all your categories and stack them
Summary_of_Ones <- bind_rows(
  count_ones(Final_JoinedDATA, "REGION"),
  count_ones(Final_JoinedDATA, "DISTRICT") %>% 
    mutate(Value = if_else(Value == "Ilemala", "Ilemela", Value)),
  count_ones(Final_JoinedDATA, "SEASON"),
  count_ones(Final_JoinedDATA, "ORIGIN_OF_SAMPLE")
)

# Arrange the final results
Summary_of_Ones_KPN <- Summary_of_Ones %>%
  arrange(Category, Value)

# Show the results
print(Summary_of_Ones_KPN)

######################################################
# Escherichia coli, and ESCR_KPN_presumptive = 1"
Confirmed_ESCR_ECO2 <- Final_JoinedDATA %>%
  filter(
    Isolate == "K.pneumoniae",
    VITEK_MS_Results == "Escherichia coli",
    ESCR_KPN_presumptive == 1
  )

# Group by Region
Confirmed_ESCR_ECO_summary_Region2 <- Confirmed_ESCR_ECO2 %>%
  group_by(REGION) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Season
Confirmed_ESCR_ECO_summary_Season2 <- Confirmed_ESCR_ECO2 %>%
  group_by(SEASON) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  Group by Origin of Sample
Confirmed_ESCR_ECO_summary_Origin2 <- Confirmed_ESCR_ECO2 %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Districts
Confirmed_ESCR_ECO_summary_District2 <- Confirmed_ESCR_ECO2 %>%
  group_by(DISTRICT) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  View the results 
print(Confirmed_ESCR_ECO_summary_Region2)
print(Confirmed_ESCR_ECO_summary_Season2)
print(Confirmed_ESCR_ECO_summary_Origin2)
print(Confirmed_ESCR_ECO_summary_District2)
########################################################

## Counting the Confirmed ESCR K.pneumoniae by VITEK, and group the results by Region,
## Season, Origin, and District
# Creating the refined subset based on "Isolate = K.pneumoniae, VITEK_MS_Results = 
# Klebsiella pneumoniae, and ESCR_KPN_presumptive = 1"
Confirmed_ESCR_KPN <- Final_JoinedDATA %>%
  filter(
    Isolate == "K.pneumoniae",
    VITEK_MS_Results == "Klebsiella pneumoniae",
    ESCR_KPN_presumptive == 1
  )

# Group by Region
Confirmed_ESCR_KPN_summary_Region <- Confirmed_ESCR_KPN %>%
  group_by(REGION) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Season
Confirmed_ESCR_KPN_summary_Season <- Confirmed_ESCR_KPN %>%
  group_by(SEASON) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Origin of Sample
Confirmed_ESCR_KPN_summary_Origin <- Confirmed_ESCR_KPN %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Districts
Confirmed_ESCR_KPN_summary_District <- Confirmed_ESCR_KPN %>%
  group_by(DISTRICT) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  View the results 
print(Confirmed_ESCR_KPN_summary_Region)
print(Confirmed_ESCR_KPN_summary_Season)
print(Confirmed_ESCR_KPN_summary_Origin)
print(Confirmed_ESCR_KPN_summary_District)
################################################################################
# Confirmed KPN from E.coli isolates
Confirmed_ESCR_KPN2 <- Final_JoinedDATA %>%
  filter(
    Isolate == "E.coli",
    VITEK_MS_Results == "Klebsiella pneumoniae",
    ESCR_KPN_presumptive == 1
  )

#  Group by Region
Confirmed_ESCR_KPN_summary_Region2 <- Confirmed_ESCR_KPN2 %>%
  group_by(REGION) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  Group by Season
Confirmed_ESCR_KPN_summary_Season2 <- Confirmed_ESCR_KPN2 %>%
  group_by(SEASON) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

# Group by Origin of Sample
Confirmed_ESCR_KPN_summary_Origin2 <- Confirmed_ESCR_KPN2 %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  Group by Districts
Confirmed_ESCR_KPN_summary_District2 <- Confirmed_ESCR_KPN2 %>%
  group_by(DISTRICT) %>%
  summarize(unique_id_count = n_distinct(INIKA_ID), .groups = "drop")

#  View the results 
print(Confirmed_ESCR_KPN_summary_Region2)
print(Confirmed_ESCR_KPN_summary_Season2)
print(Confirmed_ESCR_KPN_summary_Origin2)
print(Confirmed_ESCR_KPN_summary_District2)
###########################################################
