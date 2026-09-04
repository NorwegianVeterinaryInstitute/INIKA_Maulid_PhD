
# Note Madelaine found that something has gone wrong by the joining here! We need to go trough each step carefully to see what is happening!
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)


# Importing the Demographic_cleaned_data


UniqueData <- readRDS("data/CLEANED_DATA/Allunique.rds")
names(UniqueData)

# Selecting relevant columns
UniqueData <-UniqueData %>%
  select(-"SEASON.y", - "REGION.y", -"DISTRICT.y", -"INIKA_ID.y", 
         -"SEASON.x.x", -"SAMPLE FROM.y", -"NAME OF SCHOOL/HF.y", 
         -"REGION.x.x", -"DISTRICT.x.x", -"AR_RESIDUES_SAMPLE.y", 
         -"COLONY MORPHOLOGY ON C3GR.y", -"COLONY MORPHOLOGY ON CARBA.y", 
         -"COLONY MORPHOLOGY ON XLD.y", -"COLONY MORPHOLOGY ON BGA.y", 
         -"TSI Media_Slope.y", -"TSI Media_Butt.y", -"TSI Media_Gas.y", 
         -"TSIMedia_H2S.y", -"SIM Media_H2S.y", -"SIM Media Indole.y", 
         -"SIM Media Motility.y", -"CITRATE.y", -"UREASE.y", 
         -"Isolate.y", -"Freezing_ID_No.y", -"INIKA_prefix.y", 
         -"Isolate_suffix.y", -"Isolate_ID_repeat_count.y", 
         -"Isolate_ID_row_count.y", -"INIKA_ID.x.x", -"REGION.y.y", 
         -"DISTRICT.y.y", -"SEASON.y.y", -"ORIGIN_OF_SAMPLE.y", 
         -"INIKA_ID.y.y", -"SAMPLE FROM.x")
# Rename the columns
UniqueData <-UniqueData %>%
  rename(INIKA_ID = `INIKA_ID.x`, SEASON = `SEASON.x`,
         REGION = `REGION.x`, DISTRICT = `DISTRICT.x`,
         ORIGIN_OF_SAMPLE = `ORIGIN_OF_SAMPLE.x`,
         `COLONY MORPHOLOGY ON C3GR` = `COLONY MORPHOLOGY ON C3GR.x`,
         `COLONY MORPHOLOGY ON CARBA` = `COLONY MORPHOLOGY ON CARBA.x`,
         `COLONY MORPHOLOGY ON XLD` = `COLONY MORPHOLOGY ON XLD.x`,
         `COLONY MORPHOLOGY ON BGA` = `COLONY MORPHOLOGY ON BGA.x`,
         `TSI Media_Slope` = `TSI Media_Slope.x`,
         `TSI Media_Butt` = `TSI Media_Butt.x`,
         `TSI Media_Gas` = `TSI Media_Gas.x`,
         TSIMedia_H2S = `TSIMedia_H2S.x`, `SIM Media_H2S` = `SIM Media_H2S.x`,
         `SIM Media Indole` = `SIM Media Indole.x`,
         `SIM Media Motility` = `SIM Media Motility.x`,
         CITRATE = `CITRATE.x`,UREASE = `UREASE.x`,
         Isolate = `Isolate.x`,Freezing_ID_No = `Freezing_ID_No.x`,
         INIKA_prefix = `INIKA_prefix.x`,
         Isolate_suffix = `Isolate_suffix.x`,
         Isolate_ID_repeat_count = `Isolate_ID_repeat_count.x`,
         Isolate_ID_row_count = `Isolate_ID_row_count.x`,
         AR_RESIDUES_SAMPLE = `AR_RESIDUES_SAMPLE.x`,
         `NAME OF SCHOOL/HF` = `NAME OF SCHOOL/HF.x`)


## 23.12.25 Madelaine Notes that we have changed the criteria for the ESBL
##as this needs to be according to both the correct protocol as well as the criteria below - see script 4 lines 88 to 106!
## Select the true ESBL  
UniqueData <- UniqueData %>%
  mutate(
    ESBL = if_else(ESBL_Selection >= 5, 1, 0, missing = 0)
  ) 

# Save the file
write.csv(UniqueData, "data/CLEANED_DATA/UniqueData.csv")

long_data <- UniqueData %>%
  pivot_longer(
    cols = c(
      `COLONY MORPHOLOGY ON C3GR`,
      `COLONY MORPHOLOGY ON CARBA`, 
      `COLONY MORPHOLOGY ON XLD`, 
      `COLONY MORPHOLOGY ON BGA`, 
      CITRATE, `TSI Media_Slope`, 
      `TSI Media_Butt`, 
      `TSI Media_Gas`, TSIMedia_H2S, 
      TSIMedia_H2S, `SIM Media Indole`,
      `SIM Media Motility`,
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

## Save the file_ Maulid needs to complete. 
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
##MMJ 28/01/2026
# Create the summary table
tvla_verification_table <- UniqueData %>%
  # Filter for only the three isolates of interest
  filter(Isolate %in% c("E.coli", "K.pneumoniae", "S.typhimurium")) %>%
  group_by(Isolate) %>%
  summarise(
    # Count how many have a TVLA_ID
    `Tested at TVLA` = sum(!is.na(TVLA_ID), na.rm = TRUE),
    
    # Apply specific confirmation criteria for each isolate
    `Confirmed as per Criteria` = sum(
      !is.na(TVLA_ID) & (
        (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
          (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae") |
          (Isolate == "S.typhimurium" & VITEK_MS_Results != "S.typhimurium")
      ), 
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# View the result
print(tvla_verification_table)
################################################################################
# 1. Define the criteria function to keep the code clean
calculate_counts <- function(UniqueData) {
  UniqueData %>%
    summarise(
      `Tested at TVLA` = sum(!is.na(TVLA_ID), na.rm = TRUE),
      `Confirmed as per Criteria` = sum(
        !is.na(TVLA_ID) & (
          (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
            (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae") |
            (Isolate == "S.typhimurium" & VITEK_MS_Results != "S.typhimurium")
        ), 
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

# 2. Filter base data for target isolates
filtered_data <- UniqueData %>%
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
## MMJ from HERE down
## Calculating the Isolation rate (some of the isolates being grouped in "others")
# 1. Calculate the total number of unique samples (the denominator)
total_unique_samples <- UniqueData %>% 
  summarise(n = n_distinct(INIKA_ID)) %>% 
  pull(n)

# 2. Extract the "Others" list for the key
# We identify isolates that appear in < 2% of the total unique samples
others_list <- UniqueData %>%
  filter(Isolate != "No growth") %>%
  distinct(INIKA_ID, Isolate) %>%
  count(Isolate) %>%
  mutate(rate = (n / total_unique_samples) * 100) %>%
  filter(rate < 2) %>%
  pull(Isolate) %>%
  paste(collapse = ", ")

# 3. Prepare data for the Bar Chart
plot_data_bar <- UniqueData %>%
  filter(Isolate != "No growth") %>% 
  distinct(INIKA_ID, Isolate) %>% 
  group_by(Isolate) %>%
  summarise(unique_id_count = n(), .groups = "drop") %>%
  mutate(
    isolation_rate = (unique_id_count / total_unique_samples) * 100,
    Isolate_Grouped = if_else(isolation_rate < 2, "Others", Isolate)
  ) %>%
  group_by(Isolate_Grouped) %>%
  summarise(isolation_rate = sum(isolation_rate), .groups = "drop") %>%
  arrange(desc(isolation_rate))

# 4. Generate the Bar Chart
bar_chart <- ggplot(plot_data_bar, aes(x = reorder(Isolate_Grouped, isolation_rate), y = isolation_rate)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = paste0(round(isolation_rate, 1), "%")), 
            hjust = -0.2, 
            size = 3.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "plain", hjust = 0.5),
    # Formatting the caption (The "Others" Key) to wrap text
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0, margin = margin(t = 15))
  ) +
  labs(
    title = "Prevalence Rate of Pathogen Isolation",
    subtitle = paste0("Total Unique Samples (N = ", total_unique_samples, ")"),
    x = "Isolate Name",
    y = "Prevalence Rate (%)",
    # We use str_wrap to ensure the list of isolates doesn't run off the page
    caption = str_wrap(paste0("Note: 'Others' includes: ", others_list), width = 90)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Show the chart
print(bar_chart)

# 5. Save the results
# Save the chart as a high-resolution PNG
ggsave(
  filename = "Isolate_Prevalence_BarChart.png", 
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
# 1. Calculate the total number of unique samples (the denominator)
total_unique_samples <- UniqueData %>% 
  summarise(n = n_distinct(INIKA_ID)) %>% 
  pull(n)

# 2. Prepare the full data (No "Others" grouping)
plot_data_full <- UniqueData %>%
  filter(Isolate != "No growth") %>% 
  # Ensure an isolate is only counted once per ID
  distinct(INIKA_ID, Isolate) %>% 
  group_by(Isolate) %>%
  summarise(unique_id_count = n(), .groups = "drop") %>%
  mutate(
    isolation_rate = (unique_id_count / total_unique_samples) * 100
  ) %>%
  # Sort from highest prevalence to lowest
  arrange(desc(isolation_rate))

# 3. Generate the Bar Chart
full_bar_chart <- ggplot(plot_data_full, aes(x = reorder(Isolate, isolation_rate), y = isolation_rate)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", width = 0.8) +
  # Flip coordinates to handle long lists of names easily
  coord_flip() +
  # Add percentage labels at the end of each bar
  geom_text(aes(label = paste0(round(isolation_rate, 1), "%")), 
            hjust = -0.2, 
            size = 3, 
            fontface = "plain") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    # Increase plot height dynamically if needed when viewing
    axis.text.y = element_text(size = 9), 
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Prevalence Rate of All Identified Pathogens",
    subtitle = paste0("Denominator: Total Unique Samples (N = ", total_unique_samples, ")"),
    x = "Isolate Name",
    y = "Prevalence Rate (%)"
  ) +
  # Ensure the x-axis (now top axis due to flip) has room for labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Show the chart
print(full_bar_chart)

# 4. Save the results
# If you have many isolates, increase the 'height' parameter to prevent squishing
ggsave(
  filename = "All_Isolates_Prevalence.png", 
  plot = full_bar_chart,
  width = 10, 
  height = 12,  # Increased height to accommodate every isolate name
  dpi = 300, 
  bg = "white"
)
################################################################################
#Calculating the proportions distribution of isolates

# 1. Identify "Others" and Calculate Total Isolates 
# We create a base summary first to avoid repeating code
library(ggrepel)
base_summary <- UniqueData %>%
  filter(Isolate != "No growth") %>% 
  distinct(INIKA_ID, Isolate) %>% 
  count(Isolate, name = "count")

# Calculate Total Isolates (The sum of all unique occurrences)
total_isolates_val <- sum(base_summary$count)

# Identify which isolates fall into the "Others" category
others_list <- base_summary %>%
  mutate(temp_pct = (count / total_isolates_val) * 100) %>%
  filter(temp_pct < 2) %>%
  pull(Isolate) %>%
  paste(collapse = ", ")

# 2. Prepare the plot data 
plot_data <- base_summary %>%
  mutate(temp_pct = (count / total_isolates_val) * 100) %>%
  mutate(Isolate_Grouped = if_else(temp_pct < 2, "Others", Isolate)) %>%
  group_by(Isolate_Grouped) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    percentage = round(count / sum(count) * 100, 1),
    label = paste0(percentage, "%")
  ) %>%
  arrange(desc(count))

#  3. Generate the Pie Chart 
ggplot(plot_data, aes(x = 1, y = count, fill = reorder(Isolate_Grouped, -count))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text_repel(
    aes(label = label),
    position = position_stack(vjust = 0.5), 
    size = 3.5,
    color = "black",
    fontface = "plain",
    segment.size = 0.3,
    min.segment.length = 0, 
    box.padding = 0.5,       
    point.padding = 0.2,     
    show.legend = FALSE
  ) +
  theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5), # Centered subtitle
    plot.caption = element_text(size = 8, color = "grey30", hjust = 0, margin = margin(t = 15))
  ) +
  labs(
    title = "Proportional Distribution of Isolates",
    # Added Total Isolates to the subtitle
    subtitle = paste0("Total Number of Isolates identified: ", total_isolates_val),
    fill = "Isolate Name",
    caption = str_wrap(paste0("Note: 'Others' includes: ", others_list), width = 90)
  ) +
  scale_fill_brewer(palette = "Pastel1")

# 4. Save the Chart 
ggsave(
  filename = "Isolate_Distribution_Chart.png", 
  width = 10, 
  height = 8, 
  dpi = 300, 
  bg = "white"
)
################################################################################
# 1. Generate the Summary Table
summary_table <- UniqueData %>%
  filter(Isolate != "No growth") %>% 
  distinct(REGION, INIKA_ID, Isolate) %>% 
  count(REGION, Isolate, name = "count") %>%
  # Identify "Others" based on the 2% global threshold
  group_by(Isolate) %>%
  mutate(global_total = sum(count)) %>%
  ungroup() %>%
  mutate(Isolate_Grouped = if_else((global_total / sum(count)) * 100 < 2, "Others", Isolate)) %>%
  # Group by Region and the new Grouped Isolate
  group_by(REGION, Isolate_Grouped) %>%
  summarise(Total_Count = sum(count), .groups = "drop_last") %>%
  # Calculate Percentage within each Region
  mutate(
    Percentage = round((Total_Count / sum(Total_Count)) * 100, 1)
  ) %>%
  arrange(REGION, desc(Total_Count))

# View the table
print(summary_table)

# Optional: Export to CSV
# write.csv(summary_table, "Isolate_Summary_By_Region.csv", row.names = FALSE)

################################################################################
# 1. Create the Detailed Regional Table
regional_summary <- UniqueData %>%
  filter(Isolate != "No growth") %>%
  # Ensure unique patient-isolate combinations per region
  distinct(REGION, INIKA_ID, Isolate) %>% 
  # Count occurrences of each Isolate within each Region
  count(REGION, Isolate, name = "Count") %>%
  group_by(REGION) %>%
  # Calculate percentage based on the total isolates in THAT region
  mutate(
    Percentage = round((Count / sum(Count)) * 100, 1)
  ) %>%
  ungroup() %>%
  # Sort by Region and then the most frequent isolates
  arrange(REGION, desc(Count))

# 2. Display the result
print(regional_summary)
write_tsv(regional_summary,"Results/regional_summary.tsv")
################################################################################
# 1. Create the Detailed Seasonal Table
# Replace "Season" with your actual column name if it differs (e.g., "Month" or "Quarter")
seasonal_summary <- UniqueData %>%
  filter(Isolate != "No growth") %>%
  # Ensure unique patient-isolate combinations per Season
  distinct(SEASON, INIKA_ID, Isolate) %>% 
  # Count occurrences of each Isolate within each Season
  count(SEASON, Isolate, name = "Count") %>%
  group_by(SEASON) %>%
  # Calculate percentage based on the total isolates in THAT season
  mutate(
    Percentage = round((Count / sum(Count)) * 100, 1)
  ) %>%
  ungroup() %>%
  # Sort by Season and then by the most frequent isolates
  arrange(SEASON, desc(Count))

# 2. Display the result
print(seasonal_summary)
write_tsv(seasonal_summary,"Results/seasonal_summary.tsv")
################################################################################
# 1. Create the Detailed Summary by Sample Origin
Origin_of_sample_summary <- UniqueData %>%
  filter(Isolate != "No growth") %>%
  # Ensure unique patient-isolate combinations per Origin
  # Replace 'Origin' with your actual column name
  distinct(ORIGIN_OF_SAMPLE, INIKA_ID, Isolate) %>% 
  # Count occurrences of each Isolate within each Origin type
  count(ORIGIN_OF_SAMPLE, Isolate, name = "Count") %>%
  group_by(ORIGIN_OF_SAMPLE) %>%
  # Calculate percentage based on the total isolates in THAT specific origin
  mutate(
    Percentage = round((Count / sum(Count)) * 100, 1)
  ) %>%
  ungroup() %>%
  # Sort by Origin and then the highest frequency
  arrange(ORIGIN_OF_SAMPLE, desc(Count))

# 2. Display the result
print(Origin_of_sample_summary)
write_tsv(Origin_of_sample_summary,"Results/Origin_of_sample_summary.tsv")
################################################################################


## Confirmed Isolates as Presumptive ESCR
# Total unique IDs that have any result in the VITEK column
total_tested_vitek <- UniqueData %>%
  filter(!is.na(VITEK_MS_Results) & VITEK_MS_Results != "") %>%
  summarise(n = n_distinct(INIKA_ID)) %>%
  pull(n)

# Count of confirmed E. coli
e_coli_count <- UniqueData %>%
  filter(VITEK_MS_Results == "Escherichia coli") %>%
  nrow()

# Count of confirmed K. pneumoniae
k_pneumo_count <- UniqueData %>%
  filter(VITEK_MS_Results == "Klebsiella pneumoniae") %>%
  nrow()

# 2. Construct the summary table
vitek_analysis_table <- tribble(
  ~Metric,                                     ~Count,
  "Total unique IDs tested by VITEK MS",       total_tested_vitek,
  "Isolates confirmed as Escherichia coli",    e_coli_count,
  "Isolates confirmed as Klebsiella pneumoniae", k_pneumo_count,
  "Total Confirmed (E. coli + K. pneumoniae)", e_coli_count + k_pneumo_count
)

# 3. Add a percentage column relative to the total tested
vitek_analysis_table <- vitek_analysis_table %>%
  mutate(
    Percentage = round((Count / total_tested_vitek) * 100, 1),
    # Format for presentation
    Percentage = if_else(Metric == "Total unique IDs tested by VITEK MS", "100%", paste0(Percentage, "%"))
  )

# View the table
print(vitek_analysis_table)

# 4. Save the table to a CSV file
write_csv(vitek_analysis_table, "VITEK_Summary_Report.csv")
################################################################################
## Here we count the confirmed one which are true ESCR

# 1. Perform the calculations with the triple-condition (Isolate + VITEK + CTX)
# Total unique IDs that have any result in the VITEK column (Denominator)
total_tested_vitek <- UniqueData %>%
  filter(!is.na(VITEK_MS_Results) & VITEK_MS_Results != "") %>%
  summarise(n = n_distinct(INIKA_ID)) %>%
  pull(n)

# Count E. coli matching VITEK result AND CTX_ED5 < 21
e_coli_ctx_count <- UniqueData %>%
  filter(Isolate == "E.coli" & 
           VITEK_MS_Results == "Escherichia coli" & 
           as.numeric(as.character(CTX_ED5)) < 21) %>%
  nrow()

# Count K. pneumoniae matching VITEK result AND CTX_ED5 < 22
k_pneumo_ctx_count <- UniqueData %>%
  filter(Isolate == "K.pneumoniae" & 
           VITEK_MS_Results == "Klebsiella pneumoniae" & 
           as.numeric(as.character(CTX_ED5)) < 21) %>%
  nrow()

# 2. Construct the summary table
ctx_summary_table <- tribble(
  ~Metric,                                                  ~Count,
  "Total unique IDs tested by VITEK MS",                    total_tested_vitek,
  "E. coli (Confirmed & CTX_ED5 < 21)",                        e_coli_ctx_count,
  "K. pneumoniae (Confirmed & CTX_ED5 < 21)",                  k_pneumo_ctx_count,
  "Total Target Pathogens (Confirmed & CTX_ED5 < 21)",         e_coli_ctx_count + k_pneumo_ctx_count
)

# 3. Add a percentage column relative to the total VITEK tested
ctx_summary_table <- ctx_summary_table %>%
  mutate(
    Percentage = round((Count / total_tested_vitek) * 100, 1),
    Percentage = if_else(Metric == "Total unique IDs tested by VITEK MS", "100%", paste0(Percentage, "%"))
  )

# View the results
print(ctx_summary_table)

# 4. Save the results
write_csv(ctx_summary_table, "CTX_Filtered_VITEK_Report.csv")
################################################################################
# 1. Pre-process and filter for VITEK-Confirmed cases only
# This creates the "Denominator" pool based strictly on VITEK confirmation
base_confirmed <- UniqueData %>%
  mutate(
    CTX_num = as.numeric(as.character(CTX_ED5)),
    CRO_num = as.numeric(as.character(CRO_ED30))
  ) %>%
  filter(
    (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
      (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae")
  )

# 2. Define the core calculation function
calc_validation_confirmed <- function(df, grouping_var) {
  df %>%
    group_by(across(all_of(grouping_var)), Isolate) %>%
    summarise(
      # Denominator: Cases already confirmed by VITEK
      `VITEK_MS_Results Confirmed (Denominator)` = n(),
      
      # Numerator: Of those confirmed, how many meet the CTX/CRO criteria
      `Met Antibiotic Criteria` = sum(
        (Isolate == "E.coli" & CTX_num < 22 & CRO_num < 23) |
          (Isolate == "K.pneumoniae" & CTX_num < 21 & CRO_num < 23),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      `Criteria Match Rate (%)` = round((`Met Antibiotic Criteria` / `VITEK_MS_Results Confirmed (Denominator)`) * 100, 1)
    )
}

# 3. Generate the Three Separate Tables
summary_by_region <- calc_validation_confirmed(base_confirmed, "REGION")
summary_by_season <- calc_validation_confirmed(base_confirmed, "SEASON")
summary_by_origin <- calc_validation_confirmed(base_confirmed, "ORIGIN_OF_SAMPLE")

# 4. Print Results
print("--- AGGREGATE BY REGION (VITEK Denominator) ---")
print(summary_by_region)

print("--- AGGREGATE BY SEASON (VITEK Denominator) ---")
print(summary_by_season)

print("--- AGGREGATE BY ORIGIN (VITEK Denominator) ---")
print(summary_by_origin)
################################################################################