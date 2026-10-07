## Final joined data

library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)
library(flextable)
library(officer)
library(survey)

#It seem I have to include something here
Sys.getlocale()
Sys.setlocale("LC_CTYPE", "Norwegian_Norway.utf8")

# Importing the Demographic_cleaned_data
# Import the file using file.path()
Final_JoinedDATA <- read_csv ("data/CLEANED_DATA/joined_data_16.3.26.csv")
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

# We filter out the columns where Biochemical and TVLA results are the same for E.coli
ECO_correct_Identified<-grouped_summary%>%
  filter(Isolate =="E.coli", VITEK_MS_Results =="Escherichia coli") ## 130 observation

## Save the file as tsv and rds. 
write_tsv(ECO_correct_Identified,"data/CLEANED_DATA/ECO_correct_Identified.tsv")

saveRDS(ECO_correct_Identified,"data/CLEANED_DATA/ECO_correct_Identified.rds")

## We filter out the columns where Biochemical and TVLA results are not the same for E.coli
ECO_NOT_correct_Identified_Bio<-grouped_summary%>%
  filter(Isolate =="E.coli", VITEK_MS_Results!="Escherichia coli") 

write_tsv(ECO_NOT_correct_Identified_Bio,"data/CLEANED_DATA/ECO_NOT_correct_Identified_Bio.tsv")

saveRDS(ECO_NOT_correct_Identified_Bio,"data/CLEANED_DATA/ECO_NOT_correct_Identified_Bio.rds")

## We filter out the columns where Isolate is not E.coli, but turned E.coli by VITEK
ECO_NOT_correct_Identified_VITEK<-grouped_summary%>%
  filter(Isolate !="E.coli", VITEK_MS_Results =="Escherichia coli") # 21 

write_tsv(ECO_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/ECO_NOT_correct_Identified_VITEK.tsv")

saveRDS(ECO_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/ECO_NOT_correct_Identified_VITEK.rds")

## We filter out columns where Biochemical and TVLA results are not the same for K.pneumoniae
Kleb_correct_Identified<-grouped_summary%>%
  filter(Isolate =="K.pneumoniae", VITEK_MS_Results =="Klebsiella pneumoniae") # 52 

write_tsv(Kleb_correct_Identified,"data/CLEANED_DATA/Kleb_correct_Identified.tsv")

saveRDS(Kleb_correct_Identified,"data/CLEANED_DATA/Kleb_correct_Identified.rds") 

## We filter out columns where Biochemical revealed Isolate as K.pneumoniae but VITEK did not confirm as K.pneumoniae
Kleb_NOT_correct_Identified_Bio<-grouped_summary%>%
  filter(Isolate =="K.pneumoniae", VITEK_MS_Results!="Klebsiella pneumoniae") # 31

write_tsv(Kleb_NOT_correct_Identified_Bio,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_Bio.tsv")

saveRDS(Kleb_NOT_correct_Identified_Bio,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_Bio.rds")

## We filter out columns where Biochemical did not revealed Isolate as K.pneumoniae but VITEK confirmed as K.pneumoniae
Kleb_NOT_correct_Identified_VITEK<-grouped_summary%>%
  filter(Isolate !="K.pneumoniae", VITEK_MS_Results=="Klebsiella pneumoniae") # 8 

write_tsv(Kleb_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_VITEK.tsv")

saveRDS(Kleb_NOT_correct_Identified_VITEK,"data/CLEANED_DATA/Kleb_NOT_correct_Identified_VITEK.rds")
################################################################################
# 1.10.26 # We create a combined table for publication
# 1. Define function to calculate counts and formatted percentages
calculate_counts <- function(df, group_var) {
  df %>%
    group_by(across(all_of(group_var)), Isolate) %>%
    summarise(
      Tested = sum(!is.na(TVLA_ID), na.rm = TRUE),
      Confirmed_Count = sum(
        !is.na(TVLA_ID) & (
          (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
            (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae") |
            (Isolate == "S.typhimurium" & VITEK_MS_Results == "S.typhimurium")
        ), 
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Confirmation_Pct = ifelse(Tested > 0, (Confirmed_Count / Tested) * 100, 0),
      `Confirmed as per Criteria (%)` = sprintf("%d (%.1f%%)", Confirmed_Count, Confirmation_Pct)
    )
}

# 2. Filter base data for target isolates
filtered_data <- Final_JoinedDATA %>%
  filter(Isolate %in% c("E.coli", "K.pneumoniae", "S.typhimurium"))

# 3. Calculate metrics across each category
origin_counts <- calculate_counts(filtered_data, "ORIGIN_OF_SAMPLE") %>%
  rename(Subgroup = ORIGIN_OF_SAMPLE) %>%
  mutate(Category = "Origin of Sample")

region_counts <- calculate_counts(filtered_data, "REGION") %>%
  rename(Subgroup = REGION) %>%
  mutate(Category = "Region")

season_counts <- calculate_counts(filtered_data, "SEASON") %>%
  rename(Subgroup = SEASON) %>%
  mutate(Category = "Season")

# 4. Combine data and set standardized species names
combined_table_data <- bind_rows(origin_counts, region_counts, season_counts) %>%
  mutate(
    Isolate = case_when(
      Isolate == "E.coli" ~ "Escherichia coli",
      Isolate == "K.pneumoniae" ~ "Klebsiella pneumoniae",
      Isolate == "S.typhimurium" ~ "Salmonella Typhimurium",
      TRUE ~ Isolate
    )
  ) %>%
  select(
    Isolate,
    Category,
    `Origin / Region / Season` = Subgroup,
    `Tested at TVLA` = Tested,
    `Confirmed as per Criteria (%)`
  ) %>%
  arrange(Isolate, Category, `Origin / Region / Season`)

# 5. Convert to grouped data frame so 'Isolate' acts as a single header banner
grouped_df <- as_grouped_data(x = combined_table_data, groups = c("Isolate"))

# 6. Build the publication-ready flextable with merged Category cells
formatted_flextable <- grouped_df %>%
  flextable() %>%
  theme_booktabs() %>%
  # Merge vertically so Category name (e.g. "Origin of Sample") appears only once
  merge_v(j = "Category") %>%
  valign(j = "Category", val = "top") %>%
  # Style the Isolate section header banner rows
  italic(i = ~ !is.na(Isolate), j = 1, part = "body") %>%
  bold(i = ~ !is.na(Isolate), part = "body") %>%
  bg(i = ~ !is.na(Isolate), bg = "#EAEAEA", part = "body") %>%
  # Header & Alignment formatting
  bold(part = "header") %>%
  bold(j = "Category", part = "body") %>%
  align(j = 1:2, align = "left", part = "all") %>%
  align(j = 3:4, align = "center", part = "all") %>%
  autofit() %>%
  set_caption(caption = "Table 1: TVLA sample confirmation counts and rates (%) categorized by origin, region, and season across target isolates.")

# 7. Save Word Document in Results directory
if (!dir.exists("Results")) {
  dir.create("Results")
}

doc <- read_docx() %>%
  body_add_flextable(formatted_flextable)

print(doc, target = "Results/Publication_Table_Merged_Categories.docx")

cat("Table successfully saved to 'Results/Publication_Table_Merged_Categories.docx'\n")
####################################################
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
# 17/9/2026
#Computing the adjusted prevalence for E.coli

# ----------------------------------------------------------------------
# Directory & Data Preprocessing
# ----------------------------------------------------------------------
if (!dir.exists("Results")) {
  dir.create("Results")
}

###################################################################

# Helper function to remove .0 decimals (shows integer if decimal is 0, else 1 decimal place)
fmt_num <- function(x) {
  sub("\\.0$", "", sprintf("%.1f", x))
}

# Construct initial dataset with raw count data
df_raw <- tibble(
  Variables = c(
    "Overall",
    "Origin of Sample", "Origin of Sample",
    "Region", "Region",
    "Season", "Season",
    "District", "District", "District", "District", "District", "District"
  ),
  Values = c(
    "Total Population",
    "Adult outpatient", "Schoolchildren",
    "Kilimanjaro", "Mwanza",
    "Dry", "Wet/Rainy",
    "Hai", "Moshi Rural", "Moshi Urban", "Ilemela", "Magu", "Nyamagana"
  ),
  N  = c(2204, 840, 1364, 1098, 1106, 1108, 1096, 366, 365, 367, 365, 374, 367),
  A  = c(972,  430,  542,  451,  521,  444,  528, 148, 140, 163, 193, 195, 133),
  B  = c(163,   66,   97,   45,  118,   82,   81,  16,  15,  14,  28,  49,  41),
  MA = c(193,   85,  108,   89,  104,   85,  108,  29,  28,  32,  38,  39,  27),
  CA = c(129,   63,   66,   64,   65,   59,   70,  24,  15,  25,  24,  24,  24),
  MB = c(153,   61,   92,   43,  110,   78,   75,  16,  14,  13,  26,  44,  40),
  CB = c(23,     9,   14,    4,   19,   11,   12,   1,   2,   1,   4,  11,   4)
)

# Chi-square and p-values by Variable category
chi_tbl <- tibble(
  Variables = c("Overall", "Origin of Sample", "Region", "Season", "District"),
  Chi2      = c("-", "45.537", "0.422", "3.101", "25.045"),
  p_value   = c("-", "<0.001", "0.516", "0.078", "<0.001")
)

# Compute whole integer adjusted n, prevalence %, and 95% CI (formatted without trailing .0)
Final_JoinedDATA <- df_raw %>%
  mutate(
    `ĥ (n)`  = round((CA / MA) * A + (CB / MB) * B),
    p_hat     = `ĥ (n)` / N,
    p_pct_num = p_hat * 100,
    `ĥ (%)`  = fmt_num(p_pct_num),
    se        = sqrt(p_hat * (1 - p_hat) / N),
    ci_low    = pmax(0, (p_hat - 1.96 * se) * 100),
    ci_high   = pmin(100, (p_hat + 1.96 * se) * 100),
    `95% CI`  = paste0(fmt_num(ci_low), " - ", fmt_num(ci_high))
  ) %>%
  left_join(chi_tbl, by = "Variables") %>%
  select(Variables, Values, N, A, B, MA, CA, MB, CB, `ĥ (n)`, `ĥ (%)`, `95% CI`, Chi2, p_value)

# Save raw dataset files to Results folder
write_csv(Final_JoinedDATA, file = "Results/adjusted_prevalence_summary.csv")
write_tsv(Final_JoinedDATA, file = "Results/adjusted_prevalence_summary.tsv")
saveRDS(Final_JoinedDATA, file = "Results/adjusted_prevalence_summary.rds")

# Format table for Word document output (Variables, Chi2, and p-value appear ONCE per category)
df_word <- Final_JoinedDATA %>%
  mutate(
    Variables = if_else(duplicated(Variables), "", Variables),
    Chi2      = if_else(Variables == "", "", Chi2),
    p_value   = if_else(Variables == "", "", p_value)
  )

# Build publication-ready flextable
ft <- flextable(df_word) %>%
  set_header_labels(
    Variables = "Variables",
    Values    = "Values",
    N         = "N",
    A         = "A",
    B         = "B",
    MA        = "MA",
    CA        = "CA",
    MB        = "MB",
    CB        = "CB",
    `ĥ (n)`   = "ĥ (n)",
    `ĥ (%)`   = "ĥ (%)",
    `95% CI`  = "95% CI",
    Chi2      = "χ²",
    p_value   = "p-value"
  ) %>%
  add_footer_lines(
    values = c(
      "Abbreviations: N = Total study population analyzed.",
      "3GCR = 3rd Generation Cephalosporin-Resistant.",
      "A = Presumptive 3GCR E. coli isolates; B = Presumptive 3GCR K. pneumoniae isolates.",
      "MA = Selected/tested presumptive 3GCR E. coli isolates on VITEK MS.",
      "CA = Confirmed true E. coli isolates from group MA.",
      "MB = Selected/tested presumptive 3GCR K. pneumoniae isolates on VITEK MS.",
      "CB = Confirmed true E. coli isolates from group MB.",
      "ĥ (n) = Adjusted number of confirmed E. coli (rounded to whole integer).",
      "ĥ (%) = Adjusted prevalence percentage; 95% CI = 95% Confidence Interval for prevalence percentage (expressed as numbers without trailing .0 for whole values).",
      "χ² = Chi-square test statistic; p-value = Statistical significance level."
    )
  ) %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%
  align(j = c("Variables", "Values"), align = "left", part = "all") %>%
  autofit()

# Save Word document (.docx)
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "Results/adjusted_prevalence_summary.docx")

cat("Files generated successfully in the 'Results' folder!\n")
################################################################################
#Computing the adjusted prevalence for K.pneumoniae

# Helper function to strip trailing .0 from rounded single-decimal values
fmt_num <- function(x) {
  sub("\\.0$", "", sprintf("%.1f", x))
}

#  Raw Dataset Construction for K. pneumoniae
df_kp_raw <- tibble(
  Variables = c(
    "Overall",
    "Origin of Sample", "Origin of Sample",
    "Region", "Region",
    "Season", "Season",
    "District", "District", "District", "District", "District", "District"
  ),
  Values = c(
    "Total Population",
    "Adult outpatient", "Schoolchildren",
    "Kilimanjaro", "Mwanza",
    "Dry", "Wet/Rainy",
    "Hai", "Moshi Rural", "Moshi Urban", # Kilimanjaro districts (alphabetical)
    "Ilemela", "Magu", "Nyamagana"      # Mwanza districts (alphabetical)
  ),
  N  = c(2204, 840, 1364, 1098, 1106, 1108, 1096, 366, 365, 367, 365, 374, 367),
  A  = c(163,   66,   97,   45,  118,   82,   81,  16,  15,  14,  28,  49,  41), # Presumptive 3GCR K. pneumoniae
  B  = c(972,  430,  542,  451,  521,  444,  528, 148, 140, 163, 193, 195, 133), # Presumptive 3GCR E. coli
  MA = c(153,   61,   92,   43,  110,   78,   75,  16,  14,  13,  26,  44,  40), # Tested presumptive K. pneumoniae
  CA = c(52,     21,   31,   9,   43,   25,   27,   0,   1,   8,   4,  17,  22), # Confirmed K. pneumoniae from MA
  MB = c(193,   85,  108,   89,  104,   85,  108,  29,  28,  32,  38,  39,  27), # Tested presumptive E. coli
  CB = c(8,      4,    4,    0,    8,    3,    5,   0,   0,   0,   2,   2,   4)  # Confirmed K. pneumoniae from MB
)

# Chi-square statistics and p-values for K. pneumoniae
chi_tbl_kp <- tibble(
  Variables = c("Overall", "Origin of Sample", "Region", "Season", "District"),
  Chi2      = c("-", "32.031", "1.762", "1.206", "92.821"),
  p_value   = c("-", "<0.001", "0.184", "0.272", "<0.001")
)

# Compute whole integer adjusted n, prevalence %, and 95% CIs
Final_KP_DATA <- df_kp_raw %>%
  mutate(
    `ĥ (n)`  = round((CA / MA) * A + (CB / MB) * B),
    p_hat     = `ĥ (n)` / N,
    p_pct_num = p_hat * 100,
    `ĥ (%)`  = fmt_num(p_pct_num),
    se        = sqrt(p_hat * (1 - p_hat) / N),
    ci_low    = pmax(0, (p_hat - 1.96 * se) * 100),
    ci_high   = pmin(100, (p_hat + 1.96 * se) * 100),
    `95% CI`  = paste0(fmt_num(ci_low), " - ", fmt_num(ci_high))
  ) %>%
  left_join(chi_tbl_kp, by = "Variables") %>%
  select(Variables, Values, N, A, B, MA, CA, MB, CB, `ĥ (n)`, `ĥ (%)`, `95% CI`, Chi2, p_value)

# Export files to Results directory (.csv, .tsv, .rds)
write_csv(Final_KP_DATA, file = "Results/adjusted_prevalence_kp_summary.csv")
write_tsv(Final_KP_DATA, file = "Results/adjusted_prevalence_kp_summary.tsv")
saveRDS(Final_KP_DATA, file = "Results/adjusted_prevalence_kp_summary.rds")

# Format table for Word document (blank duplicate Variable, Chi2, and p-value labels)
df_kp_word <- Final_KP_DATA %>%
  mutate(
    Variables = if_else(duplicated(Variables), "", Variables),
    Chi2      = if_else(Variables == "", "", Chi2),
    p_value   = if_else(Variables == "", "", p_value)
  )

# Create publication-ready flextable
ft_kp <- flextable(df_kp_word) %>%
  set_header_labels(
    Variables = "Variables",
    Values    = "Values",
    N         = "N",
    A         = "A",
    B         = "B",
    MA        = "MA",
    CA        = "CA",
    MB        = "MB",
    CB        = "CB",
    `ĥ (n)`   = "ĥ (n)",
    `ĥ (%)`   = "ĥ (%)",
    `95% CI`  = "95% CI",
    Chi2      = "χ²",
    p_value   = "p-value"
  ) %>%
  add_footer_lines(
    values = c(
      "Abbreviations:",
      "3GCR = 3rd Generation Cephalosporin-Resistant.",
      "N = Total study population analyzed.",
      "A = Presumptive 3GCR K. pneumoniae isolates; B = Presumptive 3GCR E. coli isolates.",
      "MA = Selected presumptive 3GCR K. pneumoniae isolates tested on VITEK MS.",
      "CA = Confirmed true K. pneumoniae isolates from group MA.",
      "MB = Selected presumptive 3GCR E. coli isolates tested on VITEK MS.",
      "CB = Confirmed true K. pneumoniae isolates from group MB.",
      "ĥ (n) = Adjusted number of confirmed K. pneumoniae cases (rounded to whole integer).",
      "ĥ (%) = Adjusted prevalence percentage; 95% CI = 95% Confidence Interval for adjusted prevalence (expressed as lower - upper bound numbers).",
      "χ² = Chi-square test statistic; p-value = Statistical significance level (p < 0.05 indicates significant proportion difference)."
    )
  ) %>%
  theme_vanilla() %>%
  align(align = "center", part = "all") %>%
  align(j = c("Variables", "Values"), align = "left", part = "all") %>%
  autofit()

# Save Word document (.docx) to Results folder
doc_kp <- read_docx()
doc_kp <- body_add_flextable(doc_kp, value = ft_kp)
print(doc_kp, target = "Results/adjusted_prevalence_kp_summary.docx")
cat("K. pneumoniae adjusted prevalence files generated successfully in 'Results' folder!\n")
######################################################
# 22/9/2026
# Checking Isolates with CRO < 23, as 3rdGCR
# Filter criteria and calculate individual + aggregate counts

# Define and guarantee directory creation
if (!dir.exists("Results")) {
  dir.create("Results")
}

# 30.09.26 Here you have forgotten that the dataset Final_JoinedDATA you refer to have changed - you have been overwriting it! 
# So either you need to change the names of the datasets you create above or we do just reimport it again below!
Final_JoinedDATA <- read_csv ("data/CLEANED_DATA/Final_JoinedDATA.csv")

spec(Final_JoinedDATA)
names(Final_JoinedDATA) 


# Clean data and define base CGR3 dataset
base_cgr3 <- Final_JoinedDATA %>%
  mutate(across(where(is.character), str_trim)) %>%
  filter(
    str_detect(PROTOCOL, regex("CGR3", ignore_case = TRUE)),
    str_detect(VITEK_MS_Results, regex("coli|pneumoniae", ignore_case = TRUE))
  ) %>%
  mutate(
    Confirmed_Species = case_when(
      str_detect(VITEK_MS_Results, regex("coli", ignore_case = TRUE)) ~ "E. coli",
      str_detect(VITEK_MS_Results, regex("pneumoniae", ignore_case = TRUE)) ~ "K. pneumoniae",
      TRUE ~ VITEK_MS_Results
    ),
    Initial_Isolate = case_when(
      str_detect(Isolate, regex("coli", ignore_case = TRUE)) ~ "E. coli",
      str_detect(Isolate, regex("pneumoniae", ignore_case = TRUE)) ~ "K. pneumoniae",
      TRUE ~ Isolate
    ),
    Is_Reclassified = Initial_Isolate != Confirmed_Species
  )

# Total unique TVLA_ID tested by species
tested_summary <- base_cgr3 %>%
  group_by(Confirmed_Species) %>%
  summarise(Tested_TVLA_ID = n_distinct(TVLA_ID), .groups = "drop")

# Confirmed cases (CRO_ED30 < 23) + Reclassified Count
confirmed_summary <- base_cgr3 %>%
  filter(CRO_ED30 < 23) %>%
  group_by(Confirmed_Species) %>%
  summarise(
    Confirmed_Isolates = n(),
    Reclassified_Count = sum(Is_Reclassified, na.rm = TRUE),
    .groups = "drop"
  )

# Combine into final summary table with Total Aggregate row
final_table <- tested_summary %>%
  left_join(confirmed_summary, by = "Confirmed_Species") %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
  bind_rows(
    summarise(
      .,
      Confirmed_Species = "Total Aggregate",
      Tested_TVLA_ID = n_distinct(base_cgr3$TVLA_ID),
      Confirmed_Isolates = sum(Confirmed_Isolates),
      Reclassified_Count = sum(Reclassified_Count)
    )
  )

# Format flextable
ft <- final_table %>%
  flextable() %>%
  set_header_labels(
    Confirmed_Species = "Confirmed Species",
    Tested_TVLA_ID = "Tested Samples (n)",
    Confirmed_Isolates = "Confirmed Isolates (n)",
    Reclassified_Count = "Reclassified Isolates (n)"
  ) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  italic(i = ~ Confirmed_Species != "Total Aggregate", j = "Confirmed_Species") %>%
  bold(part = "header") %>%
  bold(i = ~ Confirmed_Species == "Total Aggregate", part = "body") %>%
  align(j = "Confirmed_Species", align = "left", part = "all") %>%
  align(j = c("Tested_TVLA_ID", "Confirmed_Isolates", "Reclassified_Count"), align = "center", part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1), part = "header") %>%
  hline_bottom(border = fp_border(color = "black", width = 1.5), part = "body") %>%
  padding(padding = 4, part = "all") %>%
  autofit()

# Save Word Document safely
doc <- read_docx() %>%
  body_add_par("Table 1. Summary of Isolate Confirmation and Reclassification under CGR3 Protocol", style = "table title") %>%
  body_add_flextable(ft) %>%
  body_add_par("Note: Confirmed isolates reflect CRO_ED30 < 23 mm cutoff confirmed by VITEK MS.", style = "Normal")

print(doc, target = "Results/Summary of Isolate Confirmation and Reclassification.docx")

write.csv(final_table, "Results/Summary_of_Isolate_Confirmation_and_Reclassification.csv", row.names = FALSE)

cat("Summary of Isolate Confirmation and Reclassification files generated successfully in 'Results' folder!\n")
################################################################################
