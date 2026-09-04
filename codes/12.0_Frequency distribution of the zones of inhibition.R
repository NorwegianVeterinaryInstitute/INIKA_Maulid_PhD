library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)
library(tidyr)
library(readr)

joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
spec(joined_data)

joined_datam <- joined_data %>%
  select(INIKA_ID,REGION,DISTRICT,SEASON, ORIGIN_OF_SAMPLE,Isolate,VITEK_MS_Results, ESBL_Selection,ESBL,
         AMX_ED10,AZM_ED15,CRO_ED30,CIP_ED5,DOX_ED30,FLR_ED30,
         GEN_ED10,MEM_ED10,OXY_ED30,POL_ED300,SXT_ED1_2,
         CTX_ED5,CTC_ED30)
  
original_spelling <- c("AMX_ED10", 
                       "AZM_ED15", 
                       "CRO_ED30", 
                       "CIP_ED5",
                       "DOX_ED30", 
                       "FLR_ED30",
                       "GEN_ED10", 
                       "MEM_ED10", 
                       "OXY_ED30",
                       "POL_ED300", 
                       "SXT_ED1_2",
                       "CTX_ED5",
                       "CTC_ED30")

# correct names 
tested_antibiotics  <- c("Amoxicillin", 
                         "Azithromycin",
                         "Ceftriaxone",
                         "Ciprofloxacin",
                         "Doxyciline",
                         "Florfenicol",
                         "Gentamicin", 
                         "Meropenem", 
                         "Oxytetracyline",
                         "Polymixin B (PB)",
                         "Sulfamethoxazole/Trimethoprim", 
                         "Cefotaxime",
                         "Cefotaxime + Clavulanic acid")


# checking same length
length(original_spelling) == length(tested_antibiotics)                         

# Changing the names - copy the joined table
joined_datamm <- joined_datam

# selecting the names of the columns and changing with the corresponding new names
names(joined_datamm)[names(joined_datamm) %in% original_spelling] <- tested_antibiotics
names(joined_datamm)


# Filter the data to include ONLY the isolates that match the following criteria:
# 1. Isolate is "E.coli" OR "K.pneumoniae"
# 2. VITEK_MS_Results is "Escherichia coli"
# 3. ESBL is 1
Ecoli_filtered_for_zones <- joined_datamm %>%
  filter(
    (Isolate == "E.coli" | Isolate == "K.pneumoniae") &
      VITEK_MS_Results == "Escherichia coli" &
      ESBL == 1
  )

Total_Matching_Isolates <- nrow(Ecoli_filtered_for_zones)

# Output the total number of isolates used in the analysis
cat("--- Filtering Status ---\n")
cat("Total isolates used for zone analysis (ESBL-positive and matched criteria):", Total_Matching_Isolates, "\n")


# --- 1. Raw Zone of Inhibition Distribution Analysis ---

# a. Transforming to long format
Ecoli_all_long_filtered <-
  Ecoli_filtered_for_zones %>%
  # Convert zone columns to numeric
  mutate(across(all_of(tested_antibiotics), as.numeric)) %>%
  pivot_longer(
    cols = all_of(tested_antibiotics),
    names_to = "substance",
    values_to = "value" # Zone of inhibition in mm
  ) %>%
  # Remove NA values (isolates not tested for that substance)
  filter(!is.na(value))

# b. Counting unique zone values (Frequency)
Ecoli_df_count <- Ecoli_all_long_filtered %>%
  count(substance, value)

# c. Count to wide format 
# Unique zone values (e.g., 6, 7, 8... mm) become column headers
Ecoli_df_count_wide <-
  Ecoli_df_count %>%
  arrange(value) %>% 
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  ) 

# d. Percent table calculation
# Calculates the percentage of total isolates tested that resulted in that specific zone value
Ecoli_percent_table_filtered <- 
  Ecoli_df_count_wide %>%
  # Computes percentages row-wise (within each antibiotic)
  rowwise() %>% 
  mutate(across(-substance, ~ . / sum(c_across(-substance)) * 100)) %>%
  ungroup()

# Display the final distribution table
cat("\n--- Raw Zone of Inhibition Percentage Distribution (Filtered Isolates) ---\n")
print(Ecoli_percent_table_filtered)

# Save the file
write_tsv(Ecoli_percent_table_filtered, "Results/131_Obs_Ecoli_antibiotics_percent_table.tsv")
################################################################################
# 1. Filtering Logic
Ecoli_filtered_for_zones <- joined_datamm %>%
  filter(
    Isolate == "E.coli" & 
      VITEK_MS_Results == "Escherichia coli" & 
      Cefotaxime < 21
  )

# 2. Data Cleaning and Long Format
Ecoli_long <- Ecoli_filtered_for_zones %>%
  mutate(across(all_of(tested_antibiotics), as.numeric)) %>%
  pivot_longer(
    cols = all_of(tested_antibiotics),
    names_to = "substance",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(substance = str_to_title(str_trim(substance)))

# 3. Calculate Summary: Total Isolates and a single 95%CI column
Ecoli_summary <- Ecoli_long %>%
  group_by(substance) %>%
  summarise(
    Total_Isolates = n(),
    # Calculating the 95% percentile interval for the zones
    L = round(quantile(value, 0.025, na.rm = TRUE), 1),
    U = round(quantile(value, 0.975, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(
    # Clean formatting: remove .0 from the CI bounds
    L_str = sub("\\.0$", "", sprintf("%.1f", L)),
    U_str = sub("\\.0$", "", sprintf("%.1f", U)),
    `95%CI` = paste0(L_str, "-", U_str)
  ) %>%
  select(substance, Total_Isolates, `95%CI`)

# 4. Calculate Zone Distribution (Percentages with clean formatting)
Ecoli_dist_wide <- Ecoli_long %>%
  count(substance, value) %>%
  group_by(substance) %>%
  mutate(pct = round((n / sum(n)) * 100, 1)) %>%
  ungroup() %>%
  mutate(
    # Clean formatting: remove .0 from percentages
    pct_str = sub("\\.0$", "", sprintf("%.1f", pct)),
    value = as.numeric(value)
  ) %>%
  select(substance, value, pct_str) %>%
  arrange(value) %>%
  pivot_wider(
    names_from = value,
    values_from = pct_str,
    values_fill = "0"
  )

# 5. Final Join: One row per antibiotic
Ecoli_final_report <- Ecoli_summary %>%
  left_join(Ecoli_dist_wide, by = "substance")

# Output results
cat("\n--- Final Table: Single 95%CI Column & Clean Formatting ---\n")
print(Ecoli_final_report)

# Save as Excel (Best for reports)
write_xlsx(Ecoli_final_report, "Results/20_12_25Ecoli_final_report.xlsx")
################################################################################
# ## Klebsiella ----
# ### Transforming to long format ----
# Klebsiella_all_long <-
#   Klebsiella_all %>%
#   mutate_at(vars(all_of(tested_antibiotics)), as.numeric) %>%
#   pivot_longer(cols = all_of(tested_antibiotics),
#                names_to = "substance",
#                values_to = "value")
# 
# ### counting ----
# Klebsiella_df_count <- 
#   Klebsiella_all_long %>%
#   count(substance, value)
# 
# Klebsiella_df_count
# 
# ### count to wide format ----
# Klebsiella_df_count_wide <-
#   Klebsiella_df_count %>%
#   # sorting the values in increasing order 
#   arrange(value) %>% 
#   pivot_wider(
#     names_from = value,
#     values_from = n,
#     values_fill = 0) 
# 
# Klebsiella_df_count_wide  
# 
# ### percent table ----
# Klebsiella_percent_table <- 
#   Klebsiella_df_count_wide %>%
#   # Does computation by row instead of column (default)
#   rowwise() %>%
#   mutate(across(-substance, ~ . / sum(c_across(-substance)) * 100))
# 
# Klebsiella_percent_table
# 
# ### Export tables ----
# write_tsv(Klebsiella_df_count_wide, "Results/Klebsiella_antibiotics_df_count_wide.tsv")
# write_tsv(Klebsiella_percent_table, "Results/Klebsiella_antibiotics_percent_table.tsv")
################################################################################
## Here is the computation for Klebsiella pneumoniae

# selecting the names of the columns and changing with the corresponding new names
names(joined_datamm)[names(joined_datamm) %in% original_spelling] <- tested_antibiotics
names(joined_datamm)

# Filter the data to include ONLY the isolates that match the following criteria:
# 1. Isolate is "E.coli" OR "K.pneumoniae"
# 2. VITEK_MS_Results is "Klebsiella pneumoniae"
# 3. ESBL is 1
Kpn_filtered_for_zones <- joined_datamm %>%
  filter(
    (Isolate.x == "E.coli" | Isolate.x == "K.pneumoniae") &
      VITEK_MS_Results == "Klebsiella pneumoniae" &
      ESBL == 1
  )

Total_Matching_Isolates <- nrow(Kpn_filtered_for_zones)

# Output the total number of isolates used in the analysis
cat("--- Filtering Status ---\n")
cat("Total isolates used for zone analysis (ESBL-positive and matched criteria):", Total_Matching_Isolates, "\n")


# --- 1. Raw Zone of Inhibition Distribution Analysis ---

# a. Transforming to long format
Kpn_all_long_filtered <-
  Kpn_filtered_for_zones %>%
  # Convert zone columns to numeric
  mutate(across(all_of(tested_antibiotics), as.numeric)) %>%
  pivot_longer(
    cols = all_of(tested_antibiotics),
    names_to = "substance",
    values_to = "value" # Zone of inhibition in mm
  ) %>%
  # Remove NA values (isolates not tested for that substance)
  filter(!is.na(value))

# b. Counting unique zone values (Frequency)
Kpn_df_count <- Kpn_all_long_filtered %>%
  count(substance, value)

# c. Count to wide format 
# Unique zone values (e.g., 6, 7, 8... mm) become column headers
Kpn_df_count_wide <-
  Kpn_df_count %>%
  arrange(value) %>% 
  pivot_wider(
    names_from = value,
    values_from = n,
    values_fill = 0
  ) 

# d. Percent table calculation
# Calculates the percentage of total isolates tested that resulted in that specific zone value
Kpn_percent_table_filtered <- 
  Kpn_df_count_wide %>%
  # Computes percentages row-wise (within each antibiotic)
  rowwise() %>% 
  mutate(across(-substance, ~ . / sum(c_across(-substance)) * 100)) %>%
  ungroup()

# Display the final distribution table
cat("\n--- Raw Zone of Inhibition Percentage Distribution (Filtered Isolates) ---\n")
print(Kpn_percent_table_filtered)

# Save the file
write_tsv(Kpn_percent_table_filtered, "Results/54_Obs_Kpn_antibiotics_percent_table.tsv")
################################################################################