################################################################################
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)
library(tidyr)
library(readr)

# Importing the file
joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
spec(joined_data)


# --------------------------------------------------------------------------
DATASET_COLUMN_NAMES <- list(
  id = "INIKA_ID.x", 
  age = "Age_yrs",
  region = "REGION.x", 
  district = "DISTRICT.x", 
  origin = "ORIGIN_OF_SAMPLE.x", 
  gender = "GENDER",
  season = "SEASON.x"
)
# --------------------------------------------------------------------------

# Assign column names to easy-to-use variables
ID_COL <- DATASET_COLUMN_NAMES$id
AGE_COL <- DATASET_COLUMN_NAMES$age
REGION_COL <- DATASET_COLUMN_NAMES$region
DISTRICT_COL <- DATASET_COLUMN_NAMES$district
ORIGIN_COL <- DATASET_COLUMN_NAMES$origin
GENDER_COL <- DATASET_COLUMN_NAMES$gender
SEASON_COL <- DATASET_COLUMN_NAMES$season

# Define the final variables for the pivot table analysis
variables_to_analyze <- c(REGION_COL, DISTRICT_COL, ORIGIN_COL, GENDER_COL, "Age_group", SEASON_COL)


# --- PART 1: Robust Type Conversion and Pre processing ---

joined_data_processed <- joined_data %>%
  mutate(
    # Explicitly convert all categorical and ID columns to CHARACTER for safety
    !!ID_COL := as.character(!!sym(ID_COL)),
    !!REGION_COL := as.character(!!sym(REGION_COL)),
    !!DISTRICT_COL := as.character(!!sym(DISTRICT_COL)),
    !!ORIGIN_COL := as.character(!!sym(ORIGIN_COL)),
    !!GENDER_COL := as.character(!!sym(GENDER_COL)),
    !!SEASON_COL := as.character(!!sym(SEASON_COL)),
    
    # 1. Group the 'Age' variable
    Age_group = cut(
      !!sym(AGE_COL),
      # Irregular breaks specified: 10-17, 18-27, 28-37, etc.
      breaks = c(10, 18, 28, 38, 48, 58, 68, Inf), 
      right = FALSE, 
      labels = c("10 - 17", "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 +"),
      exclude.lowest = TRUE 
    )
  )

# --- PART 2: Frequency and Percentage Calculation (Based on Unique ID) ---

# STEP 1: Calculate the total number of unique IDs for the percentage denominator
total_unique_ids <- joined_data_processed %>%
  summarise(Total_Unique = n_distinct(!!sym(ID_COL), na.rm = TRUE)) %>% 
  pull(Total_Unique)

cat(paste("Total Unique IDs (for percentage denominator):", total_unique_ids, "\n\n"))

# Function to calculate frequency and percentage for a single variable
calculate_summary_by_id <- function(data, var_name, total_ids, id_col) {
  data %>%
    group_by(!!sym(var_name)) %>%
    summarise(
      Frequency = n_distinct(!!sym(id_col), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Variable = var_name,
      Percentage = (Frequency / total_ids) * 100
    ) %>%
    rename(Value = !!sym(var_name)) %>%
    select(Variable, Value, Frequency, Percentage)
}

# STEP 2: Loop through all variables, apply the function, and create the base table
final_pivot_table_base <- map_dfr(variables_to_analyze, 
                                  ~ calculate_summary_by_id(joined_data_processed, .x, total_unique_ids, ID_COL)) %>%
  # Clean up variable names for final output
  mutate(
    Variable = case_match(
      Variable,
      ORIGIN_COL ~ "Origin of the sample",
      "Age_group" ~ "Age group",
      "REGION_COL" ~ "REGION",
      "DISTRICT_COL" ~ "DISTRICT",
      "SEASON_COL" ~ "SEASON",
      .default = Variable
    )
  )

# --- PART 3: Calculate and Insert Median/IQR Row ---

# Calculate Q1, Median, and Q3 for the entire numerical 'Age' variable
age_stats <- joined_data_processed %>%
  summarise(
    Median = median(!!sym(AGE_COL), na.rm = TRUE),
    Q1 = quantile(!!sym(AGE_COL), 0.25, na.rm = TRUE),
    Q3 = quantile(!!sym(AGE_COL), 0.75, na.rm = TRUE)
  )

# 1. Create the IQR summary row
iqr_row <- tibble(
  Variable = "Age group (Median [Q1-Q3])",
  Value = as.character(paste0(round(age_stats$Median, 1), " (", 
                              round(age_stats$Q1, 1), " - ", 
                              round(age_stats$Q3, 1), ")")),
  Frequency = NA_real_, 
  Percentage = NA_real_ 
) 

# 2. Prepare the base table by enforcing character/double types just before binding
final_pivot_table_base_ready <- final_pivot_table_base %>%
  mutate(
    Percentage = round(Percentage, 2),
    Value = as.character(Value),      # Enforce Character
    Frequency = as.double(Frequency)  # Enforce Double/Numeric
  )

# 3. Combine the tables
Demographic_final_pivot_table <- final_pivot_table_base_ready %>%
  filter(Variable != "Age group") %>%
  bind_rows(
    iqr_row, 
    filter(final_pivot_table_base_ready, Variable == "Age group")
  ) %>%
  # Reorder the final table for presentation
  arrange(match(Variable, 
                c(variables_to_analyze %>% 
                    replace(which(. == "Age_group"), "Age group (Median [Q1-Q3])"))
  )
  )

# --- PART 4: Display the Final Table ---

cat("\n--- Final Combined Frequency and Percentage Pivot Table ---\n")
print(Demographic_final_pivot_table)
# Save the table as CSV file for reporting

 write_csv(Demographic_final_pivot_table, "Results/Demographic_final_pivot_table.csv")
write_tsv(Demographic_final_pivot_table, "Results/Demographic_final_pivot_table.tsv") 
##################################################################################