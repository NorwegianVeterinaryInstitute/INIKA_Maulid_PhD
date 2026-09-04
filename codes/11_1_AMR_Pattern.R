#############################################################################
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(Hmisc)
library(magrittr)

# NOT WORKED AS I EPECTED #####

# # Importing the file 
# 
# joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
# spec(joined_data)

# 
# antibiotic_map_final <- c(
#   "Amoxicillin" = "AMX_ED10",
#   "Azithromycin" = "AZM_ED15",
#   "Ceftriaxone" = "CRO_ED30",
#   "Ciprofloxacin" = "CIP_ED5",
#   "Doxycycline" = "DOX_ED30",
#   "Florfenicol" = "FLR_ED30",
#   "Gentamicin" = "GEN_ED10",
#   "Meropenem" = "MEM_ED10",
#   "Oxytetracycline" = "OXY_ED30",
#   "Polymyxin B (PB)" = "POL_ED300",
#   "Sulphonamides/Trimethoprim" = "SXT_ED1_2",
#   "Cefotaxime" = "CTX_ED5",
#   "Cefotaxime+clavulanic acid" = "CTC_ED30"
# )
# abx_names <- names(antibiotic_map_final)
# 
# # A. Clean up ALL column names using trimws() to remove trailing spaces (The Fix!)
# names(joined_data) <- trimws(names(joined_data))
# 
# # B. Rename the columns using base R for stability
# old_names_to_find <- unname(antibiotic_map_final)
# new_names_to_assign <- names(antibiotic_map_final)
# col_indices <- match(old_names_to_find, names(joined_data))
# current_names <- names(joined_data)
# current_names[col_indices] <- new_names_to_assign
# names(joined_data) <- current_names
# 
# # C. Robust Data Type Coercion (Fixes the "list object" error)
# # Function to convert column contents safely
# convert_to_double_safe <- function(x) {
#   x %>%
#     as.character() %>%
#     trimws() %>%
#     as.numeric()
# }
# 
# # Apply the conversion to all antibiotic columns using base R's lapply
# joined_data[abx_names] <- lapply(joined_data[abx_names], convert_to_double_safe)
# 
# # D. Filter for E. coli
# ecoli_data <- joined_data %>%
#   filter(VITEK_MS_Results == "Escherichia coli")
# 
# # ----Part 2: ESC Resistance E.coli calculation ----
# # Calculate the unique counts
# denominator_n <- ecoli_data %>%
#   distinct(INIKA_ID.x) %>%
#   nrow()
# 
# numerator_n <- ecoli_data %>%
#   filter(ESBL == 1) %>%
#   distinct(INIKA_ID.x) %>%
#   nrow()
# 
# # Calculate the proportion and 95% CI (Agresti-Coull method)
# prop_esc_resistance <- data.frame(
#   Proportion = numerator_n / denominator_n,
#   Numerator = numerator_n,
#   Denominator = denominator_n
# ) %>%
#   mutate(
#     CI_result = list(Hmisc::binconf(Numerator, Denominator, method = "wilson")),
#     Lower_CI = map_dbl(CI_result, ~ .x[2]),
#     Upper_CI = map_dbl(CI_result, ~ .x[3]),
#     `95% CI` = paste0("(", round(Lower_CI * 100, 2), "%, ", round(Upper_CI * 100, 2), "%)")
#   ) %>%
#   select(Proportion, `95% CI`, Numerator, Denominator)
# 
# print("### 2. True ESC Resistance E. coli Proportion (Unique INIKA_ID.x)")
# print(prop_esc_resistance)
# 
# # ---- Part 3: Breakpoint -----
# 
# # Inserting breakpoints (CLSI/EUCAST) in mm
# breakpoints <- tribble(
#   ~Antibiotic, ~S_BP, ~R_BP, # S_BP: Susceptible breakpoint (>=), R_BP: Resistant breakpoint (<=)
#   "Amoxicillin", NA, NA,
#   "Azithromycin", 13, 12,
#   "Ceftriaxone", 23, 22,
#   "Ciprofloxacin", 26, 25,
#   "Doxycycline", 14, 13,
#   "Florfenicol", 18, 17,
#   "Gentamicin", 18, 17,
#   "Meropenem", 23, 22,
#   "Oxytetracycline", 15, 14,
#   "Polymyxin B (PB)", NA, NA,
#   "Sulphonamides/Trimethoprim", 16, 15,
#   "Cefotaxime", 26, 25,
#   "Cefotaxime+clavulanic acid", NA, NA
# )
# 
# # Function to apply breakpoints and categorize
# categorize_susceptibility <- function(zone_diameter, S_BP, R_BP) {
#   case_when(
#     is.na(zone_diameter) ~ NA_character_,
#     zone_diameter >= S_BP ~ "S", # Susceptible
#     zone_diameter <= R_BP ~ "R", # Resistant
#     TRUE ~ "I" # Intermediate
#   )
# }
# ecoli_categorized <- ecoli_data
# 
# 
# ## HERE 6/10/2025----
# for (abx in abx_names) {
#   bp_row <- breakpoints %>% filter(Antibiotic == abx)
#   S_val <- bp_row$S_BP[1]
#   R_val <- bp_row$R_BP[1]
#   
#   # Create a new column for the S/I/R category
#   ecoli_categorized <- ecoli_categorized %>%
#     mutate(
#       !!paste0(abx, "_Category") := categorize_susceptibility(.data[[abx]], S_val, R_val)
#     )
# }
# 
# #---- Part 4: Computing pivot table ----
# category_cols <- paste0(abx_names, "_Category")
# group_vars <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x")
# 
# final_pivot_table <- ecoli_categorized %>%
#   select(all_of(group_vars), all_of(category_cols)) %>%
#   
#   # Pivot the data longer
#   pivot_longer(
#     cols = all_of(category_cols),
#     names_to = "Antibiotic",
#     values_to = "Category"
#   ) %>%
#   mutate(
#     Antibiotic = str_remove(Antibiotic, "_Category$")
#   ) %>%
#   filter(!is.na(Category)) %>%
#   drop_na(all_of(group_vars)) %>%
#   
#   # Aggregate: Count Frequency and Total Tested
#   group_by(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE.x, Antibiotic) %>%
#   mutate(`Total tested` = n()) %>%
#   group_by(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE.x, Antibiotic, Category, `Total tested`) %>%
#   summarise(
#     Frequency = n(),
#     .groups = 'drop_last'
#   ) %>%
#   
#   # Calculate Percentage and 95% CI (Wilson method used here too)
#   mutate(Percentage = (Frequency / `Total tested`) * 100) %>%
#   rowwise() %>%
#   mutate(
#     ci_result = list(Hmisc::binconf(Frequency, `Total tested`, method = "wilson")),
#     `95% CI` = paste0("(", round(ci_result[2] * 100, 2), "%, ", round(ci_result[3] * 100, 2), "%)")
#   ) %>%
#   ungroup() %>%
#   select(-ci_result) %>%
#   
#   # Format to the requested Variable/Varue structure
#   arrange(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE.x, Antibiotic, Category) %>%
#   mutate(
#     Region_Display = if_else(!duplicated(REGION.x), REGION.x, NA_character_),
#     Season_Display = if_else(!duplicated(paste(REGION.x, SEASON.x)), SEASON.x, NA_character_),
#     Origin_Display = if_else(!duplicated(paste(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE.x)), ORIGIN_OF_SAMPLE.x, NA_character_)
#   ) %>%
#   pivot_longer(
#     cols = c(Region_Display, Season_Display, Origin_Display),
#     names_to = "Variable_Temp",
#     values_to = "Varue"
#   ) %>%
#   filter(!is.na(Varue)) %>%
#   mutate(
#     Variable = case_when(
#       Variable_Temp == "Region_Display" ~ "Region",
#       Variable_Temp == "Season_Display" ~ "Season",
#       Variable_Temp == "Origin_Display" ~ "Origin of Sample" # Fixed in prior response
#     )
#   ) %>%
#   # Merge the pivot back to the original frequency data
#   group_by(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE.x) %>%
#   fill(Variable, Varue, .direction = "down") %>%
#   ungroup() %>%
#   
#   # Final selection and cleaning
#   mutate(Percentage = round(Percentage, 2)) %>% # <--- THE FIX for 'object 'Percentage' not found'
#   select(
#     Variable, Varue, Antibiotic, Category, Frequency, `Total tested`,
#     Percentage, `95% CI`
#   ) %>%
#   distinct() %>%
#   arrange(Variable, Varue, Antibiotic, Category)
# 
# print("---")
# print("### Grouped Susceptibility Pivot Table (Fixed `select()` error)")
# print(final_pivot_table, n = 50)
