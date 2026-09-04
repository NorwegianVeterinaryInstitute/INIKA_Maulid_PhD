# k y---
  title: "Descriptive analysis of AMR data"
author: "Madelaine Norström"
date: "2025-10-14"
output:
  pdf_document: default
html_document: default
---
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)



# First you will need to calculate the occurrence (Prevalence) of the different ESC resistant bacterial isolates, also where you have #not find #anything!
# 
# How many resistant isolates have you find in total/ per season/ per region ( depending on at what level you think you will report!)
# 
# If not many samples you can not break down this to to small strata!
#   
#   Define the criteria for what is an ESC- resistant E.coli/ Carbapenemase resistant E.coli etc.
# 
# If 0 it is 0. #Then how many unique samples can be included as the total number of samples. #For example if sampled from a total of 430 farms and #no Carbapenemase resistant isolates were found this will be a special case. #You should need to report the 95% Confidence interval for each result as well. #The 95% confidence intervals were calculated using the exact binomial test in R version ( Check which version you run!) Copyright (C) 2025 (The R Foundation for Statistical Computing Platform).#
# 
# Below is the function you run:
  
 #{r}

binom.test(0,430, p = 0.5, "two.sided", conf.level = 0.95)
binom.test(0,430, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100


# This will be reported as follows: None of the samples were positive for Carbapeneamse resistant E.coli [95% CI: 0.0 – 1.1]
# 
# If you have found 10 ESC resistant Klebsiella from the unique farms you will calculate accordingly:

binom.test(10,430, p = 0.5, "two.sided", conf.level = 0.95)
binom.test(10,430, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100

# This (example) will be reported as follows: 2.3% [95% CI: 1.1 – 4.2] of the samples were positive for ESC-resistant Klebsiella.
# 
# For the % Resistances to different antibiotics you will only use the total number of the isolates you have within each category, and of course it might be confusing if you first detected more ESC resistant E.coli from the selective plates and then you have less ESC resistant E.coli according to the criteria <22mm for Cefotaxime. But it should be the last criteria that will decide how many you have! We may want to change the criteria for defining the isolates as resistant after checking the distributions. This is therefore also very important.I think this should be discussed.
# 
# Example: But for now if you have 110 ESC resistant isolates of for example E.coli according only to be growing on the selective plate and where the VITEK has also confirmed that it is an E.coli, but from these only 90 isolates are resistant to Cefotaxime. Then the total number of ESC resistant E.coli will be 90. If then 10 of these are resistant to for a substance "X", you calculate the % Resistance to Substance X as follows: 10/90, but you also want to include the 95% confidence interval so then you do as shown here:

binom.test(10,90, p = 0.5, "two.sided", conf.level = 0.95)
binom.test(10,90, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100

# This means that the % resistance are 11.1% but it can be anywhere between 5.5% to 19.5%.
# The example below use the original data from Maulid, so you need to change and read in your correct dataset as well as change the variable names to reflect your own dataset: like ID and ORGANISM to your throughout the code, so read carefully and change accordingly.

joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")

# Renaming the column
joined_data <- joined_data %>%
  rename(INIKA_ID = INIKA_ID.x,
         REGION = REGION.x,
         SEASON = SEASON.x,
         ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
         Isolate = Isolate.x)

## Importing the E.coli_ECOFF_EUCAST break point file 
ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx")

EPI_CUTOFF <- ECOFF_EUCAST_BREAK_POINT


calculate_ecoli_esbl_susceptibility <- function(joined_data, EPI_CUTOFF) {
  
  # 1. Select relevant columns and rename drug columns for clarity
  # FIX: Combine select and rename for clarity, and use clean column names.
  Data <- joined_data %>%
    select(
      INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
      Amoxicillin = AMX_ED10,
      Azithromycin = AZM_ED15,
      Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5,
      Doxycycline = DOX_ED30,
      Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10,
      Meropenem = MEM_ED10,
      Oxytetracycline = OXY_ED30,
      # Renaming to a syntactically valid name:
      Polymyxin_B_PB = POL_ED300,
      # Renaming to a syntactically valid name:
      Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
      Cefotaxime = CTX_ED5,
      # Renaming to a syntactically valid name:
      Cefotaxime_ClavulanicAcid = CTC_ED30
    )
  
  # 2. Pivot the drug columns into long format
  # The cols = -c(...) correctly excludes the listed identifier/grouping columns.
  AMR_Data <- Data %>%
    pivot_longer(
      # The listed columns are the ones NOT to pivot (ID and Grouping Vars)
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    )
  
  # 3. Join with the cutoff table
  # This step requires 'EPI_CUTOFF' to have an 'Antimicrobial_substance' column
  # matching the drug names defined in step 1.
  AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
  
  # Define the list of grouping variables for the final summary
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
  
  # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
  
  Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
    
    # Filter for E. coli AND true ESBL (ESBL == 1)
    filter(VITEK_MS_Results == "Escherichia coli",
           ESBL == 1) %>%
    
    # Clean the S and R cutoff columns and convert all key values to numeric
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
      
      Measured_Zone = as.numeric(str_trim(value)),
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    
    # Implement the breakpoint logic for zone diameters (in mm)
    mutate(
      Categorical_Result = case_when(
        Measured_Zone >= S_Cutoff ~ "S", # Susceptible (S): Zone diameter >= S breakpoint
        Measured_Zone <= R_Cutoff ~ "R", # Resistant (R): Zone diameter <= R breakpoint
        Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I", # Intermediate (I)
        TRUE ~ NA_character_ # Cases where data is missing or doesn't fit
      )
    ) %>%
    
    # Ensure only one result is counted per isolate for each antimicrobial test
    distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- Step 5: Define the Summary Function (Unchanged) ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    Total_Count <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')
    
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
      
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var
      ) %>%
      
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Categorical_Result,
        Frequency = N,
        Total_Tests,
        Percentage
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 6: Apply the function to all grouping variables and combine results (Unchanged) ---
  ECO_ESBL_Susceptibility_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  return(ECO_ESBL_Susceptibility_Summary)
}

# Execution:
ECO_ESBL_Susceptibility_Summary <- calculate_ecoli_esbl_susceptibility(joined_data, EPI_CUTOFF)

# Save the file
write_tsv(ECO_ESBL_Susceptibility_Summary, "Results/135_Obs_ECO_ESBL_Susceptibility_Summary-9-11-2025.tsv")
######################################################################
#AMR for the Confirmed ESC Resistance E.coli



calculate_ecoli_esbl_susceptibility <- function(joined_data, EPI_CUTOFF) {
  
  # 1. Select relevant columns and rename drug columns for clarity
  # NOTE: The 'distinct(INIKA_ID, .keep_all = TRUE)' line is intentionally omitted here
  Data <- joined_data %>%
    select(
      INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
      Amoxicillin = AMX_ED10, Azithromycin = AZM_ED15, Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5, Doxycycline = DOX_ED30, Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10, Meropenem = MEM_ED10, Oxytetracycline = OXY_ED30,
      Polymyxin_B_PB = POL_ED300, Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
      Cefotaxime = CTX_ED5, Cefotaxime_ClavulanicAcid = CTC_ED30
    )
  
  # 2. Pivot the drug columns into long format
  AMR_Data <- Data %>%
    pivot_longer(
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    )
  
  # 3. Join with the cutoff table
  AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
  
  # Define the list of grouping variables for the final summary
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
  
  # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
  
  Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
    
    # PRECISE FILTER IMPLEMENTATION:
    filter(
      # Condition 1: Isolate is E.coli OR K.pneumoniae
      Isolate %in% c("E.coli", "K.pneumoniae"),
      # Condition 2: VITEK_MS_Results MUST be Escherichia coli
      VITEK_MS_Results == "Escherichia coli",
      # Condition 3: ESBL must be 1
      ESBL == 1
    ) %>%
    
    # Clean the S and R cutoff columns and convert all key values to numeric
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
      Measured_Zone = as.numeric(str_trim(value)),
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    
    # Implement the breakpoint logic for zone diameters (in mm)
    mutate(
      Categorical_Result = case_when(
        Measured_Zone >= S_Cutoff ~ "S",
        Measured_Zone <= R_Cutoff ~ "R",
        Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    
    # Ensure only one result is counted per isolate for each antimicrobial test
    # This maintains data integrity for the prevalence calculation denominator.
    distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- Step 5: Define the Summary Function (CI Calculation for Resistance 'R') ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    # 5a. Total count (denominator)
    Total_Count <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')
    
    # 5b. Count the S, I, R results
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      
      # 5c. Join the total counts back
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
      
      # Filter to Resistance only for CI calculation
      filter(Categorical_Result == "R") %>%
      
      # 5d. Calculate Percentage and Confidence Interval
      rowwise() %>%
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var,
        # Calculate 95% Confidence Interval for Resistance percentage using prop.test
        ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
        CI_Lower = round(ci_result[1] * 100, 1),
        CI_Upper = round(ci_result[2] * 100, 1)
      ) %>%
      ungroup() %>%
      select(-ci_result) %>% # Remove the temporary list column
      
      # 5e. Select and arrange columns for final output
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Resistance_Frequency = N,
        Total_Tests,
        Resistance_Percentage = Percentage,
        CI_95_Lower = CI_Lower,
        CI_95_Upper = CI_Upper
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 6: Apply the function to all grouping variables and combine results ---
  ECO_ESBL_Susceptibility_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  return(ECO_ESBL_Susceptibility_Summary)
}

# --- Function Execution ---
ECO_AMR_Pattern <- calculate_ecoli_esbl_susceptibility(joined_data, EPI_CUTOFF)

write_tsv(ECO_AMR_Pattern, "Results/131_Obs_ECO_AMR_Pattern.tsv")
###############################################################################
# Statistical analyses

# Load required libraries
# NOTE: The 'tidyverse' package is assumed to be loaded for 'dplyr', 'tidyr', and 'rlang' functions.
# library(tidyverse)
# library(dplyr)
 library(tidyr)
 library(stringr)

calculate_ecoli_esbl_susceptibility_modified <- function(joined_data, EPI_CUTOFF) {
  
  # ... (Steps 1, 2, 3, and 4 are unchanged) ...
  
  # 1. Select relevant columns and rename drug columns for clarity
  Data <- joined_data %>%
    select(
      INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
      Amoxicillin = AMX_ED10, Azithromycin = AZM_ED15, Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5, Doxycycline = DOX_ED30, Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10, Meropenem = MEM_ED10, Oxytetracycline = OXY_ED30,
      Polymyxin_B_PB = POL_ED300, Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
      Cefotaxime = CTX_ED5, Cefotaxime_ClavulanicAcid = CTC_ED30
    )
  
  # 2. Pivot the drug columns into long format
  AMR_Data <- Data %>%
    pivot_longer(
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, DISTRICT, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    )
  
  # 3. Join with the cutoff table
  AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
  
  # Define the list of grouping variables for the final summary
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE", "DISTRICT")
  
  # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
  
  Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
    
    # PRECISE FILTER IMPLEMENTATION:
    filter(
      # Condition 1: Isolate is E.coli OR K.pneumoniae
      Isolate %in% c("E.coli", "K.pneumoniae"),
      # Condition 2: VITEK_MS_Results MUST be Escherichia coli
      VITEK_MS_Results == "Escherichia coli",
      # Condition 3: ESBL must be 1
      ESBL == 1
    ) %>%
    
    # Clean the S and R cutoff columns and convert all key values to numeric
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
      Measured_Zone = as.numeric(str_trim(value)),
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    
    # Implement the breakpoint logic for zone diameters (in mm)
    mutate(
      Categorical_Result = case_when(
        Measured_Zone >= S_Cutoff ~ "S",
        Measured_Zone <= R_Cutoff ~ "R",
        Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    
    # Ensure only one result is counted per isolate for each antimicrobial test
    distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- NEW STEP: Define the Chi-Squared Comparison Function (for heterogeneity) ---
  
  calculate_chi_square_comparison <- function(data, group_var) {
    
    # 1. Group by drug and the grouping variable
    Chi_Square_Results <- data %>%
      group_by(Antimicrobial_substance) %>%
      
      # 2. Perform the test for each drug
      summarise(
        Chi_Squared_P_value = {
          
          # Create a contingency table (Grouping Variable vs. Resistance (R vs Not R))
          contingency_table <- table(
            Group_Var = !!sym(group_var),
            Is_Resistant = factor(
              ifelse(Categorical_Result == "R", "R", "Not_R"),
              levels = c("R", "Not_R")
            )
          )
          
          # Check if table is valid (at least 2x2 with non-zero marginals)
          if(min(dim(contingency_table)) > 1 && all(colSums(contingency_table) > 0) && all(rowSums(contingency_table) > 0)) {
            
            # --- FIX APPLIED HERE: Suppress the warning when calculating expected counts ---
            expected_counts <- suppressWarnings(chisq.test(contingency_table))$expected
            
            # Use Fisher's Exact test if any expected cell count is < 5, otherwise use Chi-squared
            chi_test <- tryCatch({
              if(any(expected_counts < 5)) {
                fisher.test(contingency_table)$p.value
              } else {
                chisq.test(contingency_table)$p.value
              }
            }, error = function(e) NA_real_) 
            
          } else {
            NA_real_ # Not enough data or no variance to perform a comparison
          }
          
        }, .groups = 'drop'
      ) %>%
      
      # Add the grouping variable name for clarity
      mutate(Grouping_Variable = group_var)
    
    return(Chi_Square_Results)
  }
  
  # --- Step 5: Define the Summary Function (CI Calculation for Resistance 'R') ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    # 5a. Total count (denominator)
    Total_Count <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')
    
    # 5b. Count the S, I, R results
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      
      # 5c. Join the total counts back
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
      
      # Filter to Resistance only for CI calculation
      filter(Categorical_Result == "R") %>%
      
      # 5d. Calculate Percentage and Confidence Interval
      rowwise() %>%
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var,
        # Calculate 95% Confidence Interval for Resistance percentage using prop.test
        ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
        CI_Lower = round(ci_result[1] * 100, 1),
        CI_Upper = round(ci_result[2] * 100, 1) 
      ) %>%
      ungroup() %>%
      select(-ci_result) %>% # Remove the temporary list column
      
      # 5e. Select and arrange columns for final output
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Resistance_Frequency = N,
        Total_Tests,
        Resistance_Percentage = Percentage,
        CI_95_Lower = CI_Lower,
        CI_95_Upper = CI_Upper 
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 6: Apply the functions and combine results ---
  
  # 6a. Calculate the summary of resistance percentages
  ECO_ESBL_Susceptibility_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  # 6b. Calculate the comparison p-values
  Comparison_P_values <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_chi_square_comparison(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  # 6c. Join the p-values back to the summary table
  Final_Summary <- left_join(
    ECO_ESBL_Susceptibility_Summary,
    Comparison_P_values,
    by = c("Grouping_Variable", "Antimicrobial_substance")
  ) %>%
    # Re-order columns for a clearer output
    select(
      Grouping_Variable, Grouping_Value, Antimicrobial_substance,
      Resistance_Frequency, Total_Tests, Resistance_Percentage,
      CI_95_Lower, CI_95_Upper, Chi_Squared_P_value
    )
  
  return(Final_Summary)
}

# --- Function Execution ---
ECO_AMR_Pattern_pValue <- calculate_ecoli_esbl_susceptibility_modified(joined_data, EPI_CUTOFF)

# Save the file
write_tsv(ECO_AMR_Pattern_pValue, "Results/p_Values_ECO_AMR_Pattern_pValue.tsv")
################################################################################
# Aggregated results for distribution zone of inhibition for ECO 

calculate_ecoli_esbl_susceptibility_modified <- function(joined_data, EPI_CUTOFF) {
  
  # 1. Select relevant columns and rename drug columns for clarity
  Data <- joined_data %>%
    select(
      INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
      Amoxicillin = AMX_ED10, Azithromycin = AZM_ED15, Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5, Doxycycline = DOX_ED30, Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10, Meropenem = MEM_ED10, Oxytetracycline = OXY_ED30,
      Polymyxin_B_PB = POL_ED300, Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
      Cefotaxime = CTX_ED5, Cefotaxime_ClavulanicAcid = CTC_ED30
    )
  
  # 2. Pivot the drug columns into long format
  AMR_Data <- Data %>%
    pivot_longer(
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    )
  
  # 3. Join with the cutoff table
  AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
  
  # Define the list of grouping variables for the group-wise summary
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
  
  # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
  
  Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
    
    # PRECISE FILTER IMPLEMENTATION:
    filter(
      Isolate %in% c("E.coli", "K.pneumoniae"),
      VITEK_MS_Results == "Escherichia coli",
      ESBL == 1
    ) %>%
    
    # Clean the S and R cutoff columns and convert all key values to numeric
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
      Measured_Zone = as.numeric(str_trim(value)),
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    
    # Implement the breakpoint logic for zone diameters (in mm)
    mutate(
      Categorical_Result = case_when(
        Measured_Zone >= S_Cutoff ~ "S",
        Measured_Zone <= R_Cutoff ~ "R",
        Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    
    # Ensure only one result is counted per isolate for each antimicrobial test
    distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- NEW STEP: Define the Chi-Squared Comparison Function (for heterogeneity) ---
  
  calculate_chi_square_comparison <- function(data, group_var) {
    
    Chi_Square_Results <- data %>%
      group_by(Antimicrobial_substance) %>%
      
      summarise(
        Chi_Squared_P_value = {
          contingency_table <- table(
            Group_Var = !!sym(group_var),
            Is_Resistant = factor(
              ifelse(Categorical_Result == "R", "R", "Not_R"),
              levels = c("R", "Not_R")
            )
          )
          
          if(min(dim(contingency_table)) > 1 && all(colSums(contingency_table) > 0) && all(rowSums(contingency_table) > 0)) {
            
            # Suppress the Chi-squared warning when calculating expected counts for the check
            expected_counts <- suppressWarnings(chisq.test(contingency_table))$expected
            
            chi_test <- tryCatch({
              if(any(expected_counts < 5)) {
                fisher.test(contingency_table)$p.value
              } else {
                chisq.test(contingency_table)$p.value
              }
            }, error = function(e) NA_real_)
            
          } else {
            NA_real_
          }
          
        }, .groups = 'drop'
      ) %>%
      
      mutate(Grouping_Variable = group_var)
    
    return(Chi_Square_Results)
  }
  
  # --- Step 5a: Summary Function for Group-Wise Results ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    Total_Count <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')
    
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
      
      filter(Categorical_Result == "R") %>%
      
      rowwise() %>%
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var,
        ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
        CI_Lower = round(ci_result[1] * 100, 1),
        CI_Upper = round(ci_result[2] * 100, 1)
      ) %>%
      ungroup() %>%
      select(-ci_result) %>%
      
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Resistance_Frequency = N,
        Total_Tests,
        Resistance_Percentage = Percentage,
        CI_95_Lower = CI_Lower,
        CI_95_Upper = CI_Upper
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 5b: NEW Summary Function for OVERALL Results ---
  
  calculate_overall_susceptibility_summary <- function(data) {
    
    Total_Count <- data %>%
      group_by(Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')
    
    Summary_Table <- data %>%
      group_by(Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      left_join(Total_Count, by = "Antimicrobial_substance") %>%
      filter(Categorical_Result == "R") %>%
      
      rowwise() %>%
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        # Use "OVERALL" as both the variable name and the value
        Grouping_Variable = "OVERALL", 
        ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
        CI_Lower = round(ci_result[1] * 100, 1),
        CI_Upper = round(ci_result[2] * 100, 1)
      ) %>%
      ungroup() %>%
      select(-ci_result) %>%
      
      select(
        Grouping_Variable,
        Grouping_Value = Grouping_Variable, # Assign "OVERALL" as the Grouping_Value
        Antimicrobial_substance,
        Resistance_Frequency = N,
        Total_Tests,
        Resistance_Percentage = Percentage,
        CI_95_Lower = CI_Lower,
        CI_95_Upper = CI_Upper
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 6: Apply the functions and combine results ---
  
  # 6a-i. Calculate the OVERALL summary first
  Overall_Summary <- calculate_overall_susceptibility_summary(Ecoli_ESBL_Susceptibility)
  
  # 6a-ii. Calculate the group-wise summary
  Group_Wise_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  # 6a-iii. Combine Overall and Group-Wise Summaries
  ECO_ESBL_Susceptibility_Summary <- bind_rows(Overall_Summary, Group_Wise_Summary)
  
  # 6b. Calculate the comparison p-values (ONLY for group-wise data)
  Comparison_P_values <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_chi_square_comparison(Ecoli_ESBL_Susceptibility, g_var)
    })
  )
  
  # 6c. Join the p-values back to the summary table
  Final_Summary <- left_join(
    ECO_ESBL_Susceptibility_Summary,
    # P-values are NA for the "OVERALL" rows, which is correct since no comparison is made
    Comparison_P_values,
    by = c("Grouping_Variable", "Antimicrobial_substance")
  ) %>%
    # Re-order columns for a clearer output
    select(
      Grouping_Variable, Grouping_Value, Antimicrobial_substance,
      Resistance_Frequency, Total_Tests, Resistance_Percentage,
      CI_95_Lower, CI_95_Upper, Chi_Squared_P_value
    )
  
  return(Final_Summary)
}

# --- Function Execution ---
 ECO_AMR_Pattern_Overall_Pvalue <- calculate_ecoli_esbl_susceptibility_modified(joined_data, EPI_CUTOFF)

 # Save the file
 write_tsv(ECO_AMR_Pattern_Overall_Pvalue,"Results/Overall_Pvalue_ECO_AMR_Pattern.tsv")
################################################################################

 ## Aggregated results for distribution zone of inhibition for KPN
 
 # Importing the K.pneumoniae_ECOFF_EUCAST break point file
 
 ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx", sheet = 2)
 
 
 EPI_CUTOFF<-ECOFF_EUCAST_BREAK_POINT
 
 calculate_Kpn_esbl_susceptibility_modified <- function(joined_data, EPI_CUTOFF) {
   
   # 1. Select relevant columns and rename drug columns for clarity
   Data <- joined_data %>%
     select(
       INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
       Amoxicillin = AMX_ED10, Azithromycin = AZM_ED15, Ceftriaxone = CRO_ED30,
       Ciprofloxacin = CIP_ED5, Doxycycline = DOX_ED30, Florfenicol = FLR_ED30,
       Gentamicin = GEN_ED10, Meropenem = MEM_ED10, Oxytetracycline = OXY_ED30,
       Polymyxin_B_PB = POL_ED300, Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
       Cefotaxime = CTX_ED5, Cefotaxime_ClavulanicAcid = CTC_ED30
     )
   
   # 2. Pivot the drug columns into long format
   AMR_Data <- Data %>%
     pivot_longer(
       cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                 ORIGIN_OF_SAMPLE, ESBL),
       names_to = "Antimicrobial_substance",
       values_to = "value"
     )
   
   # 3. Join with the cutoff table
   AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
   
   # Define the list of grouping variables for the group-wise summary
   group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
   
   # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
   
   Kpn_ESBL_Susceptibility <- AMR_Cutoff %>%
     
     # PRECISE FILTER IMPLEMENTATION:
     filter(
       Isolate %in% c("E.coli", "K.pneumoniae"),
       VITEK_MS_Results == "Klebsiella pneumoniae",
       ESBL == 1
     ) %>%
     
     # Clean the S and R cutoff columns and convert all key values to numeric
     mutate(
       S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
       R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
       Measured_Zone = as.numeric(str_trim(value)),
       S_Cutoff = as.numeric(S_Cleaned),
       R_Cutoff = as.numeric(R_Cleaned)
     ) %>%
     
     # Implement the breakpoint logic for zone diameters (in mm)
     mutate(
       Categorical_Result = case_when(
         Measured_Zone >= S_Cutoff ~ "S",
         Measured_Zone <= R_Cutoff ~ "R",
         Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
         TRUE ~ NA_character_
       )
     ) %>%
     
     # Ensure only one result is counted per isolate for each antimicrobial test
     distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
     
     # Remove tests where the result could not be determined
     filter(!is.na(Categorical_Result))
   
   
   # --- NEW STEP: Define the Chi-Squared Comparison Function (for heterogeneity) ---
   
   calculate_chi_square_comparison <- function(data, group_var) {
     
     Chi_Square_Results <- data %>%
       group_by(Antimicrobial_substance) %>%
       
       summarise(
         Chi_Squared_P_value = {
           contingency_table <- table(
             Group_Var = !!sym(group_var),
             Is_Resistant = factor(
               ifelse(Categorical_Result == "R", "R", "Not_R"),
               levels = c("R", "Not_R")
             )
           )
           
           if(min(dim(contingency_table)) > 1 && all(colSums(contingency_table) > 0) && all(rowSums(contingency_table) > 0)) {
             
             # Suppress the Chi-squared warning when calculating expected counts for the check
             expected_counts <- suppressWarnings(chisq.test(contingency_table))$expected
             
             chi_test <- tryCatch({
               if(any(expected_counts < 5)) {
                 fisher.test(contingency_table)$p.value
               } else {
                 chisq.test(contingency_table)$p.value
               }
             }, error = function(e) NA_real_)
             
           } else {
             NA_real_
           }
           
         }, .groups = 'drop'
       ) %>%
       
       mutate(Grouping_Variable = group_var)
     
     return(Chi_Square_Results)
   }
   
   # --- Step 5a: Summary Function for Group-Wise Results ---
   
   calculate_susceptibility_summary <- function(data, group_var) {
     
     Total_Count <- data %>%
       group_by(!!sym(group_var), Antimicrobial_substance) %>%
       summarise(Total_Tests = n(), .groups = 'drop')
     
     Summary_Table <- data %>%
       group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
       summarise(N = n(), .groups = 'drop_last') %>%
       
       left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
       
       filter(Categorical_Result == "R") %>%
       
       rowwise() %>%
       mutate(
         Percentage = round((N / Total_Tests) * 100, 1),
         Grouping_Variable = group_var,
         ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
         CI_Lower = round(ci_result[1] * 100, 1),
         CI_Upper = round(ci_result[2] * 100, 1)
       ) %>%
       ungroup() %>%
       select(-ci_result) %>%
       
       select(
         Grouping_Variable,
         Grouping_Value = !!sym(group_var),
         Antimicrobial_substance,
         Resistance_Frequency = N,
         Total_Tests,
         Resistance_Percentage = Percentage,
         CI_95_Lower = CI_Lower,
         CI_95_Upper = CI_Upper
       )
     
     return(as.data.frame(Summary_Table))
   }
   
   # --- Step 5b: NEW Summary Function for OVERALL Results ---
   
   calculate_overall_susceptibility_summary <- function(data) {
     
     Total_Count <- data %>%
       group_by(Antimicrobial_substance) %>%
       summarise(Total_Tests = n(), .groups = 'drop')
     
     Summary_Table <- data %>%
       group_by(Antimicrobial_substance, Categorical_Result) %>%
       summarise(N = n(), .groups = 'drop_last') %>%
       left_join(Total_Count, by = "Antimicrobial_substance") %>%
       filter(Categorical_Result == "R") %>%
       
       rowwise() %>%
       mutate(
         Percentage = round((N / Total_Tests) * 100, 1),
         # Use "OVERALL" as both the variable name and the value
         Grouping_Variable = "OVERALL", 
         ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
         CI_Lower = round(ci_result[1] * 100, 1),
         CI_Upper = round(ci_result[2] * 100, 1)
       ) %>%
       ungroup() %>%
       select(-ci_result) %>%
       
       select(
         Grouping_Variable,
         Grouping_Value = Grouping_Variable, # Assign "OVERALL" as the Grouping_Value
         Antimicrobial_substance,
         Resistance_Frequency = N,
         Total_Tests,
         Resistance_Percentage = Percentage,
         CI_95_Lower = CI_Lower,
         CI_95_Upper = CI_Upper
       )
     
     return(as.data.frame(Summary_Table))
   }
   
   # --- Step 6: Apply the functions and combine results ---
   
   # 6a-i. Calculate the OVERALL summary first
   Overall_Summary <- calculate_overall_susceptibility_summary(Kpn_ESBL_Susceptibility)
   
   # 6a-ii. Calculate the group-wise summary
   Group_Wise_Summary <- bind_rows(
     lapply(group_vars_list, function(g_var) {
       calculate_susceptibility_summary(Kpn_ESBL_Susceptibility, g_var)
     })
   )
   
   # 6a-iii. Combine Overall and Group-Wise Summaries
   Kpn_ESBL_Susceptibility_Summary <- bind_rows(Overall_Summary, Group_Wise_Summary)
   
   # 6b. Calculate the comparison p-values (ONLY for group-wise data)
   Comparison_P_values <- bind_rows(
     lapply(group_vars_list, function(g_var) {
       calculate_chi_square_comparison(Kpn_ESBL_Susceptibility, g_var)
     })
   )
   
   # 6c. Join the p-values back to the summary table
   Final_Summary <- left_join(
     Kpn_ESBL_Susceptibility_Summary,
     # P-values are NA for the "OVERALL" rows, which is correct since no comparison is made
     Comparison_P_values,
     by = c("Grouping_Variable", "Antimicrobial_substance")
   ) %>%
     # Re-order columns for a clearer output
     select(
       Grouping_Variable, Grouping_Value, Antimicrobial_substance,
       Resistance_Frequency, Total_Tests, Resistance_Percentage,
       CI_95_Lower, CI_95_Upper, Chi_Squared_P_value
     )
   
   return(Final_Summary)
 }
 
 # --- Function Execution ---
 Kpn_AMR_Pattern_Overall_Pvalue <- calculate_Kpn_esbl_susceptibility_modified(joined_data, EPI_CUTOFF)
 
 # Save the file
 write_tsv(Kpn_AMR_Pattern_Overall_Pvalue,"Results/Overall_Pvalue_Kpn_AMR_Pattern.tsv")
 
###############################################################################
 
# AMR Frequencies for the confirmed ESC Resistance Klebsiella pneumoniae
 library(DescTools)
  ## Importing the K.pneumoniae_ECOFF_EUCAST break point file
  
  ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx", sheet = 2)
  
  
  EPI_CUTOFF<-ECOFF_EUCAST_BREAK_POINT
  
  # 
  calculate_Kpn_esbl_susceptibility <- function(joined_data, EPI_CUTOFF) {
    
    # 1. Select relevant columns and rename drug columns for clarity
    # NOTE: The 'distinct(INIKA_ID, .keep_all = TRUE)' line is intentionally omitted here
    Data <- joined_data %>%
      select(
        INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate, VITEK_MS_Results, ESBL,
        Amoxicillin = AMX_ED10, Azithromycin = AZM_ED15, Ceftriaxone = CRO_ED30,
        Ciprofloxacin = CIP_ED5, Doxycycline = DOX_ED30, Florfenicol = FLR_ED30,
        Gentamicin = GEN_ED10, Meropenem = MEM_ED10, Oxytetracycline = OXY_ED30,
        Polymyxin_B_PB = POL_ED300, Sulfamethoxazole_Trimethoprim = SXT_ED1_2,
        Cefotaxime = CTX_ED5, Cefotaxime_ClavulanicAcid = CTC_ED30
      )
    
    # 2. Pivot the drug columns into long format
    AMR_Data <- Data %>%
      pivot_longer(
        cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                  ORIGIN_OF_SAMPLE, ESBL),
        names_to = "Antimicrobial_substance",
        values_to = "value"
      )
    
    # 3. Join with the cutoff table
    AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
    
    # Define the list of grouping variables for the final summary
    group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
    
    # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
    
    Kpn_ESBL_Susceptibility <- AMR_Cutoff %>%
      
      # PRECISE FILTER IMPLEMENTATION:
      filter(
        # Condition 1: Isolate is E.coli OR K.pneumoniae
        Isolate %in% c("E.coli", "K.pneumoniae"),
        # Condition 2: VITEK_MS_Results MUST be Klebsiella pneumoniae
        VITEK_MS_Results == "Klebsiella pneumoniae",
        # Condition 3: ESBL must be 1
        ESBL == 1
      ) %>%
      
      # Clean the S and R cutoff columns and convert all key values to numeric
      mutate(
        S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
        R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
        Measured_Zone = as.numeric(str_trim(value)),
        S_Cutoff = as.numeric(S_Cleaned),
        R_Cutoff = as.numeric(R_Cleaned)
      ) %>%
      
      # Implement the breakpoint logic for zone diameters (in mm)
      mutate(
        Categorical_Result = case_when(
          Measured_Zone >= S_Cutoff ~ "S",
          Measured_Zone <= R_Cutoff ~ "R",
          Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
          TRUE ~ NA_character_
        )
      ) %>%
      
      # Ensure only one result is counted per isolate for each antimicrobial test
      # This maintains data integrity for the prevalence calculation denominator.
      distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
      
      # Remove tests where the result could not be determined
      filter(!is.na(Categorical_Result))
    
    
    # --- Step 5: Define the Summary Function (CI Calculation for Resistance 'R') ---
    
    calculate_susceptibility_summary <- function(data, group_var) {
      
      # 5a. Total count (denominator)
      Total_Count <- data %>%
        group_by(!!sym(group_var), Antimicrobial_substance) %>%
        summarise(Total_Tests = n(), .groups = 'drop')
      
      # 5b. Count the S, I, R results
      Summary_Table <- data %>%
        group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
        summarise(N = n(), .groups = 'drop_last') %>%
        
        # 5c. Join the total counts back
        left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
        
        # Filter to Resistance only for CI calculation
        filter(Categorical_Result == "R") %>%
        
        # 5d. Calculate Percentage and Confidence Interval
        rowwise() %>%
        mutate(
          Percentage = round((N / Total_Tests) * 100, 1),
          Grouping_Variable = group_var,
          # Calculate 95% Confidence Interval for Resistance percentage using prop.test
          ci_result = if(Total_Tests > 0) list(prop.test(N, Total_Tests, conf.level = 0.95)$conf.int) else list(c(NA, NA)),
          CI_Lower = round(ci_result[1] * 100, 1),
          CI_Upper = round(ci_result[2] * 100, 1)
        ) %>%
        ungroup() %>%
        select(-ci_result) %>% # Remove the temporary list column
        
        # 5e. Select and arrange columns for final output
        select(
          Grouping_Variable,
          Grouping_Value = !!sym(group_var),
          Antimicrobial_substance,
          Resistance_Frequency = N,
          Total_Tests,
          Resistance_Percentage = Percentage,
          CI_95_Lower = CI_Lower,
          CI_95_Upper = CI_Upper
        )
      
      return(as.data.frame(Summary_Table))
    }
    
    # --- Step 6: Apply the function to all grouping variables and combine results ---
    Kpn_ESBL_Susceptibility_Summary <- bind_rows(
      lapply(group_vars_list, function(g_var) {
        calculate_susceptibility_summary(Kpn_ESBL_Susceptibility, g_var)
      })
    )
    
    return(Kpn_ESBL_Susceptibility_Summary)
  }
  
  # --- Function Execution ---
  Kpn_AMR_Pattern <- calculate_Kpn_esbl_susceptibility(joined_data, EPI_CUTOFF)
  # Save file
  write_tsv(Kpn_AMR_Pattern, "Results/54_Obs_Kpn_AMR_Pattern.tsv")
  
################################################################################
################################################################################  
 
  ## Distribution of Zone diameter in millimeters (mm) 
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map
  # The original column names (the zone diameter values)
  zone_diameter_cols <- c(
    "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
    "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
    "SXT_ED1_2", "CTX_ED5", "CTC_ED30"
  )
  
  rename_map <- c(
    "AMX_ED10" = "Amoxicillin", 
    "AZM_ED15" = "Azithromycin", 
    "CRO_ED30" = "Ceftriaxone", 
    "CIP_ED5" = "Ciprofloxacin", 
    "DOX_ED30" = "Doxycycline", 
    "FLR_ED30" = "Florfenicol",
    "GEN_ED10" = "Gentamicin", 
    "MEM_ED10" = "Meropenem", 
    "OXY_ED30" = "Oxytetracycline", 
    "POL_ED300" = "Popymyxin_B(PB)",
    "SXT_ED1_2" = "Sulfamethoxazole/Trimethoprim",
    "CTX_ED5" = "Cefotaxime",
    "CTC_ED30" = "Cefotaxime/clavulanic acid"
  )
  
  # 2. Calculate the Distribution Statistics
  Zone_Distribution_Table <- joined_data %>% 
    # Select the unique ID column along with the zone diameter columns
    select(INIKA_ID, all_of(zone_diameter_cols)) %>% ## 420 Total tested!!
    
    # Convert data from wide to long format
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance", 
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # Apply descriptive renaming
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance, !!!rename_map)
    ) %>%
    
    # Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
    # Calculate statistics, 95% CI, and round to 1 decimal place
    summarise(
      # N_Tested is the sample size (N) since all IDs are unique here
      N_Tested = n(), 
      
      # Calculate Mean and SD directly rounded
      Mean_Zone_mm = round(mean(Zone_Diameter_mm), 1),
      SD_Zone_mm = round(sd(Zone_Diameter_mm), 1),
      
      # Calculate 95% CI Lower Bound (Lower CI = Mean - MOE)
      CI_95_Lower = round(
        mean(Zone_Diameter_mm) - 
          qt(0.975, df = n() - 1) * (sd(Zone_Diameter_mm) / sqrt(n())), 
        1
      ),
      
      # Calculate 95% CI Upper Bound (Upper CI = Mean + MOE)
      CI_95_Upper = round(
        mean(Zone_Diameter_mm) + 
          qt(0.975, df = n() - 1) * (sd(Zone_Diameter_mm) / sqrt(n())), 
        1
      ),
      
      # Other descriptive statistics rounded
      Median_Zone_mm = round(median(Zone_Diameter_mm), 1),
      Min_Zone_mm = round(min(Zone_Diameter_mm), 1),
      Max_Zone_mm = round(max(Zone_Diameter_mm), 1),
      .groups = 'drop'
    ) %>%
    
    # Arrange by the number tested or mean zone (optional)
    arrange(desc(N_Tested), Mean_Zone_mm)
  # Print the resulting table
  print("--- Distribution of Zone Diameters (mm) per Antimicrobial Substance ---")
  print(Zone_Distribution_Table)
  
  # Save the file
  
  write_tsv(Zone_Distribution_Table, "Results/420_Obs_95%CI_Zone_Distribution_Table-10-11-2025.tsv")
  
  #####################################################################
  ## Here is for those UNIQUE with no Duplicates
  
  # 2. Calculate the Distribution Statistics
  Zone_Distribution_Table <- MMJ %>%  ## 
    # Select the unique ID column along with the zone diameter columns
    select(INIKA_ID, all_of(zone_diameter_cols)) %>% ## 399 Total tested!!
    
    # Convert data from wide to long format
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance", 
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # Apply descriptive renaming
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance, !!!rename_map)
    ) %>%
    
    # Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
  # Calculate distribution statistics including 95% CI
  summarise(
    N_Tested = n_distinct(INIKA_ID), 
    
    # Calculate Mean and SD first
    Mean_Zone_mm_raw = mean(Zone_Diameter_mm),
    SD_Zone_mm_raw = sd(Zone_Diameter_mm),
    N = n(), # Total number of observations (including duplicates from same ID if any)
    
    # Calculate Standard Error (SE)
    SE = SD_Zone_mm_raw / sqrt(N),
    
    # Calculate 95% Margin of Error (MOE) using the t-distribution
    # Use qt(0.975, df = N - 1) for the t-score (critical value)
    MOE = qt(0.975, df = N - 1) * SE,
    
    # --- Final Rounded Results (1 Decimal Place) ---
    Mean_Zone_mm = round(Mean_Zone_mm_raw, 1),
    SD_Zone_mm = round(SD_Zone_mm_raw, 1),
    
    # 95% CI Lower Bound
    CI_95_Lower = round(Mean_Zone_mm_raw - MOE, 1),
    
    # 95% CI Upper Bound
    CI_95_Upper = round(Mean_Zone_mm_raw + MOE, 1),
    
    Median_Zone_mm = round(median(Zone_Diameter_mm), 1),
    Min_Zone_mm = round(min(Zone_Diameter_mm), 1),
    Max_Zone_mm = round(max(Zone_Diameter_mm), 1),
    .groups = 'drop'
  ) %>%
    arrange(desc(N_Tested), Mean_Zone_mm)
  # Print the resulting table
  print("--- Distribution of Zone Diameters (mm) per Antimicrobial Substance ---")
  print(Zone_Distribution_Table)
  
  # Save the file
  
  write_tsv(Zone_Distribution_Table, "Results/399_Obs_95%CI_Zone_Distribution_Table-10-11-2025.tsv")
######################################################################## 
  ## Here is the computation for all E.coli and K.pneumoniae from Isolate column 
  #are Escherichia coli from the VITEK_MS_Results column and are ESBL
  
  Zone_Distribution_Table_Filtered <- joined_data %>%
    # 1. Select all necessary columns, including those for filtering
    select(INIKA_ID, VITEK_MS_Results, ESBL, Isolate, all_of(zone_diameter_cols)) %>% 
    
    # 2. FILTER the rows based on the Isolate Type and ESBL status
    filter(
      # Condition: The isolate must be ESBL-positive (1) AND meet the species/VITEK criteria
      ESBL == 1 & (
        # Case 1: Isolate is E.coli AND VITEK confirms Escherichia coli
        (Isolate == "E.coli" & VITEK_MS_Results == "Escherichia coli") |
          
          # Case 2: Isolate is K.pneumoniae (no VITEK check explicitly requested for K.pneumoniae)
          (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Escherichia coli")
      )
    ) %>% ## 131 Observations
    
    # 3. Convert data from wide to long format
    # Filtering columns (VITEK_MS_Results, ESBL, Isolate) are automatically dropped from the pivot
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance", 
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # 4. Apply descriptive renaming
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance, !!!rename_map)
    ) %>%
    
    # 5. Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # 6. Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
    # 7. Calculate distribution statistics, 95% CI, and round to 1 decimal place
    summarise(
      N_Tested = n_distinct(INIKA_ID), 
      
      # Calculate Mean and SD first
      Mean_Zone_mm_raw = mean(Zone_Diameter_mm),
      SD_Zone_mm_raw = sd(Zone_Diameter_mm),
      N = n(), 
      
      # Calculate Standard Error (SE)
      SE = SD_Zone_mm_raw / sqrt(N),
      
      # Calculate 95% Margin of Error (MOE)
      MOE = qt(0.975, df = N - 1) * SE,
      
      # --- Final Rounded Results (1 Decimal Place) ---
      Mean_Zone_mm = round(Mean_Zone_mm_raw, 1),
      SD_Zone_mm = round(SD_Zone_mm_raw, 1),
      
      # 95% CI Lower Bound
      CI_95_Lower = round(Mean_Zone_mm_raw - MOE, 1),
      
      # 95% CI Upper Bound
      CI_95_Upper = round(Mean_Zone_mm_raw + MOE, 1),
      
      Median_Zone_mm = round(median(Zone_Diameter_mm), 1),
      Min_Zone_mm = round(min(Zone_Diameter_mm), 1),
      Max_Zone_mm = round(max(Zone_Diameter_mm), 1),
      .groups = 'drop'
    ) %>%
    
    # 8. Arrange by the number tested or mean zone (optional)
    arrange(desc(N_Tested), Mean_Zone_mm)
  
  # Print the resulting table
  print("--- Distribution of Zone Diameters (mm) per Antimicrobial Substance ---")
  print(Zone_Distribution_Table_Filtered)
  
  # Save the file
  
  write_tsv(Zone_Distribution_Table_Filtered, "Results/ECO_95%CI_Zone_Distribution_Table-10-11-2025.tsv")
  
###########################################################################  
  ## Here is the computation for all E.coli and K.pneumoniae from Isolate column 
  #are Klebsiella pneumoniae from the VITEK_MS_Results column and are ESBL
  
  Zone_Distribution_Table_Filtered <- joined_data %>%
    # 1. Select all necessary columns, including those for filtering
    select(INIKA_ID, VITEK_MS_Results, ESBL, Isolate, all_of(zone_diameter_cols)) %>% 
    
    # 2. FILTER the rows based on the Isolate Type and ESBL status
    filter(
      # Condition: The isolate must be ESBL-positive (1) AND meet the species/VITEK criteria
      ESBL == 1 & (
        # Case 1: Isolate is E.coli AND VITEK confirms Klebsiella pneumoniae
        (Isolate == "E.coli" & VITEK_MS_Results == "Klebsiella pneumoniae") |
          
          # Case 2: Isolate is K.pneumoniae (no VITEK check explicitly requested for K.pneumoniae)
          (Isolate == "K.pneumoniae" & VITEK_MS_Results == "Klebsiella pneumoniae")
      )
    ) %>% ## 131 Observations
    
    # 3. Convert data from wide to long format
    # Filtering columns (VITEK_MS_Results, ESBL, Isolate) are automatically dropped from the pivot
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance", 
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # 4. Apply descriptive renaming
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance, !!!rename_map)
    ) %>%
    
    # 5. Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # 6. Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
    # 7. Calculate distribution statistics, 95% CI, and round to 1 decimal place
    summarise(
      N_Tested = n_distinct(INIKA_ID), 
      
      # Calculate Mean and SD first
      Mean_Zone_mm_raw = mean(Zone_Diameter_mm),
      SD_Zone_mm_raw = sd(Zone_Diameter_mm),
      N = n(), 
      
      # Calculate Standard Error (SE)
      SE = SD_Zone_mm_raw / sqrt(N),
      
      # Calculate 95% Margin of Error (MOE)
      MOE = qt(0.975, df = N - 1) * SE,
      
      # --- Final Rounded Results (1 Decimal Place) ---
      Mean_Zone_mm = round(Mean_Zone_mm_raw, 1),
      SD_Zone_mm = round(SD_Zone_mm_raw, 1),
      
      # 95% CI Lower Bound
      CI_95_Lower = round(Mean_Zone_mm_raw - MOE, 1),
      
      # 95% CI Upper Bound
      CI_95_Upper = round(Mean_Zone_mm_raw + MOE, 1),
      
      Median_Zone_mm = round(median(Zone_Diameter_mm), 1),
      Min_Zone_mm = round(min(Zone_Diameter_mm), 1),
      Max_Zone_mm = round(max(Zone_Diameter_mm), 1),
      .groups = 'drop'
    ) %>%
    
    # 8. Arrange by the number tested or mean zone (optional)
    arrange(desc(N_Tested), Mean_Zone_mm)
  
  # Print the resulting table
  print("--- Distribution of Zone Diameters (mm) per Antimicrobial Substance ---")
  print(Zone_Distribution_Table_Filtered)
  
  # Save the file
  
  write_tsv(Zone_Distribution_Table_Filtered, "Results/KPN_95%CI_Zone_Distribution_Table-10-11-2025.tsv")
############################################################################  
   
  ## Distribution of Zone of inhibition diameter in millimeters (mm) for each Antibiotic tested to E.coli 
   # Load necessary libraries
  
  library(ggplot2) # For plotting the distributions
  library(cowplot) # For plotting the distributions
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map (using your previous definitions)
  # 1. Define the Antimicrobial Column Names and the Renaming Map
  # (No changes needed here)
  zone_diameter_cols <- c(
    "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
    "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
    "SXT_ED1_2", "CTX_ED5", "CTC_ED30"
  )
  
  rename_map <- c(
    "AMX_ED10" = "Amoxicillin", 
    "AZM_ED15" = "Azithromycin", 
    "CRO_ED30" = "Ceftriaxone", 
    "CIP_ED5" = "Ciprofloxacin", 
    "DOX_ED30" = "Doxycycline", 
    "FLR_ED30" = "Florfenicol",
    "GEN_ED10" = "Gentamicin", 
    "MEM_ED10" = "Meropenem", 
    "OXY_ED30" = "Oxytetracycline", 
    "POL_ED300" = "Polymyxin_B(PB)",
    "SXT_ED1_2" = "Sulfamethoxazole/Trimethoprim",
    "CTX_ED5" = "Cefotaxime",
    "CTC_ED30" = "Cefotaxime/ClavulanicAcid"
  )
  
  # 2. Prepare the Filtered and Long-Format Data (E. coli population)
  # Ecoli_Zone_Data_Long <- joined_data %>%  ## 119 observed
  #   filter(
  #     Isolate == "E.coli", 
  #     VITEK_MS_Results == "Escherichia coli",
  #     ESBL == "1"
  #   ) %>%
    # 2. Prepare the Filtered and Long-Format Data (E. coli population)
    # Ecoli_Zone_Data_Long <- joined_data %>%  ## 129 observed
    # filter(
    #   Isolate == "E.coli", 
    #   VITEK_MS_Results == "Escherichia coli"
    # ) %>%
  # 2. Prepare the Filtered and Long-Format Data (E. coli population)
  Ecoli_Zone_Data_Long <- joined_data %>%  ## 135 observed
    filter(
      VITEK_MS_Results == "Escherichia coli",
      ESBL == "1"
    ) %>%
    select(all_of(zone_diameter_cols)) %>%
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance_Code",
      values_to = "Zone_Diameter_mm"
    ) %>%
    filter(!is.na(Zone_Diameter_mm)) %>%
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance_Code, !!!rename_map)
    ) %>%
    select(-Antimicrobial_Substance_Code) 
  
  # 3. Generate Histograms for Each Antibiotic
  antibiotics_list <- unique(Ecoli_Zone_Data_Long$Antimicrobial_Substance)
  plot_list <- list()
  
  for (abx in antibiotics_list) {
    # Filter data for the current antibiotic
    abx_data <- Ecoli_Zone_Data_Long %>% 
      filter(Antimicrobial_Substance == abx)
    
    # Create the histogram
    p <- ggplot(abx_data, aes(x = Zone_Diameter_mm)) +
      geom_histogram(
        binwidth = 1, # Set bin size to 1 mm for clear distribution
        fill = "#4C78A8", 
        color = "white",
        boundary = 0
      ) +
      labs(
        title = paste("Zone Diameter Distribution for", abx),
        subtitle = paste0("E.coli Isolates (N = ", nrow(abx_data), ")"),
        x = "Zone Diameter (mm)",
        y = "Frequency (Count)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
    
    # Store the plot in the list
    plot_list[[abx]] <- p
    
    # OPTIONAL: Uncomment the line below to display each plot immediately
    # print(p)
  }
  
  # 4. Display or Save the Plots
  
  # You can display all plots saved in the list:
  # If you have many, you might want to view them one by one.
  # For example, to view Ciprofloxacin:
  # print(plot_list[["Ciprofloxacin"]])
  
  # Or, if you use the 'cowplot' or 'patchwork' packages, you can combine them:
  # library(cowplot)
  # plot_grid(plotlist = plot_list, ncol = 3)
  
  print(plot_list[["Ciprofloxacin"]])
  print(plot_list[["Amoxicillin"]])
  print(plot_list[["Polymyxin_B(PB)"]])
  
  ## Saving the file
  # After installing the package: 
  install.packages("patchwork")
  library(patchwork)
  combined_plot_freq <- patchwork::wrap_plots(plot_list, ncol = 3)
  
  ggsave(
    filename = paste0(output_dir, "Results/Figures", "Combined_Frequency_Distribution.tiff"),
      plot = combined_plot_freq,
      device = "tiff",
      dpi = 600,
      width = 12,
      height = 10,
      units = "in"
  )

  ###########################################################################
# Percentage of distribution
  # Load necessary library (assuming you haven't already)
  library(tidyverse) 
  library(ggplot2) 
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map
  zone_diameter_cols <- c(
    "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
    "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
    "SXT_ED1_2", "CTX_ED5", "CTC_ED30"
  )
  
  rename_map <- c(
    "AMX_ED10" = "Amoxicillin", 
    "AZM_ED15" = "Azithromycin", 
    "CRO_ED30" = "Ceftriaxone", 
    "CIP_ED5" = "Ciprofloxacin", 
    "DOX_ED30" = "Doxycycline", 
    "FLR_ED30" = "Florfenicol",
    "GEN_ED10" = "Gentamicin", 
    "MEM_ED10" = "Meropenem", 
    "OXY_ED30" = "Oxytetracycline", 
    "POL_ED300" = "Polymyxin_B(PB)",
    "SXT_ED1_2" = "Sulfamethoxazole/Trimethoprim",
    "CTX_ED5" = "Cefotaxime",
    "CTC_ED30" = "Cefotaxime/ClavulanicAcid"
  )
  
  # 2. Prepare the Filtered and Long-Format Data (E. coli population)
  Ecoli_Zone_Data_Long <- joined_data %>% 
    filter(
      VITEK_MS_Results == "Escherichia coli",
      ESBL == "1"
    ) %>%
    select(all_of(zone_diameter_cols)) %>%
    pivot_longer(
      cols = all_of(zone_diameter_cols), 
      names_to = "Antimicrobial_Substance_Code",
      values_to = "Zone_Diameter_mm"
    ) %>%
    filter(!is.na(Zone_Diameter_mm)) %>%
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance_Code, !!!rename_map)
    ) %>%
    select(-Antimicrobial_Substance_Code) 
  
  # 3. Generate Histograms for Each Antibiotic (Percentage Distribution)
  antibiotics_list <- unique(Ecoli_Zone_Data_Long$Antimicrobial_Substance)
  plot_list <- list()
  
  # 🌟 FIX APPLIED HERE: Define the variable and create the folder 🌟
  output_dir <- "Results/Figures/Percentage_Distribution_SVGs"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (abx in antibiotics_list) {
    # Filter data for the current antibiotic
    abx_data <- Ecoli_Zone_Data_Long %>% 
      filter(Antimicrobial_Substance == abx)
    
    # Calculate the total number of tested isolates for the subtitle
    N_total <- nrow(abx_data) 
    
    # Create the histogram
    p <- ggplot(abx_data, aes(x = Zone_Diameter_mm)) +
      # KEY CHANGE 1: Use 'y = after_stat(count / sum(count) * 100)' to calculate the percentage
      geom_histogram(
        aes(y = after_stat(count / sum(count) * 100)),
        binwidth = 1, 
        fill = "#4C78A8", 
        color = "white",
        boundary = 0
      ) +
      labs(
        title = paste("Zone Diameter Distribution for", abx),
        subtitle = paste0("E.coli Isolates (N = ", N_total, ")"),
        x = "Zone Diameter (mm)",
        # KEY CHANGE 2: Update Y-axis label
        y = "Relative Frequency (%)" 
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold")) +
      # KEY CHANGE 3: Add continuous scale for Y-axis with percentage labels
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
    
    install.packages("systemfonts")
    library(systemfonts)
    # Store the plot in the list
    plot_list[[abx]] <- p
    # SAVE INDIVIDUAL PLOT AS SVG 
    clean_abx_name <- gsub("/", "_", abx) 
    # 'output_dir' is now defined and the path works
    file_name <- file.path(output_dir, paste0(clean_abx_name, "_Zone_Distribution.svg"))
    
    ggsave(
      filename = file_name,
      plot = p,
      device = "svg",
      width = 6,
      height = 4.5,
      units = "in"
    )
    # print(p) # Uncomment to display plots in the loop
  }
  
  1
  # 4. Display or Save the Plots
  print(plot_list[["Ciprofloxacin"]])
  print(plot_list[["Meropenem"]])
  
  # Combine all plots in the list into a single layout (e.g., 3 columns)
  combined_plot <- wrap_plots(plot_list, ncol = 3) 
  
  # 4. Combine All Plots and Save as Single SVG
  combined_plot <- wrap_plots(plot_list, ncol = 3) 
  
  # This ggsave call will now find 'output_dir'
  ggsave(
    filename = file.path(output_dir, "Combined_Percentage_Distribution.svg"),
    plot = combined_plot,
    device = "svg",
    width = 12,
    height = 10,
    units = "in"
  )
################################################################################  
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map
  zone_diameter_cols <- c(
    "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
    "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
    "SXT_ED1_2", "CTX_ED5", "CTC_ED30"
  )
  
  rename_map <- c(
    "AMX_ED10" = "Amoxicillin", 
    "AZM_ED15" = "Azithromycin", 
    "CRO_ED30" = "Ceftriaxone", 
    "CIP_ED5" = "Ciprofloxacin", 
    "DOX_ED30" = "Doxycycline", 
    "FLR_ED30" = "Florfenicol",
    "GEN_ED10" = "Gentamicin", 
    "MEM_ED10" = "Meropenem", 
    "OXY_ED30" = "Oxytetracycline", 
    "POL_ED300" = "Polymyxin B (colistin resistance)",
    "SXT_ED1_2" = "Sulphonamides/Trimethoprim",
    "CTX_ED5" = "Cefotaxime",
    "CTC_ED30" = "Cefotaxime + clavulanic acid"
  )
  
  # 2. Filter, pivot, rename, and calculate summary statistics (distribution)
  Ecoli_Zone_Distribution_Table <- joined_data %>%
    # Filter for the target population: E.coli (Isolate.x) confirmed as Escherichia coli (VITEK_MS_Results)
    filter(
      Isolate.x == "E.coli", 
      VITEK_MS_Results == "Escherichia coli"
    ) %>%
    # Select only the relevant zone diameter columns
    select(all_of(zone_diameter_cols)) %>%
    
    # Convert data from wide to long format
    pivot_longer(
      cols = everything(),
      names_to = "Antimicrobial_Substance_Code",
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # Rename the substance codes to descriptive names
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance_Code, !!!rename_map)
    ) %>%
    select(-Antimicrobial_Substance_Code) %>% # Clean up the code column
    
    # Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
    # Calculate distribution statistics (N, Mean, SD, Median, Min, Max)
    summarise(
      N_Tested = n(),
      Mean_Zone_mm = mean(Zone_Diameter_mm),
      SD_Zone_mm = sd(Zone_Diameter_mm),
      Median_Zone_mm = median(Zone_Diameter_mm),
      Min_Zone_mm = min(Zone_Diameter_mm),
      Max_Zone_mm = max(Zone_Diameter_mm),
      .groups = 'drop'
    ) %>%
    
    # Arrange by the number tested for clarity
    arrange(desc(N_Tested), Antimicrobial_Substance)
  
  # 3. Print the resulting table to your console
  print("--- Zone Diameter Distribution Table for Confirmed E. coli Isolates ---")
  print(Ecoli_Zone_Distribution_Table)
  # Save table
  write_tsv(Ecoli_Zone_Distribution_Table, "Results/Ecoli_Zone_Distribution_Table.tsv")
#######################################################################
  # 1. Ensure the ggplot2 library is loaded
  library(ggplot2)
  
  # 2. Use ggsave() to export your plot with high-quality settings
  
  ggsave(
    filename = "Ecoli_Ciprofloxacin_Distribution.tiff", # Use .tiff or .pdf for publication
    plot = p,                                           # The plot object you want to save
    device = "tiff",                                    # Explicitly set the device type (optional for TIFF/PDF)
    width = 6,                                          # Width of the final image in the specified units
    height = 4,                                         # Height of the final image in the specified units
    units = "in",                                       # Units can be "in", "cm", or "mm"
    dpi = 300                                           # Dots Per Inch (DPI). Use 300-600 for publication quality.
  )
  # Assuming 'plot_list' is a list of your 13 antibiotic histograms
  
  for (name in names(plot_list)) {
    # Create a clean filename from the antibiotic name
    Ecoli_Zone_Distribution <- paste0("Ecoli_Zone_Distribution_", gsub(" ", "_", name), ".tiff")
    
    # Export the plot
    ggsave(
      filename = file_name,
      plot = plot_list[[name]],
      device = "tiff",
      width = 7, 
      height = 5,
      units = "in",
      dpi = 300 
    )
  }
  
  # Save file
  saveRDS(Ecoli_Zone_Distribution, file = "Ecoli_Zone_Distribution.rds")
  # This code will save 13 separate TIFF files, one for each antibiotic. 
 ############################################################################# 
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  
  joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
  spec(joined_data)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  
  # --- 1. Corrected Data Preparation (Ensuring ESBL is included) ---
  
  # CRITICAL FIX: Add the ESBL column to the Data selection
  Data <- joined_data %>%
    select(INIKA_ID.x,REGION.x,SEASON.x,ORIGIN_OF_SAMPLE.x,Isolate.x,
           VITEK_MS_Results, ESBL, # <--- ESBL ADDED HERE
           AMX_ED10,AZM_ED15,CRO_ED30,
           CIP_ED5,DOX_ED30,FLR_ED30,GEN_ED10,MEM_ED10,
           OXY_ED30,POL_ED300,SXT_ED1_2,CTX_ED5,CTC_ED30)%>%
    rename(Amoxicillin="AMX_ED10", Azithromyzin="AZM_ED15",
           Ceftriaxone="CRO_ED30", Ciprofloxacin="CIP_ED5",
           Doxycycline="DOX_ED30" , Florfenicol="FLR_ED30",
           Gentamicin="GEN_ED10", Meropenem="MEM_ED10",
           Oxytetracycline="OXY_ED30",
           "Polymyxin B (PB) "="POL_ED300",
           "Sulphonamides/Trimethoprim"="SXT_ED1_2",
           Cefotaxime ="CTX_ED5",
           "Cefotaxime+clavulanic acid"="CTC_ED30")
  
  AMR_Data<-pivot_longer(Data, cols=-c(INIKA_ID.x,
                                       VITEK_MS_Results,
                                       Isolate.x,REGION.x,SEASON.x,
                                       ORIGIN_OF_SAMPLE.x, ESBL), # <--- ESBL ADDED HERE
                         names_to = "Antimicrobial_substance",
                         values_to="value")
  
  # Assuming EPI_CUTOFF is available (from your previous code block)
  # AMR_Cutoff<-left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
  # 
  # 
  # # --- 2. Corrected Ecoli_Susceptibility Preparation ---
  # 
  # Ecoli_Susceptibility_Stricter <- AMR_Cutoff %>%
  #   
  #   # Filter only for Escherichia coli in VITEK_MS_Results
  #   filter(VITEK_MS_Results == "Escherichia coli") %>%
  #   
  #   # Ensure the measurement, breakpoints, and ESBL are numeric for comparison
  #   mutate(
  #     S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
  #     R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
  #     Measured_Zone = as.numeric(str_trim(value)),
  #     S_Cutoff = as.numeric(S_Cleaned),
  #     R_Cutoff = as.numeric(R_Cleaned),
  #     ESBL_Numeric = as.numeric(ESBL) # Convert ESBL to numeric
  #   ) %>%
  #   
  #   # Implement the CLSI/EUCAST breakpoint logic (R, I, S)
  #   mutate(
  #     Categorical_Result = case_when(
  #       Measured_Zone >= S_Cutoff ~ "S",
  #       Measured_Zone <= R_Cutoff ~ "R",
  #       TRUE ~ NA_character_
  #     )
  #   ) %>%
  #   
  #   # Ensure only one result per isolate-antimicrobial test
  #   distinct(INIKA_ID.x, Antimicrobial_substance, .keep_all = TRUE) %>%
  #   filter(!is.na(Categorical_Result))
  # 
  # 
  # # --- 3. REPLACING THE Summary Function with Strict Logic ---
  # 
  # # Denominator: All unique E. coli from VITEK_MS_Results
  # # Numerator: Unique isolates where Isolate.x=E.coli AND VITEK=E.coli AND ESBL=1 AND Result=R
  # 
  # calculate_stricter_amr_summary <- function(data, group_var) {
  #   
  #   # Denominator: Total count of UNIQUE E. coli isolates identified by VITEK
  #   Total_Ecoli_Isolates <- data %>%
  #     distinct(INIKA_ID.x, !!sym(group_var)) %>% # Unique isolate IDs per group
  #     group_by(!!sym(group_var)) %>%
  #     summarise(Denominator_Total = n(), .groups = 'drop')
  #   
  #   # Numerator: Count of unique isolates meeting the strict criteria AND being resistant (R)
  #   Summary_Table <- data %>%
  #     # Filter for the strict numerator criteria:
  #     filter(
  #       Isolate.x == "Escherichia coli" &
  #         # VITEK_MS_Results == "Escherichia coli" (already filtered in Step 2) &
  #         ESBL_Numeric == 1 &
  #         Categorical_Result == "R" # Must be resistant based on breakpoints
  #     ) %>%
  #     # Count the unique isolates meeting these criteria for each group/substance
  #     group_by(!!sym(group_var), Antimicrobial_substance) %>%
  #     summarise(
  #       Numerator_Resistant = n_distinct(INIKA_ID.x), # Count unique isolates
  #       .groups = 'drop_last'
  #     ) %>%
  #     
  #     # Join the denominator (Total E.coli Isolates) back
  #     left_join(Total_Ecoli_Isolates, by = group_var) %>%
  #     
  #     # Calculate Percentage (Resistant / Denominator_Total)
  #     mutate(
  #       Percentage = round((Numerator_Resistant / Denominator_Total) * 100, 1),
  #       Grouping_Variable = group_var
  #     ) %>%
  #     
  #     # Select and arrange columns for final output
  #     select(
  #       Grouping_Variable,
  #       Grouping_Value = !!sym(group_var),
  #       Antimicrobial_substance,
  #       Frequency = Numerator_Resistant,
  #       Total_Tests = Denominator_Total,
  #       Percentage
  #     )
  #   
  #   return(as.data.frame(Summary_Table))
  # }
  # 
  # # --- 4. Apply the function to all grouping variables and combine results ---
  # group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x")
  # 
  # # Use the NEW function here!
  # ECO_Stricter_AMR_Summary <- bind_rows(
  #   lapply(group_vars_list, function(g_var) {
  #     calculate_stricter_amr_summary(Ecoli_Susceptibility_Stricter, g_var)
  #   })
  # )
  
  # Display the final summary table (the pivot-longer output)
  print(ECO_Stricter_AMR_Summary)