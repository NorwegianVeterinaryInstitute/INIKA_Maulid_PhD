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
library(writexl)  # For saving to Excel


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

# Importing the file 

joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
spec(joined_data)

# Renaming of the columns
#joined_data <- joined_data %>% 
 # rename(INIKA_ID = INIKA_ID.x,
         # REGION = REGION.x,
         # SEASON = SEASON.x,
         # ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
         # Isolate = Isolate.x)

## Importing the E.coli_ECOFF_EUCAST break point file
#ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx")

ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_E.coli_K.pneumoniae.xlsx")
EPI_CUTOFF<-ECOFF_EUCAST_BREAK_POINT

#calculate_ecoli_esbl_susceptibility <- function(joined_data, EPI_CUTOFF) {
  
 
  
#calculate_esbl_ecoli_resistance <- function(joined_data, EPI_CUTOFF) {
  
  # Define the list of grouping variables (Only defined ONCE here)
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
  
  # --- Step 1: Select, Rename, Filter, and Ensure Uniqueness (INIKA_ID) ---
  # library(tidyverse)
  # library(stringr)
  
  # Assuming 'joined_data' and 'EPI_CUTOFF' are loaded in your environment
  
  calculate_esbl_ecoli_resistance <- function(joined_data, EPI_CUTOFF) {
    
    # --- Step 1: Select, Rename, Filter, and Ensure Uniqueness (INIKA_ID) ---
    Data <- joined_data %>%
      select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
             VITEK_MS_Results, ESBL, AMX_ED10, AZM_ED15, CRO_ED30,
             CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
             OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
      rename(
        Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
        Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
        Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
        Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
        Oxytetracycline = "OXY_ED30",
        "Popymyxin_B(PB)" = "POL_ED300",
        "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
        Cefotaxime = "CTX_ED5",
        "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
      ) %>%
      
      # Filter for E. coli criteria (either Isolate.x or VITEK_MS_Results) AND ESBL == 1
      # Note: K.pneumoniae condition in Isolate.x is added based on the prompt's condition one.
      filter((Isolate %in% c("E.coli", "K.pneumoniae") | VITEK_MS_Results == "Escherichia coli") &
               ESBL == 1) %>%
      
      # Ensure only unique INIKA_ID are considered for the count (denominator)
      distinct(INIKA_ID, .keep_all = TRUE)
    
    # --- Step 2: Pivot the drug columns into long format ---
    AMR_Data <- pivot_longer(
      Data,
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value" 
    )
    
    # --- Step 3: Join with the cutoff table ---
    AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
    
    # Define the list of grouping variables for the final summary
    group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
    
    # --- Step 4: Clean and Classify Susceptibility (S/I/R) ---
    Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
      
      # Clean the S and R cutoff columns from non-numeric characters (e.g., '<', '>')
      # and convert all key values to numeric
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
          # Susceptible (S): Zone diameter >= S breakpoint
          Measured_Zone >= S_Cutoff ~ "S",
          # Resistant (R): Zone diameter <= R breakpoint
          Measured_Zone <= R_Cutoff ~ "R",
          # Intermediate (I): Zone diameter is between the R and S breakpoints
          Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
          TRUE ~ NA_character_ # Cases where data is missing or doesn't fit
        )
      ) %>%
      
      # Ensure only one result is counted per isolate for each antimicrobial test
      distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
      
      # Remove tests where the result could not be determined
      filter(!is.na(Categorical_Result))
    
    
    # --- Step 5: Define the Summary Function (Only Resistance with CI) ---
    
    calculate_resistance_ci_summary <- function(data, group_var) {
      
      # 1. Calculate Total Tests (Denominator) per group/antimicrobial
      Total_Count <- data %>%
        group_by(!!sym(group_var), Antimicrobial_substance) %>%
        summarise(Total_Tests = n(), .groups = 'drop')
      
      # 2. Count Resistance (Numerator)
      Resistance_Count <- data %>%
        filter(Categorical_Result == "R") %>% # Only count Resistant isolates
        group_by(!!sym(group_var), Antimicrobial_substance) %>%
        summarise(N_Resistant = n(), .groups = 'drop')
      
      # 3. Join counts and calculate percentage/CI
      Summary_Table <- Total_Count %>%
        left_join(Resistance_Count, by = c(group_var, "Antimicrobial_substance")) %>%
        # Replace NA (0 Resistance) with 0
        mutate(N_Resistant = replace_na(N_Resistant, 0)) %>%
        
        # Calculate Percentage and 95% CI using prop.test
        rowwise() %>%
        mutate(
          # Run prop.test only if Total_Tests > 0
          Prop_Test = ifelse(Total_Tests > 0, 
                             list(prop.test(x = N_Resistant, n = Total_Tests)), 
                             list(NA)),
          
          Percentage = round((N_Resistant / Total_Tests) * 100, 1),
          
          # Extract CI from prop.test result
          CI_Lower = ifelse(Total_Tests > 0, 
                            round(Prop_Test$conf.int[1] * 100, 1), NA),
          CI_Upper = ifelse(Total_Tests > 0, 
                            round(Prop_Test$conf.int[2] * 100, 1), NA),
          
          # Format the CI column
          `95%_CI` = paste0("(", CI_Lower, " - ", CI_Upper, ")"),
          
          Grouping_Variable = group_var
        ) %>%
        ungroup() %>%
        
        # Select and arrange columns for final output
        select(
          Grouping_Variable,
          Grouping_Value = !!sym(group_var),
          Antimicrobial_substance,
          Categorical_Result = R, # Explicitly label as 'R'
          Frequency = N_Resistant,
          Total_Tests,
          Percentage,
          `95%_CI`
        )
      
      return(as.data.frame(Summary_Table))
    }
    
    # --- Step 6: Apply the function to all grouping variables and combine results ---
    ECO_ESBL_Resistance_Summary <- bind_rows(
      lapply(group_vars_list, function(g_var) {
        calculate_resistance_ci_summary(Ecoli_ESBL_Susceptibility, g_var)
      })
    )
    
    return(ECO_ESBL_Resistance_Summary)
  }
 ###############################################################################
  
  # --- 1. Define the Processing Function ---
  process_ecoli_data <- function(joined_data, EPI_CUTOFF) {
    
    # Step 1: Select, Rename, and Filter based on your specific criteria
    Data <- joined_data %>%
      select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
             VITEK_MS_Results, ESBL, AMX_ED10, AZM_ED15, CRO_ED30,
             CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
             OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
      rename(
        Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
        Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
        Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
        Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
        Oxytetracycline = "OXY_ED30",
        "Popymyxin_B(PB)" = "POL_ED300",
        "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
        Cefotaxime = "CTX_ED5",
        "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
      ) %>%
      # Filter: E.coli in Isolate AND Escherichia coli in VITEK AND CTX_ED5 < 21
      filter(
        Isolate %in% c("E.coli", "K.pneumoniae"),
        VITEK_MS_Results == "Escherichia coli"
         & Ceftriaxone <23, # Cefotaxime < 22 # Note: We use the renamed column "Cefotaxime" here
      ) %>%
      distinct(INIKA_ID, .keep_all = TRUE)
    
    # Step 2: Pivot to long format
    AMR_Data <- pivot_longer(
      Data,
      cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                ORIGIN_OF_SAMPLE, ESBL),
      names_to = "Antimicrobial_substance",
      values_to = "value" 
    )
    
    # Step 3: Join with Cutoff
    AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
    
    # Step 4: Classify Susceptibility
    Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
      mutate(
        S_Cleaned = as.numeric(str_replace_all(S, "[^0-9.]", "")),
        R_Cleaned = as.numeric(str_replace_all(R, "[^0-9.]", "")),
        Measured_Zone = as.numeric(str_trim(value))
      ) %>%
      mutate(
        Categorical_Result = case_when(
          Measured_Zone >= S_Cleaned ~ "S",
          Measured_Zone <= R_Cleaned ~ "R",
          Measured_Zone > R_Cleaned & Measured_Zone < S_Cleaned ~ "I",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Categorical_Result))
    
    # Step 5: Summary Logic
    group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
    
    calculate_summary <- function(df, group_var) {
      df %>%
        group_by(!!sym(group_var), Antimicrobial_substance) %>%
        summarise(
          Total_Tests = n(),
          N_Resistant = sum(Categorical_Result == "R"),
          .groups = 'drop'
        ) %>%
        rowwise() %>%
        mutate(
          Prop_Test = list(prop.test(N_Resistant, Total_Tests)),
          Percentage = round((N_Resistant / Total_Tests) * 100, 1),
          CI_Lower = round(Prop_Test$conf.int[1] * 100, 1),
          CI_Upper = round(Prop_Test$conf.int[2] * 100, 1),
          `95%_CI` = paste0("(", CI_Lower, " - ", CI_Upper, ")"),
          Grouping_Variable = group_var
        ) %>%
        select(Grouping_Variable, Grouping_Value = !!sym(group_var), 
               Antimicrobial_substance, N_Resistant, Total_Tests, Percentage, `95%_CI`)
    }
    
    # Combine results
    Final_Summary <- map_dfr(group_vars_list, ~calculate_summary(Ecoli_ESBL_Susceptibility, .x))
    
    return(Final_Summary)
  }
  
  # --- 2. Run the code ---
   result_table <- process_ecoli_data(joined_data, EPI_CUTOFF)
  
  # --- 3. Save the results ---
  
  # Option A: Save as CSV (Best for simple data sharing)
   write_csv(result_table, "Results/28_01_26Ecoli_ESBL_Resistance_Summary.csv")
  
  # Option B: Save as Excel (Best for reports)
   write_xlsx(result_table, "Results/19_12_25Ecoli_ESBL_Resistance_Summary.xlsx")
################################################################################
   # library(dplyr)
   # library(tidyr)
   # library(stringr)
   # library(purrr)
   # library(readr)
   # library(writexl)
   # library(flextable)
   # library(officer)
   # 
   # process_ecoli_data <- function(joined_data, EPI_CUTOFF) {
   #   
   #   # Ensure export folder exists
   #   if (!dir.exists("Results")) {
   #     dir.create("Results")
   #   }
   #   
   #   # --- Step 1: Define Antibiotic Class Mapping ---
   #   drug_class_map <- tibble::tribble(
   #     ~Antimicrobial_substance,         ~Antibiotic_Class,
   #     "Amoxicillin",                     "Penicillins",
   #     "Ceftriaxone",                     "3rd Gen Cephalosporins",
   #     "Cefotaxime",                      "3rd Gen Cephalosporins",
   #     "Cefotaxime/ClavulanicAcid",      "Beta-lactam / Inhibitor combinations",
   #     "Meropenem",                      "Carbapenems",
   #     "Azithromycin",                    "Macrolides",
   #     "Ciprofloxacin",                  "Fluoroquinolones",
   #     "Gentamicin",                      "Aminoglycosides",
   #     "Doxycycline",                     "Tetracyclines",
   #     "Oxytetracycline",                "Tetracyclines",
   #     "Florfenicol",                    "Phenicols",
   #     "Popymyxin_B(PB)",                "Polymyxins",
   #     "Sulfamethoxazole/Trimethoprim",  "Folate Pathway Inhibitors"
   #   )
   #   
   #   # --- Step 2: Strict Filtering of EPI_CUTOFF by EPI_ECOFF column ---
   #   # Retains ONLY antibiotics that have a valid, non-missing ECOFF value in EPI_ECOFF
   #   valid_cutoff <- EPI_CUTOFF %>%
   #     filter(!is.na(Antimicrobial_substance) & str_trim(as.character(Antimicrobial_substance)) != "")
   #   
   #   if ("EPI_ECOFF" %in% names(valid_cutoff)) {
   #     valid_cutoff <- valid_cutoff %>%
   #       filter(
   #         !is.na(EPI_ECOFF) & 
   #           str_trim(as.character(EPI_ECOFF)) != "" & 
   #           !str_trim(as.character(EPI_ECOFF)) %in% c("-", "ND", "NA", "N/A")
   #       )
   #   } else if ("ECOFF" %in% names(valid_cutoff)) {
   #     valid_cutoff <- valid_cutoff %>%
   #       filter(
   #         !is.na(ECOFF) & 
   #           str_trim(as.character(ECOFF)) != "" & 
   #           !str_trim(as.character(ECOFF)) %in% c("-", "ND", "NA", "N/A")
   #       )
   #   }
   #   
   #   # --- Step 3: Select, Rename, and Filter Base Data ---
   #   Data <- joined_data %>%
   #     select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
   #            VITEK_MS_Results, ESBL, PROTOCOL, AMX_ED10, AZM_ED15, CRO_ED30,
   #            CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
   #            OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
   #     rename(
   #       Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
   #       Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
   #       Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
   #       Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
   #       Oxytetracycline = "OXY_ED30",
   #       "Popymyxin_B(PB)" = "POL_ED300",
   #       "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
   #       Cefotaxime = "CTX_ED5",
   #       "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
   #     ) %>%
   #     filter(
   #       Isolate %in% c("E.coli", "K.pneumoniae"),
   #       VITEK_MS_Results == "Escherichia coli",
   #       Ceftriaxone < 23
   #     ) %>%
   #     distinct(INIKA_ID, .keep_all = TRUE)
   #   
   #   # --- Step 4: Pivot to Long Format and Inner Join Valid ECOFF Cutoffs ---
   #   AMR_Data <- pivot_longer(
   #     Data,
   #     cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
   #               ORIGIN_OF_SAMPLE, ESBL, PROTOCOL),
   #     names_to = "Antimicrobial_substance",
   #     values_to = "value" 
   #   )
   #   
   #   # Inner join guarantees that ONLY antibiotics with valid EPI_ECOFF values remain
   #   AMR_Cutoff <- inner_join(AMR_Data, valid_cutoff, by = "Antimicrobial_substance")
   #   
   #   # --- Step 5: Classify Resistance ---
   #   Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
   #     filter(!is.na(value) & str_trim(as.character(value)) != "") %>%
   #     mutate(
   #       S_Cleaned = as.numeric(str_replace_all(S, "[^0-9.]", "")),
   #       R_Cleaned = as.numeric(str_replace_all(R, "[^0-9.]", "")),
   #       Measured_Zone = as.numeric(str_trim(value)),
   #       
   #       # Clean and normalize text variables for robust string matching
   #       ORIGIN_OF_SAMPLE_CLEAN = str_squish(tolower(as.character(ORIGIN_OF_SAMPLE))),
   #       REGION_CLEAN = str_squish(tolower(as.character(REGION))),
   #       SEASON_CLEAN = str_squish(tolower(as.character(SEASON)))
   #     ) %>%
   #     mutate(
   #       Categorical_Result = case_when(
   #         Measured_Zone >= S_Cleaned ~ "S",
   #         Measured_Zone <= R_Cleaned ~ "R",
   #         Measured_Zone > R_Cleaned & Measured_Zone < S_Cleaned ~ "I",
   #         TRUE ~ NA_character_
   #       )
   #     )
   #   
   #   # --- Helper Function for Formatting Percentages ---
   #   fmt_pct <- function(val) {
   #     if (is.null(val) || length(val) == 0 || is.na(val)) return("-")
   #     if (val %% 1 == 0) {
   #       return(paste0(as.integer(val), "%"))
   #     } else {
   #       return(paste0(sprintf("%.1f", val), "%"))
   #     }
   #   }
   #   
   #   # --- Step 6: Compute Overall Resistance (WITH 95% CI) ---
   #   overall_res <- Ecoli_ESBL_Susceptibility %>%
   #     group_by(Antimicrobial_substance) %>%
   #     summarise(
   #       N_Total = n(),
   #       n_Resistant = sum(Categorical_Result == "R", na.rm = TRUE),
   #       .groups = "drop"
   #     ) %>%
   #     rowwise() %>%
   #     mutate(
   #       Raw_Pct = if_else(N_Total > 0, (n_Resistant / N_Total) * 100, NA_real_),
   #       Pct_Str = fmt_pct(if_else(!is.na(Raw_Pct), round(Raw_Pct, 1), NA_real_)),
   #       
   #       ci_obj = list(if (N_Total > 0) tryCatch(prop.test(n_Resistant, N_Total, conf.level = 0.95)$conf.int, error = function(e) c(NA_real_, NA_real_)) else c(NA_real_, NA_real_)),
   #       CI_Low = if_else(!is.null(ci_obj[[1]]), round(ci_obj[[1]][1] * 100, 1), NA_real_),
   #       CI_High = if_else(!is.null(ci_obj[[1]]), round(ci_obj[[1]][2] * 100, 1), NA_real_),
   #       
   #       CI_Low_Str = if_else(!is.na(CI_Low), if_else(CI_Low %% 1 == 0, as.character(as.integer(CI_Low)), sprintf("%.1f", CI_Low)), "-"),
   #       CI_High_Str = if_else(!is.na(CI_High), if_else(CI_High %% 1 == 0, as.character(as.integer(CI_High)), sprintf("%.1f", CI_High)), "-"),
   #       
   #       `Overall Resistance [%, n/N (95% CI)]` = if_else(
   #         N_Total == 0 | is.na(CI_Low), 
   #         "-", 
   #         paste0(Pct_Str, ", ", n_Resistant, "/", N_Total, " (", CI_Low_Str, "-", CI_High_Str, ")")
   #       )
   #     ) %>%
   #     ungroup() %>%
   #     select(Antimicrobial_substance, `Total Tested Sample (N)` = N_Total, `Overall Resistance [%, n/N (95% CI)]`)
   #   
   #   # --- Flexible Helper Function for Stratified Resistance ---
   #   calc_stratified_res <- function(df, group_col, pattern) {
   #     df %>%
   #       filter(str_detect(!!sym(group_col), pattern)) %>%
   #       group_by(Antimicrobial_substance) %>%
   #       summarise(
   #         N_Sub = n(),
   #         n_Resistant = sum(Categorical_Result == "R", na.rm = TRUE),
   #         .groups = "drop"
   #       ) %>%
   #       rowwise() %>%
   #       mutate(
   #         Raw_Pct = if_else(N_Sub > 0, (n_Resistant / N_Sub) * 100, NA_real_),
   #         Pct_Str = fmt_pct(if_else(!is.na(Raw_Pct), round(Raw_Pct, 1), NA_real_)),
   #         Formatted_Val = if_else(N_Sub == 0 | is.na(Raw_Pct), "-", paste0(Pct_Str, ", ", n_Resistant, "/", N_Sub))
   #       ) %>%
   #       ungroup() %>%
   #       select(Antimicrobial_substance, Formatted_Val)
   #   }
   #   
   #   # --- Step 7: Stratifications ---
   #   origin_adult <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "ORIGIN_OF_SAMPLE_CLEAN", "adult|outpatient") %>%
   #     rename(`Adult Outpatient [%, n/N]` = Formatted_Val)
   #   
   #   origin_school <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "ORIGIN_OF_SAMPLE_CLEAN", "school") %>%
   #     rename(`Schoolchildren [%, n/N]` = Formatted_Val)
   #   
   #   region_kili <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "REGION_CLEAN", "kili") %>%
   #     rename(`Kilimanjaro [%, n/N]` = Formatted_Val)
   #   
   #   region_mwanza <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "REGION_CLEAN", "mwanza") %>%
   #     rename(`Mwanza [%, n/N]` = Formatted_Val)
   #   
   #   season_dry <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "SEASON_CLEAN", "dry") %>%
   #     rename(`Dry Season [%, n/N]` = Formatted_Val)
   #   
   #   season_wet <- calc_stratified_res(Ecoli_ESBL_Susceptibility, "SEASON_CLEAN", "wet") %>%
   #     rename(`Wet Season [%, n/N]` = Formatted_Val)
   #   
   #   # --- Step 8: Combine Columns, Reorder (Class First, Antibiotic Second) ---
   #   Final_Table <- drug_class_map %>%
   #     inner_join(valid_cutoff %>% select(Antimicrobial_substance) %>% distinct(), by = "Antimicrobial_substance") %>%
   #     left_join(overall_res, by = "Antimicrobial_substance") %>%
   #     left_join(origin_adult, by = "Antimicrobial_substance") %>%
   #     left_join(origin_school, by = "Antimicrobial_substance") %>%
   #     left_join(region_kili, by = "Antimicrobial_substance") %>%
   #     left_join(region_mwanza, by = "Antimicrobial_substance") %>%
   #     left_join(season_dry, by = "Antimicrobial_substance") %>%
   #     left_join(season_wet, by = "Antimicrobial_substance") %>%
   #     rename(
   #       `Antibiotic Class` = Antibiotic_Class,
   #       `Antimicrobial Tested` = Antimicrobial_substance
   #     ) %>%
   #     select(
   #       `Antibiotic Class`,
   #       `Antimicrobial Tested`,
   #       everything()
   #     ) %>%
   #     arrange(`Antibiotic Class`, `Antimicrobial Tested`)
   #   
   #   Final_Table[is.na(Final_Table)] <- "-"
   #   
   #   # --- Step 9: Enhanced Visibility Flextable Formatting ---
   #   std_border <- fp_border(color = "#B0C4DE", width = 1)
   #   
   #   ft <- flextable(Final_Table) %>%
   #     theme_vanilla() %>%
   #     fontsize(size = 10, part = "body") %>%
   #     fontsize(size = 11, part = "header") %>%
   #     padding(padding.top = 5, padding.bottom = 5, padding.left = 5, padding.right = 5, part = "all") %>%
   #     align(j = 1:2, align = "left", part = "all") %>%
   #     align(j = 3:ncol(Final_Table), align = "center", part = "all") %>%
   #     bold(part = "header") %>%
   #     bold(j = 1, part = "body") %>%
   #     bg(part = "header", bg = "#1F497D") %>%
   #     color(part = "header", color = "white") %>%
   #     bg(i = seq(2, nrow(Final_Table), by = 2), bg = "#F2F4F8", part = "body") %>%
   #     border_inner_h(border = std_border, part = "body") %>%
   #     border_inner_v(border = std_border, part = "all") %>%
   #     autofit() %>%
   #     fit_to_width(max_width = 11)
   #   
   #   doc <- read_docx() %>%
   #     body_end_section_landscape() %>%
   #     body_add_par("Table: Overall and Stratified AMR Resistance Summary for Confirmed ESC Resistant E. coli", style = "heading 1") %>%
   #     body_add_flextable(ft) %>%
   #     body_end_section_landscape()
   #   
   #   word_path <- "Results/Ecoli_ESBL_Resistance_Stratified_Summary.docx"
   #   print(doc, target = word_path)
   #   
   #   # Export CSV and Excel
   #   write_csv(Final_Table, "Results/Ecoli_ESBL_Resistance_Stratified_Summary.csv")
   #   write_xlsx(Final_Table, "Results/Ecoli_ESBL_Resistance_Stratified_Summary.xlsx")
   #   
   #   cat("\nResults successfully saved to Results/ folder in Word, CSV, and Excel formats.\n")
   #   
   #   return(Final_Table)
   # }
   # 
   # # --- Function Execution ---
   # result_table <- process_ecoli_data(joined_data, EPI_CUTOFF)
##########################################################   
   
   
   # Example of how to call the function (assuming you have joined_data and EPI_CUTOFF)
  # ECO_ESBL_Resistance_Results <- calculate_esbl_ecoli_resistance(joined_data, EPI_CUTOFF)
  
  # return(ECO_ESBL_Resistance_Results)
  # Example of how to call the function (assuming you have joined_data and EPI_CUTOFF)
  ECO_ESBL_Resistance_Results <- calculate_esbl_ecoli_resistance(joined_data, EPI_CUTOFF)
  
  # return(ECO_ESBL_Resistance_Results)
  
  
  
  
  
  
  
  
#   # 1. Select relevant columns and rename drug columns for clarity
#   Data <- joined_data %>%
#     select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
#            VITEK_MS_Results, ESBL, AMX_ED10, AZM_ED15, CRO_ED30,
#            CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
#            OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
#     rename(
#       Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
#       Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
#       Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
#       Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
#       Oxytetracycline = "OXY_ED30",
#       "Popymyxin_B(PB)" = "POL_ED300",
#       "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
#       Cefotaxime = "CTX_ED5",
#       "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
#     )
#   
#   # 2. Pivot the drug columns into long format
#   AMR_Data <- pivot_longer(
#     Data,
#     cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
#               ORIGIN_OF_SAMPLE, ESBL),
#     names_to = "Antimicrobial_substance",
#     values_to = "value" 
#   )
#   
#   # 3. Join with the cutoff table
#   AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
#   
#   # Define the list of grouping variables for the final summary
#   group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
#   
#   # --- Step 4: Filter, Clean, and Classify Susceptibility (S/I/R) ---
#   
#   Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
#     
#     # Filter for E. coli AND true ESBL (ESBL == 1)
#     filter(VITEK_MS_Results == "Escherichia coli",
#            ESBL == 1) %>%
#     
#     # Clean the S and R cutoff columns from non-numeric characters (e.g., '<', '>')
#     # and convert all key values to numeric
#     mutate(
#       S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
#       R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
#       
#       Measured_Zone = as.numeric(str_trim(value)),
#       S_Cutoff = as.numeric(S_Cleaned),
#       R_Cutoff = as.numeric(R_Cleaned)
#     ) %>%
#     
#     # Implement the breakpoint logic for zone diameters (in mm)
#     mutate(
#       Categorical_Result = case_when(
#         # Susceptible (S): Zone diameter >= S breakpoint
#         Measured_Zone >= S_Cutoff ~ "S",
#         # Resistant (R): Zone diameter <= R breakpoint
#         Measured_Zone <= R_Cutoff ~ "R",
#         # Intermediate (I): Zone diameter is between the R and S breakpoints
#         Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
#         TRUE ~ NA_character_ # Cases where data is missing or doesn't fit
#       )
#     ) %>%
#     
#     # Ensure only one result is counted per isolate for each antimicrobial test
#     distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
#     
#     # Remove tests where the result could not be determined
#     filter(!is.na(Categorical_Result))
#   
#   
#   # --- Step 5: Define the Summary Function ---
#   
#   calculate_susceptibility_summary <- function(data, group_var) {
#     
#     # Total count (denominator) for the current grouping variable and antimicrobial
#     Total_Count <- data %>%
#       group_by(!!sym(group_var), Antimicrobial_substance) %>%
#       summarise(Total_Tests = n(), .groups = 'drop')
#     
#     # Count the S, I, R results
#     Summary_Table <- data %>%
#       group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
#       summarise(N = n(), .groups = 'drop_last') %>%
#       
#       # Join the total counts back
#       left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
#       
#       # Calculate Percentage (rounded to 1 decimal place)
#       mutate(
#         Percentage = round((N / Total_Tests) * 100, 1),
#         Grouping_Variable = group_var
#       ) %>%
#       
#       # Select and arrange columns for final output
#       select(
#         Grouping_Variable,
#         Grouping_Value = !!sym(group_var),
#         Antimicrobial_substance,
#         Categorical_Result,
#         Frequency = N,
#         Total_Tests,
#         Percentage
#       )
#     
#     return(as.data.frame(Summary_Table))
#   }
#   
#   # --- Step 6: Apply the function to all grouping variables and combine results ---
#   ECO_ESBL_Susceptibility_Summary <- bind_rows(
#     lapply(group_vars_list, function(g_var) {
#       calculate_susceptibility_summary(Ecoli_ESBL_Susceptibility, g_var)
#     })
#   )
#   
#   return(ECO_ESBL_Susceptibility_Summary)
# }
# Save the file
write_tsv(ECO_ESBL_Susceptibility_Summary, "Results/ECO_ESBL_Susceptibility_Summary-9-11-2025.tsv")
######################################################################
  
# AMR Frequencies for Klebsiella pneumoniae
 
  ## Importing the K.pneumoniae_ECOFF_EUCAST break point file
  
  ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx", sheet = 2)
  
  
  EPI_CUTOFF<-ECOFF_EUCAST_BREAK_POINT
  
  
  process_kpn_esbl_susceptibility <- function(joined_data, EPI_CUTOFF) {
    
    # Step 1: Select, Rename, and Filter
    Filtered_Wide_Data <- joined_data %>%
      select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
             VITEK_MS_Results, ESBL, AMX_ED10, AZM_ED15, CRO_ED30,
             CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
             OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
      rename(
        Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
        Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
        Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
        Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
        Oxytetracycline = "OXY_ED30",
        "Popymyxin_B(PB)" = "POL_ED300",
        "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
        Cefotaxime = "CTX_ED5",
        "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
      ) %>%
      filter(
        Isolate == "K.pneumoniae",
        VITEK_MS_Results == "Klebsiella pneumoniae",
        Cefotaxime < 21 & Ceftriaxone < 23,
      ) %>%
      distinct(INIKA_ID, .keep_all = TRUE)
    
    # Step 2: Pivot to long format
    AMR_Long <- Filtered_Wide_Data %>%
      pivot_longer(
        cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
                  ORIGIN_OF_SAMPLE, ESBL),
        names_to = "Antimicrobial_substance",
        values_to = "value" 
      )
    
    # Step 3: Join with Cutoff
    AMR_Cutoff <- left_join(AMR_Long, EPI_CUTOFF, by = "Antimicrobial_substance")
    
    # Step 4: Classify Susceptibility
    Kpn_Classified <- AMR_Cutoff %>%
      mutate(
        S_Cleaned = as.numeric(str_replace_all(S, "[^0-9.]", "")),
        R_Cleaned = as.numeric(str_replace_all(R, "[^0-9.]", "")),
        Measured_Zone = as.numeric(str_trim(value))
      ) %>%
      mutate(
        Categorical_Result = case_when(
          Measured_Zone >= S_Cleaned ~ "S",
          Measured_Zone <= R_Cleaned ~ "R",
          Measured_Zone > R_Cleaned & Measured_Zone < S_Cleaned ~ "I",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Categorical_Result))
    
    # Step 5: Summary Calculation
    group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
    
    calculate_summary <- function(df, group_var) {
      df %>%
        group_by(!!sym(group_var), Antimicrobial_substance) %>%
        summarise(
          Total_Tests = n(),
          N_Resistant = sum(Categorical_Result == "R"),
          .groups = 'drop'
        ) %>%
        rowwise() %>%
        mutate(
          # FIXED: Using binom.test instead of prop.test to avoid Chi-sq warnings
          # Added a check for Total_Tests > 0 to prevent errors
          BT = if(Total_Tests > 0) list(binom.test(N_Resistant, Total_Tests)) else list(NULL),
          
          Percentage = round((N_Resistant / Total_Tests) * 100, 1),
          
          CI_Lower = if(!is.null(BT)) round(BT$conf.int[1] * 100, 1) else NA,
          CI_Upper = if(!is.null(BT)) round(BT$conf.int[2] * 100, 1) else NA,
          
          `95%_CI` = paste0("(", CI_Lower, " - ", CI_Upper, ")"),
          Grouping_Variable = group_var
        ) %>%
        select(Grouping_Variable, Grouping_Value = !!sym(group_var), 
               Antimicrobial_substance, N_Resistant, Total_Tests, Percentage, `95%_CI`)
    }
    
    Final_Summary <- map_dfr(group_vars_list, ~calculate_summary(Kpn_Classified, .x))
    
    return(as.data.frame(Final_Summary))
  }
  
  # --- 2. Run and Save ---
   result_table <- process_kpn_esbl_susceptibility(joined_data, EPI_CUTOFF)
  
  # --- 3. Save the results ---
  
  # Option A: Save as CSV (Best for simple data sharing)
  write_csv(result_table, "Results/19_12_25Kpn_ESBL_Resistance_Summary.csv")
  
  # Option B: Save as Excel (Best for reports)
  write_xlsx(result_table, "Results/19_12_25Kpn_ESBL_Resistance_Summary.xlsx")
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
  Zone_Distribution_Table <- joined_data %>%  ## 
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
    select(INIKA_ID,REGION,SEASON,ORIGIN_OF_SAMPLE,Isolate,
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
  
  AMR_Data<-pivot_longer(Data, cols=-c(INIKA_ID,
                                       VITEK_MS_Results,
                                       Isolate,REGION,SEASON,
                                       ORIGIN_OF_SAMPLE, ESBL), # <--- ESBL ADDED HERE
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
  #######################################################
  # 22/9/2026
 #  # AMR prevalence and zone of inhibition distribution
 #  # Function to calculate
 #  binom.test(0,430, p = 0.5, "two.sided", conf.level = 0.95)
 #  binom.test(0,430, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100
 #  
 #  
 #  # This will be reported as follows: None of the samples were positive for Carbapeneamse resistant E.coli [95% CI: 0.0 – 1.1]
 #  # 
 #  # If you have found 10 ESC resistant Klebsiella from the unique farms you will calculate accordingly:
 #  
 #  binom.test(10,430, p = 0.5, "two.sided", conf.level = 0.95)
 #  binom.test(10,430, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100
 #  
 #  # This (example) will be reported as follows: 2.3% [95% CI: 1.1 – 4.2] of the samples were positive for ESC-resistant Klebsiella.
 #  # 
 #  # For the % Resistances to different antibiotics you will only use the total number of the isolates you have within each category, and of course it might be confusing if you first detected more ESC resistant E.coli from the selective plates and then you have less ESC resistant E.coli according to the criteria <22mm for Cefotaxime. But it should be the last criteria that will decide how many you have! We may want to change the criteria for defining the isolates as resistant after checking the distributions. This is therefore also very important.I think this should be discussed.
 #  # 
 #  # Example: But for now if you have 110 ESC resistant isolates of for example E.coli according only to be growing on the selective plate and where the VITEK has also confirmed that it is an E.coli, but from these only 90 isolates are resistant to Cefotaxime. Then the total number of ESC resistant E.coli will be 90. If then 10 of these are resistant to for a substance "X", you calculate the % Resistance to Substance X as follows: 10/90, but you also want to include the 95% confidence interval so then you do as shown here:
 #  
 #  binom.test(10,90, p = 0.5, "two.sided", conf.level = 0.95)
 #  binom.test(10,90, p = 0.5, "two.sided", conf.level = 0.95)$conf.int*100
 #  #library(here)
 #  # Importing the file
 # Final_JoinedDATA <- read_csv("data/CLEANED_DATA/Final_JoinedDATA.csv")
 # names(Final_JoinedDATA)
 # 
 # ## Importing the E.coli_ECOFF_EUCAST
 # ECOFF_E_coli<- read_excel("data/ECOFF_E.coli_K.pneumoniae.xlsx",
 #                                         sheet = "E.coli_ECOFF_EUCAST")
 # 
 # 
 # View(ECOFF_E_coli)
 # 
 # EPI_CUTOFF<-ECOFF_E_coli
 # 
 # # Define the list of grouping variables (Only defined ONCE here)
 # group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
 # 
 # # Helper function to dynamically format numbers (1 decimal place, drop .0 for integers)
 # fmt_num <- function(x) {
 #   if (is.na(x)) return("-")
 #   rounded <- round(x, 1)
 #   if (rounded %% 1 == 0) {
 #     return(as.character(as.integer(rounded)))
 #   } else {
 #     return(sprintf("%.1f", rounded))
 #   }
 # }
 # 
 # # Vectorized wrapper for fmt_num
 # fmt_num_vec <- function(vec) {
 #   sapply(vec, fmt_num)
 # }
 # 
 # # ==============================================================================
 # # 1. DEFINE MAIN PROCESSING FUNCTION
 # # ==============================================================================
 # 
 # process_amr_summary <- function(Final_JoinedDATA, EPI_CUTOFF) {
 #   
 #   # --- Step 1: Select, Rename, Filter, and Ensure Uniqueness ---
 #   Data <- Final_JoinedDATA %>%
 #     select(
 #       INIKA_ID, PROTOCOL, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
 #       VITEK_MS_Results, ESBL_Presumptivefinal, AMX_ED10, AZM_ED15, CRO_ED30,
 #       CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
 #       OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30
 #     ) %>%
 #     rename(
 #       Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
 #       Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
 #       Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
 #       Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
 #       Oxytetracycline = "OXY_ED30",
 #       "Polymyxin_B(PB)" = "POL_ED300",
 #       "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
 #       Cefotaxime = "CTX_ED5",
 #       "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
 #     ) %>%
 #     filter(
 #       Ceftriaxone <= 23,
 #       ESBL_Presumptivefinal == 1,
 #       str_detect(PROTOCOL, regex("CGR3", ignore_case = TRUE)),
 #       (
 #         str_detect(Isolate, regex("coli|pneumoniae", ignore_case = TRUE)) |
 #           str_detect(VITEK_MS_Results, regex("Escherichia coli|Klebsiella pneumoniae", ignore_case = TRUE))
 #       )
 #     ) %>%
 #     distinct(INIKA_ID, .keep_all = TRUE)
 #   
 #   # --- Step 2: Pivot drug columns into long format ---
 #   AMR_Data <- pivot_longer(
 #     Data,
 #     cols = -c(INIKA_ID, PROTOCOL, VITEK_MS_Results, Isolate, REGION, SEASON,
 #               ORIGIN_OF_SAMPLE, ESBL_Presumptivefinal),
 #     names_to = "Antimicrobial_substance",
 #     values_to = "value"
 #   )
 #   
 #   # --- Step 3: Join with the cutoff table ---
 #   AMR_Cutoff <- left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")
 #   
 #   # --- Step 4: Clean and Classify Susceptibility (S/I/R) ---
 #   Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
 #     mutate(
 #       S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
 #       R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
 #       Measured_Zone = as.numeric(str_trim(value)),
 #       S_Cutoff = as.numeric(S_Cleaned),
 #       R_Cutoff = as.numeric(R_Cleaned),
 #       Categorical_Result = case_when(
 #         Measured_Zone >= S_Cutoff ~ "S",
 #         Measured_Zone <= R_Cutoff ~ "R",
 #         Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
 #         TRUE ~ NA_character_
 #       )
 #     ) %>%
 #     distinct(INIKA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
 #     filter(!is.na(Categorical_Result))
 #   
 #   # --- Step 5: Internal Summary Helper Function ---
 #   calculate_resistance_ci_summary <- function(data, group_var) {
 #     
 #     data %>%
 #       group_by(Antimicrobial_substance, Grouping_Value = .data[[group_var]]) %>%
 #       summarise(
 #         Total_Tests = n(),
 #         N_Resistant = sum(Categorical_Result == "R", na.rm = TRUE),
 #         .groups = "drop"
 #       ) %>%
 #       mutate(
 #         Grouping_Variable = group_var,
 #         Prop_Raw = ifelse(Total_Tests > 0, (N_Resistant / Total_Tests) * 100, NA_real_),
 #         
 #         # Calculate Wilson 95% Confidence Intervals safely
 #         Prop_Test = map2(N_Resistant, Total_Tests, function(x, y) {
 #           if (y > 0) {
 #             prop.test(x, y)
 #           } else {
 #             NULL
 #           }
 #         }),
 #         CI_Lower_Raw = map_dbl(Prop_Test, function(res) {
 #           if (!is.null(res)) res$conf.int[1] * 100 else NA_real_
 #         }),
 #         CI_Upper_Raw = map_dbl(Prop_Test, function(res) {
 #           if (!is.null(res)) res$conf.int[2] * 100 else NA_real_
 #         }),
 #         
 #         # Clean formatting without brackets and optional decimals
 #         Pct_Formatted = fmt_num_vec(Prop_Raw),
 #         CI_Lower_Formatted = fmt_num_vec(CI_Lower_Raw),
 #         CI_Upper_Formatted = fmt_num_vec(CI_Upper_Raw),
 #         
 #         `Resistant (n, %)` = paste0(N_Resistant, " (", Pct_Formatted, "%)"),
 #         `95%_CI` = if_else(!is.na(CI_Lower_Raw), paste0(CI_Lower_Formatted, " - ", CI_Upper_Formatted), "-")
 #       ) %>%
 #       # Order strictly: Antibiotic -> Stratifier -> Subgroup -> Metrics
 #       select(
 #         Antimicrobial_substance,
 #         Grouping_Variable,
 #         Grouping_Value,
 #         Total_Tests,
 #         `Resistant (n, %)`,
 #         `95%_CI`
 #       )
 #   }
 #   
 #   # --- Step 6: Map across demographic variables ---
 #   group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
 #   
 #   map_dfr(group_vars_list, ~ calculate_resistance_ci_summary(Ecoli_ESBL_Susceptibility, .x)) %>%
 #     arrange(Antimicrobial_substance, Grouping_Variable, Grouping_Value)
 # }
 # 
 # # ==============================================================================
 # # 2. RUN PIPELINE & BUILD PUBLICATION TABLE
 # # ==============================================================================
 # 
 # # Execute pipeline
 # summary_df <- process_amr_summary(Final_JoinedDATA, EPI_CUTOFF)
 # 
 # # Format flextable (APA / ICMJE standard)
 # ft <- summary_df %>%
 #   flextable() %>%
 #   set_header_labels(
 #     Antimicrobial_substance = "Antimicrobial Agent",
 #     Grouping_Variable = "Stratifier",
 #     Grouping_Value = "Subgroup",
 #     Total_Tests = "Tested (N)",
 #     `Resistant (n, %)` = "Resistant (n, %)",
 #     `95%_CI` = "95% CI"
 #   ) %>%
 #   font(fontname = "Times New Roman", part = "all") %>%
 #   fontsize(size = 9, part = "all") %>%
 #   bold(part = "header") %>%
 #   align(j = c("Total_Tests", "Resistant (n, %)", "95%_CI"), align = "center", part = "all") %>%
 #   align(j = c("Antimicrobial_substance", "Grouping_Variable", "Grouping_Value"), align = "left", part = "all") %>%
 #   border_remove() %>%
 #   hline_top(border = fp_border(color = "black", width = 1.5), part = "header") %>%
 #   hline_bottom(border = fp_border(color = "black", width = 1.0), part = "header") %>%
 #   hline_bottom(border = fp_border(color = "black", width = 1.5), part = "body") %>%
 #   padding(padding = 3, part = "all") %>%
 #   merge_v(j = c("Antimicrobial_substance", "Grouping_Variable", "Grouping_Value")) %>% # <--- MOVED HERE
 #   autofit()
 # # Build Word document structure
 # doc <- read_docx() %>%
 #   body_add_par(
 #     "Table 1. Antimicrobial Resistance Profiles and 95% Confidence Intervals Stratified by Demographic Factors",
 #     style = "Normal"
 #   ) %>%
 #   body_add_flextable(ft) %>%
 #   body_add_par(
 #     "Note: 95% Confidence Intervals (95% CI) computed using Pearson's chi-squared test for proportions (prop.test). Numbers reflect 1 decimal place where applicable, ignoring zeros for whole numbers.",
 #     style = "Normal"
 #   )
 # 
 # # ==============================================================================
 # # 3. DIRECT FILE SAVING TO "Results" FOLDER
 # # ==============================================================================
 # 
 # # Guarantee folder exists
 # if (!dir.exists("Results")) {
 #   dir.create("Results")
 # }
 # 
 # # Save target Word Document directly to Results directory
 # print(doc, target = "Results/ESBL_Antimicrobial_Resistance_Summary.docx")
 # 
 # cat("\nSUCCESS: Updated publication-formatted table saved to 'Results/ESBL_Antimicrobial_Resistance_Summary.docx'!\n")
 # ####################################################
 # Trial2
 # # Helper function to dynamically format numbers (1 decimal place, drop .0 for integers)
 # fmt_num <- function(x) {
 #   if (is.na(x)) return("-")
 #   rounded <- round(x, 1)
 #   if (rounded %% 1 == 0) {
 #     return(as.character(as.integer(rounded)))
 #   } else {
 #     return(sprintf("%.1f", rounded))
 #   }
 # }
 # 
 # # Vectorized wrapper for fmt_num
 # fmt_num_vec <- function(vec) {
 #   sapply(vec, fmt_num)
 # }
 # 
 # # ==============================================================================
 # # 1. DEFINE MAIN PROCESSING FUNCTION
 # # ==============================================================================
 # 
 # process_amr_summary <- function(Final_JoinedDATA, EPI_CUTOFF) {
 #   
 #   # --- Step 1: Base CGR3 Filtering & Species Confirmation ---
 #   base_cgr3 <- Final_JoinedDATA %>%
 #     mutate(across(where(is.character), str_trim)) %>%
 #     filter(
 #       str_detect(PROTOCOL, regex("CGR3", ignore_case = TRUE)),
 #       str_detect(VITEK_MS_Results, regex("coli|pneumoniae", ignore_case = TRUE))
 #     ) %>%
 #     mutate(
 #       Confirmed_Species = case_when(
 #         str_detect(VITEK_MS_Results, regex("coli", ignore_case = TRUE)) ~ "E. coli",
 #         str_detect(VITEK_MS_Results, regex("pneumoniae", ignore_case = TRUE)) ~ "K. pneumoniae",
 #         TRUE ~ VITEK_MS_Results
 #       )
 #     )
 #   
 #   # --- Step 2: Calculate Total Confirmed Baseline Count per Species ---
 #   confirmed_summary <- base_cgr3 %>%
 #     filter(CRO_ED30 < 23) %>%
 #     group_by(Confirmed_Species) %>%
 #     summarise(
 #       Confirmed_Isolates = n_distinct(TVLA_ID),
 #       .groups = "drop"
 #     )
 #   
 #   # --- Step 3: Select, Rename, Filter Confirmed Unique Cases ---
 #   Data <- base_cgr3 %>%
 #     select(
 #       TVLA_ID, PROTOCOL, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,
 #       VITEK_MS_Results, Confirmed_Species, ESBL_Presumptivefinal, 
 #       AMX_ED10, AZM_ED15, CRO_ED30, CIP_ED5, DOX_ED30, FLR_ED30, 
 #       GEN_ED10, MEM_ED10, OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30
 #     ) %>%
 #     rename(
 #       Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
 #       Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
 #       Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
 #       Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
 #       Oxytetracycline = "OXY_ED30",
 #       "Polymyxin_B(PB)" = "POL_ED300",
 #       "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
 #       Cefotaxime = "CTX_ED5",
 #       "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
 #     ) %>%
 #     filter(
 #       Ceftriaxone < 23, # Confirmed resistance cutoff
 #       ESBL_Presumptivefinal == 1
 #     ) %>%
 #     distinct(TVLA_ID, .keep_all = TRUE)
 #   
 #   # --- Step 4: Pivot drug columns into long format ---
 #   AMR_Data <- Data %>%
 #     pivot_longer(
 #       cols = -c(TVLA_ID, PROTOCOL, VITEK_MS_Results, Isolate, Confirmed_Species,
 #                 REGION, SEASON, ORIGIN_OF_SAMPLE, ESBL_Presumptivefinal),
 #       names_to = "Antimicrobial_substance",
 #       values_to = "value"
 #     )
 #   
 #   # Align EPI_CUTOFF column name if needed
 #   if ("Antibiotic_tested" %in% names(EPI_CUTOFF)) {
 #     EPI_CUTOFF <- rename(EPI_CUTOFF, Antimicrobial_substance = "Antibiotic_tested")
 #   }
 #   
 #   # --- Step 5: Join Cutoffs & Filter ONLY Agents with Valid Cutoff Values ---
 #   AMR_Cutoff <- inner_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance") %>%
 #     filter(!is.na(S) & !is.na(R)) # Strictly keep agents with defined cutoffs
 #   
 #   # --- Step 6: Categorize Results & Drop Missing/Un-assayed Measurements ---
 #   Confirmed_Susceptibility <- AMR_Cutoff %>%
 #     mutate(
 #       S_Cleaned = str_replace_all(S, "[^0-9.]", ""),
 #       R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
 #       Measured_Zone = as.numeric(str_trim(value)),
 #       S_Cutoff = as.numeric(S_Cleaned),
 #       R_Cutoff = as.numeric(R_Cleaned),
 #       Categorical_Result = case_when(
 #         Measured_Zone >= S_Cutoff ~ "S",
 #         Measured_Zone <= R_Cutoff ~ "R",
 #         Measured_Zone > R_Cutoff & Measured_Zone < S_Cutoff ~ "I",
 #         TRUE ~ NA_character_
 #       )
 #     ) %>%
 #     distinct(TVLA_ID, Antimicrobial_substance, .keep_all = TRUE) %>%
 #     filter(!is.na(Categorical_Result)) # Drops missing test values for that drug
 #   
 #   # --- Step 7: Summary Function Using Dynamic Actual Tested Denominators ---
 #   calculate_resistance_ci_summary <- function(data, group_var, conf_sum) {
 #     
 #     data %>%
 #       group_by(Confirmed_Species, Antimicrobial_substance, Grouping_Value = .data[[group_var]]) %>%
 #       summarise(
 #         Subgroup_Tested = n_distinct(TVLA_ID), # Actual validly tested isolates in this subgroup
 #         N_Resistant = sum(Categorical_Result == "R", na.rm = TRUE),
 #         .groups = "drop"
 #       ) %>%
 #       mutate(Grouping_Variable = group_var) %>%
 #       left_join(conf_sum, by = "Confirmed_Species") %>%
 #       mutate(
 #         # Percentage relative to the TRUE tested count for that specific drug & subgroup
 #         Prop_Raw = ifelse(Subgroup_Tested > 0, (N_Resistant / Subgroup_Tested) * 100, NA_real_),
 #         
 #         # 95% CI based on actual tested denominator (Subgroup_Tested)
 #         Prop_Test = map2(N_Resistant, Subgroup_Tested, function(x, n) {
 #           if (!is.na(n) && n > 0 && x <= n) {
 #             prop.test(x, n)
 #           } else {
 #             NULL
 #           }
 #         }),
 #         CI_Lower_Raw = map_dbl(Prop_Test, function(res) {
 #           if (!is.null(res)) res$conf.int[1] * 100 else NA_real_
 #         }),
 #         CI_Upper_Raw = map_dbl(Prop_Test, function(res) {
 #           if (!is.null(res)) res$conf.int[2] * 100 else NA_real_
 #         }),
 #         
 #         # Clean formatting
 #         Pct_Formatted = fmt_num_vec(Prop_Raw),
 #         CI_Lower_Formatted = fmt_num_vec(CI_Lower_Raw),
 #         CI_Upper_Formatted = fmt_num_vec(CI_Upper_Raw),
 #         
 #         `Resistant (n, %)` = paste0(N_Resistant, " (", Pct_Formatted, "%)"),
 #         `95%_CI` = if_else(!is.na(CI_Lower_Raw), paste0(CI_Lower_Formatted, " - ", CI_Upper_Formatted), "-")
 #       ) %>%
 #       select(
 #         Confirmed_Species,
 #         Confirmed_Isolates,
 #         Antimicrobial_substance,
 #         Grouping_Variable,
 #         Grouping_Value,
 #         Subgroup_Tested,
 #         `Resistant (n, %)`,
 #         `95%_CI`
 #       )
 #   }
 #   
 #   # --- Step 8: Map across demographic variables ---
 #   group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
 #   
 #   map_dfr(
 #     group_vars_list,
 #     ~ calculate_resistance_ci_summary(Confirmed_Susceptibility, .x, confirmed_summary)
 #   ) %>%
 #     arrange(Confirmed_Species, Antimicrobial_substance, Grouping_Variable, Grouping_Value)
 # }
 # 
 # # ==============================================================================
 # # 2. RUN PIPELINE & BUILD PUBLICATION TABLE
 # # ==============================================================================
 # 
 # # Execute pipeline
 # summary_df <- process_amr_summary(Final_JoinedDATA, EPI_CUTOFF)
 # 
 # # Format flextable (APA / ICMJE standard)
 # ft <- summary_df %>%
 #   flextable() %>%
 #   set_header_labels(
 #     Confirmed_Species = "Confirmed Species",
 #     Confirmed_Isolates = "Confirmed Total (N)",
 #     Antimicrobial_substance = "Antimicrobial Agent",
 #     Grouping_Variable = "Stratifier",
 #     Grouping_Value = "Subgroup",
 #     Subgroup_Tested = "Tested (n)",
 #     `Resistant (n, %)` = "Resistant (n, %)",
 #     `95%_CI` = "95% CI"
 #   ) %>%
 #   font(fontname = "Times New Roman", part = "all") %>%
 #   fontsize(size = 9, part = "all") %>%
 #   bold(part = "header") %>%
 #   align(j = c("Confirmed_Isolates", "Subgroup_Tested", "Resistant (n, %)", "95%_CI"), align = "center", part = "all") %>%
 #   align(j = c("Confirmed_Species", "Antimicrobial_substance", "Grouping_Variable", "Grouping_Value"), align = "left", part = "all") %>%
 #   border_remove() %>%
 #   hline_top(border = fp_border(color = "black", width = 1.5), part = "header") %>%
 #   hline_bottom(border = fp_border(color = "black", width = 1.0), part = "header") %>%
 #   hline_bottom(border = fp_border(color = "black", width = 1.5), part = "body") %>%
 #   padding(padding = 3, part = "all") %>%
 #   merge_v(j = c("Confirmed_Species", "Confirmed_Isolates", "Antimicrobial_substance", "Grouping_Variable", "Grouping_Value")) %>%
 #   autofit()
 # 
 # # Build Word document structure
 # doc <- read_docx() %>%
 #   body_add_par(
 #     "Table 1. Antimicrobial Resistance Profiles and 95% Confidence Intervals Stratified by Demographic Factors",
 #     style = "Normal"
 #   ) %>%
 #   body_add_flextable(ft) %>%
 #   body_add_par(
 #     "Note: Subgroup Tested (n) represents isolates with valid test results for each specific antimicrobial agent. Percentages and 95% Confidence Intervals (95% CI) are calculated using the subgroup tested count (n) as the denominator via Pearson's chi-squared test (prop.test).",
 #     style = "Normal"
 #   )
 # 
 # # ==============================================================================
 # # 3. DIRECT FILE SAVING TO "Results" FOLDER
 # # ==============================================================================
 # 
 # if (!dir.exists("Results")) {
 #   dir.create("Results")
 # }
 # 
 # print(doc, target = "Results/ESBL_Antimicrobial_Resistance_Summary.docx")
 # 
 # cat("\nSUCCESS: Publication-formatted table saved to 'Results/ESBL_Antimicrobial_Resistance_Summary.docx'!\n")