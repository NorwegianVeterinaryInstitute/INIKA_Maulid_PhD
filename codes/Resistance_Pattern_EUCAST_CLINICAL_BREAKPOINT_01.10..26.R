#SETUP

#IMPORT DATA

#RENAME VARIABLES

#IMPORT ECOFF FILES

#E. COLI RESISTANCE

#K. PNEUMONIAE RESISTANCE

#ZONE DIAMETER DISTRIBUTION

#FIGURE GENERATION

#EXPORT RESULTS


##############################################################################
# SETUP
##############################################################################

library(tidyverse)
library(readxl)
library(writexl)
library(flextable)
library(officer)
library(ggplot2)

##############################################################################
# MASTER DATA
##############################################################################

joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")

spec(joined_data)
names(joined_data)[grep("INIKA", names(joined_data))]

joined_data <- joined_data %>%
  rename(
    INIKA_ID = INIKA_ID.x.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    ESBL = ESBL_Presumptivefinal
  )%>%
  select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,PROTOCOL,
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
           ) 

names(joined_data)



##############################################################################
# ECOFF TABLES
##############################################################################
# Imorting E.coli ECOFF
Ecoli_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx"
)

Ecoli_EPI_CUTOFF <- Ecoli_ECOFF_EUCAST_BREAK_POINT

# Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
#   "data/ECOFF_E.coli_K.pneumoniae.xlsx",
#   sheet = 2
# )
# 
# Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT
# 
# 
# ECOFF_EUCAST_BREAK_POINT <- read_excel(
#   "data/ECOFF_E.coli_K.pneumoniae.xlsx"
# )
# 
# EPI_CUTOFF <- ECOFF_EUCAST_BREAK_POINT

#Import ECOFF table with the sheet for Klebsiella
Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx",
  sheet = "K.pneumoniae_ECOFF_EUCAST"
)

Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT


##############################################################################
# Process E.coli Data
##############################################################################
Data<-joined_data

AMR_Data <- Data %>%
  pivot_longer(
    cols = -c(
      INIKA_ID,
      VITEK_MS_Results,
      PROTOCOL,
      Isolate,
      REGION,
      SEASON,
      ORIGIN_OF_SAMPLE,
      ESBL
    ),
    names_to = "Antimicrobial_substance",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

unique(Ecoli_EPI_CUTOFF$Antimicrobial_substance)
names(Data)
## The commented bellow seems to be a repetition of the above
# AMR_Data <- Data %>%
#   pivot_longer(
#     cols = -c(
#       INIKA_ID,
#       VITEK_MS_Results,
#       PROTOCOL,
#       Isolate,
#       REGION,
#       SEASON,
#       ORIGIN_OF_SAMPLE,
#       ESBL
#     ),
#     names_to = "Antimicrobial_substance",
#     values_to = "value"
#   ) %>%
#   filter(!is.na(value))
  AMR_Data <- AMR_Data %>%
  mutate(
    Antimicrobial_substance = recode(
      Antimicrobial_substance,
      "Popymyxin_B(PB)" = "Polymyxin B",
      "Sulfamethoxazole/Trimethoprim" = "Sulfamethoxazole_Trimethoprim",
      "Cefotaxime/ClavulanicAcid" = "Cefotaxime _Clavulanic Acid"
    )
  )

unique(Ecoli_EPI_CUTOFF$Antimicrobial_substance)
names(Data)

AMR_Cutoff <- AMR_Data %>%
  left_join(Ecoli_EPI_CUTOFF, by = "Antimicrobial_substance")

AMR_Cutoff %>%
  summarise(
    Missing_S = sum(is.na(S)),
    Missing_R = sum(is.na(R))
  )

# Checking out what is the true number of ESCR E.coli!
Test_N_ESCR<-AMR_Cutoff%>%
  filter(PROTOCOL=="CGR3",
         VITEK_MS_Results=="Escherichia coli",
         Antimicrobial_substance=="Ceftriaxone",
         value <23)

# There are 147 isolates fulfilling these criteria!
Test_N_ESCR %>%
  summarise(
    N_rows = n(),
    N_unique_isolates = n_distinct(INIKA_ID)
  )

Test_N_ESCR %>%
  count(INIKA_ID, sort = TRUE) %>%
  filter(n > 1)

CGR3_Data <- joined_data %>%
  filter(PROTOCOL == "CGR3")

joined_data %>%
  count(PROTOCOL, sort = TRUE)

CGR3_Data %>%
  summarise(
    Total_CGR3 = n_distinct(INIKA_ID)
  )

Total_Ecoli_CGR3 <- AMR_Cutoff %>%
  filter(
    PROTOCOL == "CGR3",
    VITEK_MS_Results == "Escherichia coli",
    Antimicrobial_substance == "Ceftriaxone"
  ) %>%
  distinct(INIKA_ID) # 151 observations

ESCR_ECO_confirmed <- AMR_Cutoff %>%
  filter(
    PROTOCOL == "CGR3",
    VITEK_MS_Results == "Escherichia coli",
    Antimicrobial_substance == "Ceftriaxone",
    value < 23
  )%>%
  distinct(INIKA_ID) # 146 Observations


ESCR_ECO_confirmed <- joined_data %>%
  filter(
    PROTOCOL == "CGR3",
    VITEK_MS_Results == "Escherichia coli",
    Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE) # 146 Observations


 
names(ESCR_ECO_confirmed)

tested_antibiotics <- names(ESCR_ECO_confirmed)[
  names(ESCR_ECO_confirmed) %in% unique(Ecoli_EPI_CUTOFF$Antimicrobial_substance)
]

tested_antibiotics

ESCR_long <- ESCR_ECO_confirmed %>%
  pivot_longer(
    cols = all_of(tested_antibiotics),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm))

names(joined_data)

ESCR_ECO_confirmed <- joined_data %>%
  filter(
    PROTOCOL == "CGR3",
    VITEK_MS_Results == "Escherichia coli",
    Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE)


N_Total <- nrow(Total_Ecoli_CGR3)
#N_ESCR <- nrow(ESCR_Ecoli) # Not found!!
N_ESCR <- nrow(ESCR_ECO_confirmed)


tested_antibiotics <- c(
  "Amoxicillin",
  "Azithromycin",
  "Ceftriaxone",
  "Ciprofloxacin",
  "Doxycycline",
  "Florfenicol",
  "Gentamicin",
  "Meropenem",
  "Oxytetracycline",
  "Popymyxin_B(PB)",
  "Sulfamethoxazole/Trimethoprim",
  "Cefotaxime",
  "Cefotaxime/ClavulanicAcid"
)

# From Maulid
drug_class_map <- tibble::tribble(
  ~Antimicrobial_substance,         ~Antibiotic_Class,
  "Amoxicillin",                     "Penicillins",
  "Ceftriaxone",                     "3rd Gen Cephalosporins",
  "Cefotaxime",                      "3rd Gen Cephalosporins",
  "Cefotaxime/ClavulanicAcid",       "Beta-lactam / Inhibitor combinations",
  "Meropenem",                       "Carbapenems",
  "Azithromycin",                    "Macrolides",
  "Ciprofloxacin",                   "Fluoroquinolones",
  "Gentamicin",                      "Aminoglycosides",
  "Doxycycline",                     "Tetracyclines",
  "Oxytetracycline",                 "Tetracyclines",
  "Florfenicol",                     "Phenicols",
  "Popymyxin_B(PB)",                 "Polymyxins",
  "Sulfamethoxazole/Trimethoprim",   "Folate Pathway Inhibitors"
)

ESCR_long <- ESCR_ECO_confirmed %>%
  pivot_longer(
    cols = any_of(tested_antibiotics),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm))

 


ESCR_long <- ESCR_long %>%
  mutate(
    R_limit =
      as.numeric(str_replace_all("Antimicrobial_substance", "[^0-9.]", "")),
    is_resistant =
      if_else(mm <= R_limit, 1, 0)
  ) # R has replaced with "Antimicrobial_substance" since was not found





intersect(
  tested_antibiotics,
  names(ESCR_ECO_confirmed)
)


names(ESCR_ECO_confirmed)
  
# Filter for E. coli criteria (either Isolate or VITEK_MS_Results) AND ESBL == 1
# Note: K.pneumoniae condition in Isolate is added based on the prompt's condition one.
joined_data <- joined_data %>%
  filter((Isolate %in% c("E.coli", "K.pneumoniae") | VITEK_MS_Results == "Escherichia coli") &
                                       ESBL == 1) %>%
  
  # Ensure only unique INIKA_ID are considered for the count (denominator)
  distinct(INIKA_ID, .keep_all = TRUE) # 50 Observations

Data <- joined_data %>%
  filter(
    (Isolate %in% c("E.coli", "K.pneumoniae") |
       VITEK_MS_Results == "Escherichia coli") &
      ESBL == 1
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE) # 50 Observations

# --- Step 2: Pivot the drug columns into long format ---
AMR_Data <- AMR_Data %>%
  pivot_longer(
  cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
            ORIGIN_OF_SAMPLE, ESBL),
  names_to = "Antimicrobial_substance",
  values_to = "value",
  values_transform = list(value = as.character)
)

# --- Step 3: Join with the cutoff table ---
AMR_Cutoff <- left_join(AMR_Data, Ecoli_EPI_CUTOFF, by = "Antimicrobial_substance")

# Define the list of grouping variables for the final summary
group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")

# --- Step 4: Clean and Classify Susceptibility (S/I/R) ---
Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
  
  # Clean the S and R cutoff columns from non-numeric characters (e.g., '<', '>')
  # and convert all key values to numeric
  mutate(
    S_Cutoff      = as.numeric(na_if(str_replace_all(S, "[^0-9.]", ""), "")),
    R_Cutoff      = as.numeric(na_if(str_replace_all(R, "[^0-9.]", ""), "")),
    Measured_Zone = as.numeric(na_if(str_replace_all(value, "[^0-9.]", ""), "")) ) %>%
  
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
}
##############################################
# --- 1. Define the Processing Function ---
# ==============================================================================
# 1. MAIN DATA PROCESSING FUNCTION
# ==============================================================================
process_ecoli_data <- function(joined_data, Ecoli_EPI_CUTOFF) {
  
  # Step 1: Filter raw data and assign to 'Data'
  Data <- joined_data %>%
    filter(
      Isolate %in% c("E.coli", "K.pneumoniae"),
      VITEK_MS_Results == "Escherichia coli",
      Ceftriaxone < 23 
    ) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
  
  # Step 2: Pivot to long format
  AMR_Data <- pivot_longer(
    Data,
    cols = -c(INIKA_ID, VITEK_MS_Results, Isolate, REGION, SEASON,
              ORIGIN_OF_SAMPLE, ESBL),
    names_to = "Antimicrobial_substance",
    values_to = "value",
    values_transform = list(value = as.character)
  )
  
  # Step 3: Join with Cutoff
  AMR_Cutoff <- left_join(AMR_Data, Ecoli_EPI_CUTOFF, by = "Antimicrobial_substance")
  
  # Step 4: Classify Susceptibility
  Ecoli_ESBL_Susceptibility <- AMR_Cutoff %>%
    mutate(
      S_Cleaned     = as.numeric(na_if(str_replace_all(S, "[^0-9.]", ""), "")),
      R_Cleaned     = as.numeric(na_if(str_replace_all(R, "[^0-9.]", ""), "")),
      Measured_Zone = as.numeric(na_if(str_replace_all(value, "[^0-9.]", ""), ""))
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
  
  # Step 5: Summary Logic & Statistical Tests
  
  # Helper to compute Fisher's Exact Test p-value across categories for a given drug
  calc_stratum_pvalue <- function(df, group_var, drug) {
    sub_df <- df %>% filter(Antimicrobial_substance == drug)
    tbl <- table(sub_df[[group_var]], sub_df$Categorical_Result == "R")
    
    # Use Fisher's exact test across categories
    if (nrow(tbl) > 1 && ncol(tbl) > 0) {
      tryCatch({
        p_val <- fisher.test(tbl)$p.value
        if (p_val < 0.001) "< 0.001" else sub("\\.0$", "", as.character(round(p_val, 3)))
      }, error = function(e) "NA")
    } else {
      "NA"
    }
  }
  
  calculate_summary <- function(df, group_var = NULL) {
    df_grouped <- if (is.null(group_var)) {
      df %>% group_by(Antimicrobial_substance)
    } else {
      df %>% group_by(!!sym(group_var), Antimicrobial_substance)
    }
    
    res <- df_grouped %>%
      summarise(
        Total_Tests = n(),
        N_Resistant = sum(Categorical_Result == "R"),
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        Prop_Test    = list(binom.test(N_Resistant, Total_Tests)),
        Pct_Raw      = round((N_Resistant / Total_Tests) * 100, 1),
        CI_Lower_Raw = round(Prop_Test$conf.int[1] * 100, 1),
        CI_Upper_Raw = round(Prop_Test$conf.int[2] * 100, 1),
        
        Percentage   = sub("\\.0$", "", as.character(Pct_Raw)),
        CI_Lower     = sub("\\.0$", "", as.character(CI_Lower_Raw)),
        CI_Upper     = sub("\\.0$", "", as.character(CI_Upper_Raw)),
        
        `95%_CI`     = paste0("(", CI_Lower, " - ", CI_Upper, ")"),
        Grouping_Variable = if (is.null(group_var)) "OVERALL" else group_var,
        Grouping_Value    = if (is.null(group_var)) "Overall" else as.character(!!sym(group_var))
      ) %>%
      ungroup()
    
    # Add p-values for stratified groups (Overall does not have a p-value comparison)
    if (!is.null(group_var)) {
      p_vals <- sapply(unique(res$Antimicrobial_substance), function(d) {
        calc_stratum_pvalue(df, group_var, d)
      })
      res <- res %>% mutate(p_value = p_vals[Antimicrobial_substance])
    } else {
      res <- res %>% mutate(p_value = "—")
    }
    
    res %>% select(
      Grouping_Variable, Grouping_Value, Antimicrobial_substance, 
      Total_Tests, N_Resistant, Percentage, `95%_CI`, p_value
    )
  }
  
  group_vars_list <- c("REGION", "SEASON", "ORIGIN_OF_SAMPLE")
  
  Stratified_Summary <- map_dfr(group_vars_list, ~calculate_summary(Ecoli_ESBL_Susceptibility, .x))
  Overall_Summary    <- calculate_summary(Ecoli_ESBL_Susceptibility, group_var = NULL)
  
  Final_Summary <- bind_rows(Overall_Summary, Stratified_Summary)
  return(Final_Summary)
}

# ==============================================================================
# 2. WORD PUBLICATION EXPORT FUNCTION
# ==============================================================================
format_publication_word <- function(data, output_filepath) {
  
  clean_data <- data %>%
    # 1. Recode Variable names and Value categories
    mutate(
      Grouping_Variable = case_when(
        Grouping_Variable == "OVERALL"          ~ "Overall",
        Grouping_Variable == "REGION"           ~ "Region",
        Grouping_Variable == "SEASON"           ~ "Season",
        Grouping_Variable == "ORIGIN_OF_SAMPLE" ~ "Origin of sample",
        TRUE ~ Grouping_Variable
      ),
      Grouping_Value = case_when(
        Grouping_Value == "KILI" ~ "Kilimanjaro",
        Grouping_Value == "MWA"  ~ "Mwanza",
        Grouping_Value == "d"    ~ "Dry",
        Grouping_Value == "w"    ~ "Wet/Rainy",
        TRUE ~ Grouping_Value
      )
    ) %>%
    # 2. Select columns including p-value
    select(
      `Variable`       = Grouping_Variable,
      `Value`          = Grouping_Value,
      `Antimicrobial`  = Antimicrobial_substance,
      `Tested (N)`     = Total_Tests,
      `Resistant (n)`  = N_Resistant,
      `Resistance (%)` = Percentage,
      `95% CI`         = `95%_CI`,
      `p-value`        = p_value
    ) %>%
    # 3. Deduplicate repeating Variable, Value, AND p-value labels within groups
    mutate(
      p_value_dup = if_else(Variable == lag(Variable, default = "") & 
                              Antimicrobial == lag(Antimicrobial, default = ""), TRUE, FALSE),
      `p-value`   = if_else(p_value_dup, "", `p-value`),
      Variable    = if_else(row_number() == 1 | Variable != lag(Variable, default = ""), Variable, ""),
      Value       = if_else(row_number() == 1 | Value != lag(Value, default = ""), Value, "")
    ) %>%
    select(-p_value_dup)
  
  # 4. Build Flextable
  ft <- flextable(clean_data) %>%
    theme_vanilla() %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    align(align = "left", part = "all") %>%
    align(j = c("Tested (N)", "Resistant (n)", "Resistance (%)", "95% CI", "p-value"), align = "center", part = "all") %>%
    
    # APA Standard 3-line borders
    border_remove() %>%
    hline_top(border = fp_border(color = "black", width = 1.5), part = "header") %>%
    hline_bottom(border = fp_border(color = "black", width = 1.0), part = "header") %>%
    hline_bottom(border = fp_border(color = "black", width = 1.5), part = "body") %>%
    padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
    autofit()
  
  # 5. Write Word Document
  doc <- read_docx() %>%
    body_add_par("Table 1. Antimicrobial Resistance Prevalence across Demographics", style = "heading 1") %>%
    body_add_par("") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par("Note. CI = Confidence Interval calculated via Clopper-Pearson exact binomial test. p-values calculated using Fisher's exact test across demographic categories.", style = "Normal")
  
  print(doc, target = output_filepath)
  message("Saved Publication Word Table with p-values to: ", output_filepath)
}

# ==============================================================================
# 3. EXECUTION & EXPORT
# ==============================================================================

# Run data pipeline
result_table <- process_ecoli_data(joined_data, Ecoli_EPI_CUTOFF)

# Set output directory
output_dir <- "Results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export CSV
write_csv(result_table, file.path(output_dir, "28_01_26Ecoli_ESBL_Resistance_Summary.csv"))

# Export Excel
write_xlsx(result_table, file.path(output_dir, "19_12_25Ecoli_ESBL_Resistance_Summary.xlsx"))

# Export APA Publication Word Table
word_path <- file.path(output_dir, "AMR_Summary_Table_Publication.docx")
format_publication_word(result_table, output_filepath = word_path)
################################################################################

#Supplementary Table 3
## Computing the AMR frequency for ESCR E.coli and Distribution of Zone diameter in millimeters (mm)

#  PREPARATION 
original_spelling <- c("AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
                       "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
                       "SXT_ED1_2", "CTX_ED5", "CTC_ED30") 

# Note: Removed leading space from "Azithromycin" to ensure matching with EPI_CUTOFF
tested_antibiotics <- c("Amoxicillin", "Azithromycin", "Ceftriaxone", "Ciprofloxacin", 
                        "Doxycycline", "Florfenicol", "Gentamicin", "Meropenem", 
                        "Oxytetracycline", "Polymyxin_B(PB)", "Sulfamethoxazole/Trimethoprim", 
                        "Cefotaxime", "Cefotaxime/ClavulanicAcid")

#  DATA PROCESSING 
ESCR_confirmed_data <- joined_data %>%
  rename(any_of(setNames(original_spelling, tested_antibiotics))) %>%
  filter(
    (
      # Group 1: Standard E. coli
      grepl("E", Isolate, ignore.case = TRUE) & 
        grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & 
        grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
    ) | 
      (
        # Group 2: Metallic blue isolates identified as E. coli by VITEK
        grepl("K", Isolate, ignore.case = TRUE) & 
          grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & 
          grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
      ) & Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

message(paste("Isolates processed:", nrow(ESCR_confirmed_data)))

# CALCULATIONS USING EPI_CUTOFF 
ESCR_long <- ESCR_confirmed_data %>% 
  pivot_longer(cols = any_of(tested_antibiotics), 
               names_to = "Antimicrobial_substance", 
               values_to = "mm") %>%
  filter(!is.na(mm)) %>%
  # Join with your cutoff table (Ensure EPI_CUTOFF has 'Antimicrobial_substance' and 'R' columns)
  left_join(EPI_CUTOFF, by = "Antimicrobial_substance") %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    # Extract numeric value from 'R' column (e.g., handles "<=13" or "13")
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    # Resistant if zone is less than or equal to the R limit
    is_resistant = ifelse(mm <= R_limit, 1, 0)
  )

Summary_Resistance <- ESCR_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    R_count = sum(is_resistant, na.rm = TRUE),
    Prop = R_count / N,
    # Binomial test for 95% CI
    Lower = binom.test(R_count, N)$conf.int[1] * 100,
    Upper = binom.test(R_count, N)$conf.int[2] * 100,
    `Resistance % [95% CI]` = sprintf("%.1f%% [%.1f-%.1f]", Prop * 100, Lower, Upper),
    .groups = "drop"
  )

# Frequency distribution of mm values
Percentage_Wide <- ESCR_long %>%
  count(Antimicrobial_substance, mm) %>%
  group_by(Antimicrobial_substance) %>%
  mutate(Percentage = (n / sum(n)) * 100) %>%
  ungroup() %>%
  pivot_wider(id_cols = Antimicrobial_substance, names_from = mm, 
              values_from = Percentage, values_fill = 0)

# Sort mm columns numerically
mm_cols <- as.character(sort(as.numeric(names(Percentage_Wide)[-1])))

# Combine Resistance Summary with mm distribution
Final_Formatted <- Summary_Resistance %>% 
  select(Antimicrobial_substance, N, `Resistance % [95% CI]`) %>%
  left_join(Percentage_Wide, by = "Antimicrobial_substance") %>%
  select(Antimicrobial_substance, N, `Resistance % [95% CI]`, all_of(mm_cols)) %>%
  mutate(across(all_of(mm_cols), function(x) {
    res <- round(x, 1)
    txt <- format(res, nsmall = 0)
    ifelse(res == 0, "-", as.character(res))
  }))

# CREATE & STYLE TABLE 
# We use 'autofit' to let the table expand to its content
ESCR_ECO_AMR <- Final_Formatted %>%
  flextable() %>%
  set_table_properties(layout = "autofit", width = 1) %>% 
  
  theme_booktabs() %>%
  set_caption(caption = "Table: Resistance Percentages and ZOI Frequency Distribution for ESCR E.coli Isolates (N=149).") %>%
  
  # Font size 6.5 is the "sweet spot" for many columns on one landscape page
  fontsize(size = 6.5, part = "all") %>%
  bold(part = "header") %>%
  
  # Set padding to zero to save every millimeter of horizontal space
  padding(padding.top = 1, padding.bottom = 1, 
          padding.left = 0, padding.right = 0, part = "all") %>%
  
  # Align text
  align(align = "center", j = 2:ncol(Final_Formatted), part = "all") %>%
  bg(j = 3, bg = "#F2F2F2", part = "body") %>%
  
  # Important: Fix the header rotation if names are too long
  valign(valign = "bottom", part = "header")

#  EXPORT 
if(!dir.exists("Results")) dir.create("Results")

# Create the document using the 'officer' workflow
doc <- read_docx() %>%
  # Use body_add_par to add a little space before the table if needed
  body_add_flextable(value = ESCR_ECO_AMR) %>%
  body_end_section_landscape()

print(doc, target = "Results/ESCR_ECO_AMR_Frequency_Distribution_Table.docx")
save_as_docx( ESCR_ECO_AMR, path = "Results/ESCR_ECO_AMR_Frequency_Distribution_Table.docx")
################################################################################  
# Supplementary Table 4
## Computing the ESCR K.pneumoniae frequency Distribution of Zone diameter in millimeters (mm)

#  1. DATA PREP 
ESCR_Kpn_confirmed <- joined_data %>%
  rename(any_of(setNames(original_spelling, tested_antibiotics))) %>%
  filter(
    # 1. Species Confirmation via VITEK
    trimws(VITEK_MS_Results) == "Klebsiella pneumoniae",
    
    # 2. Specific Isolate + Morphology logic
    (
      (grepl("K", Isolate, ignore.case = TRUE) & 
         `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") |
        (grepl("E", Isolate, ignore.case = TRUE) & 
           `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish"))
    ) & Ceftriaxone < 23
  ) %>%
  arrange(INIKA_ID, Isolate) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

message(paste("K. pneumoniae Isolates processed:", nrow(ESCR_Kpn_confirmed)))

# 2. CALCULATIONS (Resistance Only) 
ESCR_Kpn_long <- ESCR_Kpn_confirmed %>%
  pivot_longer(cols = any_of(tested_antibiotics), 
               names_to = "Antimicrobial_substance", 
               values_to = "mm") %>%
  filter(!is.na(mm)) %>%
  # Join with Kpn-specific Cutoffs
  left_join(Kpn_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    # Resistant if zone is <= R limit
    is_resistant = ifelse(mm <= R_limit, 1, 0)
  )

Summary_Resistance <- ESCR_Kpn_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    R_count = sum(is_resistant, na.rm = TRUE),
    Prop = (R_count / N) * 100,
    # Binomial exact CI
    Lower = binom.test(R_count, N)$conf.int[1] * 100,
    Upper = binom.test(R_count, N)$conf.int[2] * 100,
    `Resistance % [95% CI]` = sprintf("%.1f%% [%.1f-%.1f]", Prop, Lower, Upper),
    .groups = "drop"
  )

# Frequency distribution (Percentage)
Percentage_Wide <- ESCR_Kpn_long %>%
  count(Antimicrobial_substance, mm) %>%
  group_by(Antimicrobial_substance) %>%
  mutate(Percentage = (n / sum(n)) * 100) %>%
  ungroup() %>%
  pivot_wider(id_cols = Antimicrobial_substance, names_from = mm, 
              values_from = Percentage, values_fill = 0)

mm_cols <- as.character(sort(as.numeric(names(Percentage_Wide)[-1])))

# Combine Resistance Summary with mm distribution
Final_Kpn_Formatted <- Summary_Resistance %>% 
  select(Antimicrobial_substance, N, `Resistance % [95% CI]`) %>%
  left_join(Percentage_Wide, by = "Antimicrobial_substance") %>%
  select(Antimicrobial_substance, N, `Resistance % [95% CI]`, all_of(mm_cols)) %>%
  mutate(across(all_of(mm_cols), function(x) {
    res <- round(x, 1)
    # Clean zeros to "-" for readability
    ifelse(res == 0, "-", as.character(res))
  }))

#  3. CREATE COMPACT TABLE 
final_ESCR_Kpn_table <- Final_Kpn_Formatted %>%
  flextable() %>%
  set_caption(caption = paste0("Table: Resistance Percentages and ZOI Frequency Distribution for ESCR K. pneumoniae (N=58", nrow(ESCR_Kpn_confirmed), ").")) %>%
  theme_vanilla() %>%
  fontsize(size = 7, part = "all") %>%
  padding(padding.top = 1, padding.bottom = 1, padding.left = 0, padding.right = 0, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", j = 2:ncol(Final_Kpn_Formatted), part = "all") %>%
  set_table_properties(layout = "autofit", width = 1) %>%
  bg(j = 3, bg = "#F9F9F9", part = "body") # Light highlight on Resistance column

#  4. EXPORT 
if(!dir.exists("Results")) dir.create("Results")

doc <- read_docx() %>%
  body_add_flextable(value = final_ESCR_Kpn_table) %>%
  body_end_section_landscape()

print(doc, target = "Results/Kpn_Resistance_Frequency_Table.docx")
save_as_docx(final_ESCR_Kpn_table, path = "Results/final_ESCR_Kpn_table.docx")

################################################################################   
# Supplementary Table 5
#Computing AMR for ESCR E.coli
process_ecoli_amr_analysis <- function(joined_data, EPI_CUTOFF) {
  
  # Internal Reference Data 
  abx_classes <- tribble(
    ~Antimicrobial_substance,         ~Class,
    "Amoxicillin",                    "Penicillins",
    "Azithromycin",                   "Macrolides",
    "Ceftriaxone",                    "Cephalosporins (3rd Gen)",
    "Cefotaxime",                     "Cephalosporins (3rd Gen)",
    "Cefotaxime/ClavulanicAcid",      "β-lactam/Inhibitor",
    "Ciprofloxacin",                  "Fluoroquinolones",
    "Doxycycline",                    "Tetracyclines",
    "Oxytetracycline",                "Tetracyclines",
    "Gentamicin",                     "Aminoglycosides",
    "Meropenem",                      "Carbapenems",
    "Florfenicol",                    "Phenicols",
    "Polymyxin_B(PB)",                "Polymyxins",
    "Sulfamethoxazole/Trimethoprim",  "Sulfonamides"
  )
  
  rename_map <- c(
    Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
    Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
    Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
    Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
    Oxytetracycline = "OXY_ED30", "Polymyxin_B(PB)" = "POL_ED300",
    "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
    Cefotaxime = "CTX_ED5", "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
  )
  
  # Data Cleaning & Filtering 
  # Clean the column names of the input data first to remove hidden spaces
  colnames(joined_data) <- trimws(colnames(joined_data))
  
  # Check if the columns exist before trying to select them
  missing_cols <- setdiff(rename_map, colnames(joined_data))
  
  if(length(missing_cols) > 0) {
    stop(paste("The following columns are missing from joined_data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  cleaned_data <- joined_data |>
    select(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE, 
           `COLONY MORPHOLOGY ON C3GR`, Isolate, VITEK_MS_Results, 
           ESCR_KPN_presumptive, any_of(rename_map)) |>
    rename(any_of(rename_map)) |>
    filter(
      # 1. Species Confirmation via VITEK (Primary Filter)
      trimws(VITEK_MS_Results) == "Escherichia coli",
      
      # 2. Specific Isolate + Morphology logic
      (
        # Group A: E. coli that are Pinkish/Reddish (but VITEK says Kpn)
        (grepl("E", Isolate, ignore.case = TRUE) & 
           `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish") |
          
          # Group B: K. pneumoniae that are Metallic blue
           (grepl("K", Isolate, ignore.case = TRUE) & 
              `COLONY MORPHOLOGY ON C3GR` == "Metallic blue"))
      ),
      
      # 3. Third-generation cephalosporin susceptibility breakpoint (< 23 mm)
      Ceftriaxone < 23
    ) |>
    arrange(INIKA_ID, Isolate) |>
    # Ensure unique participants are counted once
    distinct(INIKA_ID, .keep_all = TRUE)
  
  message(paste("Isolates processed:", nrow(cleaned_data)))
  
  # Pivot & Susceptibility Classification 
  kpn_long <- cleaned_data |>
    pivot_longer(
      cols = any_of(names(rename_map)),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) |>
    left_join(EPI_CUTOFF, by = "Antimicrobial_substance") |>
    mutate(
      across(c(S, R), ~as.numeric(str_remove_all(.x, "[^0-9.]")), .names = "{.col}_num"),
      Measured_Zone = as.numeric(Measured_Zone),
      Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(Result))
  
  # Summary Helper Function 
  get_summary <- function(df, var_name) {
    group_cols <- if(var_name == "Overall") "Antimicrobial_substance" else c(var_name, "Antimicrobial_substance")
    
    df |>
      group_by(across(all_of(group_cols))) |>
      summarise(
        n_res = sum(Result == "R"),
        total = n(),
        .groups = "drop"
      ) |>
      mutate(
        Grouping_Var = var_name,
        Grouping_Val = if(var_name == "Overall") "Overall" else as.character(.data[[var_name]])
      )
  }
  
  # Generate and Format Final Table
  group_vars <- c("Overall", "REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE")
  
  final_table <- map_dfr(group_vars, ~get_summary(ecoli_long, .x)) |>
    rowwise() |>
    mutate(
      perc = round((n_res / total) * 100, 1),
      ci = if(Grouping_Var == "Overall") {
        bt <- binom.test(n_res, total)
        paste0(" (", round(bt$conf.int[1]*100, 1), "–", round(bt$conf.int[2]*100, 1), ")")
      } else "",
      cell_text = paste0(perc, "% (", n_res, "/", total, ")", ci)
    ) |>
    ungroup() |>
    select(Antimicrobial_substance, Grouping_Val, cell_text) |>
    pivot_wider(names_from = Grouping_Val, values_from = cell_text) |>
    left_join(abx_classes, by = "Antimicrobial_substance") |>
    select(Class, Antimicrobial_substance, Overall, everything()) |>
    arrange(Class, Antimicrobial_substance)
  
  # Export to Word (Landscape & Fitted) 
  if(!dir.exists("Results")) dir.create("Results")
  
  # Create the flextable
  ft <- flextable(final_table) |>
    theme_booktabs() |>
    fontsize(size = 8, part = "all") |>        # Small font to fit many columns
    padding(padding = 1, part = "all") |>      # Tighten cell space
    merge_v(j = ~ Class) |>                    # Group by Class vertically
    bold(part = "header") |>
    set_table_properties(layout = "autofit")
  
  # Define Landscape & Narrow Margins
  sect_props <- prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5),
    type = "continuous"
  )
  
  # Save as Word
  save_as_docx(ft, path = "Results/ecoli_AMR_Table_Landscape.docx", pr_section = sect_props)
  
  # Also save Excel as backup
  write_xlsx(final_table, "Results/ecoli_AMR_Table.xlsx")
  
  message("Success: Table saved to Results folder as Word (Landscape) and Excel.")
  return(final_table)
}

#  Execute 
ecoli_amr_final_results <- process_ecoli_amr_analysis(joined_data, EPI_CUTOFF)
#######################################################
# Supplementary Table 6
# Computing AMR for ESCR K.pneumoniae
process_kpn_amr_analysis <- function(joined_data, Kpn_EPI_CUTOFF) {
  
  # Internal Reference Data 
  abx_classes <- tribble(
    ~Antimicrobial_substance,         ~Class,
    "Amoxicillin",                    "Penicillins",
    "Azithromycin",                   "Macrolides",
    "Ceftriaxone",                    "Cephalosporins (3rd Gen)",
    "Cefotaxime",                     "Cephalosporins (3rd Gen)",
    "Cefotaxime/ClavulanicAcid",      "β-lactam/Inhibitor",
    "Ciprofloxacin",                  "Fluoroquinolones",
    "Doxycycline",                    "Tetracyclines",
    "Oxytetracycline",                "Tetracyclines",
    "Gentamicin",                     "Aminoglycosides",
    "Meropenem",                      "Carbapenems",
    "Florfenicol",                    "Phenicols",
    "Polymyxin_B(PB)",                "Polymyxins",
    "Sulfamethoxazole/Trimethoprim",  "Sulfonamides"
  )
  
  rename_map <- c(
    Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
    Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
    Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
    Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
    Oxytetracycline = "OXY_ED30", "Polymyxin_B(PB)" = "POL_ED300",
    "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
    Cefotaxime = "CTX_ED5", "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
  )
  
  # Data Cleaning & Filtering 
  # Clean the column names of the input data first to remove hidden spaces
  colnames(joined_data) <- trimws(colnames(joined_data))
  
  # Check if the columns exist before trying to select them
  missing_cols <- setdiff(rename_map, colnames(joined_data))
  
  if(length(missing_cols) > 0) {
    stop(paste("The following columns are missing from joined_data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  cleaned_data <- joined_data |>
    select(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE, 
           `COLONY MORPHOLOGY ON C3GR`, Isolate, VITEK_MS_Results, 
           ESCR_KPN_presumptive, any_of(rename_map)) |>
    rename(any_of(rename_map)) |>
    filter(
      # 1. Species Confirmation via VITEK (Primary Filter)
      trimws(VITEK_MS_Results) == "Klebsiella pneumoniae",
      
      # 2. Specific Isolate + Morphology logic
      (
        # Group A: K. pneumoniae that are Metallic blue
        (grepl("K", Isolate, ignore.case = TRUE) & 
           `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") |
          
          # Group B: E. coli that are Pinkish/Reddish (but VITEK says Kpn)
          (grepl("E", Isolate, ignore.case = TRUE) & 
             `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish"))
      ),
      
      # 3. Third-generation cephalosporin susceptibility breakpoint (< 23 mm)
      Ceftriaxone < 23
    ) |>
    arrange(INIKA_ID, Isolate) |>
    # Ensure unique participants are counted once
    distinct(INIKA_ID, .keep_all = TRUE)
  
  message(paste("Isolates processed:", nrow(cleaned_data)))
  
  # Pivot & Susceptibility Classification 
  kpn_long <- cleaned_data |>
    pivot_longer(
      cols = any_of(names(rename_map)),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) |>
    left_join(Kpn_EPI_CUTOFF, by = "Antimicrobial_substance") |>
    mutate(
      across(c(S, R), ~as.numeric(str_remove_all(.x, "[^0-9.]")), .names = "{.col}_num"),
      Measured_Zone = as.numeric(Measured_Zone),
      Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(Result))
  
  # Summary Helper Function 
  get_summary <- function(df, var_name) {
    group_cols <- if(var_name == "Overall") "Antimicrobial_substance" else c(var_name, "Antimicrobial_substance")
    
    df |>
      group_by(across(all_of(group_cols))) |>
      summarise(
        n_res = sum(Result == "R"),
        total = n(),
        .groups = "drop"
      ) |>
      mutate(
        Grouping_Var = var_name,
        Grouping_Val = if(var_name == "Overall") "Overall" else as.character(.data[[var_name]])
      )
  }
  
  # Generate and Format Final Table
  group_vars <- c("Overall", "REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE")
  
  final_table <- map_dfr(group_vars, ~get_summary(kpn_long, .x)) |>
    rowwise() |>
    mutate(
      perc = round((n_res / total) * 100, 1),
      ci = if(Grouping_Var == "Overall") {
        bt <- binom.test(n_res, total)
        paste0(" (", round(bt$conf.int[1]*100, 1), "–", round(bt$conf.int[2]*100, 1), ")")
      } else "",
      cell_text = paste0(perc, "% (", n_res, "/", total, ")", ci)
    ) |>
    ungroup() |>
    select(Antimicrobial_substance, Grouping_Val, cell_text) |>
    pivot_wider(names_from = Grouping_Val, values_from = cell_text) |>
    left_join(abx_classes, by = "Antimicrobial_substance") |>
    select(Class, Antimicrobial_substance, Overall, everything()) |>
    arrange(Class, Antimicrobial_substance)
  
  # Export to Word (Landscape & Fitted) 
  if(!dir.exists("Results")) dir.create("Results")
  
  # Create the flextable
  ft <- flextable(final_table) |>
    theme_booktabs() |>
    fontsize(size = 8, part = "all") |>        # Small font to fit many columns
    padding(padding = 1, part = "all") |>      # Tighten cell space
    merge_v(j = ~ Class) |>                    # Group by Class vertically
    bold(part = "header") |>
    set_table_properties(layout = "autofit")
  
  # Define Landscape & Narrow Margins
  sect_props <- prop_section(
    page_size = page_size(orient = "landscape"),
    page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5),
    type = "continuous"
  )
  
  # Save as Word
  save_as_docx(ft, path = "Results/Kpn_AMR_Table_Landscape.docx", pr_section = sect_props)
  
  # Also save Excel as backup
  write_xlsx(final_table, "Results/Kpn_AMR_Table.xlsx")
  
  message("Success: Table saved to Results folder as Word (Landscape) and Excel.")
  return(final_table)
}

#  Execute 
kpn_amr_final_results <- process_kpn_amr_analysis(joined_data, Kpn_EPI_CUTOFF)
################################################################
# Supplementary Table 7
# Computing the MDR for ESCR E.coli
# DEFINE LOOKUP TABLES 
abx_classes <- tibble(
  Antimicrobial_substance = c(
    "Amoxicillin", "Azithromycin", "Ceftriaxone", "Cefotaxime", 
    "Cefotaxime/ClavulanicAcid", "Ciprofloxacin", "Doxycycline", 
    "Oxytetracycline", "Gentamicin", "Meropenem", "Florfenicol", 
    "Polymyxin_B(PB)", "Sulfamethoxazole/Trimethoprim"
  ),
  Class = c(
    "Penicillins", "Macrolides", "Cephalosporins (3rd Gen)", "Cephalosporins (3rd Gen)", 
    "β-lactam/Inhibitor", "Fluoroquinolones", "Tetracyclines", 
    "Tetracyclines", "Aminoglycosides", "Carbapenems", "Phenicols", 
    "Polymyxins", "Sulfonamides"
  )
)

# DEFINE THE PROCESSING FUNCTION 
process_ecoli_MDR_results <- function(joined_data, EPI_CUTOFF, abx_classes) {
  
  # Clean, Filter, and De-duplicate
  Data_Clean <- joined_data %>%
    filter(
      (
        (
          # Group 1: Standard E. coli
          grepl("E", Isolate, ignore.case = TRUE) & 
            grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & 
            grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
        ) | 
          (
            # Group 2: Metallic blue isolates that VITEK says are actually E. coli
            grepl("K", Isolate, ignore.case = TRUE) & 
              grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & 
              grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
          )
      ) & 
        # Ceftriaxone resistance breakpoint criterion (< 23 mm)
        CRO_ED30 < 23
    ) %>%
    # Use any_of for safety during the select
    select(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE, `COLONY MORPHOLOGY ON C3GR`, Isolate,
           VITEK_MS_Results, ESCR_ECO_presumptive, AMX_ED10, AZM_ED15, CRO_ED30,
           CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
           OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
    rename(
      Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
      Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
      Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
      Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
      Oxytetracycline = "OXY_ED30", "Polymyxin_B(PB)" = "POL_ED300",
      "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
      Cefotaxime = "CTX_ED5", "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
    ) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
  
  # Console message to track N
  message(paste("Isolates processed for MDR analysis:", nrow(Data_Clean)))
  
  # Categorize (SIR)
  Ecoli_SIR <- Data_Clean %>%
    pivot_longer(
      cols = any_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) %>%
    left_join(EPI_CUTOFF, by = "Antimicrobial_substance") %>%
    left_join(abx_classes, by = "Antimicrobial_substance") %>%
    mutate(
      Measured_Zone = as.numeric(str_trim(Measured_Zone)),
      S_num = as.numeric(str_replace_all(S, "[^0-9.]", "")),
      R_num = as.numeric(str_replace_all(R, "[^0-9.]", "")),
      Categorical_Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ "S" 
      )
    )
  
  #  MDR Calculation (Isolate Level)
  MDR_Data <- Ecoli_SIR %>%
    group_by(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE, Class) %>%
    summarise(Class_Resistant = any(Categorical_Result %in% c("R", "I")), .groups = "drop_last") %>%
    group_by(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE) %>%
    summarise(
      Resistant_Classes_Count = sum(Class_Resistant, na.rm = TRUE),
      MDR_Status = ifelse(Resistant_Classes_Count >= 3, "MDR", "Non-MDR"),
      .groups = "drop"
    )
  
  return(MDR_Data)
}

#  RUN ANALYSIS 
ecoli_mdr_results <- process_ecoli_MDR_results(joined_data, EPI_CUTOFF, abx_classes)

# PREPARE PUBLICATION TABLE 
calc_stats <- function(df) {
  n_total <- nrow(df)
  n_mdr <- sum(df$MDR_Status == "MDR")
  ci <- binconf(n_mdr, n_total, method = "wilson")
  tibble(N = n_total, MDR_n = n_mdr, Percent = ci[1]*100, Lower = ci[2]*100, Upper = ci[3]*100)
}

overall <- ecoli_mdr_results %>% calc_stats() %>% mutate(Variable = "Overall", Category = "Total Population")
grouped <- ecoli_mdr_results %>%
  pivot_longer(cols = c(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE), names_to = "Variable", values_to = "Category") %>%
  group_by(Variable, Category) %>%
  group_modify(~ calc_stats(.x)) %>%
  ungroup()

final_ecoli_mdr_table <- bind_rows(overall, grouped) %>%
  mutate(
    Variable = recode(Variable, "REGION.x" = "Region", "SEASON.x" = "Season", "ORIGIN_OF_SAMPLE" = "Origin"),
    `MDR n (%)` = sprintf("%d (%.1f%%)", MDR_n, Percent),
    `95% CI` = sprintf("[%.1f - %.1f]", Lower, Upper)
  ) %>%
  select(Variable, Category, N, `MDR n (%)`, `95% CI`)

# SAVE TO WORD 
doc_table <- final_ecoli_mdr_table %>%
  as_grouped_data(groups = "Variable") %>%
  flextable() %>%
  set_header_labels(N = "Total N", `MDR n (%)` = "MDR n (%)", `95% CI` = "95% CI (Wilson)") %>%
  bold(part = "header") %>%
  autofit()

save_as_docx(doc_table, path = "Ecoli_MDR_Final_Table.docx")

save_as_docx(doc_table, path = "Results/Ecoli_MDR_Final_Table.docx")
##################################################################
# Supplementary Table 8
# Computing the MDR for ESCR K.pneumoniae

process_kpn_mdr_analysis <- function(joined_data, Kpn_EPI_CUTOFF) {
  
  # Internal Reference Data 
  abx_classes <- tribble(
    ~Antimicrobial_substance,         ~Class,
    "Amoxicillin",                    "Penicillins",
    "Azithromycin",                   "Macrolides",
    "Ceftriaxone",                    "Cephalosporins (3rd Gen)",
    "Cefotaxime",                     "Cephalosporins (3rd Gen)",
    "Cefotaxime/ClavulanicAcid",      "β-lactam/Inhibitor",
    "Ciprofloxacin",                  "Fluoroquinolones",
    "Doxycycline",                    "Tetracyclines",
    "Oxytetracycline",                "Tetracyclines",
    "Gentamicin",                     "Aminoglycosides",
    "Meropenem",                      "Carbapenems",
    "Florfenicol",                    "Phenicols",
    "Polymyxin_B(PB)",                "Polymyxins",
    "Sulfamethoxazole/Trimethoprim",  "Sulfonamides"
  )
  
  rename_map <- c(
    Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
    Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
    Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
    Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
    Oxytetracycline = "OXY_ED30", "Polymyxin_B(PB)" = "POL_ED300",
    "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
    Cefotaxime = "CTX_ED5", "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
  )
  
  # Data Cleaning 
  cleaned_data <- joined_data %>%
    rename(any_of(rename_map)) %>%
    filter(
      # 1. Species Confirmation via VITEK (Primary Filter)
      trimws(VITEK_MS_Results) == "Klebsiella pneumoniae",
      
      # 2. Specific Isolate + Morphology logic
      (
        # Group A: K. pneumoniae that are Metallic blue
        (grepl("K", Isolate, ignore.case = TRUE) & 
           `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") |
          
          # Group B: E. coli that are Pinkish/Reddish (but VITEK says Kpn)
          (grepl("E", Isolate, ignore.case = TRUE) & 
             `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish"))
      ) & Ceftriaxone < 23
    ) |>
    arrange(INIKA_ID, Isolate) |>
    # Ensure unique participants are counted once
    distinct(INIKA_ID, .keep_all = TRUE)
  # Pivot & Classify 
  kpn_long <- cleaned_data %>%
    pivot_longer(cols = any_of(names(rename_map)), 
                 names_to = "Antimicrobial_substance", values_to = "Measured_Zone") %>%
    left_join(Kpn_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
    left_join(abx_classes, by = "Antimicrobial_substance") %>%
    mutate(
      Measured_Zone = as.numeric(Measured_Zone),
      S_num = as.numeric(str_remove_all(S, "[^0-9.]")),
      R_num = as.numeric(str_remove_all(R, "[^0-9.]")),
      Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ "S"
      )
    )
  
  # MDR Calculation 
  Kpn_mdr_summary <- kpn_long %>%
    group_by(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE, Class) %>%
    summarise(Resistant_in_Class = any(Result %in% c("R", "I")), .groups = "drop") %>%
    group_by(INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE) %>%
    summarise(
      N_Classes_Resistant = sum(Resistant_in_Class),
      MDR_Status = if_else(N_Classes_Resistant >= 3, "MDR", "Non-MDR"),
      .groups = "drop"
    )
  
  return(Kpn_mdr_summary)
}

# Execute Analysis
kpn_mdr_results <- process_kpn_mdr_analysis(joined_data, Kpn_EPI_CUTOFF)

# PREPARE PUBLICATION TABLE 
calc_stats <- function(df) {
  n_total <- nrow(df)
  n_mdr <- sum(df$MDR_Status == "MDR")
  ci <- binconf(n_mdr, n_total, method = "wilson")
  tibble(N = n_total, MDR_n = n_mdr, Percent = ci[1]*100, Lower = ci[2]*100, Upper = ci[3]*100)
}

overall <- kpn_mdr_results %>% calc_stats() %>% mutate(Variable = "Overall", Category = "Total Population")
grouped <- kpn_mdr_results %>%
  pivot_longer(cols = c(REGION.x, SEASON.x, ORIGIN_OF_SAMPLE), names_to = "Variable", values_to = "Category") %>%
  group_by(Variable, Category) %>%
  group_modify(~ calc_stats(.x)) %>%
  ungroup()

final_Kpn_mdr_table <- bind_rows(overall, grouped) %>%
  mutate(
    Variable = recode(Variable, "REGION.x" = "Region", "SEASON.x" = "Season", "ORIGIN_OF_SAMPLE" = "Origin"),
    `MDR n (%)` = sprintf("%d (%.1f%%)", MDR_n, Percent),
    `95% CI` = sprintf("[%.1f - %.1f]", Lower, Upper)
  ) %>%
  select(Variable, Category, N, `MDR n (%)`, `95% CI`)

# STEP 5: SAVE TO WORD 
doc_table <- final_Kpn_mdr_table %>%
  as_grouped_data(groups = "Variable") %>%
  flextable() %>%
  set_header_labels(N = "Total N", `MDR n (%)` = "MDR n (%)", `95% CI` = "95% CI (Wilson)") %>%
  bold(part = "header") %>%
  autofit()

save_as_docx(doc_table, path = "Kpn_MDR_Final_Table.docx")
save_as_docx(doc_table, path = "Results/Kpn_MDR_Final_Table.docx")
##################################################################
#Figure 1
#Drawing a histogram for ESCR E.coli AMR frequency
# 1. DATA FILTERING & DEDUPLICATION
ESCR_ECO_confirmed <- joined_data %>%
  rename(any_of(setNames(original_spelling, tested_antibiotics))) %>%
  filter(
    # Species Confirmation via VITEK
    trimws(VITEK_MS_Results) == "Escherichia coli",
    # Specific Isolate + Morphology logic
    (
      (grepl("K", Isolate, ignore.case = TRUE) & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") |
        (grepl("E", Isolate, ignore.case = TRUE) & `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish"))
    ) & `Ceftriaxone` < 23
  ) %>%
  arrange(INIKA_ID, Isolate) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

# 2. LONG TRANSFORM & RESISTANCE CALCULATION
ESCR_long <- ESCR_ECO_confirmed %>% 
  pivot_longer(
    cols = any_of(tested_antibiotics), 
    names_to = "Antimicrobial_substance", 
    values_to = "mm"
  ) %>%
  filter(!is.na(mm)) %>%
  mutate(
    # String cleaning before numeric conversion
    mm = as.numeric(str_trim(as.character(mm)))
  ) %>%
  left_join(EPI_CUTOFF, by = "Antimicrobial_substance") %>%
  mutate(
    R_limit = as.numeric(str_replace_all(as.character(R), "[^0-9.]", "")),
    is_resistant = if_else(mm <= R_limit, 1, 0, missing = 0)
  )

# 3. GROUP SUMMARY & EXACT 95% CIs
Summary_Resistance <- ESCR_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    R_count = sum(is_resistant, na.rm = TRUE),
    Prop = R_count / N,
    .groups = "drop"
  ) %>%
  filter(Prop > 0) %>%
  # Safe vectorized confidence interval calculation using purrr
  mutate(
    ci = map2(R_count, N, ~ binom.test(.x, .y)$conf.int * 100),
    Lower = map_dbl(ci, 1),
    Upper = map_dbl(ci, 2)
  ) %>%
  select(-ci)

# 4. GENERATE REFINED PLOT
AMR_CI_Plot <- ggplot(
  Summary_Resistance, 
  aes(
    x = reorder(Antimicrobial_substance, Prop), 
    y = Prop * 100, 
    fill = Antimicrobial_substance
  )
) + 
  geom_col(color = "black", alpha = 0.9, width = 0.7, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper), 
    width = 0.25, color = "gray20", linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", Prop * 100),
      hjust = if_else(Antimicrobial_substance == "Florfenicol", -0.5, 1.2),
      color = if_else(Antimicrobial_substance == "Florfenicol", "black", "white")
    ), 
    size = 3.5, fontface = "bold"
  ) +
  scale_color_identity() +
  coord_flip() + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Resistance Percentage (%)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold")
  )

# 5. SAVE AND EXPORT TO WORD
if (!dir.exists("Results")) dir.create("Results")

ggsave("Results/Resistance_Refined_Plot.png", plot = AMR_CI_Plot, width = 9, height = 6, dpi = 300)

doc <- read_docx() %>%
  body_add_img(src = "Results/Resistance_Refined_Plot.png", width = 6.5, height = 4.3)

print(doc, target = "Results/ESCR_AMR_Refined_Histogram.docx")


################################################################################   
# Figure 2
# Drawing a histogram for ESCR K.pneumoniae AMR frequency
# DATA PREP  
ESCR_Kpn_confirmed <- joined_data %>%
  rename(any_of(setNames(original_spelling, tested_antibiotics))) %>%
  filter(
    # Species Confirmation via VITEK
    trimws(VITEK_MS_Results) == "Klebsiella pneumoniae",
    # Specific Isolate + Morphology logic
    (
      (grepl("K", Isolate, ignore.case = TRUE) & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") |
        (grepl("E", Isolate, ignore.case = TRUE) & `COLONY MORPHOLOGY ON C3GR` %in% c("Pink", "Pinkish", "Reddish", "Pinkish/Reddish"))
    ) & Ceftriaxone < 23
  ) %>%
  arrange(INIKA_ID, Isolate) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

#  CALCULATIONS 
ESCR_Kpn_long <- ESCR_Kpn_confirmed %>%
  pivot_longer(cols = any_of(tested_antibiotics), 
               names_to = "Antimicrobial_substance", 
               values_to = "mm") %>%
  filter(!is.na(mm)) %>%
  left_join(Kpn_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    is_resistant = ifelse(mm <= R_limit, 1, 0)
  )

Summary_Kpn <- ESCR_Kpn_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    R_count = sum(is_resistant, na.rm = TRUE),
    Prop = (R_count / N) * 100,
    Lower = binom.test(R_count, N)$conf.int[1] * 100,
    Upper = binom.test(R_count, N)$conf.int[2] * 100,
    .groups = "drop"
  ) %>%
  # Remove antibiotics with 0% resistance
  filter(Prop > 0)

# GENERATE PLOT 
# Dynamic N for the title
total_n_kpn <- nrow(ESCR_Kpn_confirmed)

Kpn_AMR_Plot <- ggplot(Summary_Kpn, 
                       aes(x = reorder(Antimicrobial_substance, Prop), 
                           y = Prop, 
                           fill = Antimicrobial_substance)) + 
  geom_col(color = "black", alpha = 0.9, width = 0.7, show.legend = FALSE) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                width = 0.25, color = "gray20", linewidth = 0.7) +
  # Conditional positioning: Florfenicol outside/black, others inside/white
  geom_text(aes(
    label = sprintf("%.1f%%", Prop),
    hjust = ifelse(Antimicrobial_substance == "Florfenicol", -0.5, 1.2),
    color = ifelse(Antimicrobial_substance == "Florfenicol", "black", "white")
  ), size = 3.5, fontface = "bold") +
  scale_color_identity() +
  coord_flip() + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = paste0("Resistance Percentage (%)"),
    x = " ",
    y = " "
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold")
  )

# EXPORT 
if(!dir.exists("Results")) dir.create("Results")
ggsave("Results/Kpn_Resistance_Histogram.png", plot = Kpn_AMR_Plot, width = 9, height = 6, dpi = 300)

doc_kpn <- read_docx() %>%
  body_add_img(src = "Results/Kpn_Resistance_Histogram.png", width = 6.5, height = 4.3)

print(doc_kpn, target = "Results/ESCR_Kpn_AMR_Histogram.docx")
################################################################################   
