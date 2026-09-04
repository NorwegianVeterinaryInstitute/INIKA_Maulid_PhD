## Factors associated with ESBL carriage
# By MMJ

library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(leaflet)
library(gtsummary)
library(flextable)
library(broom.helpers)
library(logistf)
library(read_csv2)



joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")
names(joined_data)
# Maulid: I would not include the criteria for ESBL_Selection for this paper:
#I put a # before filter(ESBL_Selection >= 5).
#You also should include the CRE as positives:
# Then you have defined all others as  "Negative/No Growth", but which are these samples?
#This should have been those that were NO growth from the biochemical results or where you dound the Salmonella typhiumirum which was afterwards neagtive all of them.

#  DATA PREPARATION
df_final <- joined_data %>%
  # Apply the ESBL screening threshold
  #filter(ESBL_Selection >= 5) %>%
  mutate(
    # Create the text label for species identification
    Confirmed_Species = case_when(
      (str_detect(Isolate, "E.coli|E. coli") & VITEK_MS_Results == "Escherichia coli" & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish") ~ "ESBL E. coli",
      (str_detect(Isolate, "K.pneumoniae") & VITEK_MS_Results == "Escherichia coli" & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") ~ "ESBL E. coli",
      (str_detect(Isolate, "K.pneumoniae") & VITEK_MS_Results == "Klebsiella pneumoniae" & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") ~ "ESBL K. pneumoniae",
      (str_detect(Isolate, "E.coli|E. coli") & VITEK_MS_Results == "Klebsiella pneumoniae" & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish") ~ "ESBL K. pneumoniae",
      TRUE ~ "Negative/No Growth"
    ),
    # Create numeric Outcome (1 = Positive, 0 = Negative) for the regression
    Outcome_ESC_Bact = if_else(Confirmed_Species == "Negative/No Growth", 0, 1)
  ) %>%
  # Select and Rename columns
  select(
    Outcome_ESC_Bact,
    GENDER,
    REGION = REGION.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    Heard_AMR = "Have you ever heard about AMR?",
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    Animals_Treated = "Are your animals being treated for any diseases?",
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?"
  ) %>%
  # Convert to factors to ensure levels (Yes/No, etc.) display sequentially
  mutate(across(-Outcome_ESC_Bact, ~as.factor(.x))) %>%
  # Create a factor version of the outcome for table display
  mutate(Outcome_Factor = factor(Outcome_ESC_Bact, levels = c(1, 0), labels = c("ESBL Positive", "Negative")))

# DEFINE LABELS
var_labels <- list(
  GENDER ~ "Gender",
  REGION ~ "Region",
  ORIGIN_OF_SAMPLE ~ "Source of Sample",
  Heard_AMR ~ "Heard about AMR",
  Prescription_After_Lab ~ "Prescription after Lab Result",
  Owns_Farm_Animals ~ "Owns Farm Animals",
  Animals_Treated ~ "Animals Treated for Disease",
  Withdrawal_Time ~ "Withdrawal Time Knowledge",
  Has_Toilet ~ "Availability of Toilet",
  Wash_Hands_Toilet ~ "Hand washing after toilet visit"
)

# RUN THE MODEL
firth_model <- logistf(
  Outcome_ESC_Bact ~ GENDER + REGION + ORIGIN_OF_SAMPLE + Heard_AMR + 
    Prescription_After_Lab + Owns_Farm_Animals + Animals_Treated + 
    Withdrawal_Time + Has_Toilet + Wash_Hands_Toilet,
  data = df_final,
  control = logistf.control(maxit = 400)
)

# CREATE THE SEQUENTIAL TABLES
# Descriptive Table ( I would rather call thos the univariate table?)
tab_counts <- df_final %>%
  select(Outcome_Factor, GENDER, REGION, ORIGIN_OF_SAMPLE, Heard_AMR, 
         Prescription_After_Lab, Owns_Farm_Animals, Animals_Treated, 
         Withdrawal_Time, Has_Toilet, Wash_Hands_Toilet) %>%
  tbl_summary(
    by = Outcome_Factor,
    label = var_labels,
    # Every variable to show its sub-responses (Yes/No)
    type = all_categorical() ~ "categorical", 
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Total**")

# Regression Table
tab_reg <- firth_model %>%
  tbl_regression(
    exponentiate = TRUE,
    label = var_labels,
    tidy_fun = broom.helpers::tidy_parameters,
    # This ensures "Yes" and "No" appear in the regression column too
    add_estimate_to_reference_rows = TRUE 
  )

# MERGE AND FINAL FORMATTING 
final_table <- tbl_merge(
  tbls = list(tab_counts, tab_reg),
  tab_spanner = c("**Counts (ESBL Status)**", "**Firth Logistic Regression**")
) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_header(label = "**Variable & Response**") %>%
  # Adding the table caption here
  modify_caption("**Table 1. Sociodemographic and Behavioral Factors Associated with Confirmed ESBL Carriage**")

# EXPORT TO "Results" FOLDER
if (!dir.exists("Results")) {
  dir.create("Results")
}

final_table %>%
  as_flex_table() %>%
  flextable::autofit() %>%
  flextable::save_as_docx(path = "Results/Factors_Associated_with_ESBL_Sequential.docx")

message("Table with caption successfully saved to 'Results/Factors_Associated_with_ESBL_Sequential.docx'.")

################################################################################

# 1. DATA PREPARATION (Maintaining logic from Code 1)
df_final <- joined_data %>%
  distinct(INIKA_ID.x.x, .keep_all = TRUE) %>%
  mutate(
    Confirmed_Status = case_when(
      (str_detect(Isolate, "E.coli|E. coli") & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish" & VITEK_MS_Results == "Escherichia coli") ~ "ESCR Positive",
      (str_detect(Isolate, "K.pneumoniae") & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue" & VITEK_MS_Results == "Escherichia coli" & CTX_ED5 < 22) ~ "ESCR Positive",
      (str_detect(Isolate, "K.pneumoniae") & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue" & VITEK_MS_Results == "Klebsiella pneumoniae") ~ "ESCR Positive",
      (str_detect(Isolate, "E.coli|E. coli") & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish" & VITEK_MS_Results == "Klebsiella pneumoniae" & CTX_ED5 < 21) ~ "ESCR Positive",
      TRUE ~ "Negative"
    ),
    Outcome_ESC_Bact = if_else(Confirmed_Status == "ESCR Positive", 1, 0)
  ) %>%
  select(
    Outcome_ESC_Bact, 
    GENDER, 
    REGION = REGION.x, 
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    Heard_AMR = "Have you ever heard about AMR?",
   Antibiotics_use_lifetime = "Have you or your children used any antibiotics at any time?" ,
   Antibiotics_use_last_six_month = "Have you or your children used any antibiotics in the past six months?" ,
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
   Taking_care_animals = "If Yes, do you participate in taking care of animals i.e. cleaning, Feeding?" ,
    Animals_Treated = "Are your animals being treated for any diseases?",
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
   Owns_poultry_farm = "Do you have a poultry farm at your home?",
   Taking_care_poultry = "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?",
   Poultry_Treated = "Are your Poultry being treated by veterinary doctor?" ,
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?"
  ) %>%
  mutate(across(-Outcome_ESC_Bact, ~as.factor(.x))) %>%
  mutate(Outcome_Factor = factor(Outcome_ESC_Bact, levels = c(1, 0), labels = c("ESCR Positive", "Negative")))

# 2. DEFINE LABELS & VARIABLES
var_labels <- list(
  GENDER ~ "Gender", 
  REGION ~ "Region", 
  ORIGIN_OF_SAMPLE ~ "Source of Sample",
  Heard_AMR ~ "Heard about AMR",
  Antibiotics_use_lifetime ~ "Antibiotics use in lifetime",
  Antibiotics_use_last_six_month ~ "Antibiotics use in last six months",
  Prescription_After_Lab ~ "Prescription after Lab Result",
  Owns_Farm_Animals ~ "Owns Farm Animals",
  Taking_care_animals ~ "Taking care of animals",
  Animals_Treated ~ "Animals Treated for Disease",
  Withdrawal_Time ~ "Withdrawal Time Knowledge",
  Owns_poultry_farm ~ "Owns Poultry farm",
  Taking_care_poultry ~ "Taking care of poultry",
  Poultry_Treated ~ "Poultry treated for disease",
  Has_Toilet ~ "Availability of Toilet",
  Wash_Hands_Toilet ~ "Hand washing after toilet visit"
)

analysis_vars <- c("GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR",
                   "Antibiotics_use_lifetime","Antibiotics_use_last_six_month",
                   "Prescription_After_Lab", "Owns_Farm_Animals","Taking_care_animals",
                   "Animals_Treated","Owns_poultry_farm","Taking_care_poultry",
                   "Poultry_Treated","Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet")

df_analysis <- df_final %>%
  select(Outcome_Factor, Outcome_ESC_Bact, all_of(analysis_vars))

# 3. BASE DESCRIPTIVE TABLE (Updated to show all response levels)
tab_counts <- df_analysis %>%
  select(Outcome_Factor, all_of(analysis_vars)) %>%
  tbl_summary(
    by = Outcome_Factor,
    label = var_labels,
    # This ensures "Yes/No" levels both appear in the Label column
    type = all_categorical() ~ "categorical", 
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Total**")

# 4. TABLE 1: DESCRIPTIVE + UNIVARIATE
tab_univariate <- df_analysis %>%
  select(Outcome_ESC_Bact, all_of(analysis_vars)) %>%
  tbl_uvregression(
    method = logistf::logistf,
    y = Outcome_ESC_Bact,
    exponentiate = TRUE,
    label = var_labels,
    add_estimate_to_reference_rows = TRUE, # Ensures alignment with tab_counts rows
    method.args = list(control = logistf::logistf.control(maxit = 400))
  ) %>%
  modify_header(estimate ~ "**OR (95% CI)**")

table1_uni <- tbl_merge(
  tbls = list(tab_counts, tab_univariate),
  tab_spanner = c("**Descriptive Statistics**", "**Univariate Analysis**")
) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_header(label = "**Variable & Response**") %>%
  modify_caption("**Table 1. Descriptive and Univariate Analysis of Factors Associated with ESCR**")

# 5. TABLE 2: DESCRIPTIVE + MULTIVARIATE
multivariate_model <- logistf::logistf(
  Outcome_ESC_Bact ~ ., 
  data = df_analysis %>% select(-Outcome_Factor),
  control = logistf::logistf.control(maxit = 400)
)

tab_multivariate <- multivariate_model %>%
  tbl_regression(
    exponentiate = TRUE,
    label = var_labels,
    add_estimate_to_reference_rows = TRUE, # Ensures alignment with tab_counts rows
    tidy_fun = broom.helpers::tidy_parameters
  ) %>%
  modify_header(estimate ~ "**aOR (95% CI)**")

table2_multi <- tbl_merge(
  tbls = list(tab_counts, tab_multivariate),
  tab_spanner = c("**Descriptive Statistics**", "**Multivariate Analysis**")
) %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_header(label = "**Variable & Response**") %>%
  modify_caption("**Table 2. Descriptive and Multivariate Analysis of Factors Associated with ESCR**")

# 6. EXPORT
if (!dir.exists("Results")) dir.create("Results")

# Export Table 1
table1_uni %>%
  as_flex_table() %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit() %>%
  flextable::save_as_docx(path = "Results/Table1_Descriptive_Univariate.docx")

# Export Table 2
table2_multi %>%
  as_flex_table() %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit() %>%
  flextable::save_as_docx(path = "Results/Table2_Descriptive_Multivariate.docx")

message("Tables generated successfully with detailed response levels.")
#################################################################################
