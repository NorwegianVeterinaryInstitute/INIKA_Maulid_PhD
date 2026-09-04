install.packages("logistf")
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(leaflet)
library(maps)
library(car)
library(pscl)
library(pROC)
library(ResourceSelection)
library(forcats)
library(broom)
library(rcompanion)
library(logistf)
library(officer)
library(flextable)

# Importing the data set

joined_data_16_3_26<-read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")
names(joined_data_16_3_26)

# I define the cases as all those that have been identified by VITEK and ESC=1 as Cases=1), then all presumptive salmonella or no growth are defined as controls
# Then we remove all those that could be cases but are either found to be another spp. or not identified by the VITEK confirmation-very conservative.
#I did this by checking that the resulting file handles the data as I want it to be handled so I only selected essential variables, afterwards I just do the same without selection
Original_joined <- joined_data_16_3_26 %>%
  select(INIKA_ID.x.x, PROTOCOL, Isolate, Isolate_ID, ESC, VITEK_MS_Results, Isolate_NVI) %>%
mutate(
  Case = case_when(
    ESC == 1 | Isolate_NVI %in% c("E.coli", "Klebsiella pneumoniae")  ~ "1",
    Isolate == "No growth" ~ "0",
    Isolate %in% c("S.typhimurium", "S.paratyphi A", "S.typhi") ~ "0",
    
    !is.na(VITEK_MS_Results) ~ "EX",
    is.na(ESC) ~ "EX",
    TRUE ~ as.character(ESC)
  )
)%>%
  filter(!Case=="EX")

Original_joined <- joined_data_16_3_26 %>%
  select(INIKA_ID.x.x, PROTOCOL, Isolate, Isolate_ID, ESC, VITEK_MS_Results, Isolate_NVI) %>%
  mutate(
    Case = case_when(
      ESC == 1 | Isolate_NVI %in% c("E.coli", "Klebsiella pneumoniae") ~ "1",
      
      Isolate == "No growth" ~ "0",
      Isolate %in% c("S.typhimurium", "S.paratyphi A", "S.typhi") ~ "0",
      
      # ONLY non-NA AND not the two species → EX
      !is.na(VITEK_MS_Results) & 
        !VITEK_MS_Results %in% c("E.coli", "Klebsiella pneumoniae") ~ "EX",
      
      is.na(ESC) ~ "EX",
      
      TRUE ~ as.character(ESC)
    )
  ) %>%
  filter(Case != "EX")

Original_joined <- joined_data_16_3_26 %>%
  mutate(
    Case = case_when(
      ESC == 1 | Isolate_NVI %in% c("E.coli", "Klebsiella pneumoniae") ~ "1",
      
      Isolate == "No growth" ~ "0",
      Isolate %in% c("S.typhimurium", "S.paratyphi A", "S.typhi") ~ "0",
      
      # ONLY non-NA AND not the two species → EX
      !is.na(VITEK_MS_Results) & 
        !VITEK_MS_Results %in% c("E.coli", "Klebsiella pneumoniae") ~ "EX",
      
      is.na(ESC) ~ "EX",
      
      TRUE ~ as.character(ESC)
    )
  ) %>%
  filter(Case != "EX")%>%# We only keep the columns(variables needed)
  # Create numeric Outcome (1 = Positive, 0 = Negative) for the regression
  mutate(Case=as.numeric(Case))%>%
  # Select and Rename columns
  select(
    Case,
    GENDER,
    REGION = REGION.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    Heard_AMR = "Have you ever heard about AMR?",
    Used_AB_EVER="Have you or your children used any antibiotics at any time?",                                     
    Used_AB_SIXM= "Have you or your children used any antibiotics in the past six months?",
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    OWNS_Poultry = "Do you have a poultry farm at your home?" , 
    Animals_Treated = "Are your animals being treated for any diseases?",
    Caretaking_Animals = "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?",  
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?"
  ) %>%
  # Convert to factors to ensure levels (Yes/No, etc.) display sequentially
  mutate(across(-Case, ~as.factor(.x))) %>%
  # Create a factor version of the outcome for table display
  mutate(Outcome_Factor = factor(Case, levels = c(1, 0), labels = c("ESCR Positive", "Negative")))

# Before doing anything get the frequencies of the answers from all the Independent variables


vars <- c("GENDER","REGION", "ORIGIN_OF_SAMPLE","Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
          "Prescription_After_Lab","Owns_Farm_Animals","OWNS_Poultry","Animals_Treated",
          "Caretaking_Animals","Withdrawal_Time","Has_Toilet","Wash_Hands_Toilet"
          )

result <- lapply(vars, function(v) {
  Original_joined %>%
    mutate(
      var = .data[[v]],
      var = ifelse(is.na(var), "NA", as.character(var)),
      Case = ifelse(is.na(Case), "NA", as.character(Case))
    ) %>%
    count(variable = v, value = var, Case)
}) %>%
  bind_rows()

result

##################################################
# Table 1. Frequencies in %
# Clean the data set by filtering out "Missing" entries for either controls or cases
data_model_clean <- Original_joined %>% 
  filter(
    Used_AB_SIXM != "Missing",
    Owns_Farm_Animals != "Missing"
  )

# Define the variables for frequency calculation
vars <- c("GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
          "Prescription_After_Lab", "Owns_Farm_Animals", "OWNS_Poultry", "Animals_Treated",
          "Caretaking_Animals", "Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet")

# Compute the frequencies using the newly cleaned data (data_model_clean)
result <- lapply(vars, function(v) {
  data_model_clean %>% 
    mutate(
      var = .data[[v]],
      var = ifelse(is.na(var), "NA", as.character(var)),
      Case = ifelse(is.na(Case), "NA", as.character(Case))
    ) %>%
    count(variable = v, value = var, Case)
}) %>%
  bind_rows()

# Re-code 0/1 to Control/Case, pivot wider, and calculate percentages
table_ready <- result %>%
  mutate(Case = case_when(
    Case == "0" ~ "Control",
    Case == "1" ~ "Case",
    TRUE        ~ as.character(Case) # Keeps "NA" or any other value intact
  )) %>%
  pivot_wider(
    names_from = Case, 
    values_from = n, 
    values_fill = 0
  ) %>%
  # Ensure both Control and Case columns exist even if data is thin
  { if(!"Control" %in% names(.)) mutate(., Control = 0) else . } %>%
  { if(!"Case" %in% names(.)) mutate(., Case = 0) else . } %>%
  # Group by Variable to calculate percentages within each specific factor
  group_by(variable) %>%
  mutate(
    Control_Pct = (Control / sum(Control)) * 100,
    Case_Pct    = (Case / sum(Case)) * 100
  ) %>%
  ungroup() %>%
  mutate(Total = Control + Case) %>%
  rename(Variable = variable, Stratum = value) %>%
  # Format the text to combine counts and percentages (e.g., "150 (45.2%)")
  mutate(
    Control_fmt = sprintf("%s (%.1f%%)", format(Control, big.mark = ","), Control_Pct),
    Case_fmt    = sprintf("%s (%.1f%%)", format(Case, big.mark = ","), Case_Pct)
  ) %>%
  # Rearrange for final display
  select(Variable, Stratum, Control_fmt, Case_fmt, Total)

# Create a professional, publishable Word Table 
ft <- flextable(table_ready) %>%
  # Set professional column headers
  set_header_labels(
    Variable = "Variable",
    Stratum = "Stratum",
    Control_fmt = "Control (N, %)",
    Case_fmt = "Case (N, %)",
    Total = "Total (N)"
  ) %>%
  # Merge variable name rows vertically for clean presentation
  merge_v(j = ~ Variable) %>% 
  # Apply formal academic theme (black/grey lines, clean layout)
  theme_vanilla() %>%
  # Formatting headers and fonts
  bold(part = "header") %>%
  color(color = "#1F4E78", part = "header") %>% # Deep Navy professional header text
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  # Set alignment (Left-align labels, Right-align metrics)
  align(j = c("Variable", "Stratum"), align = "left", part = "all") %>%
  align(j = c("Control_fmt", "Case_fmt", "Total"), align = "right", part = "all") %>%
  # Add thousands separators for the raw Total column
  colformat_int(j = "Total", big.mark = ",") %>%
  # Autofit column widths safely
  autofit()

# Initialize Word Document and export to Results folder
dir.create("Results", showWarnings = FALSE) # Ensures folder exists

doc <- read_docx() %>%
  body_add_par("Statistical Analysis Report: Case-Control Crosstabulation", style = "heading 1") %>%
  body_add_par("Generated automatically from consolidated R workspace metrics.", style = "Normal") %>%
  body_add_par("", style = "Normal") %>% 
  body_add_flextable(ft) %>%
  body_add_par("Note: Missing values were treated as explicit 'NA' components. Rows were omitted if 'Used_AB_SIXM' or 'Owns_Farm_Animals' data was recorded as 'Missing'. Percentages reflect proportions within column totals per variable.", style = "Table Caption")

# Save file
print(doc, target = "Results/analysis_results1.docx")
#############################################################
# Now you can see how many that have answered and there are some NAs that needs to be handled.
# Below, we categorize them as missing and then we can see if they are associated with the outcome or not! This requires quite a lot of sensitivity testing
# Other option would be to imputate them or just remove all observations where any variables have missing content- but then we loose a lot of power of the analyses and it could also introduce bias
#########################################################
#MMJ 27/5/2026 
#MN 03.07.2026 Note you need to use the data_model_clean
#data_model <- Original_joined %>%
data_model <- data_model_clean %>%
  mutate
    Used_AB_SIXM           = fct_na_value_to_level(as.factor(Used_AB_SIXM), "Missing"),
    Prescription_After_Lab = fct_na_value_to_level(as.factor(Prescription_After_Lab), "Missing"),
    Owns_Farm_Animals      = fct_na_value_to_level(as.factor(Owns_Farm_Animals), "Missing"),
    OWNS_Poultry           = fct_na_value_to_level(as.factor(OWNS_Poultry), "Missing"),
    Animals_Treated        = fct_na_value_to_level(as.factor(Animals_Treated), "Missing"),
    Caretaking_Animals     = fct_na_value_to_level(as.factor(Caretaking_Animals), "Missing"),
    Withdrawal_Time        = fct_na_value_to_level(as.factor(Withdrawal_Time), "Missing"),
    Case                   = as.factor(Case)
  )
####################################################################
# Better to use Case as numeric! avoids factor confusion- what is the reference level 

data_model <- data_model %>%
  mutate(Case = as.numeric(as.character(Case)))

model_1 <- glm(Case ~ REGION,
               data = data_model,
               family = binomial)

exp(cbind(OR = coef(model_1), confint(model_1)))

summary(model_1)
anova(model_1, test = "Chisq")
# Region is not significant

model_2 <- glm(Case ~ ORIGIN_OF_SAMPLE,
               data = data_model,
               family = binomial)

exp(cbind(OR = coef(model_2), confint(model_2)))

summary(model_2)
# ORIGIN_OF_SAMPLE significant

model_3 <- glm(Case ~ Heard_AMR,
               data = data_model,
               family = binomial)

exp(cbind(OR = coef(model_3), confint(model_3)))

summary(model_3)
# this is at least <0.20 so keep 

####

#MN 03.07.26_ You should not include the variables that are conditional ("Animals_Treated", Prescription_After_Lab, Withdrawal_Time )   in the analyses:
# Variables to include
vars <- c(
  "GENDER",
  "REGION",
  "ORIGIN_OF_SAMPLE",
  "Heard_AMR",
  "Used_AB_EVER",
  "Used_AB_SIXM",
  "Owns_Farm_Animals",
  "OWNS_Poultry",
  "Caretaking_Animals",
  "Has_Toilet",
  "Wash_Hands_Toilet"
)

# Ensure Case is numeric (0/1)
data_model <- data_model %>%
  mutate(Case = as.numeric(as.character(Case)))

# Function to run one model
run_univ <- function(var) {
  formula <- as.formula(paste("Case ~", var))
  
  model <- glm(formula, data = data_model, family = binomial)
  
  broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%   # remove intercept
    mutate(
      Variable = var,
      OR = exp(estimate)
    ) %>%
    select(Variable, term, estimate, OR, p.value)
}

# Run all models
univ_results <- lapply(vars, run_univ) %>%
  bind_rows()

# View results
univ_results

# We should write this table out.

write.csv(univ_results,"Results/Univariate_results.csv")

#################################
# MMJ 3/6/2026 Creating a publishable table (Table 2)

# vars <- c(
#   "GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
#   "Prescription_After_Lab", "Owns_Farm_Animals", "OWNS_Poultry", "Animals_Treated",
#   "Caretaking_Animals", "Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet"
# )
# 
# # Ensure Case is numeric (0/1)
# data_model <- data_model %>%
#   mutate(Case = as.numeric(as.character(Case)))
# 
# # Updated Function to run model and extract 95% CI
# run_univ <- function(var) {
#   formula <- as.formula(paste("Case ~", var))
#   model <- glm(formula, data = data_model, family = binomial)
#   
#   # Set conf.int = TRUE to calculate 95% Confidence Intervals
#   broom::tidy(model, conf.int = TRUE) %>%
#     filter(term != "(Intercept)") %>%   # remove intercept
#     mutate(
#       Variable = var,
#       OR = exp(estimate),
#       Lower_CI = exp(conf.low),         # Exponentiate lower bound
#       Upper_CI = exp(conf.high)         # Exponentiate upper bound
#     ) %>%
#     select(Variable, term, estimate, OR, Lower_CI, Upper_CI, p.value)
# }
# 
# # Run all models with CIs
# univ_results <- lapply(vars, run_univ) %>%
#   bind_rows()
# 
# # Format the table for publication (combining OR and CI)
# table_ready <- univ_results %>%
#   mutate(
#     # Rounds and formats as: OR (Lower, Upper) e.g., 1.45 (1.12, 1.89)
#     OR_CI = sprintf("%.2f (%.2f, %.2f)", OR, Lower_CI, Upper_CI),
#     estimate = round(estimate, 3)
#   ) %>%
#   select(Variable, term, estimate, OR_CI, p.value)
# 
# # Create a professional Word Table (Flextable)
# ft_regression <- flextable(table_ready) %>%
#   # Merge variable labels vertically for a clean presentation
#   merge_v(j = ~ Variable) %>%
#   # Use a formal academic layout theme
#   theme_vanilla() %>%
#   # Rename column headers to clear, publication-grade syntax
#   set_header_labels(
#     Variable = "Variable",
#     term = "Term / Group Level",
#     estimate = "Log-Odds Estimate",
#     OR_CI = "Odds Ratio (95% CI)",
#     p.value = "p-value"
#   ) %>%
#   # Formatting text styling and fonts
#   bold(part = "header") %>%
#   color(color = "#1F4E78", part = "header") %>% # Deep Navy Blue
#   font(fontname = "Calibri", part = "all") %>%
#   fontsize(size = 10, part = "all") %>%
#   # Alignment: Text labels left-aligned, statistical metrics right-aligned
#   align(j = c("Variable", "term"), align = "left", part = "all") %>%
#   align(j = c("estimate", "OR_CI", "p.value"), align = "right", part = "all") %>%
#   # Precision formatting for floats
#   colformat_double(j = c("estimate"), digits = 3) %>%
#   colformat_double(j = c("p.value"), digits = 3) %>%
#   # Highlight significant findings: bold p-values if less than 0.05
#   bold(i = ~ p.value < 0.05, j = ~ p.value) %>%
#   autofit()
# # Initialize Word Document and export to Results folder
# dir.create("Results", showWarnings = FALSE) # Assures folder exists
# 
# doc_regression <- read_docx() %>%
#   body_add_par("Statistical Analysis Report: Univariate Logistic Regression Models", style = "heading 1") %>%
#   body_add_par("Generated automatically from consolidated R workspace metrics.", style = "Normal") %>%
#   body_add_par("", style = "Normal") %>% # Fixed paragraph element spacer
#   body_add_flextable(ft_regression) %>%
#   body_add_par("Note: Bold values highlight statistical significance at the p < 0.05 tier. The baseline reference level intercepts are omitted from the layout.", style = "Table Caption")
# 
# # Save file 
# print(doc_regression, target = "Results/univariate_regression_results.docx")
########################################################
vars <- c(
  "GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
  "Prescription_After_Lab", "Owns_Farm_Animals", "OWNS_Poultry", "Animals_Treated",
  "Caretaking_Animals", "Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet"
)

# Ensure Case is numeric (0/1) and predictors are factors to reliably pull reference levels
data_model <- data_model %>%
  mutate(
    Case = as.numeric(as.character(Case)),
    across(all_of(vars), as.factor)
  )

# Function to run model, extract 95% CI, and include baseline reference levels
run_univ <- function(var) {
  formula <- as.formula(paste("Case ~", var))
  model <- glm(formula, data = data_model, family = binomial)
  
  #  Extract the baseline reference level for this factor
  ref_level <- levels(data_model[[var]])[1]
  
  # Create a reference row tibble
  ref_row <- tibble(
    Variable = var,
    term = paste0(var, ref_level),
    estimate = 0,
    OR = 1,
    Lower_CI = 1,
    Upper_CI = 1,
    p.value = NA_real_
  )
  
  #  Extract model results and combine with the reference row
  broom::tidy(model, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%   # remove mathematical intercept
    mutate(
      Variable = var,
      OR = exp(estimate),
      Lower_CI = exp(conf.low),          # Exponentiate lower bound
      Upper_CI = exp(conf.high)          # Exponentiate upper bound
    ) %>%
    select(Variable, term, estimate, OR, Lower_CI, Upper_CI, p.value) %>%
    # Bind reference row at the top of this variable's block
    bind_rows(ref_row, .) 
}

# Run all models with CIs and Reference levels included
univ_results <- lapply(vars, run_univ) %>%
  bind_rows()

# Clean up term labels for presentation (e.g., "GENDERMale" becomes "Male")
univ_results <- univ_results %>%
  mutate(
    term = mapply(function(v, t) gsub(paste0("^", v), "", t), Variable, term)
  )

# Format the table for publication (combining OR and CI, adding stars)
table_ready <- univ_results %>%
  mutate(
    # Format stars based on significance tier
    star = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    # If it's a reference group, show just 1.00 (omitting the word "Reference")
    OR_CI = ifelse(OR == 1 & is.na(p.value), 
                   "1.00", 
                   sprintf("%.2f (%.2f, %.2f)", OR, Lower_CI, Upper_CI)),
    # Format p-value with its stars, keeping text clean if NA
    p_fmt = ifelse(is.na(p.value), "—", sprintf("%.3f%s", p.value, star)),
    estimate_fmt = ifelse(is.na(p.value), "0.000", sprintf("%.3f", estimate))
  ) %>%
  select(Variable, term, estimate_fmt, OR_CI, p_fmt, p.value)

# Create a professional Word Table (Flextable)
ft_regression <- flextable(table_ready) %>%
  # Merge variable labels vertically for a clean presentation
  merge_v(j = ~ Variable) %>%
  # Use a formal academic layout theme
  theme_vanilla() %>%
  # Rename column headers to clear, publication-grade syntax
  set_header_labels(
    Variable = "Variable",
    term = "Term / Group Level",
    estimate_fmt = "Log-Odds Estimate",
    OR_CI = "Odds Ratio (95% CI)",
    p_fmt = "p-value"
  ) %>%
  # Formatting headers and fonts
  bold(part = "header") %>%
  color(color = "#1F4E78", part = "header") %>% # Deep Navy Blue
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  # Alignment: Text labels left-aligned, statistical metrics right-aligned
  align(j = c("Variable", "term"), align = "left", part = "all") %>%
  align(j = c("estimate_fmt", "OR_CI", "p_fmt"), align = "right", part = "all") %>%
  # Hide the raw p.value helper column from the final layout
  void(j = ~ p.value, part = "all") %>%
  autofit()

# Initialize Word Document and export to Results folder
dir.create("Results", showWarnings = FALSE) # Assures folder exists

doc_regression <- read_docx() %>%
  body_add_par("Statistical Analysis Report: Univariate Logistic Regression Models", style = "heading 1") %>%
  body_add_par("Generated automatically from consolidated R workspace metrics.", style = "Normal") %>%
  body_add_par("", style = "Normal") %>% 
  body_add_flextable(ft_regression) %>%
  body_add_par("Note: Significance tiers are indicated by stars (*p < 0.05, **p < 0.01, ***p < 0.001). First level listed for each variable (OR = 1.00) represents the baseline reference group.", style = "Table Caption")

# Save file 
print(doc_regression, target = "Results/univariate_regression_results1.docx")
###########################################################
#Check for dependencies before doing the multivariate analyses
# Conditional variables not included_ becuase they are clarly not independent and the have also many missing informations, "Withdrawal_Time","Prescription_After_Lab","Animals_Treated","Caretaking_Animals",
vars <- c(
  "GENDER","REGION","ORIGIN_OF_SAMPLE","Heard_AMR",
  "Used_AB_EVER","Used_AB_SIXM",
  "Owns_Farm_Animals","OWNS_Poultry",
  "Has_Toilet","Wash_Hands_Toilet"
)

check_assoc <- function(v1, v2) {
  tab <- table(data_model[[v1]], data_model[[v2]])
  
 # test <- chisq.test(tab) # MMJ 27/5:This gave a warning message:Chi-squared approximation may be incorrect; to rectify it I used the following function
  
  test <- chisq.test(tab, simulate.p.value = TRUE, B = 2000)
  cv <- cramerV(tab)
  
  data.frame(
    Var1 = v1,
    Var2 = v2,
    p_value = test$p.value,
    Cramers_V = cv
  )
}

# all pairwise combinations
pairs <- combn(vars, 2, simplify = FALSE)

assoc_results <- lapply(pairs, function(x) check_assoc(x[1], x[2])) %>%
  bind_rows()

assoc_results 

############################################
# Create a publishable table 3
# Pairwise Chi-Squared Tests & Cramer's V Strength of Association

# Ensure the Results folder exists
if(!dir.exists("Results")) dir.create("Results")

# Format the data frame for presentation
table_data <- assoc_results %>%
  mutate(
    # Format p-value to standard journal guidelines
    p_value_fmt = ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value)),
    Cramers_V_fmt = sprintf("%.4f", Cramers_V),
    # Add a formal interpretation column based on standard thresholds
    Interpretation = case_when(
      p_value < 0.05 & Cramers_V >= 0.30 ~ "Significant (Strong)",
      p_value < 0.05 & Cramers_V >= 0.10 ~ "Significant (Moderate)",
      p_value < 0.05 ~ "Significant (Weak)",
      TRUE ~ "Not Significant"
    )
  ) %>%
  select(Var1, Var2, p_value_fmt, Cramers_V_fmt, Interpretation)

# Build a high-quality publication-ready table layout
ft <- flextable(table_data) %>%
  set_header_labels(
    Var1 = "Variable 1",
    Var2 = "Variable 2",
    p_value_fmt = "p-value",
    Cramers_V_fmt = "Cramer's V",
    Interpretation = "Interpretation"
  ) %>%
  # Apply formal styles (Professional font size, alignment)
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(j = c("Var1", "Var2", "Interpretation"), align = "left", part = "all") %>%
  align(j = c("p_value_fmt", "Cramers_V_fmt"), align = "right", part = "all") %>%
  
  # Explicitly named padding parameters for flextable
  padding(padding.top = 5, padding.bottom = 5, 
          padding.left = 6, padding.right = 6, part = "all") %>%
  
  # Theme handling alternating row backgrounds properly
  theme_zebra(
    odd_header = "#003366", even_header = "#003366",
    odd_body = "#F2F5F8", even_body = "#FFFFFF"
  ) %>%
  
  # Design visual accents: Deep corporate blue headers and white header text
  color(color = "white", part = "header") %>%
  bold(part = "header") %>%
  
  # Visual highlights: Bold cells that achieve a statistically significant p-value
  bold(i = ~ as.numeric(gsub("<", "", p_value_fmt)) < 0.05 | p_value_fmt == "<0.001", 
       j = c("p_value_fmt", "Cramers_V_fmt", "Interpretation")) %>%
  
  # Set professional column width constraints
  autofit()

# Define custom professional text formatting for our titles
title_style <- fp_text(font.size = 22, bold = TRUE, color = "#003366", font.family = "Calibri")
subtitle_style <- fp_text(font.size = 12, italic = TRUE, color = "#4682B4", font.family = "Calibri")

# Initialize the formal Word Document structure
doc <- read_docx() %>%
  # CORRECTED: Custom styled Main Title & Subtitle using fpar + ftext
  body_add_fpar(fpar(ftext("Statistical Association Analysis Report", prop = title_style))) %>%
  body_add_fpar(fpar(ftext("Pairwise Chi-Squared Tests & Cramer's V Strength of Association", prop = subtitle_style))) %>%
  body_add_par("", style = "Normal") %>% # Empty spacer paragraph
  
  # Executive brief text block (Using recognized heading 2 style)
  body_add_par("1. Methodology & Overview", style = "heading 2") %>%
  body_add_par(
    "This document summarizes the pairwise associations evaluated across the ten core categorical variables. Pearson's Chi-Squared test of independence was performed using Monte Carlo simulated p-values (B = 2,000 simulations) to correct for cell dispersion issues and structural sparsity. Cramer's V metrics measure the strength of relationship effect size.",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>% # Spacer
  
  # Dynamic Data Table Insertion
  body_add_par("2. Pairwise Association Table", style = "heading 2") %>%
  body_add_par("Table 1: Matrix of associations, significance scores, and descriptive effect sizes.", style = "Normal") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>% # Spacer
  body_add_par("* Note: Rows with significant relationships (p < 0.05) emphasize analytical targets in bold typeface.", style = "Normal")

# Save file 
print(doc, target = "Results/Pairwise_Associations_Report.docx")
###############################################
## Cleaner table 3

# Ensure the output directory exists
if(!dir.exists("Results")) dir.create("Results")

# 1. Format the data frame for manuscript presentation
table_data <- assoc_results %>%
  mutate(
    # Determine association strength stars based on significance and effect size
    star = case_when(
      p_value < 0.05 & Cramers_V >= 0.30 ~ "***", # Strong
      p_value < 0.05 & Cramers_V >= 0.10 ~ "**",  # Moderate
      p_value < 0.05                     ~ "*",   # Weak
      TRUE                               ~ ""     # Not Significant
    ),
    # Merge Cramer's V and its significance star into a single column
    Cramers_V_fmt = sprintf("%.3f%s", Cramers_V, star),
    
    # Standardize p-values to 3 decimal places per standard journal guidelines
    p_value_fmt = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  ) %>%
  select(Var1, Var2, p_value_fmt, Cramers_V_fmt)

# 2. Build a high-quality, APA-style manuscript table
ft <- flextable(table_data) %>%
  set_header_labels(
    Var1 = "Variable 1",
    Var2 = "Variable 2",
    p_value_fmt = "p-value",
    Cramers_V_fmt = "Cramer's V"
  ) %>%
  theme_booktabs() %>% # Clean horizontal lines, no vertical lines
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(j = c("Var1", "Var2"), align = "left", part = "all") %>%
  align(j = c("p_value_fmt", "Cramers_V_fmt"), align = "right", part = "all") %>%
  padding(padding.top = 4, padding.bottom = 4, 
          padding.left = 6, padding.right = 6, part = "all") %>%
  color(color = "#333333", part = "header") %>%
  bold(part = "header") %>%
  autofit()

# Define text properties for document headers
title_style    <- fp_text(font.size = 18, bold = TRUE, color = "#1F4E78", font.family = "Calibri")
subtitle_style <- fp_text(font.size = 11, italic = TRUE, color = "#595959", font.family = "Calibri")

# 3. Initialize and construct the Word Document
doc <- read_docx() %>%
  body_add_fpar(fpar(ftext("Statistical Association Analysis Report", prop = title_style))) %>%
  body_add_fpar(fpar(ftext("Pairwise Chi-Squared Tests & Cramer's V Strength of Association", prop = subtitle_style))) %>%
  body_add_par("", style = "Normal") %>% 
  
  body_add_par("1. Methodology & Overview", style = "heading 2") %>%
  body_add_par(
    "This document summarizes the pairwise associations evaluated across the categorical variables. Pearson's Chi-Squared test of independence was performed using Monte Carlo simulated p-values (B = 2,000 simulations) to correct for cell dispersion issues and structural sparsity. Cramer's V metrics measure the strength of relationship effect size.",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>% 
  
  body_add_par("2. Pairwise Association Table", style = "heading 2") %>%
  body_add_par("Table 1: Matrix of pairwise associations, significance scores, and descriptive effect sizes.", style = "Table Caption") %>%
  body_add_flextable(ft) %>%
  body_add_par("", style = "Normal") %>% 
  
  # RECTIFIED: Changed style from 'Table Footnote' to 'Table Caption' to prevent the crash
  body_add_par(
    "Note: Association strength tiers are indicated via Cramer's V coefficients using asterisks based on statistical significance (p < 0.05): * = Weak relationship (Cramer's V < 0.10), ** = Moderate relationship (Cramer's V 0.10 to < 0.30), and *** = Strong relationship (Cramer's V ≥ 0.30).", 
    style = "Table Caption"
  )

# Save the document
print(doc, target = "Results/Pairwise_Associations_Report1.docx")
############################################################

# Now I do dare to do it automatically...

# Variables to test
vars <- c(
  "GENDER","REGION","ORIGIN_OF_SAMPLE","Heard_AMR",
  "Used_AB_EVER","Used_AB_SIXM",
  "Owns_Farm_Animals","OWNS_Poultry",
 "Has_Toilet","Wash_Hands_Toilet"
)

# Ensure Case is numeric
data_model <- data_model %>%
  mutate(Case = as.numeric(as.character(Case)))

# --- UNIVARIATE PIPELINE ---
run_univ <- function(var) {
  formula <- as.formula(paste("Case ~", var))
  model <- glm(formula, data = data_model, family = binomial)
  
  broom::tidy(model) %>%
    filter(term != "(Intercept)") %>%
    mutate(Variable = var) %>%
    select(Variable, term, estimate, p.value)
}

univ_results <- lapply(vars, run_univ) %>%
  bind_rows()
# Write out as a csv or table for the paper.
##################################################
# Creating Table 3 for publication

# Create Results folder if it doesn't exist
if(!dir.exists("Results")) { dir.create("Results") }

# Variables 
vars <- c(
  "GENDER","REGION","ORIGIN_OF_SAMPLE","Heard_AMR",
  "Used_AB_EVER","Used_AB_SIXM",
  "Owns_Farm_Animals","OWNS_Poultry",
  "Has_Toilet","Wash_Hands_Toilet"
)

data_model <- data_model %>%
  mutate(Case = as.numeric(as.character(Case)))

#  UNIVARIATE PIPELINE 
run_univ = function(var) {
  formula <- as.formula(paste("Case ~", var))
  model <- glm(formula, data = data_model, family = binomial)
  
  # Calculate Odds Ratios and 95% CI alongside estimate
  broom::tidy(model, exponentiate = FALSE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Variable = var,
      Odds_Ratio = exp(estimate),
      CI_Lower = exp(conf.low),
      CI_Upper = exp(conf.high)
    ) %>%
    select(Variable, term, estimate, std.error, Odds_Ratio, CI_Lower, CI_Upper, p.value)
}

# Run pipeline across variables
univ_results <- lapply(vars, run_univ) %>% bind_rows()

#  EXPORT TO CSV 
write_csv(univ_results, "Results/univariate_analysis_results.csv")


#  EXPORT TO PUBLICATION-READY WORD TABLE 
# Clean terms for professional appearance
table_data <- univ_results %>%
  mutate(
    # Save a numeric raw p-value column specifically for flextable filtering
    raw_p = p.value, 
    term = mapply(function(v, t) gsub(v, "", t), Variable, term),
    `Odds Ratio (95% CI)` = sprintf("%.2f (%.2f–%.2f)", Odds_Ratio, CI_Lower, CI_Upper),
    `p-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    estimate = round(estimate, 3)
  ) %>%
  select(Variable, `Term/Level` = term, `Log Odds` = estimate, `Odds Ratio (95% CI)`, `p-value`, raw_p)

# Build APA format flextable with bolding condition
ft <- flextable(table_data, col_keys = c("Variable", "Term/Level", "Log Odds", "Odds Ratio (95% CI)", "p-value")) %>%
  merge_v(j = ~ Variable) %>%     # Merges repeating variable rows for a clean look
  valign(j = ~ Variable, valign = "top") %>%
  theme_booktabs() %>%            # Applies thin top/bottom publication lines
  
  # Highlight significant findings: bold p-values if less than 0.05
  bold(i = ~ raw_p < 0.05, j = ~ `p-value`) %>%
  
  autofit() %>%
  set_caption("Table 1: Univariate Logistic Regression Analysis of Case Status Predictors") %>%
  add_footer_lines("Note: Models built via generalized linear model (binomial family, logit link). Bold p-values indicate significance at p < 0.05.")

# Save table to a Word file
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "Results/univariate_analysis_table.docx")
####################################################

# --- SELECT VARIABLES (p < 0.20) ---
vars_selected <- univ_results %>%
  filter(p.value < 0.20) %>%
  pull(Variable) %>%
  unique()

vars_selected
# !There are six variables selected out of which I rather like to keep those that makes more biologigal sense...
# I will first offer one at the time to the model
##############!!!!!!!!!!!!!!!_ NOTE All below has to be run agian and interpreteded once more as the dataset was changed becuase I had to change the criteria for classifying the Cases! We need to check that the dataset used are as we like it to be!

##################################################

model_Mv1 <- glm(Case ~ Heard_AMR+Used_AB_SIXM,
                         data = data_model,
                         family = binomial)

exp(cbind(OR = coef(model_Mv1), confint(model_Mv1)))

summary(model_Mv1)

model_Mv2 <- glm(Case ~ Heard_AMR+Used_AB_SIXM + Owns_Farm_Animals,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv2), confint(model_Mv2)))

summary(model_Mv2)

model_Mv3 <- glm(Case ~ Heard_AMR+Used_AB_SIXM + Animals_Treated ,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv3), confint(model_Mv3)))

summary(model_Mv3)
# This is not better than the Mv2

model_Mv4 <- glm(Case ~ Heard_AMR+Used_AB_SIXM + ORIGIN_OF_SAMPLE ,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv4), confint(model_Mv4)))

summary(model_Mv4)

#Has Toilet NS- to be removed

model_Mv4 <- glm(Case ~ Heard_AMR+Used_AB_SIXM + ORIGIN_OF_SAMPLE+Used_AB_EVER,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv4), confint(model_Mv4)))

summary(model_Mv4)


AIC(model_Mv4, model_Mv3, model_Mv2, model_Mv1)


prob <- predict(model_Mv4, type = "response")
roc_curve <- roc(data_model$Case, prob)

auc(roc_curve)

roc1 <- roc(data_model$Case, predict(model_Mv4, type="response"))
roc2 <- roc(data_model$Case, predict(model_Mv2, type="response"))

auc(roc1)

vif(model_Mv4)


exp(cbind(
  OR = coef(model_Mv4),
  confint(model_Mv4)
))


auc(roc2)


# This is to check what happens if we remove records where the content is missing

table(data_model$Used_AB_SIXM, data_model$Case)
ls()


data_model_clean <- data_model %>%
  filter(
    Used_AB_SIXM != "Missing",
    Owns_Farm_Animals != "Missing"
  )

# Here we need to model the forward selection process once more from MV1 to MV4 or 5 as long as something is significant in the model it should be kept, otherwise excluded.
model_Mv4_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE+Used_AB_EVER,
  data = data_model_clean,
  family = binomial
)

summary(model_Mv4_clean)

model_Mv3_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE+Used_AB_EVER+Has_Toilet,
  data = data_model_clean,
  family = binomial
)

summary(model_Mv3_clean)
# 29.05.26_ This model seems to be the best fit.

AIC(model_Mv4_clean, model_Mv3_clean)

table(data_model_clean$Case)

sum(is.na(prob))
sum(is.na(data_model_clean$Case))



prob <- predict(model_Mv3_clean, type = "response")
roc_curve <- roc(data_model_clean$Case, prob)

auc(roc_curve)

roc1 <- roc(data_model_clean$Case, predict(model_Mv3_clean, type="response"))
roc2 <- roc(data_model_clean$Case, predict(model_Mv4_clean, type="response"))

auc(roc1)

vif(model_Mv3)


exp(cbind(
  OR = coef(model_Mv4),
  confint(model_Mv4)
))

exp(cbind(
  OR = coef(model_Mv3_clean),
  confint(model_Mv3_clean)
))


# We would like to see the Risks so need to change the levels: You need to do this for each of the variables included in the final model

data_model_clean$Heard_AMR <- relevel(data_model_clean$Heard_AMR, ref = "Yes")
data_model_clean$Used_AB_EVER <- relevel(data_model_clean$Used_AB_EVER, ref = "Yes")
data_model_clean$Has_Toilet <- relevel(data_model_clean$Has_Toilet, ref = "Yes")
data_model_clean$ORIGIN_OF_SAMPLE <- relevel(data_model_clean$ORIGIN_OF_SAMPLE, ref = "Schoolchildren")


model_Mv3_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER + Has_Toilet,
  data = data_model_clean,
  family = binomial
)


summary(model_Mv3_clean)

exp(cbind(
  OR = coef(model_Mv3_clean),
  confint(model_Mv3_clean)
))

# When checking the 95% CI of the OR for the model_Mv3_clean the Used_AB_EVER- had a confodence interval from below 1 to slightly over 1 , which makes me wonder whether it should be removed from the model as it will not be a real risk factor
# Try below without including

model_Mv5_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE +  Has_Toilet,
  data = data_model_clean,
  family = binomial
)


summary(model_Mv5_clean)
#Now the Predictor Has_Toilet changes the coefficients and are no more significant!
#  Something is happening here...

exp(cbind(
  OR = coef(model_Mv5_clean),
  confint(model_Mv5_clean)
))


vif(model_Mv3_clean)
##################################
#After we had the meeting I just checked what happened if I remove one of the variables see above in green.

#I therefore do the full_model including all relevat variables and then exclude one by one
full_model <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER+ Used_AB_SIXM,
  data = data_model_clean,
  family = binomial
)
summary(full_model)

reduced_model <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER,
  data = data_model_clean,
  family = binomial
)
summary(reduced_model)

exp(cbind(
  OR = coef(reduced_model),
  confint(reduced_model)
))

####
#Clean up
#################################
#CLEAN VERSION
library(pROC)

model_Mv3_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER + Has_Toilet+ Used_AB_SIXM,
  data = data_model_clean,
  family = binomial
)

summary(model_Mv3_clean)

# Compare models (only if both exist)
AIC(model_Mv4_clean, model_Mv3_clean)

# Predictions
prob <- predict(model_Mv3_clean, type = "response")

# ROC + AUC
roc_curve <- roc(data_model_clean$Case, prob)
auc(roc_curve)

plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")  # reference line
legend("bottomright", legend = paste("AUC =", round(auc(roc_curve), 3)))
``
######################################################
# Creating a multivariate table (Table 4) and figure 1 for publication

# Create Results folder if it doesn't exist
if(!dir.exists("Results")) { dir.create("Results") }

# MODEL DEFINITION & PREDICTIONS 
model_Mv3_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER + Has_Toilet + Used_AB_SIXM,
  data = data_model_clean,
  family = binomial
)

# Get AIC (Handles comparison safely if model_Mv4_clean exists)
aic_val <- round(AIC(model_Mv3_clean), 2)
if(exists("model_Mv4_clean")) {
  aic_comp <- AIC(model_Mv4_clean, model_Mv3_clean)
  print(aic_comp)
}

# Predictions & AUC
prob <- predict(model_Mv3_clean, type = "response")
roc_curve <- roc(data_model_clean$Case, prob)
auc_val <- round(auc(roc_curve), 3)

# Save ROC Plot directly to the Results folder
png("Results/multivariate_roc_curve.png", width = 6, height = 6, units = "in", res = 300)
plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = paste("AUC =", auc_val), bty = "n")
dev.off()


#  TIDY RESULTS & EXTRACT METRICS 
# Exponentiate the coefficients to get Adjusted Odds Ratios (aOR)
multiv_results <- broom::tidy(model_Mv3_clean, exponentiate = FALSE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    aOR = exp(estimate),
    CI_Lower = exp(conf.low),
    CI_Upper = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, aOR, CI_Lower, CI_Upper, p.value)


# EXPORT TO CSV 
write_csv(multiv_results, "Results/multivariate_analysis_results.csv")


#  EXPORT TO PUBLICATION-READY WORD TABLE 
# Map variables to clean labels for the final manuscript table
 
table_data <- multiv_results %>%
  mutate(
    Variable = case_when(
      grepl("Heard_AMR", term) ~ "Heard of AMR",
      grepl("ORIGIN_OF_SAMPLE", term) ~ "Origin of Sample",
      grepl("Used_AB_EVER", term) ~ "Used Antibiotics Ever",
      grepl("Has_Toilet", term) ~ "Has Toilet",
      grepl("Used_AB_SIXM", term) ~ "Used Antibiotics (Past 6M)",
      TRUE ~ term
    ),
    `Term/Level` = gsub("Heard_AMR|ORIGIN_OF_SAMPLE|Used_AB_EVER|Has_Toilet|Used_AB_SIXM", "", term),
    `Term/Level` = ifelse(`Term/Level` == "" | `Term/Level` == term, "Yes/Variant", `Term/Level`),
    `Adjusted OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", aOR, CI_Lower, CI_Upper),
    `p-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    estimate = round(estimate, 3),
    raw_p = p.value 
  ) %>%
  select(Variable, `Term/Level`, `Log Odds` = estimate, `Adjusted OR (95% CI)`, `p-value`, raw_p)

# Build the styled APA/Journal Table 
ft <- flextable(table_data, col_keys = c("Variable", "Term/Level", "Log Odds", "Adjusted OR (95% CI)", "p-value")) %>%
  merge_v(j = ~ Variable) %>%     
  valign(j = ~ Variable, valign = "top") %>%
  theme_booktabs() %>%            
  bold(i = ~ raw_p < 0.05, j = ~ `p-value`) %>%
  autofit() %>%
  set_caption("Table 2: Multivariable Logistic Regression Model of Positive Case Predictors") %>%
  add_footer_lines(paste0("Note: Model fit using a generalized linear model (binomial family, logit link). ",
                          "Model Fit: AIC = ", aic_val, "; Discrimination: AUC = ", auc_val, ". ",
                          "Bold p-values indicate statistical significance at p < 0.05."))


# Save Table Only
doc_table <- read_docx() %>%
  body_add_flextable(ft)

print(doc_table, target = "Results/multivariate_analysis_table.docx")


#  Save ROC Figure Only
doc_figure <- read_docx() %>%
  # Label the Figure title at the top of the page
  body_add_par("Figure 1", style = "Normal") %>%
  body_add_par("Receiver Operating Characteristic (ROC) Curve of the Multivariate Predictive Model", style = "Normal") %>%
  
  # Embed the high-resolution image generated earlier
  body_add_img(
    src = "Results/multivariate_roc_curve.png", 
    width = 5.0,    
    height = 5.0,   
    pos = "after"
  ) %>%
  
  # Add the figure legend below the chart area
  body_add_par(paste0("Note: The blue solid line represents model performance (AUC = ", auc_val, 
                      "). The gray dashed line represents the reference line of chance (AUC = 0.50)."), 
               style = "Normal")

print(doc_figure, target = "Results/multivariate_analysis_figure.docx")
##################################################
# Table 4. Multivariable Logistic regression analysis and Fig.1
if(!dir.exists("Results")) { dir.create("Results") }

# Define model variables for factor treatment and reference parsing
model_vars <- c("Heard_AMR", "ORIGIN_OF_SAMPLE", "Used_AB_EVER", "Has_Toilet", "Used_AB_SIXM")

# Ensure predictors are factors in the underlying dataset to grab reference baselines
data_model_clean <- data_model_clean %>%
  mutate(across(all_of(model_vars), as.factor))

# --- MODEL DEFINITION & PREDICTIONS ---
model_Mv3_clean <- glm(
  Case ~ Heard_AMR + ORIGIN_OF_SAMPLE + Used_AB_EVER + Has_Toilet + Used_AB_SIXM,
  data = data_model_clean,
  family = binomial
)

# Get AIC
aic_val <- round(AIC(model_Mv3_clean), 2)
if(exists("model_Mv4_clean")) {
  aic_comp <- AIC(model_Mv4_clean, model_Mv3_clean)
  print(aic_comp)
}

# Predictions & AUC
prob <- predict(model_Mv3_clean, type = "response")
roc_curve <- roc(data_model_clean$Case, prob)
auc_val <- round(auc(roc_curve), 3)

# Save ROC Plot directly to the Results folder
png("Results/multivariate_roc_curve.png", width = 6, height = 6, units = "in", res = 300)
plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = paste("AUC =", auc_val), bty = "n")
dev.off()


# --- TIDY RESULTS & EXTRACT METRICS ---
multiv_results <- broom::tidy(model_Mv3_clean, exponentiate = FALSE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    aOR = exp(estimate),
    CI_Lower = exp(conf.low),
    CI_Upper = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, aOR, CI_Lower, CI_Upper, p.value)

# EXPORT RAW TO CSV
write_csv(multiv_results, "Results/multivariate_analysis_results.csv")


# --- GENERATE BASELINES & RESTRUCTURE MANUSCRIPT DATA ---
# Construct reference baseline rows dynamically for each covariate factor
reference_rows <- lapply(model_vars, function(v) {
  ref_level <- levels(data_model_clean[[v]])[1]
  tibble(
    Variable = case_when(
      v == "Heard_AMR" ~ "Heard of AMR",
      v == "ORIGIN_OF_SAMPLE" ~ "Origin of Sample",
      v == "Used_AB_EVER" ~ "Used Antibiotics Ever",
      v == "Has_Toilet" ~ "Has Toilet",
      v == "Used_AB_SIXM" ~ "Used Antibiotics (Past 6M)",
      TRUE ~ v
    ),
    `Term/Level` = ref_level,
    `Log Odds` = "0.000",
    `Adjusted OR (95% CI)` = "1.00",
    `p-value` = "—"
  )
}) %>% bind_rows()

# Clean and transform model estimates
table_model_data <- multiv_results %>%
  mutate(
    Variable = case_when(
      grepl("Heard_AMR", term) ~ "Heard of AMR",
      grepl("ORIGIN_OF_SAMPLE", term) ~ "Origin of Sample",
      grepl("Used_AB_EVER", term) ~ "Used Antibiotics Ever",
      grepl("Has_Toilet", term) ~ "Has Toilet",
      grepl("Used_AB_SIXM", term) ~ "Used Antibiotics (Past 6M)",
      TRUE ~ term
    ),
    `Term/Level` = gsub("Heard_AMR|ORIGIN_OF_SAMPLE|Used_AB_EVER|Has_Toilet|Used_AB_SIXM", "", term),
    `Term/Level` = ifelse(`Term/Level` == "" | `Term/Level` == term, "Yes/Variant", `Term/Level`),
    
    # 1. Add significance star notation directly to the p-value column
    star = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    `Adjusted OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", aOR, CI_Lower, CI_Upper),
    `p-value` = ifelse(p.value < 0.001, paste0("<0.001", star), sprintf("%.3f%s", p.value, star)),
    `Log Odds` = sprintf("%.3f", estimate)
  ) %>%
  select(Variable, `Term/Level`, `Log Odds`, `Adjusted OR (95% CI)`, `p-value`)

# Combine reference baselines and values together, keeping variables grouped
table_ready <- bind_rows(reference_rows, table_model_data) %>%
  arrange(factor(Variable, levels = c("Heard of AMR", "Origin of Sample", "Used Antibiotics Ever", "Has Toilet", "Used Antibiotics (Past 6M)")))


# --- EXPORT TO PUBLICATION-READY WORD TABLE ---
ft <- flextable(table_ready) %>%
  merge_v(j = ~ Variable) %>%     
  valign(j = ~ Variable, valign = "top") %>%
  theme_booktabs() %>%            
  autofit() %>%
  set_caption("Table 2: Multivariable Logistic Regression Model of Positive Case Predictors") %>%
  add_footer_lines(paste0("Note: Model fit using a generalized linear model (binomial family, logit link). ",
                          "Model Fit: AIC = ", aic_val, "; Discrimination: AUC = ", auc_val, ". ",
                          "Significance levels: *p < 0.05, **p < 0.01, ***p < 0.001. First level listed per variable (OR = 1.00) indicates baseline reference group."))

# Save Table Only
doc_table <- read_docx() %>%
  body_add_flextable(ft)

print(doc_table, target = "Results/multivariate_analysis_table1.docx")


# --- Save ROC Figure Only ---
doc_figure <- read_docx() %>%
  body_add_par("Figure 1", style = "Normal") %>%
  body_add_par("Receiver Operating Characteristic (ROC) Curve of the Multivariate Predictive Model", style = "Normal") %>%
  body_add_img(
    src = "Results/multivariate_roc_curve.png", 
    width = 5.0,    
    height = 5.0,   
    pos = "after"
  ) %>%
  body_add_par(paste0("Note: The blue solid line represents model performance (AUC = ", auc_val, 
                      "). The gray dashed line represents the reference line of chance (AUC = 0.50)."), 
               style = "Normal")

print(doc_figure, target = "Results/multivariate_analysis_figure1.docx")

##############################################################
#Wash_Hands_Toilet ns

# if we think we want to include the withdrawal time we loose to many observations, so it is better to keep th MV2!

