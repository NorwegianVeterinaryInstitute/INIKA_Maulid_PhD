
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

# I was unable to import the above data set, due to change of directory
joined_data_16_3_26 <- read_csv("G:/.shortcut-targets-by-id/153ra6Su4lB1-ZcUTRshuT9CGPgtCu605/INIKA_OH_TZ/PhD_AREA/Maulid/INIKA_OH_TZ_ORIGINAL_DATA/NEW_MAULID/data/CLEANED_DATA/joined_data_16.3.26.csv")
names(joined_data_16_3_26)


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
  filter(Case != "EX") %>% 
  # Create numeric Outcome (1 = Positive, 0 = Negative) for the regression
  mutate(Case = as.numeric(Case)) %>%
  # Select and Rename columns
  select(
    Case,
    GENDER,
    REGION = REGION.x,
    SEASON = SEASON.x,
    DISTRICT = DISTRICT.x,
    Age_yrs,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    Heard_AMR = "Have you ever heard about AMR?",
    Used_AB_EVER = "Have you or your children used any antibiotics at any time?",                                        
    Used_AB_SIXM = "Have you or your children used any antibiotics in the past six months?",
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    OWNS_Poultry = "Do you have a poultry farm at your home?", 
    Animals_Treated = "Are your animals being treated for any diseases?",
    Caretaking_Animals = "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?",  
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?",
    OccupationalStatus = "What is your occupation and/or of your caretaker?*"
  ) %>%
  # Rename Ilemala to Ilemela before changing to factor
  mutate(DISTRICT = ifelse(DISTRICT == "Ilemala", "Ilemela", as.character(DISTRICT))) %>%
  
  # Recode OccupationalStatus into 4 major publication subgroups
  mutate(
    OccupationalStatus = str_trim(as.character(OccupationalStatus)), # Clean trailing/leading spaces
    OccupationalStatus = case_when(
      OccupationalStatus %in% c(
        "Bank Menager", "Clerk", "Construction engineer","Construction technician",
        "Electrical engineering","Engineer", 
        "Guard", "Laboratory technologist", "Mechanical engineering", 
        "Medical clinical officer", "Medical doctor", "Nurse", "Police", 
        "Soldier", "Tanesco worker", "Teacher", "Village executive officer"
      ) ~ "Formal employment",
      
      OccupationalStatus %in% c(
        "Barber", "Business man", "Business Man", "Business woman", "Business", 
        "Carpenter", "Entrepreneur", "Food vendor", "Local construction technician", 
        "Plumber", "Tailor", "Traditional healer"
      ) ~ "Informal/Self employment",
      
      OccupationalStatus %in% c(
        "Animal husbandry", "Carter", "Driver", "Fisherman", 
        "Motorcycle rider", "Peasant"
      ) ~ "Agriculture and Transportation",
      
      OccupationalStatus %in% c("Student", "Priest") ~ "Others",
      
      TRUE ~ as.character(OccupationalStatus) # Catch-all for missing or unmapped values
    )
  ) %>%
  
  # Convert to factors to ensure levels display sequentially
  mutate(across(-Case, ~as.factor(.x))) %>%
  
  # Create factor versions and derived demographic/age categories
  mutate(
    Outcome_Factor = factor(
      Case,
      levels = c(1, 0),
      labels = c("ESCR Positive", "Negative")
    ),
    Demographic = case_when(
      DISTRICT %in% c("Magu", "Hai", "Moshi Rural") ~ "Rural",
      TRUE ~ "Urban"
    ),
    Age_yrs = as.numeric(as.character(Age_yrs)),
    Age_group = case_when(
      Age_yrs < 13 ~ "Children",
      Age_yrs >= 13 & Age_yrs < 20 ~ "Teens",
      Age_yrs >= 20 & Age_yrs < 30 ~ "Young Adults",
      Age_yrs >= 30 & Age_yrs < 45 ~ "Adults",
      Age_yrs >= 45 & Age_yrs < 60 ~ "Middle Age",
      Age_yrs >= 60 ~ "Old",
      TRUE ~ "Check"
    ),
    Age_group = factor(
      Age_group, 
      levels = c("Children", "Teens", "Young Adults", "Adults", "Middle Age", "Old", "Check")
    )
  )
View(Original_joined)

# 20.06.26 Note the age_groups gives to small groups and are so correlated with the origin of the sampels so we will omit the age from any further analyses- see script "aldersfordeling_script.R"
# So it is fine to exclude Age_group from the vars, I think! Or it has to be grouped in another way! 

# Before doing anything get the frequencies of the answers from all the Independent variables


vars <- c("GENDER","REGION","SEASON",
          "DISTRICT",
          "Age_group", "ORIGIN_OF_SAMPLE","Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
          "Prescription_After_Lab","Owns_Farm_Animals","OWNS_Poultry","Animals_Treated",
          "Caretaking_Animals","Withdrawal_Time","Has_Toilet","Wash_Hands_Toilet", 
          "Demographic", "OccupationalStatus"
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
################################################################################
# Table 1. Frequencies in %
# Clean the data set by filtering out "Missing" entries for either controls or cases
data_model_clean <- Original_joined %>% 
  filter(
    Used_AB_SIXM != "Missing",
    Owns_Farm_Animals != "Missing",
    OWNS_Poultry != "Missing",
    Used_AB_EVER != "Unsure", 
    Used_AB_SIXM  != "Unsure"
  )

# Define the variables for frequency calculation
vars <- c("GENDER", "REGION", "DISTRICT", "Demographic", "SEASON", "Age_group",
          "ORIGIN_OF_SAMPLE", "Heard_AMR", "Used_AB_EVER", "Used_AB_SIXM",
          "Prescription_After_Lab", "Owns_Farm_Animals", "OWNS_Poultry", 
          "Animals_Treated","Caretaking_Animals", "Withdrawal_Time",
          "Has_Toilet", "Wash_Hands_Toilet", "OccupationalStatus")

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
  select(Variable, Stratum, Case_fmt, Control_fmt,Total)


# Create a professional, publishable Word Table 
ft <- flextable(table_ready) %>%
  # Set professional column headers
  set_header_labels(
    Variable = "Variable",
    Stratum = "Stratum",
    Case_fmt = "Case (N, %)",
    Control_fmt = "Control (N, %)",
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
  align(j = c("Case_fmt", "Control_fmt", "Total"), align = "right", part = "all") %>%
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
  body_add_par("Note: Missing values were treated as explicit 'NA' components. 
  Rows were omitted if 'Used_AB_SIXM', 'Owns_Farm_Animals' and 'Owns_Poultry' data was recorded as 'Missing'. 
               And 'Used_AB_SIXM' or 'Used_AB_Ever' data was recorded as 'Unsure'. 
               Percentages reflect proportions within column totals per variable.", style = "Table Caption")

# Save file (Fixed: Changed '%' to 'Percent')
print(doc, target = "Results/Frequencies_Percent_of_respondents_17.07.26.docx")
################################################################################

#  Creating a publishable table (Supplementary Table S2)
#Exclude the variables that are conditional ("Animals_Treated",   "Caretaking_Animals",Prescription_After_Lab, Withdrawal_Time )   in the analyses:
# Variables to include

vars <- c( "GENDER",
           "REGION",
           "DISTRICT",
           "Demographic",
           "SEASON", 
           "Age_group", 
           "ORIGIN_OF_SAMPLE",
           "Heard_AMR",
           "Used_AB_EVER",
           "Used_AB_SIXM",
           "Owns_Farm_Animals",
           "OWNS_Poultry",
           "Has_Toilet",
           "Wash_Hands_Toilet",
           "OccupationalStatus"
 
)

data_model <- data_model_clean  %>%
  mutate(
    Case = as.numeric(as.factor(Case)) - 1,
    
    #  convert numeric placeholders back to text labels, then to factor
    across(all_of(vars), ~ {
      char_var <- as.character(.x)
      mapped_var <- case_match(char_var,
                               "1" ~ "No",
                               "2" ~ "Yes",
                               "3" ~ "Unsure",
                               .default = char_var # Keep original labels if they are already text (e.g., GENDER, REGION)
      )
      
      #  Turn into a factor
      as.factor(mapped_var)
    })
  )

# Function to run model, extract 95% CI, and include baseline reference levels
run_univ <- function(var) {
  formula <- as.formula(paste("Case ~", var))
  model <- glm(formula, data = data_model, family = binomial)
  
  # Extract the baseline reference level for this factor
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
  
  # Extract model results and combine with the reference row
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
    # Format stars based on significance tier (including *† for p < 0.20)
    star = case_when(
      is.na(p.value)   ~ "",
      p.value < 0.001  ~ "***†", # p < 0.001 & eligible for multivariable
      p.value < 0.01   ~ "**†",  # p < 0.01 & eligible for multivariable
      p.value < 0.05   ~ "*†",   # p < 0.05 & eligible for multivariable
      round(p.value, 3) == 0.050 ~ "*†", 
      p.value < 0.20   ~ "*†",   # Screened for multivariable model inclusion (0.05 <= p < 0.20)
      TRUE             ~ ""
    ),
    # Format OR and CI layout cleanly
    OR_CI = ifelse(OR == 1 & is.na(p.value), 
                   "1.00", 
                   sprintf("%.2f [%.2f - %.2f]", OR, Lower_CI, Upper_CI)),
    
    # Custom formatting for small p-values to prevent displaying 0.000
    p_fmt = case_when(
      is.na(p.value)   ~ "—",
      p.value < 0.001  ~ paste0("<0.001", star),
      TRUE             ~ sprintf("%.3f%s", p.value, star)
    ),
    
    estimate_fmt = ifelse(is.na(p.value), "0.000", sprintf("%.3f", estimate))
  ) %>%
  select(Variable, term, estimate_fmt, OR_CI, p_fmt, p.value)

# Create a professional Word Table (Flextable)
ft_regression <- flextable(table_ready) %>%
  merge_v(j = ~ Variable) %>%
  theme_vanilla() %>%
  set_header_labels(
    Variable = "Variable",
    term = "Term / Group Level",
    estimate_fmt = "Log-Odds Estimate",
    OR_CI = "Odds Ratio [95% CI]",
    p_fmt = "p-value"
  ) %>%
  bold(part = "header") %>%
  color(color = "#1F4E78", part = "header") %>% 
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(j = c("Variable", "term"), align = "left", part = "all") %>%
  align(j = c("estimate_fmt", "OR_CI", "p_fmt"), align = "right", part = "all") %>%
  void(j = ~ p.value, part = "all") %>%
  autofit()

# Initialize Word Document and export to Results folder
dir.create("Results", showWarnings = FALSE) # Assures folder exists

doc_regression <- read_docx() %>%
  body_add_par("Statistical Analysis Report: Univariate Logistic Regression Models", style = "heading 1") %>%
  body_add_par("Generated automatically from consolidated R workspace metrics.", style = "Normal") %>%
  body_add_par("", style = "Normal") %>% 
  body_add_flextable(ft_regression) %>%
  body_add_par("Note: First level listed for each variable (OR = 1.00) represents the baseline reference group. Significance tiers: *p < 0.05, **p < 0.01, ***p < 0.001. † Screened for multivariable model inclusion (p < 0.20).", style = "Table Caption")
# Save file 
print(doc_regression, target = "Results/Univariate_regression_results_17.07.26.docx")
################################################################################
#Check for dependencies before doing the multivariate analyses
# Conditional variables not included_ because they are clearly not independent and the have also many missing informations, "Withdrawal_Time","Prescription_After_Lab","Animals_Treated","Caretaking_Animals",

vars <- c(
  "GENDER","REGION","ORIGIN_OF_SAMPLE","DISTRICT", "Demographic", "SEASON", "Age_group",
  "Heard_AMR","Used_AB_EVER","Used_AB_SIXM","Owns_Farm_Animals","OWNS_Poultry",
  "Has_Toilet","Wash_Hands_Toilet","OccupationalStatus" 
)
check_assoc <- function(v1, v2) {
  tab <- table(data_model[[v1]], data_model[[v2]])
  
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

################################################################################
## Supplementary Table S3
# Pairwise Chi-Squared Tests & Cramer's V Strength of Association
# Ensure the output directory exists
if(!dir.exists("Results")) dir.create("Results")

#  Format the data frame for manuscript presentation
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

#  Build a high-quality, APA-style manuscript table
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

#  Initialize and construct the Word Document
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
  body_add_par(
    "Note: Association strength tiers are indicated via Cramer's V coefficients using asterisks based on statistical significance (p < 0.05): 
    * = Weak relationship (Cramer's V < 0.10), 
    ** = Moderate relationship (Cramer's V 0.10 to < 0.30), and 
    *** = Strong relationship (Cramer's V ≥ 0.30).", 
    style = "Table Caption"
  )

# Printing out  a csv file as well:
write.csv(table_data,"Results/Pairwise_Cramers.csv")

# Save the document
print(doc, target = "Results/Pairwise Chi-Squared Tests & Cramer's V Strength of Association_17.07.26.docx")
################################################################################

# SELECT VARIABLES (p < 0.20) 
vars_selected <- univ_results %>%
  filter(p.value < 0.20) %>%
  pull(Variable) %>%
  unique()

vars_selected

# CHOOSE ONE VARIABLE AT THE TIME.

model_Mv1 <- glm(Case ~ DISTRICT + SEASON,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv1), confint(model_Mv1)))

summary(model_Mv1)

model_Mv1.1 <- glm(Case ~ SEASON + Demographic, 
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv1.1), confint(model_Mv1.1)))

summary(model_Mv1.1)

#( ! Note it is not possible to include both demographic and Distric seem district is a better predictor than the demographic!)



model_Mv2 <- glm(Case ~ DISTRICT + SEASON + ORIGIN_OF_SAMPLE,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv2), confint(model_Mv2)))

summary(model_Mv2)

model_Mv3 <- glm(Case ~ DISTRICT+ SEASON + ORIGIN_OF_SAMPLE + Heard_AMR,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv3), confint(model_Mv3)))

summary(model_Mv3)


model_Mv4 <- glm(Case ~ DISTRICT + SEASON + ORIGIN_OF_SAMPLE + Age_group + Heard_AMR,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv4), confint(model_Mv4)))

summary(model_Mv4)

model_Mv5 <- glm(Case ~ DISTRICT + SEASON + ORIGIN_OF_SAMPLE + Heard_AMR + Used_AB_SIXM,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv5), confint(model_Mv5)))

summary(model_Mv5)

# Had to correct the model_Mv6 from Mv5 in the code below! was wrong written!
# Model_Mv6 and model_mv7 are both useful!
model_Mv6 <- glm(Case ~ DISTRICT + SEASON + ORIGIN_OF_SAMPLE + Heard_AMR + Used_AB_SIXM + Has_Toilet,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv6), confint(model_Mv6)))

summary(model_Mv6)
vif(model_Mv6)

model_Mv7 <- glm(Case ~ DISTRICT + SEASON  + Heard_AMR + Used_AB_SIXM + Has_Toilet,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv7), confint(model_Mv7)))

summary(model_Mv7)

vif(model_Mv7)

model_Mv8 <- glm(Case ~ DISTRICT + SEASON  + Heard_AMR + Used_AB_SIXM + Has_Toilet+ Age_group,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv8), confint(model_Mv8)))

summary(model_Mv8)

# model_Mv6 seems to be the best since it has lowest AIC (711.4658) compared to other models
# and is the behavioral Predictors.
AIC(model_Mv8, model_Mv7, model_Mv6,model_Mv5,model_Mv3,model_Mv2,model_Mv1.1,model_Mv1 )


prob <- predict(model_Mv8, type = "response")
roc_curve <- roc(data_model$Case, prob)

auc(roc_curve)

roc1 <- roc(data_model$Case, predict(model_Mv8, type="response"))
roc2 <- roc(data_model$Case, predict(model_Mv7, type="response"))

auc(roc1)
auc(roc2)

vif(model_Mv8)
################################################################################

# From model Mv6
# Table 2. Multivariable Logistic regression analysis and Fig.1
if(!dir.exists("Results")) { dir.create("Results", showWarnings = FALSE) }

# Define model variables for factor treatment and reference parsing
model_vars <- c("Heard_AMR", "ORIGIN_OF_SAMPLE", "Used_AB_SIXM", "DISTRICT",
                "Has_Toilet", "SEASON") 

# Ensure predictors are factors in the underlying data set to grab reference baselines
data_model <- data_model %>%
  mutate(across(all_of(model_vars), as.factor))

# MODEL DEFINITION & PREDICTIONS 
model_Mv6 <- glm(Case ~ DISTRICT + SEASON + Heard_AMR + Used_AB_SIXM + Has_Toilet + ORIGIN_OF_SAMPLE,
                 data = data_model,
                 family = binomial
)

# Get AIC
aic_val <- round(AIC(model_Mv6), 2)

# Predictions & AUC calculation with 95% CI (via pROC)
prob <- predict(model_Mv6, type = "response")
roc_curve <- roc(data_model$Case, prob)

auc_val <- round(auc(roc_curve), 3)
auc_ci  <- ci.auc(roc_curve) # Calculates DeLong 95% CI

# Format string for display (e.g., "0.593 (95% CI: 0.548 – 0.638)")
auc_str <- sprintf("%.3f (95%% CI: %.3f – %.3f)", auc_val, auc_ci[1], auc_ci[3])

# Save ROC Plot directly to the Results folder
png("Results/multivariate_roc_curve.png", width = 6, height = 6, units = "in", res = 300)
plot(roc_curve, col = "blue", lwd = 2, main = "ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = paste("AUC =", auc_str), bty = "n", cex = 0.9)
dev.off()


# TIDY RESULTS & EXTRACT METRICS 
multiv_results <- broom::tidy(model_Mv6, exponentiate = FALSE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    aOR = exp(estimate),
    CI_Lower = exp(conf.low),
    CI_Upper = exp(conf.high)
  ) %>%
  select(term, estimate, std.error, aOR, CI_Lower, CI_Upper, p.value)

# EXPORT RAW TO CSV
write_csv(multiv_results, "Results/Multivariate_analysis_results.csv")


# GENERATE BASELINES & RESTRUCTURE MANUSCRIPT DATA 
# Construct reference baseline rows dynamically for each covariate factor
reference_rows <- lapply(model_vars, function(v) {
  ref_level <- levels(data_model[[v]])[1]
  tibble(
    Variable = case_when(
      v == "Heard_AMR" ~ "Heard of AMR",
      v == "ORIGIN_OF_SAMPLE" ~ "Origin of Sample",
      v == "Used_AB_SIXM"  ~ "Used antibiotics (Past 6 Months)",
      v == "DISTRICT" ~ "District",
      v == "SEASON" ~ "Season",
      v == "Has_Toilet" ~ "Has a toilet at home place",
      TRUE ~ v
    ),
    `Term/Level` = as.character(ref_level),
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
      grepl("Used_AB_SIXM", term) ~ "Used antibiotics (Past 6 Months)",
      grepl("DISTRICT", term) ~ "District",
      grepl("SEASON", term) ~ "Season",
      grepl("Has_Toilet", term) ~ "Has a toilet at home place",
      TRUE ~ term
    ),
    
    # Dynamically strips variable prefixes to leave clean stratum name
    `Term/Level` = term,
    `Term/Level` = gsub("^Heard_AMR|^ORIGIN_OF_SAMPLE|^DISTRICT|^Used_AB_SIXM|^SEASON|^Has_Toilet", "", `Term/Level`),
    
    # Add significance star notation directly to the p-value column
    star = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    `Adjusted OR (95% CI)` = sprintf("%.2f [%.2f – %.2f]", aOR, CI_Lower, CI_Upper),
    `p-value` = ifelse(p.value < 0.001, paste0("<0.001", star), sprintf("%.3f%s", p.value, star)),
    `Log Odds` = sprintf("%.3f", estimate)
  ) %>%
  select(Variable, `Term/Level`, `Log Odds`, `Adjusted OR (95% CI)`, `p-value`)

# Combine reference baselines and values together, keeping variables grouped
table_ready <- bind_rows(reference_rows, table_model_data) %>%
  arrange(factor(Variable, levels = c("Heard of AMR", "Origin of Sample", "Used antibiotics (Past 6 Months)",
                                      "District", "Season", "Has a toilet at home place"))) 


# EXPORT TO PUBLICATION-READY WORD TABLE 
ft <- flextable(table_ready) %>%
  merge_v(j = ~ Variable) %>%     
  valign(j = ~ Variable, valign = "top") %>%
  theme_booktabs() %>%            
  autofit() %>%
  set_caption("Table 2: Multivariable Logistic Regression Model of Positive Case Predictors") %>%
  add_footer_lines(paste0("Note: Model fit using a generalized linear model (binomial family, logit link). ",
                          "Model Fit: AIC = ", aic_val, "; Discrimination: AUC = ", auc_str, ". ",
                          "Significance levels: *p < 0.05, **p < 0.01, ***p < 0.001. First level listed per variable (OR = 1.00) indicates baseline reference group."))

# Save Table Only
doc_table <- read_docx() %>%
  body_add_flextable(ft)

print(doc_table, target = "Results/Multivariate_analysis_table2.docx")


# Save ROC Figure Only 
doc_figure <- read_docx() %>%
  body_add_par("Fig. 1 Receiver Operating Characteristic (ROC) curve of the multivariable predictive model", style = "Normal") %>%
  body_add_img(
    src = "Results/Multivariate_roc_curve.png", 
    width = 5.0,    
    height = 5.0,   
    pos = "after"
  ) %>%
  body_add_par(paste0("Note: The blue solid line represents overall model performance (AUC = ", auc_str, 
                      "). The gray dashed line represents the non-discrimination reference line (AUC = 0.500)."), 
               style = "Normal")

print(doc_figure, target = "Results/ROC_CURVE_figure1.docx")
################################################################################