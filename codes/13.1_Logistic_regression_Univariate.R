#1. Simulate and set up an example Dataset
# Install required packages (run once)
install.packages(c("dplyr", "broom", "pROC", "ResourceSelection"))

# Load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(broom)
library(pROC)
library(ResourceSelection)
library(flextable) 
library(officer)

joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")
names(joined_data)
#  DATA PREPARATION
df_final <- joined_data %>%
  mutate(
    # Create the text label for species identification
    Confirmed_Species = case_when(
      (str_detect(Isolate, "E.coli|E. coli") & VITEK_MS_Results == "Escherichia coli" & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish") ~ "ESCR E. coli",
      (str_detect(Isolate, "K.pneumoniae") & VITEK_MS_Results == "Escherichia coli" & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") ~ "ESCR E. coli",
      (str_detect(Isolate, "K.pneumoniae") & VITEK_MS_Results == "Klebsiella pneumoniae" & `COLONY MORPHOLOGY ON C3GR` == "Metallic blue") ~ "ESCR K. pneumoniae",
      (str_detect(Isolate, "E.coli|E. coli") & VITEK_MS_Results == "Klebsiella pneumoniae" & `COLONY MORPHOLOGY ON C3GR` == "Pinkish/Reddish") ~ "ESCR K. pneumoniae",
      (str_detect(Isolate, "E.coli|E. coli") & Isolate_NVI == "E.coli" & `COLONY MORPHOLOGY ON CARBA` == "Pinkish/Reddish") ~ "CP E.coli",
      (str_detect(Isolate, "E.coli|E. coli") & Isolate_NVI == "Klebsiella pneumoniae" & `COLONY MORPHOLOGY ON CARBA` == "Pinkish/Reddish") ~ "CP K. pneumoniae",
      TRUE ~ "Negative/No Growth"
    ),
    # Create numeric Outcome (1 = Positive, 0 = Negative) for the regression
    Outcome_ESC_CP_Bact = if_else(Confirmed_Species == "Negative/No Growth", 0, 1)
  ) %>%
  # Select and Rename columns
  select(
    Outcome_ESC_CP_Bact,
    GENDER,
    REGION = REGION.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    Heard_AMR = "Have you ever heard about AMR?",
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    Animals_Treated = "Are your animals being treated for any diseases?",
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?",
    Antibiotics_use_lifetime = "Have you or your children used any antibiotics at any time?" ,
    Antibiotics_use_last_six_month = "Have you or your children used any antibiotics in the past six months?" ,
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    Taking_care_animals = "If Yes, do you participate in taking care of animals i.e. cleaning, Feeding?" ,
    Owns_poultry_farm = "Do you have a poultry farm at your home?",
    Taking_care_poultry = "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?",
    Poultry_Treated = "Are your Poultry being treated by veterinary doctor?" ,

  ) %>%
  # Convert to factors to ensure levels (Yes/No, etc.) display sequentially
  mutate(across(-Outcome_ESC_CP_Bact, ~as.factor(.x))) %>%
  # Create a factor version of the outcome for table display
  mutate(Outcome_Factor = factor(Outcome_ESC_CP_Bact, levels = c(1, 0), labels = c("ESCR_CP Positive", "Negative")))


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



data<-df_final
# Ensure outcome is factor for interpretation if desired
data$outcome <- factor(data$Outcome_ESC_CP_Bact, levels = c(0,1))

#2. Univariate screening
vars <- c("GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR",
         "Antibiotics_use_lifetime","Antibiotics_use_last_six_month",
         "Prescription_After_Lab", "Owns_Farm_Animals","Taking_care_animals",
         "Animals_Treated","Owns_poultry_farm","Taking_care_poultry",
         "Poultry_Treated","Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet")

# univ_results <- lapply(vars, function(v) {
#   formula <- as.formula(paste("outcome ~", v))
#   model <- glm(formula, data = data, family = binomial)
#   tidy(model) %>% mutate(variable = v)
# })
# 
# univ_table <- bind_rows(univ_results)
# univ_table
# 
# candidate_vars <- univ_table %>%
#   filter(term != "(Intercept)", p.value < 0.20) %>%
#   pull(variable) %>%
#   unique()
# 
# candidate_vars
# 
# full_formula <- as.formula(
#   paste("outcome ~", paste(candidate_vars, collapse = " + "))
# )
# 
# full_model <- glm(full_formula, data = data, family = binomial)
# 
# 
# #4. Forward stepwise selection (AIC-based)
# 
# null_model <- glm(outcome ~ 1, data = data, family = binomial)
# 
# forward_model <- step(
#   null_model,
#   scope = list(lower = null_model, upper = full_model),
#   direction = "forward",
#   trace = TRUE
# )
# 
# summary(forward_model)
# 
# summary(full_model)

#MMJ 14/5/2026
# Univariate Analysis
univ_results <- lapply(vars, function(v) {
  formula <- as.formula(paste("outcome ~", v))
  model <- glm(formula, data = data, family = binomial)
  tidy(model) %>% mutate(variable = v)
})

univ_table <- bind_rows(univ_results)

# Identify Candidate Variables (p < 0.20)
candidate_vars <- univ_table %>%
  filter(term != "(Intercept)", p.value < 0.20) %>%
  pull(variable) %>%
  unique()

# Create a dataset containing only the outcome and candidate predictors
# then remove any rows with missing values (NA).
analysis_vars <- c("outcome", candidate_vars)
clean_data <- data %>%
  select(all_of(analysis_vars)) %>%
  drop_na()

# Check how many rows were kept
print(paste("Rows in original data:", nrow(data)))
print(paste("Rows in clean data:", nrow(clean_data)))


# Define the Full Formula and Model
full_formula <- as.formula(
  paste("outcome ~", paste(candidate_vars, collapse = " + "))
)

# Use clean_data for both models
full_model <- glm(full_formula, data = clean_data, family = binomial)

# Forward stepwise selection (AIC-based)
null_model <- glm(outcome ~ 1, data = clean_data, family = binomial)

forward_model <- step(
  null_model,
  scope = list(lower = null_model, upper = full_model),
  direction = "forward",
  trace = TRUE
)

# Results
summary(forward_model)
summary(full_model)


# Madelaine_Test

# First model (1 predictor)
#model_1 <- glm(REGION ~ ORIGIN_OF_SAMPLE, 
        #       data = data, 
            #   family = binomial)
model_1 <- glm(REGION ~ 1, 
                     data = data, 
                family = binomial)
summary(model_1)

# Add second predictor
# model_2 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime, 
#                data = data, 
#                family = binomial)
model_2 <- glm(REGION ~ ORIGIN_OF_SAMPLE, 
                              data = data, 
                              family = binomial)


summary(model_2)

# model_2 <- glm(REGION~ 1,  
#                data = data, 
#                family = binomial)
# summary(model_2)
# Add 3rd predictor
model_3 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime, 
                               data = data, 
                               family = binomial)
summary(model_3)

# Add 4th predictor
model_4 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime + Owns_Farm_Animals, 
               data = data, 
               family = binomial)
summary(model_4)

# Add 5th predictor
model_5 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals, 
               data = data, 
               family = binomial)
summary(model_5)

# Add 6th predictor
model_6 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals + Animals_Treated, 
               data = data, 
               family = binomial)
summary(model_6)

# Add 7th predictor
model_7 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals + Animals_Treated +
                 Owns_poultry_farm, 
               data = data, 
               family = binomial)
summary(model_7)

# Add 8th predictor
model_8 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals + Animals_Treated +
                 Owns_poultry_farm + Taking_care_poultry, 
               data = data, 
               family = binomial)
summary(model_8)
# Add 9th predictor
model_9 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals + Animals_Treated +
                 Owns_poultry_farm + Taking_care_poultry + Poultry_Treated, 
               data = data, 
               family = binomial)
summary(model_9)

# Add 10th predictor
model_10 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                 Owns_Farm_Animals + Taking_care_animals + Animals_Treated +
                 Owns_poultry_farm + Taking_care_poultry + Poultry_Treated + 
                  Has_Toilet, 
               data = data, 
               family = binomial)
summary(model_10)

# Add 11th predictor
model_11 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime +
                  Owns_Farm_Animals + Taking_care_animals + Animals_Treated +
                  Owns_poultry_farm + Taking_care_poultry + Poultry_Treated + 
                  Has_Toilet + Wash_Hands_Toilet, 
                data = data, 
                family = binomial)
summary(model_11)

# Compare models

AIC(model_1, model_2,model_3, model_4, model_5, model_6, model_7, model_8,model_9,
    model_10, model_11)

##############################
#Below the further steps to be done after having several models that can be assessed and compared- need to check you have the right names on the dataset and 
#variables ( outcome)


# 6. Goodness-of-fit (Hosmer-Lemeshow test)
# Convert outcome back to numeric for test
hl <- hoslem.test(
  as.numeric(as.character(data$outcome)),
  fitted(forward_model),
  g = 10
)

hl

#Interpretation:

#p > 0.05 → good fit
#p < 0.05 → poor fit

#7. ROC curve and AUC
roc_obj <- roc(data$outcome, fitted(forward_model))

# AUC
auc(roc_obj)

# Plot ROC
plot(roc_obj, col = "blue", main = "ROC Curve")

#8. Optional compare ROC across models
roc_full <- roc(data$outcome, fitted(full_model))
roc_forward <- roc(data$outcome, fitted(forward_model))

plot(roc_full, col = "red", main = "ROC Comparison")
lines(roc_forward, col = "blue")

legend("bottomright", legend = c("Full", "Forward"),
       col = c("red", "blue"), lwd = 2)
#9 Odds ratios (interpretation)
exp(cbind(OR = coef(forward_model), confint(forward_model)))


#Summary of workflow


#Univariate screening
#→ Filter variables (e.g., p < 0.20)


#Build full model


#Forward selection (stepwise)
#→ optimize AIC


#Model comparison
#→ AIC()


#Model evaluation

#Hosmer–Lemeshow (calibration)
#ROC / AUC (discrimination)

#Notes (important in real research)

#Avoid relying solely on stepwise methods (they can overfit)
#Always consider:

#clinical relevance
#confounding
#interaction terms


#Validate model (cross-validation or external dataset)
################################################################################
# MMJ 14/5/2026
# DATA PREPARATION
data <- df_final
data$outcome <- factor(data$Outcome_ESC_CP_Bact, levels = c(0, 1))

# Define all possible predictors for the Univariate screening
all_potential_vars <- c(
  "GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR",
  "Antibiotics_use_lifetime", "Antibiotics_use_last_six_month",
  "Prescription_After_Lab", "Owns_Farm_Animals", "Taking_care_animals",
  "Animals_Treated", "Owns_poultry_farm", "Taking_care_poultry",
  "Poultry_Treated", "Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet"
)

#  Create a Complete Case Dataset 
# This ensures AIC comparisons are valid across all models
data_complete <- data %>%
  select(outcome, all_of(all_potential_vars)) %>%
  drop_na()

print(paste("Original rows:", nrow(data)))
print(paste("Rows kept for analysis:", nrow(data_complete)))


# FORWARD SELECTION

# Univariate screening (using data_complete)
univ_results <- lapply(all_potential_vars, function(v) {
  formula <- as.formula(paste("outcome ~", v))
  model <- glm(formula, data = data_complete, family = binomial)
  tidy(model) %>% mutate(variable = v)
})

univ_table <- bind_rows(univ_results)

# Identify Candidate Variables (p < 0.20)
candidate_vars <- univ_table %>%
  filter(term != "(Intercept)", p.value < 0.20) %>%
  pull(variable) %>%
  unique()

# Stepwise Selection
full_formula <- as.formula(paste("outcome ~", paste(candidate_vars, collapse = " + ")))
full_model   <- glm(full_formula, data = data_complete, family = binomial)
null_model   <- glm(outcome ~ 1, data = data_complete, family = binomial)

forward_model <- step(
  null_model,
  scope = list(lower = null_model, upper = full_model),
  direction = "forward",
  trace = TRUE
)

summary(forward_model)


# MODEL COMPARISON 

model_1  <- glm(REGION ~ 1, data = data_complete, family = binomial)
model_2  <- update(model_1, . ~ . + ORIGIN_OF_SAMPLE)
model_3  <- update(model_2, . ~ . + Antibiotics_use_lifetime)
model_4  <- update(model_3, . ~ . + Owns_Farm_Animals)
model_5  <- update(model_4, . ~ . + Taking_care_animals)
model_6  <- update(model_5, . ~ . + Animals_Treated)
model_7  <- update(model_6, . ~ . + Owns_poultry_farm)
model_8  <- update(model_7, . ~ . + Taking_care_poultry)
model_9  <- update(model_8, . ~ . + Poultry_Treated)
model_10 <- update(model_9, . ~ . + Has_Toilet)
model_11 <- update(model_10, . ~ . + Wash_Hands_Toilet)

# Compare all models 
aic_comparison <- AIC(model_1, model_2, model_3, model_4, model_5, 
                      model_6, model_7, model_8, model_9, model_10, model_11)

print(aic_comparison)

# Goodness-of-fit (Hosmer-Lemeshow test)
# Convert outcome back to numeric for test
hl <- hoslem.test(
  as.numeric(as.character(data$outcome)),
  fitted(forward_model),
  g = 10
)

hl

#Interpretation:

#p > 0.05 → good fit
#p < 0.05 → poor fit

# ROC curve and AUC
roc_obj <- roc(data$outcome, fitted(forward_model))

# AUC
auc(roc_obj)

# Plot ROC
plot(roc_obj, col = "blue", main = "ROC Curve")

# Optional compare ROC across models
roc_full <- roc(data$outcome, fitted(full_model))
roc_forward <- roc(data$outcome, fitted(forward_model))

plot(roc_full, col = "red", main = "ROC Comparison")
lines(roc_forward, col = "blue")

legend("bottomright", legend = c("Full", "Forward"),
       col = c("red", "blue"), lwd = 2)
# Odds ratios (interpretation)
exp(cbind(OR = coef(forward_model), confint(forward_model)))
################################################################################

data <- df_final
data$outcome <- factor(data$Outcome_ESC_CP_Bact, levels = c(0, 1))

all_potential_vars <- c(
  "GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR",
  "Antibiotics_use_lifetime", "Antibiotics_use_last_six_month",
  "Prescription_After_Lab", "Owns_Farm_Animals", "Taking_care_animals",
  "Animals_Treated", "Owns_poultry_farm", "Taking_care_poultry",
  "Poultry_Treated", "Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet"
)

# Create Complete Case Dataset (Ensures valid AIC/ROC comparisons)
data_complete <- data %>%
  select(outcome, all_of(all_potential_vars)) %>%
  drop_na()

message(paste("Rows in original data:", nrow(data)))
message(paste("Rows kept for analysis:", nrow(data_complete)))

# 2. UNIVARIATE SCREENING
univ_results <- lapply(all_potential_vars, function(v) {
  formula <- as.formula(paste("outcome ~", v))
  model <- glm(formula, data = data_complete, family = binomial)
  tidy(model) %>% mutate(variable = v)
})

univ_table <- bind_rows(univ_results)
candidate_vars <- univ_table %>% 
  filter(term != "(Intercept)", p.value < 0.20) %>% 
  pull(variable) %>% unique()

# 3. FORWARD STEPWISE SELECTION
full_formula <- as.formula(paste("outcome ~", paste(candidate_vars, collapse = " + ")))
full_model   <- glm(full_formula, data = data_complete, family = binomial)
null_model   <- glm(outcome ~ 1, data = data_complete, family = binomial)

forward_model <- step(null_model, 
                      scope = list(lower = null_model, upper = full_model),
                      direction = "forward", trace = TRUE)

# 4. MODEL DIAGNOSTICS
# Hosmer-Lemeshow Test
#hl <- hoslem.test(as.numeric(as.character(data_complete$outcome)), 
             #     fitted(forward_model), g = 10)
# # Use a safer number of bins (g)
# # Often for smaller datasets, g should be: number of variables + 1
# num_bins <- 16 
# 
# hl <- tryCatch({
#   hoslem.test(as.numeric(as.character(data_complete$outcome)), 
#               fitted(forward_model), g = num_bins)
# }, error = function(e) {
#   # If it still fails, return a list with a message
#   list(p.value = NA, method = "Hosmer-Lemeshow (Failed due to small sample/bins)")
# })
# 
# print(hl)
# 
# # ROC and AUC
# roc_forward <- roc(data_complete$outcome, fitted(forward_model))
# roc_full    <- roc(data_complete$outcome, fitted(full_model))
# auc_val     <- auc(roc_forward)
# 
# # Save ROC Plot as high-res TIFF
# tiff("Results/ROC_Comparison.tiff", width = 6, height = 6, units = 'in', res = 300)
# plot(roc_full, col = "red", main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))
# lines(roc_forward, col = "blue")
# legend("bottomright", legend = c("Full Model", "Forward Model"), col = c("red", "blue"), lwd = 2)
# dev.off()
# 
# # 5. PREPARE PUBLICATION TABLE (Odds Ratios)
# final_table_data <- tidy(forward_model, conf.int = TRUE, exponentiate = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   select(term, estimate, conf.low, conf.high, p.value) %>%
#   rename(
#     Variable = term,
#     `Adjusted OR` = estimate,
#     `Lower CI (95%)` = conf.low,
#     `Upper CI (95%)` = conf.high,
#     `p-value` = p.value
#   )
# 
# # Create Professional Flextable
# ft <- flextable(final_table_data) %>%
#   colformat_double(digits = 2) %>%
#   colformat_double(j = "p-value", digits = 3) %>%
#   bold(i = ~ `p-value` < 0.05, j = "p-value") %>%
#   set_caption("Table: Multivariable Logistic Regression of Factors Associated with ESC Resistance") %>%
#   theme_booktabs() %>%
#   autofit()
# 
# # 6. EXPORT TO WORD
# doc <- read_docx() %>%
#   body_add_par("AMR Statistical Analysis Report", style = "heading 1") %>%
#   body_add_par("Multivariable Model Results", style = "heading 2") %>%
#   body_add_flextable(ft) %>%
#   body_add_par(paste("Note: Variables selected via forward stepwise AIC from a candidate pool of p < 0.20."), style = "Normal") %>%
#   body_add_break() %>%
#   body_add_par("Model Diagnostics", style = "heading 2") %>%
#   body_add_par(paste("Hosmer-Lemeshow p-value:", round(hl$p.value, 4))) %>%
#   body_add_par(paste("Area Under the Curve (AUC):", round(auc_val, 3))) %>%
#   body_add_par(paste("AIC of Final Model:", round(AIC(forward_model), 2))) %>%
#   body_add_par(paste("Number of observations included:", nrow(data_complete)))
# 
# print(doc, target = "Results/AMR_Final_Analysis.docx")
# 
# message("All results successfully saved in the 'Results' folder.")
# 

# MODEL DIAGNOSTICS
# Use a safer number of bins (g)
num_bins <- 16 

hl <- tryCatch({
  hoslem.test(as.numeric(as.character(data_complete$outcome)), 
              fitted(forward_model), g = num_bins)
}, error = function(e) {
  list(p.value = NA, method = "Hosmer-Lemeshow (Failed due to small sample/bins)")
})

print(hl)

# ROC and AUC
roc_forward <- roc(data_complete$outcome, fitted(forward_model))
roc_full    <- roc(data_complete$outcome, fitted(full_model))
auc_val     <- auc(roc_forward)

# Save ROC Plot as high-res TIFF
tiff("Results/ROC_Comparison.tiff", width = 6, height = 6, units = 'in', res = 300)
plot(roc_full, col = "red", main = paste("ROC Curve (AUC =", round(auc_val, 3), ")"))
lines(roc_forward, col = "blue")
legend("bottomright", legend = c("Full Model", "Forward Model"), col = c("red", "blue"), lwd = 2)
dev.off()


# PREPARE PUBLICATION TABLE 
# Extract all initial candidate variables from the full model
all_variables <- tidy(full_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term)

# Extract only the variables that made it into the forward stepwise model
forward_variables <- tidy(forward_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value)

# Merge them together, flag selected ones with an asterisk, and clean names
final_table_data <- all_variables %>%
  left_join(forward_variables, by = "term") %>%
  mutate(
    # Add '*' to the variable name if it was successfully retained in the forward model
    Variable = if_else(!is.na(estimate), paste0(term, " *"), term)
  ) %>%
  select(Variable, estimate, conf.low, conf.high, p.value) %>%
  rename(
    `Adjusted OR` = estimate,
    `Lower CI (95%)` = conf.low,
    `Upper CI (95%)` = conf.high,
    `p-value` = p.value
  )


# Create Professional Flextable
ft <- flextable(final_table_data) %>%
  colformat_double(digits = 2, na_str = "-") %>%
  colformat_double(j = "p-value", digits = 3, na_str = "-") %>%
  # Bold the p-value only for significant variables that remained in the model
  bold(i = ~ `p-value` < 0.05, j = "p-value") %>%
  set_caption("Table: Multivariable Logistic Regression of Factors Associated with ESC Resistance") %>%
  theme_booktabs() %>%
  autofit()


# EXPORT TO WORD
doc <- read_docx() %>%
  body_add_par("AMR Statistical Analysis Report", style = "heading 1") %>%
  body_add_par("Multivariable Model Results", style = "heading 2") %>%
  body_add_flextable(ft) %>%
  body_add_par("Note: * Variables retained in the final model via forward stepwise AIC selection. Variables marked with '-' were dropped during selection.", style = "Normal") %>%
  body_add_break() %>%
  body_add_par("Model Diagnostics", style = "heading 2") %>%
  body_add_par(paste("Hosmer-Lemeshow p-value:", round(hl$p.value, 4))) %>%
  body_add_par(paste("Area Under the Curve (AUC):", round(auc_val, 3))) %>%
  body_add_par(paste("AIC of Final Model:", round(AIC(forward_model), 2))) %>%
  body_add_par(paste("Number of observations included:", nrow(data_complete)))

print(doc, target = "Results/AMR_Final_Analysis.docx")

message("All results successfully saved in the 'Results' folder.")

