#1. Simulate and set up an example Dataset
# Install required packages (run once)
install.packages(c("dplyr", "broom", "pROC", "ResourceSelection"))

# Load libraries
library(dplyr)
library(broom)
library(pROC)
library(ResourceSelection)


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
  mutate(across(-Outcome_ESC_Bact, ~as.factor(.x))) %>%
  # Create a factor version of the outcome for table display
  mutate(Outcome_Factor = factor(Outcome_ESC_Bact, levels = c(1, 0), labels = c("ESBL Positive", "Negative")))


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
data$outcome <- factor(data$Outcome_ESC_Bact, levels = c(0,1))

#2. Univariate screening
vars <- c("GENDER", "REGION", "ORIGIN_OF_SAMPLE", "Heard_AMR",
         "Antibiotics_use_lifetime","Antibiotics_use_last_six_month",
         "Prescription_After_Lab", "Owns_Farm_Animals","Taking_care_animals",
         "Animals_Treated","Owns_poultry_farm","Taking_care_poultry",
         "Poultry_Treated","Withdrawal_Time", "Has_Toilet", "Wash_Hands_Toilet")

univ_results <- lapply(vars, function(v) {
  formula <- as.formula(paste("outcome ~", v))
  model <- glm(formula, data = data, family = binomial)
  tidy(model) %>% mutate(variable = v)
})

univ_table <- bind_rows(univ_results)
univ_table

candidate_vars <- univ_table %>%
  filter(term != "(Intercept)", p.value < 0.20) %>%
  pull(variable) %>%
  unique()

candidate_vars

full_formula <- as.formula(
  paste("outcome ~", paste(candidate_vars, collapse = " + "))
)

full_model <- glm(full_formula, data = data, family = binomial)


#4. Forward stepwise selection (AIC-based)

null_model <- glm(outcome ~ 1, data = data, family = binomial)

forward_model <- step(
  null_model,
  scope = list(lower = null_model, upper = full_model),
  direction = "forward",
  trace = TRUE
)

summary(forward_model)

summary(full_model)

# Madelaine_Test

# First model (1 predictor)
model_1 <- glm(REGION ~ ORIGIN_OF_SAMPLE, 
               data = data, 
               family = binomial)

# Add second predictor
model_2 <- glm(REGION ~ ORIGIN_OF_SAMPLE + Antibiotics_use_lifetime+........., 
               data = data, 
               family = binomial)


summary(model_2)

model_2 <- glm(REGION~ 1,  
               data = data, 
               family = binomial)
summary(model_2)

# 5.Compare models

AIC(model_1, model_2)

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







