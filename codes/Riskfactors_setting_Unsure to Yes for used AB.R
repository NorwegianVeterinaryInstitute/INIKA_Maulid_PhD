
#17.07.26 checking what happens with the model if we recategorize unsure to Yes?

data_model_clean <- Original_joined %>% 
  filter(
    Used_AB_SIXM != "Missing",
    Owns_Farm_Animals != "Missing",
    OWNS_Poultry != "Missing",
   )%>%
  mutate(Used_AB_EVER= case_when(Used_AB_EVER=="Unsure" ~ "Yes",
                                 TRUE  ~ Used_AB_EVER))%>%
mutate(Used_AB_SIXM=case_when(Used_AB_SIXM =="Unsure" ~ "Yes",
                         TRUE  ~ Used_AB_SIXM))


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


#MMJ 27/5/2026 
#MN 03.07.2026 Note you need to use the data_model_clean
data_model <- data_model_clean %>%
  mutate(
    #  Convert to factor and explicitly turn NA into a "Missing" level
    #  Convert to numeric integers safely
    Used_AB_SIXM           = as.numeric(fct_na_value_to_level(as.factor(Used_AB_SIXM), "Missing")),
    Prescription_After_Lab = as.numeric(fct_na_value_to_level(as.factor(Prescription_After_Lab), "Missing")),
    Owns_Farm_Animals      = as.numeric(fct_na_value_to_level(as.factor(Owns_Farm_Animals), "Missing")),
    OWNS_Poultry           = as.numeric(fct_na_value_to_level(as.factor(OWNS_Poultry), "Missing")),
    Animals_Treated        = as.numeric(fct_na_value_to_level(as.factor(Animals_Treated), "Missing")),
    Caretaking_Animals     = as.numeric(fct_na_value_to_level(as.factor(Caretaking_Animals), "Missing")),
    Withdrawal_Time        = as.numeric(fct_na_value_to_level(as.factor(Withdrawal_Time), "Missing")),
    
    # Target variable 'Case' converted to 0 and 1 if it is binary, or standard numeric
    Case                   = as.numeric(as.factor(Case))
  )

################################################################################
# Better to use Case as numeric! avoids factor confusion- what is the reference level 

data_model <- data_model %>%
  mutate(Case = as.numeric(as.factor(Case)) - 1
  )

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

model_4 <- glm(Case ~ Used_AB_EVER,
               data = data_model,
               family = binomial)

exp(cbind(OR = coef(model_4), confint(model_4)))

summary(model_4)

####

#MN 03.07.26_ You should not include the variables that are conditional ("Animals_Treated",   "Caretaking_Animals",Prescription_After_Lab, Withdrawal_Time )   in the analyses:
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

# These are the 4 variables that have p<0.20 from the univariate analyses to be offered in the multivariate analyses
#"ORIGIN_OF_SAMPLE","Heard_AMR","Used_AB_EVER","Used_AB_SIXM",
model_Mv1 <- glm(Case ~ Heard_AMR+Used_AB_SIXM,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv1), confint(model_Mv1)))

summary(model_Mv1)

model_Mv2 <- glm(Case ~ Heard_AMR+Used_AB_SIXM + ORIGIN_OF_SAMPLE ,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv2), confint(model_Mv2)))

summary(model_Mv2)





model_Mv4 <- glm(Case ~ Heard_AMR + ORIGIN_OF_SAMPLE,
                 data = data_model,
                 family = binomial)

exp(cbind(OR = coef(model_Mv4), confint(model_Mv4)))

summary(model_Mv4)

#  The model_Mv4 IS THE BEST MODEL ACCORDING TO OUR DATASET.
#15.07.26 On line above I wrote that the model_Mv4 is the best model! ( including Heard_AMR+ ORIGIN_OF_SAMPLE + Used_AB_EVER)
# I leave the codes below, and make a separate check for possible onteraction terms afterwards!

AIC(model_Mv4, model_Mv2, model_Mv1)


prob <- predict(model_Mv4, type = "response")
roc_curve <- roc(data_model$Case, prob)

auc(roc_curve)

roc1 <- roc(data_model$Case, predict(model_Mv4, type="response"))
roc2 <- roc(data_model$Case, predict(model_Mv2, type="response"))

auc(roc1)
auc(roc2)

vif(model_Mv4)


