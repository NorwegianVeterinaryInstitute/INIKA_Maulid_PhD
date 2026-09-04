#17.07.26 checking what happens with the model if we recategorize unsure to Yes?

data_model_clean <- Original_joined %>% 
  filter(
    Used_AB_SIXM != "Missing",
    Owns_Farm_Animals != "Missing",
    OWNS_Poultry != "Missing",
    Used_AB_EVER  != "Unsure",
    Used_AB_SIXM  != "Unsure"
    
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
print(doc, target = "Results/Frequencies_% of respondents_17.07.26.docx")


#####################
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

data_model <- data_model_clean %>%
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
data_model <- data_model_clean %>%
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
###############################################################


################################################################################
#Check for dependencies before doing the multivariate analyses
# Conditional variables not included_ because they are clearly not independent and the have also many missing informations, "Withdrawal_Time","Prescription_After_Lab","Animals_Treated","Caretaking_Animals",
vars <- c(
  "GENDER","REGION","ORIGIN_OF_SAMPLE","Heard_AMR",
  "Used_AB_EVER","Used_AB_SIXM",
  "Owns_Farm_Animals","OWNS_Poultry",
  "Has_Toilet","Wash_Hands_Toilet"
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

############################################
## Supplementary Table 2
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

# Save the document
print(doc, target = "Results/Pairwise Chi-Squared Tests & Cramer's V Strength of Association_17.07.26.docx")
############################################################
#################################
#  Creating a publishable table (Supplementary Table 1)
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

data_model <- data_model %>%
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
    # Format stars based on significance tier
    star = case_when(
      is.na(p.value)   ~ "",
      p.value < 0.001  ~ "***",
      p.value < 0.01   ~ "**",
      p.value < 0.05   ~ "*",
      round(p.value, 3) == 0.050 ~ "*",  # Explicitly give exactly 0.05 a star if desired
      TRUE             ~ ""
    ),
    # If it's a reference group, show just 1.00 (omitting the word "Reference")
    OR_CI = ifelse(OR == 1 & is.na(p.value), 
                   "1.00", 
                   sprintf("%.2f [%.2f - %.2f]", OR, Lower_CI, Upper_CI)),
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
    OR_CI = "Odds Ratio [95% CI]",
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
print(doc_regression, target = "Results/univariate_regression_results1_17.07.26.docx")

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
