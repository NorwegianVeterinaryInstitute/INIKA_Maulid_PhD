library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(leaflet)
library(maps)
library(car)
library(pscl)
library(pROC)
library(ResourceSelection)

#joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
#spec(joined_data)

#colnames(joined_data) <- trimws(colnames(joined_data))
# # Selecting and renaming the specific column
#joined_data <- joined_data %>%
#   select("INIKA_ID.x","Age_yrs","GENDER", "SEASON.x", "REGION.x", 
#          "DISTRICT.x", "ORIGIN_OF_SAMPLE.x",
#          "Isolate.x","VITEK_MS_Results", "ESBL_Selection","ESBL",
#          "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
#          "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
#          "SXT_ED1_2", "CTX_ED5", "CTC_ED30", 
#          "Latitude", "Longitude", "GPS_Precision",
#          "If school-aged children, name of the school*", 
#          "Which class/grade are you?*", 
#          "Who is your caretaker?*", 
#          "If others, mention", 
#          "What is your occupation and/or of your caretaker?*", 
#          "Have you ever heard about AMR?", 
#          "If yes, how did you get this information?", 
#          "Have you or your children used any antibiotics at any time?", 
#          "Have you or your children used any antibiotics in the past six months?", 
#          "If yes, where did you get these drugs from?", 
#          "If it was drug sellers or pharmacy; did you have a prescription from the doctor/prescriber?", 
#          "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?", 
#          "Do you have farm animals at your home place?", 
#          "If Yes, do you participate in taking care of animals i.e. cleaning, Feeding?", 
#          "Are your animals being treated for any diseases?", 
#          "Do you have a poultry farm at your home?", 
#          "If yes, are they for eggs or meat?", 
#          "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?", 
#          "Are your Poultry being treated by veterinary doctor?", 
#          "If Yes, does the veterinary doctor tell you about the withdrawal time?", 
#          "Do you have a toilet at your home place?",
#          "Do you wash your hands with soap after every toilet visit?", 
#          "AR_RESIDUES_SAMPLE.x", "COLONY MORPHOLOGY ON C3GR.x", 
#          "COLONY MORPHOLOGY ON CARBA.x", "COLONY MORPHOLOGY ON XLD.x", 
#          "COLONY MORPHOLOGY ON BGA.x", "TSI Media_Slope.x", "TSI Media_Butt.x", 
#          "TSI Media_Gas.x", "TSIMedia_H2S.x", "SIM Media_H2S.x", "SIM Media Indole.x", 
#          "SIM Media Motility.x", "CITRATE.x", "UREASE.x") %>% 
#   rename(INIKA_ID = "INIKA_ID.x", 
#          SEASON = "SEASON.x", 
#          REGION = "REGION.x",
#          DISTRICT = "DISTRICT.x", 
#          ORIGIN_OF_SAMPLE = "ORIGIN_OF_SAMPLE.x",
#          AR_RESIDUES_SAMPLE = "AR_RESIDUES_SAMPLE.x",
#          `COLONY MORPHOLOGY ON C3GR` ="COLONY MORPHOLOGY ON C3GR.x", 
#          `COLONY MORPHOLOGY ON CARBA` ="COLONY MORPHOLOGY ON CARBA.x", 
#          `COLONY MORPHOLOGY ON XLD` ="COLONY MORPHOLOGY ON XLD.x", 
#          `COLONY MORPHOLOGY ON BGA` ="COLONY MORPHOLOGY ON BGA.x",
#          `TSI Media_Slope` ="TSI Media_Slope.x",
#          `TSI Media_Butt` ="TSI Media_Butt.x",
#          `TSI Media_Gas`="TSI Media_Gas.x",
#          TSIMedia_H2S ="TSIMedia_H2S.x",
#          `SIM Media_H2S`="SIM Media_H2S.x",
#          `SIM Media Indole` ="SIM Media Indole.x",
#          `SIM Media Motility` ="SIM Media Motility.x",
#          CITRATE ="CITRATE.x",
#          UREASE ="UREASE.x",
#          Isolate ="Isolate.x") %>% 
#   mutate(DISTRICT = ifelse(
#     DISTRICT == "Ilemala", 
#     "Ilemela",        
#     DISTRICT          
#   )
#   )
# write.csv(joined_data, "data/CLEANED_DATA/joined_data.csv")
################################################################################  
# Importing the file
joined_data <- read_csv("data/CLEANED_DATA/joined_data.csv") 
spec(joined_data)
joined_data <- joined_data %>% 
  rename(
    SAC_School = "If school-aged children, name of the school*",
    Grade_Class = "Which class/grade are you?*",
    Caretaker = "Who is your caretaker?*",
    Others = "If others, mention",
    Occupation_Caretaker = "What is your occupation and/or of your caretaker?*",
    Heard_AMR = "Have you ever heard about AMR?",
    AMR_Info_Source = "If yes, how did you get this information?",
    Ever_Used_Antibiotics = "Have you or your children used any antibiotics at any time?",
    Used_Antibiotics_6M = "Have you or your children used any antibiotics in the past six months?",
    Antibiotic_Source = "If yes, where did you get these drugs from?",
    Prescription = "If it was drug sellers or pharmacy; did you have a prescription from the doctor/prescriber?",
    Prescription_After_Lab = "If was prescribed by the prescriber, was it after Laboratory results of culture and sensitivity?",
    Owns_Farm_Animals = "Do you have farm animals at your home place?",
    Participate_Animal_Care = "If Yes, do you participate in taking care of animals i.e. cleaning, Feeding?",
    Animals_Treated = "Are your animals being treated for any diseases?",
    Owns_Poultry_Farm = "Do you have a poultry farm at your home?",
    Purpose = "If yes, are they for eggs or meat?",
    Participate_Poultry_Care = "Do you participate in taking care of Poultry i.e. Cleaning/Feeding?",
    Poultry_Treated_Vet = "Are your Poultry being treated by veterinary doctor?",
    Withdrawal_Time = "If Yes, does the veterinary doctor tell you about the withdrawal time?",
    Has_Toilet = "Do you have a toilet at your home place?",
    Wash_Hands_Toilet = "Do you wash your hands with soap after every toilet visit?",
    C3GR="COLONY MORPHOLOGY ON C3GR",
    CARBA="COLONY MORPHOLOGY ON CARBA",
    XLD="COLONY MORPHOLOGY ON XLD",
    BGA="COLONY MORPHOLOGY ON BGA"
    
  ) %>% 
  select(-...1) %>% 
  mutate(
    SAC_School = case_when(
    SAC_School == "Pasua primary school" ~ "Pasua Primary School",
    SAC_School == "Mirongo primary school" ~ "Mirongo Primary School",
    SAC_School == "Pasua primary school" ~ "Pasua Primary School",
    SAC_School == "Umbwe primary school" ~ "Umbwe Primary School",
    SAC_School == "Mirongo Primary school" ~ "Mirongo Primary School",
    TRUE ~ SAC_School

    
    ))%>%
  select("INIKA_ID.x.x","Age_yrs","GENDER", "SEASON.x", "REGION.x", "DISTRICT.x", "ORIGIN_OF_SAMPLE.x","Isolate","VITEK_MS_Results", "C3GR",
         "CARBA","XLD","BGA","ESBL_Selection","SAC_School",
         "Grade_Class",
         "Caretaker",
         "Others",
         "Occupation_Caretaker",
         "Heard_AMR" ,
         "AMR_Info_Source",
         "Ever_Used_Antibiotics" ,
         "Used_Antibiotics_6M" ,
         "Antibiotic_Source",
         "Prescription" ,
         "Prescription_After_Lab" ,
         "Owns_Farm_Animals" ,
         "Participate_Animal_Care",
         "Animals_Treated" ,
         "Owns_Poultry_Farm" ,
         "Purpose" ,
         "Participate_Poultry_Care" ,
         "Poultry_Treated_Vet", 
         "Withdrawal_Time", 
         "Has_Toilet",
         "Wash_Hands_Toilet")%>%
  mutate(Outcome_ESC_Bact = (case_when(VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae") 
                            & (C3GR %in% c("Pinkish/Reddish","Pinkish" ,"Metallic blue") 
                            | CARBA %in% c("Pinkish", "Metallic blue")
                            ) ~ "1",
                            Isolate == "No growth" ~ "0",
                            TRUE ~ "Check")))

  unique(joined_data[["Isolate"]])  # this gives you all observation              
## Filtering
    filtered_data <- joined_data %>%
      filter(ORIGIN_OF_SAMPLE.x == "Schoolchildren" & is.na(SAC_School))

    # Filling name of the school to Origin of sample if is Schoolchidren     
    joined_data <- joined_data %>%
      mutate(
        SAC_School = case_when(
          # Check all three conditions simultaneously
          DISTRICT.x == "Magu" & ORIGIN_OF_SAMPLE.x == "Schoolchildren" & is.na(SAC_School) ~ "Salisima primary school",
          DISTRICT.x == "Nyamagana" & ORIGIN_OF_SAMPLE.x == "Schoolchildren" & is.na(SAC_School) ~ "Mirongo Primary School",
          # Default: For all other rows, keep the original value
          TRUE ~ SAC_School
        )
      )
    
    
    Check_1<-joined_data%>%
      filter(Outcome_ESC_Bact=="1")
    
    Check_0<-joined_data%>%
      filter(Outcome_ESC_Bact=="0")
    
    Check_Check<-joined_data%>%
      filter(Outcome_ESC_Bact=="Check")
    
    # We can exclude all observations named "Check" from the dataset for factor analyses
    
    joined_data_Factor<-joined_data%>%
     filter(Outcome_ESC_Bact!="Check")
    
    view(joined_data_Factor)
    
    # 1. Need to check that there are no missing variables within the Questions(Factors)
    # 2. What is the frequency of the different answers/ of the factors? In realtion to the Outcome_ESC_Bact ( group by)
    # Use the frequency function you already have, 
    
    library(dplyr)
    library(purrr)
    
    # Function to compute frequencies, percentages, and optional grouping by outcome
    get_combined_frequencies <- function(df, outcome_var = NULL) {
      if (!is.null(outcome_var) && outcome_var %in% names(df)) {
        # Grouped by outcome
        result <- map_dfr(
          setdiff(names(df), outcome_var),
          function(colname) {
            freq_table <- as.data.frame(table(df[[colname]], df[[outcome_var]], useNA = "ifany"))
            colnames(freq_table) <- c("Value", "Outcome", "Frequency")
            freq_table <- freq_table %>%
              group_by(Outcome) %>%
              mutate(Percent = round(Frequency / sum(Frequency) * 100, 2),
                     Column = colname)
            freq_table
          }
        )
      } else {
        # No grouping
        result <- map_dfr(
          names(df),
          function(colname) {
            freq_table <- as.data.frame(table(df[[colname]], useNA = "ifany"))
            colnames(freq_table) <- c("Value", "Frequency")
            freq_table <- mutate(freq_table,
                                 Percent = round(Frequency / sum(Frequency) * 100, 2),
                                 Column = colname)
            freq_table
          }
        )
      }
      return(result)
    }
    
    # Example usage:
  
    
    # Get combined table grouped by Outcome
    combined_table <- get_combined_frequencies(joined_data_Factor, outcome_var = "Outcome_ESC_Bact")
    
    Factors<-combined_table%>%
      filter(Column!="INIKA_ID")
    
    # View the table
    View(Factors) 
    
    #Filter the ones that are NA to check if there is an explanation or if you might have the content in your papers/ questionaires?

 ## HERE BELLOW FAILED TO COMPUTE##
# Univariate analysis

   # A. First define the outcomevariable as numeric
  
    joined_data_Factor$Outcome_ESC_Bact<- as.numeric(joined_data_Factor$Outcome_ESC_Bact) 
     #  Step 1: Define the  Outcome and Predictor Variables
    # 1. Define the outcome (dependent) variable
    
    #outcome_var <- "Outcome_ESC_Bact" 
    
    # 2. Define  list of predictor (independent) variables
    # Replace this  with the actual column names from your dataset B
    predictor_vars <- c("Heard_AMR", "Used_Antibiotics_6M", 
                        "Prescription_After_Lab", "Owns_Farm_Animals", 
                        "Participate_Animal_Care",
                        "Animals_Treated", "Owns_Poultry_Farm", 
                        "Participate_Poultry_Care",
                        "Poultry_Treated_Vet", "Withdrawal_Time", "Has_Toilet", 
                        "Wash_Hands_Toilet",
                        "ORIGIN_OF_SAMPLE", "GENDER", "REGION")
    
  #  Step 2: Compute Univariate P-values 
    
    # for doing each variable separately
   # glm(Outcome_ESC_Bact ~ Heard_AMR, data = joined_data_Factor, family = binomial)
    
    # cjecking the levels of the predictors
    
    #table(joined_data_Factor$Outcome_ESC_Bact)
    

    # Step 1: Run univariate models and extract OR, CI, and p-values
    univariate_results <- lapply(predictor_vars, function(var) {
      f <- as.formula(paste(outcome_var, "~", var))
      model <- glm(f, data = joined_data_Factor, family = binomial)
      
      # Extract coefficients
      coefs <- summary(model)$coefficients
      
      # Odds Ratio = exp(beta)
      OR <- exp(coefs[2, "Estimate"])
      
      # 95% Confidence Interval for OR
      CI <- exp(confint(model))[2, ]   # second row = predictor
      
      # p-value
      pval <- coefs[2, "Pr(>|z|)"]
      
      # Return as data frame row
      data.frame(
        Predictor = var,
        Odds_Ratio = OR,
        CI_Lower = CI[1],
        CI_Upper = CI[2],
        P_Value = pval
      )
    })
    
    # Step 2: Combine into one summary table
    univariate_summary <- do.call(rbind, univariate_results)
    
    # Step 3: Sort by p-value if desired
    univariate_summary <- univariate_summary[order(univariate_summary$P_Value), ]
    
    # View results
    print(univariate_summary)
    
    
   # 20.11.25 All above are ok: 
   # The predictors Heard_AMR, Prescription_after_Lab, ORIGIN_OF_SAMPLE, Withdrawal_Time,REGION are all potential predicptrs to be included in the final model
    # Perhaps also include the Owns_Farm_Animals,Has_Toilet
    
    # NOW mulitivariate model to be built
    
    multi_model_1<-glm(Outcome_ESC_Bact ~  Heard_AMR + Prescription_After_Lab + ORIGIN_OF_SAMPLE + Withdrawal_Time + REGION, data=joined_data_Factor, family = binomial )
    summary(multi_model_1)
    
    multi_model_2<-glm(Outcome_ESC_Bact ~  Heard_AMR + Prescription_After_Lab + ORIGIN_OF_SAMPLE + Withdrawal_Time , data=joined_data_Factor, family = binomial )
    summary(multi_model_2)
    
    multi_model_3<-glm(Outcome_ESC_Bact ~  Heard_AMR  + ORIGIN_OF_SAMPLE + Withdrawal_Time , data=joined_data_Factor, family = binomial )
    summary(multi_model_3)
    
    
    multi_model_4<-glm(Outcome_ESC_Bact ~  Heard_AMR + ORIGIN_OF_SAMPLE +  Prescription_After_Lab, data=joined_data_Factor, family = binomial )
    summary(multi_model_4)
    
      
      # Ensure the variables are treated as factors for categorical analysis
 # Convert the outcome variable
      Factors[[outcome_var]] <- as.factor(Factors[[outcome_var]])
      
 # Convert all predictor variables
      for (var in predictor_vars) {
        Factors[[var]] <- as.factor(Factors[[var]])
      }  
      
      # Perform the Chi-squared test (suitable for two categorical variables)
      # A tryCatch block handles errors if a test fails (e.g., due to zero cells)
      test_result <- tryCatch({
        test <- chisq.test(Factors[[var]], Factors[[outcome_var]])
        test$p.value
      }, error = function(e) {
        # If Chi-squared fails, return NA or consider a different test (e.g., Fisher's exact)
        # For this example, we'll return a placeholder p-value of 1
        message(paste("Warning: Chi-square test failed for", var, ". Error:", conditionMessage(e)))
        return(1)
      })
      
      # Store the results
 
    ### START HERE 13.11.2025###

## Creating new dataset with unique data that were confirmed ESBL

# 1. Identify the unique INIKA_IDs that meet the filtering criteria
target_ids <- joined_data %>%
  filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
  distinct(INIKA_ID) %>%
  pull(INIKA_ID)

# 2. Filter the original data to keep only the first observation for each target ID
#    (We need to pick one row per ID to use as a template for the recoding)
base_data_unique_id <- joined_data %>%
  filter(INIKA_ID %in% target_ids) %>%
  distinct(INIKA_ID, .keep_all = TRUE) # This keeps the first occurrence for each unique ID

# --- Scenario A: VITEK_MS_Results = "Escherichia coli" ---
# 3. Mutate the base data for Scenario A
df_scenario_A <- base_data_unique_id %>%
  mutate(
    VITEK_MS_Results = "Escherichia coli", # Set VITEK result for this scenario
    ESBL = 1,                             # Set ESBL
    Scenario = "A"                        # Identify the scenario
  )

# --- Scenario B: VITEK_MS_Results = "Klebsiella pneumoniae" ---
# 4. Mutate the base data for Scenario B
df_scenario_B <- base_data_unique_id %>%
  mutate(
    VITEK_MS_Results = "Klebsiella pneumoniae", # Set VITEK result for this scenario
    ESBL = 1,                                  # Set ESBL
    Scenario = "B"                             # Identify the scenario
  )

# --- Final Step: Combine the two scenarios vertically ---
# 5. Bind the two scenario data frames together
combined_new_df <- bind_rows(df_scenario_A, df_scenario_B)

# Verify the final count: it should be twice the number of unique target IDs
 print(nrow(combined_new_df))

######
 unique_id_count <- combined_new_df %>%
   # 1. Filter the data based on VITEK_MS_Results and ESBL_Selection
   filter(
     VITEK_MS_Results == "Escherichia coli",
     ESBL_Selection >= 5
   ) %>%
   # 2. Count the number of unique INIKA_IDs in the filtered subset
   summarise(
     Unique_INIKA_IDs = n_distinct(INIKA_ID)
   ) %>%
   pull() # Extract the count value
 
 # Print the result
 print(unique_id_count)
 








###############################################################################
# 
#   # 1. Define the base filtered and distinct data once
# library(dplyr)
# 
# # 1. Identify the unique INIKA_IDs that meet the filtering criteria
# target_ids <- df %>%
#   filter(Isolate %in% c("E.coli", "K.pneumoniae")) %>%
#   distinct(INIKA_ID) %>%
#   pull(INIKA_ID)
# 
# # 2. Filter the original data to keep only the first observation for each target ID
# #    (We need to pick one row per ID to use as a template for the recoding)
# base_data_unique_id <- df %>%
#   filter(INIKA_ID %in% target_ids) %>%
#   distinct(INIKA_ID, .keep_all = TRUE) # This keeps the first occurrence for each unique ID
# 
# # --- Scenario A: VITEK_MS_Results = "Escherichia coli" ---
# # 3. Mutate the base data for Scenario A
# df_scenario_A <- base_data_unique_id %>%
#   mutate(
#     VITEK_MS_Results = "Escherichia coli", # Set VITEK result for this scenario
#     ESBL = 1,                             # Set ESBL
#     Scenario = "A"                        # Identify the scenario
#   )

# --- Scenario B: VITEK_MS_Results = "Klebsiella pneumoniae" ---
# 4. Mutate the base data for Scenario B
df_scenario_B <- base_data_unique_id %>%
  mutate(
    VITEK_MS_Results = "Klebsiella pneumoniae", # Set VITEK result for this scenario
    ESBL = 1,                                  # Set ESBL
    Scenario = "B"                             # Identify the scenario
  )

# --- Final Step: Combine the two scenarios vertically ---
# 5. Bind the two scenario data frames together
combined_new_df <- bind_rows(df_scenario_A, df_scenario_B)

# Verify the final count: it should be twice the number of unique target IDs
# print(nrow(combined_new_df))
###
 unique_id_count <- UNIQUE_ESC %>%
   # 1. Filter the data based on VITEK_MS_Results and ESBL_Selection
   filter(
     VITEK_MS_Results == "Escherichia coli",
     ESBL_Selection >= 5
   ) %>%
   # 2. Count the number of unique INIKA_IDs in the filtered subset
   summarise(
     Unique_INIKA_IDs = n_distinct(INIKA_ID)
   ) %>%
   pull() # Extract the count value

 # Print the result
 print(unique_id_count)



################################################################################
