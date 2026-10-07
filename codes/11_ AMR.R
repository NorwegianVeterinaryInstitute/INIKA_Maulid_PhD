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


# Step 1–5: Prepare the frequency table
df_freq_wide <- joined_data %>%
  select(INIKA_ID,REGION,ORIGIN_OF_SAMPLE,Isolate,VITEK_MS_Results, ESBL_Selection, ESBL,PROTOCOL,
         AMX_ED10,AZM_ED15,CRO_ED30,CIP_ED5,DOX_ED30,FLR_ED30,
         GEN_ED10,MEM_ED10,OXY_ED30,POL_ED300,SXT_ED1_2,
         CTX_ED5,CTC_ED30)%>%

  mutate(across(-c(VITEK_MS_Results, Isolate, ESBL), as.character)) %>%
  pivot_longer(
    cols = -c(VITEK_MS_Results, Isolate, ESBL),
    names_to = "substance_name",
    values_to = "value"
  ) %>%
  group_by(VITEK_MS_Results, substance_name, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(VITEK_MS_Results, substance_name) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(VITEK_MS_Results, substance_name),
    names_from = value,
    values_from = percentage,
    values_fill = 0
  )
  

# Step 6: Reorder value columns numerically
value_cols <- setdiff(names(df_freq_wide), c("VITEK_MS_Results", "substance_name"))

# Step 7: Convert to numeric, sort, then back to character
sorted_value_cols <- value_cols[order(suppressWarnings(as.numeric(value_cols)))]

# Step 8: Reorder the columns in df_freq_wide
df_freq_wide <- df_freq_wide %>%
  select(VITEK_MS_Results, substance_name, all_of(sorted_value_cols))

# Step 9: Round all numeric columns to one decimal
df_freq_wide <- df_freq_wide %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))



## Importing the CLSI break point file
CLSI_2024_BREAK_POINT <- read_excel("data/CLSI_2024_BREAK POINT.xlsx")
# Well we can not say that the EPI_CUTOFF is equal to the CLSI breakpoint, so we either need to import a file with the EPICUTOFFS or write code for each substance with the cutoffs where we have an epidemiological cutoff, nothing wrong with the codeing just the naming...
EPI_CUTOFF<-CLSI_2024_BREAK_POINT
  

  Data <- joined_data %>%
  select(INIKA_ID.x,REGION.x,SEASON.x,ORIGIN_OF_SAMPLE.x,Isolate.x,
         VITEK_MS_Results,ESBL, AMX_ED10,AZM_ED15,CRO_ED30,
         CIP_ED5,DOX_ED30,FLR_ED30,GEN_ED10,MEM_ED10,
         OXY_ED30,POL_ED300,SXT_ED1_2,CTX_ED5,CTC_ED30)%>%
  rename(Amoxicillin="AMX_ED10", Azithromyzin="AZM_ED15", 
         Ceftriaxone="CRO_ED30", Ciprofloxacin="CIP_ED5",
         Doxycycline="DOX_ED30" , Florfenicol="FLR_ED30",
         Gentamicin="GEN_ED10", Meropenem="MEM_ED10", 
         Oxytetracycline="OXY_ED30",
         "Polymyxin B (PB) "="POL_ED300",    # Is there a space betweeen in the file?, Should be better to clean that up beforehand , you can use the package stringr, install if not already installed.
         "Sulphonamides/Trimethoprim"="SXT_ED1_2", 
         Cefotaxime ="CTX_ED5",
         "Cefotaxime+clavulanic acid"="CTC_ED30")

  AMR_Data<-pivot_longer(Data, cols=-c(INIKA_ID.x,
                                       VITEK_MS_Results, 
                                       Isolate.x,REGION.x,SEASON.x, 
                                       ORIGIN_OF_SAMPLE.x,ESBL,),
                         names_to = "Antimicrobial_substance",
                         values_to="value")
  
# Join the pivot_longer file with the EPI_CUTOFF by substance (Note this means that the in both files needs to have been cleaned so they are spelled exactly the same in both files).
# Note! do this separately for ESC-resistant E.coli, ESC-resistant Klebsiella and Salmonella as the Cutoff values will differ and it is more easy for you to do one thing at the time.

AMR_Cutoff<-left_join(AMR_Data, EPI_CUTOFF, by = "Antimicrobial_substance")


# Categorize each isolate from the AMR_Data to either 1 or 0 ( or R and S ( Resistance and Susceptible)
# Note it is not the concentration that you use for the Breakpoint/ cutoff when classifying
# You need to decide if you want to use the S or the R and here you will have a problem as well as the file you have read in have charachter values like >= or <=, meaning you need to mae a new column, based on either S or R, 
# and then determine "Cutoff" using the case_when function, for each of the values.
# You can try yourself , if problematic let me know, but if we are going to use the epidemiological cutoff values anyway you need to read in another table.
  " MN I hav estopped here just had a snapshot below, it seems ok, but you/ we need to have this classification of Phenotypic ressitance in place first."                                                           
 AMR_Data_Classified <- AMR_Cutoff %>%
     mutate(Phenotype = case_when(
     value < Concentration ~ "1",
            TRUE ~ "0"
                    ))   

# Then you can count the total number of Phenotypes=1 
# in relation to the total number of isolates per bacterial species 
# and antibiotic substance

      library(dplyr)
                                                             
  resistance_summary <-  AMR_Data_Classified%>%
 group_by(VITEK_MS_Results, Antimicrobial_substance) %>%
 summarise(
 Total = n(),
 Resistant = sum(Phenotype == 1, na.rm = TRUE),
 Resistance_Percent = round((Resistant / Total) * 100, 1)
                    )
  ###################################################################
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  
  names(AMR_Cutoff)
  
  
  # Define the list of grouping variables
  group_vars_list <- c("REGION.x","SEASON.x", "ORIGIN_OF_SAMPLE.x")
  
  # --- Step 1: Calculate the Categorical Result (S or R) ---
  
  Ecoli_Susceptibility <- AMR_Cutoff %>%
    
    # Filter only for Escherichia coli
    filter(VITEK_MS_Results == "Escherichia coli") %>%
    
    # Ensure the measurement and breakpoints are numeric for comparison
    # The measured zone diameter is in the 'value' column
    # Use str_replace_all to remove common non-numeric symbols before conversion
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""), # Removes all characters except numbers and decimal point
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
     
    
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""), # Removes all characters except numbers and decimal point
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),

      Measured_Zone = as.numeric(str_trim(value)),
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    # Ensure only one result is counted per isolate for each antimicrobial test
    distinct(INIKA_ID.x, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- Step 2: Define the Summary Function ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    # Total count (denominator) for the current grouping variable and antimicrobial
    Total_Count <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance) %>%
      summarise(Total_Tests = n(), .groups = 'drop')

    # Count the S, I, R results
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%

      # Join the total counts back
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%

      # Calculate Percentage (rounded to 1 decimal place)
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var
      ) %>%

      # Select and arrange columns for final output
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Categorical_Result,
        Frequency = N,
        Total_Tests,
        Percentage
      )

    return(as.data.frame(Summary_Table))
  }

  # --- Step 3: Apply the function to both grouping variables and combine results ---
  ECO_Full_Susceptibility_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Ecoli_Susceptibility, g_var)
    })
  )
   
    
    
 # ECO_Final_Susceptibility_Summary <- purrr::map_dfr(group_vars_list, ~calculate_susceptibility_summary(Ecoli_Susceptibility, .x))
  
  # Save file
  write_tsv(ECO_Full_Susceptibility_Summary, "Results/ECO_Full_Susceptibility_Summary.tsv")
  
######################################################################
  # AMR Frequencies for Klebsiella pneumoniae
  
  # Define the list of grouping variables
  group_vars_list <- c("REGION.x","SEASON.x", "ORIGIN_OF_SAMPLE.x")
  
  # --- Step 1: Calculate the Categorical Result (S, or R) ---
  
  Kpn_Susceptibility <- AMR_Cutoff %>%
    
    # Filter only for Escherichia coli
    filter(VITEK_MS_Results == "Klebsiella pneumoniae") %>%
    
    # Ensure the measurement and breakpoints are numeric for comparison
    # The measured zone diameter is in the 'value' column
    # Use str_replace_all to remove common non-numeric symbols before conversion
    mutate(
      S_Cleaned = str_replace_all(S, "[^0-9.]", ""), # Removes all characters except numbers and decimal point
      R_Cleaned = str_replace_all(R, "[^0-9.]", ""),
      
      Measured_Zone = as.numeric(str_trim(value)), 
      S_Cutoff = as.numeric(S_Cleaned),
      R_Cutoff = as.numeric(R_Cleaned)
    ) %>%
    
    # Implement the breakpoint logic for zone diameters (in mm)
 
    mutate(
      Categorical_Result = case_when(
        Measured_Zone >= S_Cutoff ~ "S",
        Measured_Zone <= R_Cutoff ~ "R",
        TRUE ~ NA_character_
      )
    ) %>%
    
    # Ensure only one result is counted per isolate for each antimicrobial test
    distinct(INIKA_ID.x, Antimicrobial_substance, .keep_all = TRUE) %>%
    
    # Remove tests where the result could not be determined
    filter(!is.na(Categorical_Result))
  
  
  # --- Step 2: Define the Summary Function ---
  
  calculate_susceptibility_summary <- function(data, group_var) {
    
    # Total count (denominator) for the current grouping variable and antimicrobial
    Total_Count <- data %>% 
      group_by(!!sym(group_var), Antimicrobial_substance) %>% 
      summarise(Total_Tests = n(), .groups = 'drop')
    
    # Count the S or R results
    Summary_Table <- data %>%
      group_by(!!sym(group_var), Antimicrobial_substance, Categorical_Result) %>%
      summarise(N = n(), .groups = 'drop_last') %>%
      
      # Join the total counts back
      left_join(Total_Count, by = c(group_var, "Antimicrobial_substance")) %>%
      
      # Calculate Percentage (rounded to 1 decimal place)
      mutate(
        Percentage = round((N / Total_Tests) * 100, 1),
        Grouping_Variable = group_var
      ) %>%
      
      # Select and arrange columns for final output
      select(
        Grouping_Variable,
        Grouping_Value = !!sym(group_var),
        Antimicrobial_substance,
        Categorical_Result,
        Frequency = N,
        Total_Tests,
        Percentage
      )
    
    return(as.data.frame(Summary_Table))
  }
  
  # --- Step 3: Apply the function to both grouping variables and combine results ---
  Kpn_Full_Susceptibility_Summary <- bind_rows(
    lapply(group_vars_list, function(g_var) {
      calculate_susceptibility_summary(Kpn_Susceptibility, g_var)
    })
  )
 #HERE 5/11/2025 
 # KpnFinal_Susceptibility_Summary <- purrr::map_dfr(group_vars_list, ~calculate_susceptibility_summary(Kpn_Susceptibility, .x))
  
  # Save file
  write_tsv(Kpn_Full_Susceptibility_Summary, "Results/Kpn_Full_Susceptibility_Summary.tsv")
########################################################################  
 
   # Load necessary libraries
  library(dplyr)
  library(tidyr)
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map
  # The original column names (the zone diameter values)
  zone_diameter_cols <- c(
    "AMX_ED10", "AZM_ED15", "CRO_ED30", "CIP_ED5", "DOX_ED30", 
    "FLR_ED30", "GEN_ED10", "MEM_ED10", "OXY_ED30", "POL_ED300", 
    "SXT_ED1_2", "CTX_ED5", "CTC_ED30"
  )
  
  # The named vector for renaming (Original_Column_Name = New_Descriptive_Name)
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
  
  # 2. Calculate the Distribution Statistics with renaming post-pivot
  Zone_Distribution_Table <- joined_data %>%
    # Select the specific Zone Diameter columns using the original names
    select(all_of(zone_diameter_cols)) %>%
    
    # Convert data from wide to long format
    pivot_longer(
      cols = everything(),
      names_to = "Antimicrobial_Substance", # This column now holds the codes (e.g., "AMX_ED10")
      values_to = "Zone_Diameter_mm"
    ) %>%
    
    # *** FIX APPLIED HERE: Rename the values in the 'Antimicrobial_Substance' column ***
    mutate(
      Antimicrobial_Substance = recode(Antimicrobial_Substance, !!!rename_map)
    ) %>%
    
    # Remove rows where the zone was not tested (NA)
    filter(!is.na(Zone_Diameter_mm)) %>%
    
    # Group by the descriptive substance name
    group_by(Antimicrobial_Substance) %>%
    
    # Calculate distribution statistics
    summarise(
      N_Tested = n(),
      Mean_Zone_mm = mean(Zone_Diameter_mm),
      SD_Zone_mm = sd(Zone_Diameter_mm),
      Median_Zone_mm = median(Zone_Diameter_mm),
      Min_Zone_mm = min(Zone_Diameter_mm),
      Max_Zone_mm = max(Zone_Diameter_mm),
      .groups = 'drop'
    ) %>%
    
    # Arrange by the number tested or mean zone (optional)
    arrange(desc(N_Tested), Mean_Zone_mm)
  
  # Print the resulting table
  print("--- Distribution of Zone Diameters (mm) per Antimicrobial Substance ---")
  print(Zone_Distribution_Table)
  
  # Save the file
  write_tsv(Zone_Distribution_Table, "Results/Zone_Distribution_Table.tsv")
  #####################################################################
  # Load necessary libraries
  # Load necessary libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2) # For plotting the distributions
  library(cowplot) # For plotting the distributions
  
  # 1. Define the Antimicrobial Column Names and the Renaming Map (using your previous definitions)
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
  
  # 2. Prepare the Filtered and Long-Format Data (E. coli population)
  Ecoli_Zone_Data_Long <- joined_data %>%
    # Filter for the target population
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
    select(-Antimicrobial_Substance_Code) # Remove the code column
  
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
  
###########################################################################
  
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
    select(INIKA_ID.x,REGION.x,SEASON.x,ORIGIN_OF_SAMPLE.x,Isolate.x,
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
  
  AMR_Data<-pivot_longer(Data, cols=-c(INIKA_ID.x,
                                       VITEK_MS_Results,
                                       Isolate.x,REGION.x,SEASON.x,
                                       ORIGIN_OF_SAMPLE.x, ESBL), # <--- ESBL ADDED HERE
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