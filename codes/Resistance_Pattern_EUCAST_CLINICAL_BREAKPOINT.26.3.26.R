# k y---
 # Pivot and Classify 
    Ecoli_Long <- Data_Clean %>%
      pivot_longer(
        cols = any_of(abx_classes$Antimicrobial_substance),
        names_to = "Antimicrobial_substance",
        values_to = "value"
      ) %>%
      left_join(EPI_CUTOFF, by = "Antimicrobial_substance") %>%
      mutate(
        S_num = as.numeric(str_replace_all(S, "[^0-9.]", "")),
        R_num = as.numeric(str_replace_all(R, "[^0-9.]", "")),
        Measured_Zone = as.numeric(str_trim(value)),
        Categorical_Result = case_when(
          Measured_Zone >= S_num ~ "S",
          Measured_Zone <= R_num ~ "R",
          Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Categorical_Result))
    
    # Summary Logic 
    calculate_summary <- function(df, group_var) {
      group_cols <- if(group_var == "Overall") "Antimicrobial_substance" else c(group_var, "Antimicrobial_substance")
      df %>%
        group_by(across(all_of(group_cols))) %>%
        summarise(
          Total = n(),
          Resistant = sum(Categorical_Result == "R"),
          .groups = "drop"
        ) %>%
        mutate(
          Grouping_Variable = group_var,
          Grouping_Value = if(group_var == "Overall") "Overall" else as.character(.data[[group_var]])
        )
    }
    
    group_vars <- c("Overall", "REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE")
    Combined_Summaries <- map_dfr(group_vars, ~calculate_summary(Ecoli_Long, .x))
    
    # Final Formatting 
    Final_Table <- Combined_Summaries %>%
      rowwise() %>%
      mutate(
        Percentage = round((Resistant / Total) * 100, 1),
        CI_text = if(Grouping_Variable == "Overall") {
          test <- binom.test(Resistant, Total)
          paste0(" (", round(test$conf.int[1]*100, 1), "–", round(test$conf.int[2]*100, 1), ")")
        } else { "" },
        Display_Value = paste0(Percentage, "% (", Resistant, "/", Total, ")", CI_text)
      ) %>%
      ungroup() %>%
      select(Antimicrobial_substance, Grouping_Value, Display_Value) %>%
      pivot_wider(names_from = Grouping_Value, values_from = Display_Value) %>%
      left_join(abx_classes, by = "Antimicrobial_substance") %>%
      select(Class, Antimicrobial_substance)
    title: "Descriptive analysis of AMR data"
# Below is the function you run:

################################################################################
## Supplementary Table 6
# AMR Frequencies for Klebsiella pneumoniae
 
  ## Importing the K.pneumoniae_ECOFF_EUCAST break point file
  
  #Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_EUCAST.xlsx", sheet = 2)
  
   Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel("data/ECOFF_E.coli_K.pneumoniae.xlsx", 
                                           sheet = "K.pneumoniae_ECOFF_EUCAST")
  Kpn_EPI_CUTOFF<-Kpn_ECOFF_EUCAST_BREAK_POINT
  
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
################################################################################
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
################################################################################ 
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
#############################################################################
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