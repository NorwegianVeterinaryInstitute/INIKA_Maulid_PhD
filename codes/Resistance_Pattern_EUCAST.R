

##############################################################################
# SETUP
##############################################################################
library(Hmisc)
library(tidyverse)
library(readxl)
library(writexl)
library(flextable)
library(officer)
library(ggplot2)


##############################################################################
# MASTER DATA
##############################################################################

# Importing the file 

joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")
spec(joined_data)

names(joined_data)[grep("INIKA", names(joined_data))]

# Renaming of the columns


joined_data <- joined_data %>%
  rename(
    INIKA_ID = INIKA_ID.x.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    ESBL = ESBL_Presumptivefinal
  )

##############################################################################
# ECOFF TABLES
##############################################################################
# Importing E.coli ECOFF

Ecoli_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx"
)

Ecoli_EPI_CUTOFF <- Ecoli_ECOFF_EUCAST_BREAK_POINT



#Import ECOFF table with the sheet for Klebsiella
Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx",
  sheet = "K.pneumoniae_ECOFF_EUCAST"
)
Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT


Ecoli_EPI_CUTOFF <- Ecoli_EPI_CUTOFF %>%
  mutate(
    Antimicrobial_substance = recode(
      Antimicrobial_substance,
      "Polymyxin B" = "Polymyxin_B(PB)",
      "Sulfamethoxazole_Trimethoprim" = "Sulfamethoxazole/Trimethoprim",
      "Cefotaxime _Clavulanic Acid" = "Cefotaxime/Clavulanic Acid"
    )
  )

Kpn_EPI_CUTOFF <- Kpn_EPI_CUTOFF %>%
  mutate(
    Antimicrobial_substance = recode(
      Antimicrobial_substance,
      "Polymyxin B" = "Polymyxin_B(PB)",
      "Sulfamethoxazole_Trimethoprim" = "Sulfamethoxazole/Trimethoprim",
      "Cefotaxime_Clavulanic Acid" = "Cefotaxime/Clavulanic Acid"
    )
  )
#################################
# Supplementary Table 5
#################################
# Define the list of grouping variables (Only defined ONCE here)
group_vars_list <- c("ORIGIN_OF_SAMPLE", "REGION.x", "SEASON.x")

# 1. Define Antibiotic Class Mapping 
# This ensures drugs are grouped by their mechanism of action in your table

abx_classes <- tribble(
  ~Antimicrobial_substance,          ~Class,
  "Amoxicillin",                     "Penicillins",
  "Azithromycin",                    "Macrolides",
  "Ceftriaxone",                     "Cephalosporins (3rd Gen)",
  "Cefotaxime",                      "Cephalosporins (3rd Gen)",
  "Cefotaxime/Clavulanic Acid",      "β-lactam/Inhibitor",
  "Ciprofloxacin",                   "Fluoroquinolones",
  "Doxycycline",                     "Tetracyclines",
  "Oxytetracycline",                 "Tetracyclines",
  "Gentamicin",                      "Aminoglycosides",
  "Meropenem",                       "Carbapenems",
  "Florfenicol",                     "Phenicols",
  "Polymyxin_B(PB)",                  "Polymyxins",
  "Sulfamethoxazole/Trimethoprim",   "Sulfonamides"
)

abx_cols <- abx_classes$Antimicrobial_substance


#Used later for: MDR calculation and Table grouping
#Step 3: Define ESCR E. coli
# Define the Main Processing Function 
process_ecoli_amr_results <- function(joined_data, Ecoli_EPI_CUTOFF) {
  
  # Cleaning and Filtering 
  

  
  Data_Clean <- joined_data %>%
    select(
      INIKA_ID, REGION.x, SEASON.x, ORIGIN_OF_SAMPLE,
      `COLONY MORPHOLOGY ON C3GR`, Isolate,
      VITEK_MS_Results,
      AMX_ED10, AZM_ED15, CRO_ED30,
      CIP_ED5, DOX_ED30, FLR_ED30,
      GEN_ED10, MEM_ED10, OXY_ED30,
      POL_ED300, SXT_ED1_2,
      CTX_ED5, CTC_ED30
    ) %>%
    rename(
        Amoxicillin = AMX_ED10,
        Azithromycin = AZM_ED15,
        Ceftriaxone = CRO_ED30,
        Ciprofloxacin = CIP_ED5,
        Doxycycline = DOX_ED30,
        Florfenicol = FLR_ED30,
        Gentamicin = GEN_ED10,
        Meropenem = MEM_ED10,
        Oxytetracycline = OXY_ED30,
        `Polymyxin_B(PB)` = POL_ED300,
        `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
        Cefotaxime = CTX_ED5,
        `Cefotaxime/Clavulanic Acid` = CTC_ED30) %>%
    mutate(across(all_of(abx_cols), ~ as.numeric(trimws(.)))) %>%
    filter(
      (
        (
          grepl("E", Isolate, ignore.case = TRUE) &
            grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
            grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
        ) |
          (
            grepl("K", Isolate, ignore.case = TRUE) &
              grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
              grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
          )
      ) &
        Ceftriaxone < 23
    )%>%
    distinct(INIKA_ID, .keep_all = TRUE)
  
  message(
    paste("ESCR E. coli isolates retained:",
          nrow(Data_Clean))
  )


  # Pivot and Classify 
  Ecoli_Long <- Data_Clean %>%
    pivot_longer(
      cols = any_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    ) %>%
    left_join(Ecoli_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
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
  
  group_vars <- c("Overall", "ORIGIN_OF_SAMPLE", "REGION.x", "SEASON.x")
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
    select(Class, Antimicrobial_substance, Overall, everything()) %>%
    arrange(Class, Antimicrobial_substance)
  
  return(Final_Table)
}
# Run the Process and Save Output 
final_ESCR_ECO_AMR_table <- process_ecoli_amr_results(joined_data, Ecoli_EPI_CUTOFF)

# View results
print(final_ESCR_ECO_AMR_table)

# Export to Excel 
write_xlsx(final_ESCR_ECO_AMR_table, "Results/Ecoli_ESCR_Resistance_Table.xlsx")

# Export to Word 

# Create the flextable object
ft <- flextable(final_ESCR_ECO_AMR_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(size = 9, part = "all") %>%  # Fits better on one page
  merge_v(j = ~ Class) %>%              # Cleanly groups Penicillins, etc.
  bold(part = "header")

# Define the Landscape properties
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5),
  type = "continuous"
)

# Format the flextable with space-saving logic
ft_compact <- ft %>%
  fontsize(size = 8, part = "all") %>%           # Smaller font for dense data
  padding(padding = 1, part = "all") %>%         # Reduce cell padding
  set_table_properties(layout = "fixed") %>%     # Force columns to stay within bounds
  width(width = 1.1) %>%                         # Manually set column width
  width(j = 1:2, width = 1.5)                    # Give Class and Substance more room

# Create the Results folder if it doesn't exist
if (!dir.exists("Results")) dir.create("Results")

# Save
save_as_docx(
  ft_compact, 
  path = "Results/Ecoli_AMR_Results_Final_Fit.docx",
  pr_section = sect_properties
)
################################################################################
## Supplementary Table 6
# AMR Frequencies for Klebsiella pneumoniae
 
process_kpn_amr_analysis <- function(joined_data, Kpn_EPI_CUTOFF) {
  
  # Data Cleaning & Filtering 
  # Clean the column names of the input data first to remove hidden spaces
  colnames(joined_data) <- trimws(colnames(joined_data))

  #if(length(missing_cols) > 0) {
  #  stop(paste("The following columns are missing from joined_data:", 
   #            paste(missing_cols, collapse = ", ")))
 # }
  
  
  cleaned_data <- joined_data %>%
    rename(
      Amoxicillin = AMX_ED10,
      Azithromycin = AZM_ED15,
      Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5,
      Doxycycline = DOX_ED30,
      Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10,
      Meropenem = MEM_ED10,
      Oxytetracycline = OXY_ED30,
      `Polymyxin_B(PB)` = POL_ED300,
      `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
      Cefotaxime = CTX_ED5,
      `Cefotaxime/Clavulanic Acid` = CTC_ED30) %>%
        mutate(across(all_of(abx_cols), ~ as.numeric(trimws(.)))) %>%
       
    select(
      INIKA_ID,
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE,
      `COLONY MORPHOLOGY ON C3GR`,
      Isolate,
      VITEK_MS_Results,
      ESCR_KPN_presumptive,
      any_of(abx_classes$Antimicrobial_substance)
    ) %>%
    mutate(
      across(
        all_of(abx_classes$Antimicrobial_substance),
        ~ as.numeric(trimws(.))
      )
    ) %>%
    filter(
      trimws(VITEK_MS_Results) == "Klebsiella pneumoniae",
    
      Ceftriaxone < 23
    )%>%
    arrange(INIKA_ID, Isolate) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
  
message(paste("Isolates processed:", nrow(cleaned_data)))
  
  # Pivot & Susceptibility Classification 
  kpn_long <- cleaned_data %>%
    pivot_longer(
      cols = all_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) %>%
    left_join(Kpn_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
    mutate(
      across(c(S, R), ~as.numeric(str_remove_all(.x, "[^0-9.]")), .names = "{.col}_num"),
      Measured_Zone = as.numeric(Measured_Zone),
      Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
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
  group_vars <- c("Overall", "ORIGIN_OF_SAMPLE", "REGION.x", "SEASON.x")
  
  final_tablekpn <- map_dfr(group_vars, ~get_summary(kpn_long, .x)) |>
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
  ft <- flextable(final_tablekpn) |>
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
  write_xlsx(final_tablekpn, "Results/Kpn_AMR_Table.xlsx")
  
  message("Success: Table saved to Results folder as Word (Landscape) and Excel.")
  return(final_tablekpn)
}

#  Execute 
kpn_amr_final_results <- process_kpn_amr_analysis(joined_data, Kpn_EPI_CUTOFF)
###############################################################################

##############################################################################
# Supplementary Table 7
# MDR in Confirmed ESCR E. coli
##############################################################################

  ##############################################################################
  # Supplementary Table 7
  # MDR in Confirmed ESCR E. coli
  ##############################################################################
  
  process_ecoli_MDR_results <- function(
    joined_data,
    Ecoli_EPI_CUTOFF,
    abx_classes
  ) {
    
    ###########################################################################
    # Create confirmed ESCR E. coli dataset
    ###########################################################################
    
    ESCR_confirmed_data <- joined_data %>%
      rename(
        Amoxicillin = AMX_ED10,
        Azithromycin = AZM_ED15,
        Ceftriaxone = CRO_ED30,
        Ciprofloxacin = CIP_ED5,
        Doxycycline = DOX_ED30,
        Florfenicol = FLR_ED30,
        Gentamicin = GEN_ED10,
        Meropenem = MEM_ED10,
        Oxytetracycline = OXY_ED30,
        `Polymyxin_B(PB)` = POL_ED300,
        `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
        Cefotaxime = CTX_ED5,
        `Cefotaxime/Clavulanic Acid` = CTC_ED30
      ) %>%
      mutate(
        across(
          all_of(abx_classes$Antimicrobial_substance),
          ~ as.numeric(trimws(.))
        )
      ) %>%
      filter(
        (
          (
            grepl("E", Isolate, ignore.case = TRUE) &
              grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
              grepl(
                "Pink|Red",
                `COLONY MORPHOLOGY ON C3GR`,
                ignore.case = TRUE
              )
          ) |
            (
              grepl("K", Isolate, ignore.case = TRUE) &
                grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
                grepl(
                  "Metallic blue",
                  `COLONY MORPHOLOGY ON C3GR`,
                  ignore.case = TRUE
                )
            )
        ) &
          Ceftriaxone < 23
      ) %>%
      distinct(INIKA_ID, .keep_all = TRUE)
    
    message(
      paste(
        "ESCR E. coli isolates for MDR analysis:",
        nrow(ESCR_confirmed_data)
      )
    )
    
    ###########################################################################
    # SIR classification
    ###########################################################################
    
    Ecoli_SIR <- ESCR_confirmed_data %>%
      pivot_longer(
        cols = all_of(abx_classes$Antimicrobial_substance),
        names_to = "Antimicrobial_substance",
        values_to = "Measured_Zone"
      ) %>%
      left_join(
        Ecoli_EPI_CUTOFF,
        by = "Antimicrobial_substance"
      ) %>%
      left_join(
        abx_classes,
        by = "Antimicrobial_substance"
      ) %>%
      mutate(
        Measured_Zone = as.numeric(Measured_Zone),
        S_num = as.numeric(str_replace_all(S, "[^0-9.]", "")),
        R_num = as.numeric(str_replace_all(R, "[^0-9.]", "")),
        Categorical_Result = case_when(
          Measured_Zone >= S_num ~ "S",
          Measured_Zone <= R_num ~ "R",
          Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Categorical_Result))
    
    ###########################################################################
    # MDR calculation
    ###########################################################################
    
    MDR_Data <- Ecoli_SIR %>%
      group_by(
        INIKA_ID,
        REGION.x,
        SEASON.x,
        ORIGIN_OF_SAMPLE,
        Class
      ) %>%
      summarise(
        Class_Resistant = any(
          Categorical_Result %in% c("R", "I")
        ),
        .groups = "drop"
      ) %>%
      group_by(
        INIKA_ID,
        REGION.x,
        SEASON.x,
        ORIGIN_OF_SAMPLE
      ) %>%
      summarise(
        Resistant_Classes_Count =
          sum(Class_Resistant, na.rm = TRUE),
        MDR_Status =
          ifelse(
            Resistant_Classes_Count >= 3,
            "MDR",
            "Non-MDR"
          ),
        .groups = "drop"
      )
    
    return(MDR_Data)
  }
  
  ##############################################################################
  # RUN ANALYSIS
  ##############################################################################
  
  ecoli_mdr_results <- process_ecoli_MDR_results(
    joined_data,
    Ecoli_EPI_CUTOFF,
    abx_classes
  )
  
  ##############################################################################
  # HELPER FUNCTION
  ##############################################################################
  
  calc_stats <- function(df) {
    
    n_total <- nrow(df)
    n_mdr <- sum(df$MDR_Status == "MDR")
    
    ci <- binconf(
      n_mdr,
      n_total,
      method = "wilson"
    )
    
    tibble(
      N = n_total,
      MDR_n = n_mdr,
      Percent = ci[1] * 100,
      Lower = ci[2] * 100,
      Upper = ci[3] * 100
    )
  }
  
  ##############################################################################
  # SUMMARY TABLE
  ##############################################################################
  
  overall <- ecoli_mdr_results %>%
    calc_stats() %>%
    mutate(
      Variable = "Overall",
      Category = "Total Population"
    )
  
  grouped <- ecoli_mdr_results %>%
    pivot_longer(
      cols = c(
        REGION.x,
        SEASON.x,
        ORIGIN_OF_SAMPLE
      ),
      names_to = "Variable",
      values_to = "Category"
    ) %>%
    group_by(
      Variable,
      Category
    ) %>%
    group_modify(~ calc_stats(.x)) %>%
    ungroup()
  
  final_ecoli_mdr_table <- bind_rows(
    overall,
    grouped
  ) %>%
    mutate(
      Variable = recode(
        Variable,
        "REGION.x" = "Region",
        "SEASON.x" = "Season",
        "ORIGIN_OF_SAMPLE" = "Origin"
      ),
      `MDR n (%)` =
        sprintf(
          "%d (%.1f%%)",
          MDR_n,
          Percent
        ),
      `95% CI` =
        sprintf(
          "[%.1f - %.1f]",
          Lower,
          Upper
        )
    ) %>%
    select(
      Variable,
      Category,
      N,
      `MDR n (%)`,
      `95% CI`
    )
  
  ##############################################################################
  # EXPORT
  ##############################################################################
  
  doc_table <- final_ecoli_mdr_table %>%
    as_grouped_data(groups = "Variable") %>%
    flextable() %>%
    set_header_labels(
      N = "Total N",
      `MDR n (%)` = "MDR n (%)",
      `95% CI` = "95% CI (Wilson)"
    ) %>%
    bold(part = "header") %>%
    autofit()
  
  if (!dir.exists("Results")) {
    dir.create("Results")
  }
  
  save_as_docx(
    doc_table,
    path = "Results/Ecoli_MDR_Final_Table.docx"
  )



#  RUN ANALYSIS 
ecoli_mdr_results <- process_ecoli_MDR_results(joined_data, Ecoli_EPI_CUTOFF, abx_classes)

# PREPRE PUBLICATION TABLE 
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


##############################################################################
# Supplementary Table 8
# MDR in Confirmed ESCR K. pneumoniae
##############################################################################

process_kpn_MDR_results <- function(
    joined_data,
    Kpn_EPI_CUTOFF,
    abx_classes
) {
  
  ###########################################################################
  # Create confirmed ESCR K. pneumoniae dataset
  ###########################################################################
  
  ESCR_confirmed_data <- joined_data %>%
    rename(
      Amoxicillin = AMX_ED10,
      Azithromycin = AZM_ED15,
      Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5,
      Doxycycline = DOX_ED30,
      Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10,
      Meropenem = MEM_ED10,
      Oxytetracycline = OXY_ED30,
      `Polymyxin_B(PB)` = POL_ED300,
      `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
      Cefotaxime = CTX_ED5,
      `Cefotaxime/Clavulanic Acid` = CTC_ED30
    ) %>%
    mutate(
      across(
        all_of(abx_classes$Antimicrobial_substance),
        ~ as.numeric(trimws(.))
      )
    ) %>%
    filter(
      (
        (
          grepl("K", Isolate, ignore.case = TRUE) &
            grepl("Klebsiella pneumoniae", VITEK_MS_Results, ignore.case = TRUE) &
            grepl(
              "Metallic blue",
              `COLONY MORPHOLOGY ON C3GR`,
              ignore.case = TRUE
            )
        ) |
          (
            grepl("E", Isolate, ignore.case = TRUE) &
              grepl("Klebsiella pneumoniae", VITEK_MS_Results, ignore.case = TRUE) &
              grepl(
                "Pink|Red",
                `COLONY MORPHOLOGY ON C3GR`,
                ignore.case = TRUE
              )
          )
      ) &
        Ceftriaxone < 23
    ) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
  
  message(
    paste(
      "ESCR K. pneumoniae isolates for MDR analysis:",
      nrow(ESCR_confirmed_data)
    )
  )
  
  ###########################################################################
  # SIR classification
  ###########################################################################
  
  Kpn_SIR <- ESCR_confirmed_data %>%
    pivot_longer(
      cols = all_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) %>%
    left_join(
      Kpn_EPI_CUTOFF,
      by = "Antimicrobial_substance"
    ) %>%
    left_join(
      abx_classes,
      by = "Antimicrobial_substance"
    ) %>%
    mutate(
      Measured_Zone = as.numeric(Measured_Zone),
      S_num = as.numeric(str_replace_all(S, "[^0-9.]", "")),
      R_num = as.numeric(str_replace_all(R, "[^0-9.]", "")),
      Categorical_Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Categorical_Result))
  
  ###########################################################################
  # MDR calculation
  ###########################################################################
  
  MDR_Data <- Kpn_SIR %>%
    group_by(
      INIKA_ID,
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE,
      Class
    ) %>%
    summarise(
      Class_Resistant = any(
        Categorical_Result %in% c("R", "I")
      ),
      .groups = "drop"
    ) %>%
    group_by(
      INIKA_ID,
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE
    ) %>%
    summarise(
      Resistant_Classes_Count =
        sum(Class_Resistant, na.rm = TRUE),
      MDR_Status =
        ifelse(
          Resistant_Classes_Count >= 3,
          "MDR",
          "Non-MDR"
        ),
      .groups = "drop"
    )
  
  return(MDR_Data)
}

##############################################################################
# RUN ANALYSIS
##############################################################################

kpn_mdr_results <- process_kpn_MDR_results(
  joined_data,
  Kpn_EPI_CUTOFF,
  abx_classes
)

##############################################################################
# HELPER FUNCTION
##############################################################################

calc_stats <- function(df) {
  
  n_total <- nrow(df)
  n_mdr <- sum(df$MDR_Status == "MDR")
  
  ci <- binconf(
    n_mdr,
    n_total,
    method = "wilson"
  )
  
  tibble(
    N = n_total,
    MDR_n = n_mdr,
    Percent = ci[1] * 100,
    Lower = ci[2] * 100,
    Upper = ci[3] * 100
  )
}

##############################################################################
# SUMMARY TABLE
##############################################################################

overall <- kpn_mdr_results %>%
  calc_stats() %>%
  mutate(
    Variable = "Overall",
    Category = "Total Population"
  )

grouped <- kpn_mdr_results %>%
  pivot_longer(
    cols = c(
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE
    ),
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  group_by(
    Variable,
    Category
  ) %>%
  group_modify(~ calc_stats(.x)) %>%
  ungroup()

final_kpn_mdr_table <- bind_rows(
  overall,
  grouped
) %>%
  mutate(
    Variable = recode(
      Variable,
      "REGION.x" = "Region",
      "SEASON.x" = "Season",
      "ORIGIN_OF_SAMPLE" = "Origin"
    ),
    `MDR n (%)` =
      sprintf(
        "%d (%.1f%%)",
        MDR_n,
        Percent
      ),
    `95% CI` =
      sprintf(
        "[%.1f - %.1f]",
        Lower,
        Upper
      )
  ) %>%
  select(
    Variable,
    Category,
    N,
    `MDR n (%)`,
    `95% CI`
  )

##############################################################################
# EXPORT
##############################################################################

doc_table <- final_kpn_mdr_table %>%
  as_grouped_data(groups = "Variable") %>%
  flextable() %>%
  set_header_labels(
    N = "Total N",
    `MDR n (%)` = "MDR n (%)",
    `95% CI` = "95% CI (Wilson)"
  ) %>%
  bold(part = "header") %>%
  autofit()

if (!dir.exists("Results")) {
  dir.create("Results")
}

save_as_docx(
  doc_table,
  path = "Results/Kpn_MDR_Final_Table.docx"
)


##############################################################################
##############################################################################
# Supplementary Table 3
# Distr##############################################################################
# Computing the AMR for confirmed CR Escherichia coli
##############################################################################

process_cr_ecoli_amr_results <- function(
    joined_data,
    Ecoli_EPI_CUTOFF
) {
  
  ###########################################################################
  # Clean, confirm species and filter CR E. coli
  ###########################################################################
  
  Data_Clean <- joined_data %>%
    select(
      INIKA_ID,
      ORIGIN_OF_SAMPLE,
      SEASON.x,
      REGION.x,
      DISTRICT.x,
      `COLONY MORPHOLOGY ON CARBA`,
      VITEK_MS_Results,
      Isolate,
      Isolate_NVI,
      AMX_ED10,
      AZM_ED15,
      CRO_ED30,
      CIP_ED5,
      DOX_ED30,
      FLR_ED30,
      GEN_ED10,
      MEM_ED10,
      OXY_ED30,
      POL_ED300,
      SXT_ED1_2,
      CTX_ED5,
      CTC_ED30
    ) %>%
    rename_antibiotics() %>%
    convert_antibiotics() %>%
    rename(
      SEASON = SEASON.x,
      REGION = REGION.x,
      DISTRICT = DISTRICT.x
    ) %>%
    mutate(
      is_nvi_ecoli = str_detect(
        coalesce(Isolate_NVI, ""),
        "(?i)E\\.?\\s*coli|Escherichia\\s+coli|^E$"
      ),
      
      CONFIRMED_SPECIES = case_when(
        is_nvi_ecoli ~ "CR Escherichia coli",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      CONFIRMED_SPECIES == "CR Escherichia coli",
      str_detect(
        coalesce(`COLONY MORPHOLOGY ON CARBA`, ""),
        "(?i)Red|Pinkish|Pink|Metalic blue"
      )
    ) %>%
    distinct(
      INIKA_ID,
      .keep_all = TRUE
    )
  
  message(
    sprintf(
      "Total confirmed CR E. coli unique isolates processed: %d",
      nrow(Data_Clean)
    )
  )
  
  ###########################################################################
  # Long-format susceptibility data
  ###########################################################################
  
  Ecoli_Long <- Data_Clean %>%
    pivot_longer(
      cols = any_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    ) %>%
    left_join(
      Ecoli_EPI_CUTOFF,
      by = "Antimicrobial_substance"
    ) %>%
    mutate(
      S_num = as.numeric(
        str_replace_all(S, "[^0-9.]", "")
      ),
      R_num = as.numeric(
        str_replace_all(R, "[^0-9.]", "")
      ),
      Measured_Zone = as.numeric(
        str_trim(value)
      ),
      Categorical_Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num &
          Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(Categorical_Result)
    )
  
  ###########################################################################
  # Summary helper
  ###########################################################################
  
  calculate_summary <- function(
    df,
    group_var
  ) {
    
    group_cols <-
      if (group_var == "Overall") {
        "Antimicrobial_substance"
      } else {
        c(
          group_var,
          "Antimicrobial_substance"
        )
      }
    
    df %>%
      group_by(
        across(
          all_of(group_cols)
        )
      ) %>%
      summarise(
        Total = n(),
        Resistant = sum(
          Categorical_Result == "R"
        ),
        .groups = "drop"
      ) %>%
      mutate(
        Grouping_Variable = group_var,
        Grouping_Value =
          if (group_var == "Overall") {
            "Overall"
          } else {
            as.character(
              .data[[group_var]]
            )
          }
      )
  }
  
  ###########################################################################
  # Generate summaries
  ###########################################################################
  
  group_vars <- c(
    "Overall",
    "ORIGIN_OF_SAMPLE",
    "REGION",
    "SEASON"
  )
  
  Combined_Summaries <- map_dfr(
    group_vars,
    ~ calculate_summary(
      Ecoli_Long,
      .x
    )
  )
  
  ###########################################################################
  # Final table
  ###########################################################################
  
  Final_Table <- Combined_Summaries %>%
    rowwise() %>%
    mutate(
      Percentage = round(
        (Resistant / Total) * 100,
        1
      ),
      
      CI_text =
        if (Grouping_Variable == "Overall") {
          
          test <- binom.test(
            Resistant,
            Total
          )
          
          paste0(
            " (",
            round(
              test$conf.int[1] * 100,
              1
            ),
            "–",
            round(
              test$conf.int[2] * 100,
              1
            ),
            ")"
            
          )
          
        } else {
          
          ""
          
        },
      
      Display_Value =
        paste0(
          Percentage,
          "% (",
          Resistant,
          "/",
          Total,
          ")",
          CI_text
        )
    ) %>%
    ungroup() %>%
    select(
      Antimicrobial_substance,
      Grouping_Value,
      Display_Value
    ) %>%
    pivot_wider(
      names_from = Grouping_Value,
      values_from = Display_Value
    ) %>%
    left_join(
      abx_classes,
      by = "Antimicrobial_substance"
    ) %>%
    select(
      Class,
      Antimicrobial_substance,
      Overall,
      everything()
    ) %>%
    arrange(
      Class,
      Antimicrobial_substance
    )
  
  Final_Table
}

##############################################################################
# Run analysis
##############################################################################

final_CR_ECO_AMR_table <- process_cr_ecoli_amr_results(
  joined_data,
  Ecoli_EPI_CUTOFF
)

##############################################################################
# Export Excel
##############################################################################

write_xlsx(
  final_CR_ECO_AMR_table,
  "Results/CR_Ecoli_AMR_Resistance_Table.xlsx"
)

##############################################################################
# Export Word
##############################################################################

ft_compact <- flextable(
  final_CR_ECO_AMR_table
) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(
    size = 8,
    part = "all"
  ) %>%
  padding(
    padding = 1,
    part = "all"
  ) %>%
  merge_v(
    j = ~ Class
  ) %>%
  bold(
    part = "header"
  ) %>%
  set_table_properties(
    layout = "fixed"
  ) %>%
  width(
    width = 1.1
  ) %>%
  width(
    j = 1:2,
    width = 1.5
  )

sect_properties <- prop_section(
  page_size = page_size(
    orient = "landscape"
  ),
  page_margins = page_mar(
    bottom = 0.5,
    top = 0.5,
    right = 0.5,
    left = 0.5
  ),
  type = "continuous"
)

save_as_docx(
  ft_compact,
  path = "Results/CR_Ecoli_AMR_Results_Final_Fit.docx",
  pr_section = sect_properties
)ibution of Zone Diameters and Resistance Frequency in ESCR E. coli
##############################################################################

##############################################################################
# Create confirmed ESCR E. coli dataset
##############################################################################

ESCR_confirmed_data <- joined_data %>%
  select(
    INIKA_ID,
    Isolate,
    VITEK_MS_Results,
    `COLONY MORPHOLOGY ON C3GR`,
    AMX_ED10,
    AZM_ED15,
    CRO_ED30,
    CIP_ED5,
    DOX_ED30,
    FLR_ED30,
    GEN_ED10,
    MEM_ED10,
    OXY_ED30,
    POL_ED300,
    SXT_ED1_2,
    CTX_ED5,
    CTC_ED30
  ) %>%
  rename(
    Amoxicillin = AMX_ED10,
    Azithromycin = AZM_ED15,
    Ceftriaxone = CRO_ED30,
    Ciprofloxacin = CIP_ED5,
    Doxycycline = DOX_ED30,
    Florfenicol = FLR_ED30,
    Gentamicin = GEN_ED10,
    Meropenem = MEM_ED10,
    Oxytetracycline = OXY_ED30,
    `Polymyxin_B(PB)` = POL_ED300,
    `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
    Cefotaxime = CTX_ED5,
    `Cefotaxime/Clavulanic Acid` = CTC_ED30
  ) %>%
  mutate(
    across(
      all_of(abx_classes$Antimicrobial_substance),
      ~ as.numeric(trimws(.))
    )
  ) %>%
  filter(
    (
      (
        grepl("E", Isolate, ignore.case = TRUE) &
          grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
          grepl(
            "Pink|Red",
            `COLONY MORPHOLOGY ON C3GR`,
            ignore.case = TRUE
          )
      ) |
        (
          grepl("K", Isolate, ignore.case = TRUE) &
            grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
            grepl(
              "Metallic blue",
              `COLONY MORPHOLOGY ON C3GR`,
              ignore.case = TRUE
            )
        )
    ) &
      Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

message(
  paste(
    "ESCR E. coli isolates processed:",
    nrow(ESCR_confirmed_data)
  )
)

##############################################################################
# Calculate resistance frequencies and zone diameter distributions
##############################################################################

ESCR_long <- ESCR_confirmed_data %>%
  pivot_longer(
    cols = all_of(abx_classes$Antimicrobial_substance),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm)) %>%
  left_join(
    Ecoli_EPI_CUTOFF,
    by = "Antimicrobial_substance"
  ) %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    Has_Breakpoint = !is.na(R_limit),
    is_resistant = case_when(
      is.na(R_limit) ~ NA_real_,
      mm <= R_limit ~ 1,
      mm > R_limit ~ 0
    )
  )

##############################################################################
# Resistance Summary
##############################################################################

Summary_Resistance <- ESCR_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    N_defined = sum(!is.na(is_resistant)),
    R_count = sum(is_resistant, na.rm = TRUE),
    Has_Breakpoint = any(Has_Breakpoint),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    `Resistance % [95% CI]` = {
      
      if (!Has_Breakpoint || N_defined == 0) {
        
        "ND"
        
      } else {
        
        ci <- binom.test(
          as.integer(R_count),
          as.integer(N_defined)
        )$conf.int
        
        sprintf(
          "%.1f%% [%.1f-%.1f]",
          100 * R_count / N_defined,
          ci[1] * 100,
          ci[2] * 100
        )
      }
    }
  ) %>%
  ungroup() %>%
  select(
    Antimicrobial_substance,
    N,
    `Resistance % [95% CI]`
  )

##############################################################################
# Frequency Distribution of Zone Diameters
##############################################################################

Percentage_Wide <- ESCR_long %>%
  count(
    Antimicrobial_substance,
    mm
  ) %>%
  group_by(
    Antimicrobial_substance
  ) %>%
  mutate(
    Percentage = (n / sum(n)) * 100
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = Antimicrobial_substance,
    names_from = mm,
    values_from = Percentage,
    values_fill = 0
  )

##############################################################################
# Sort mm columns numerically
##############################################################################

mm_cols <- names(Percentage_Wide)[-1]

mm_cols <- mm_cols[
  order(as.numeric(mm_cols))
]

##############################################################################
# Combine Resistance Summary with mm Distribution
##############################################################################

Final_Formatted <- Summary_Resistance %>%
  left_join(
    Percentage_Wide,
    by = "Antimicrobial_substance"
  ) %>%
  select(
    Antimicrobial_substance,
    N,
    `Resistance % [95% CI]`,
    all_of(mm_cols)
  ) %>%
  mutate(
    across(
      all_of(mm_cols),
      ~{
        res <- round(.x, 1)
        
        ifelse(
          res == 0,
          "-",
          format(
            res,
            nsmall = 1,
            trim = TRUE
          )
        )
      }
    )
  )

##############################################################################
# Create and Style Table
##############################################################################

ESCR_ECO_AMR <- Final_Formatted %>%
  flextable() %>%
  set_table_properties(
    layout = "autofit",
    width = 1
  ) %>%
  theme_booktabs() %>%
  set_caption(
    caption = paste0(
      "Table: Resistance Percentages and ZOI Frequency Distribution for ESCR E. coli isolates (N=",
      nrow(ESCR_confirmed_data),
      ")."
    )
  ) %>%
  fontsize(
    size = 6.5,
    part = "all"
  ) %>%
  bold(part = "header") %>%
  padding(
    padding.top = 1,
    padding.bottom = 1,
    padding.left = 0,
    padding.right = 0,
    part = "all"
  ) %>%
  align(
    align = "center",
    j = 2:ncol(Final_Formatted),
    part = "all"
  ) %>%
  bg(
    j = 3,
    bg = "#F2F2F2",
    part = "body"
  ) %>%
  valign(
    valign = "bottom",
    part = "header"
  ) %>%
  add_footer_lines(
    values = "ND = No ECOFF or interpretive breakpoint available."
  )

##############################################################################
# Export
##############################################################################

if (!dir.exists("Results")) {
  dir.create("Results")
}

doc <- read_docx() %>%
  body_add_flextable(
    value = ESCR_ECO_AMR
  ) %>%
  body_end_section_landscape()

print(
  doc,
  target = "Results/ESCR_ECO_AMR_Frequency_Distribution_Table.docx"
)

save_as_docx(
  ESCR_ECO_AMR,
  path = "Results/ESCR_ECO_AMR_Frequency_Distribution_Table.docx"
)


################################################################################  
# Supplementary Table 4
## Computing the ESCR K.pneumoniae frequency Distribution of Zone diameter in millimeters (mm)
##############################################################################


##############################################################################
# Create confirmed ESCR K. pneumoniae dataset
##############################################################################

ESCR_Kpn_confirmed <- joined_data %>%
  select(
    INIKA_ID,
    Isolate,
    VITEK_MS_Results,
    `COLONY MORPHOLOGY ON C3GR`,
    AMX_ED10,
    AZM_ED15,
    CRO_ED30,
    CIP_ED5,
    DOX_ED30,
    FLR_ED30,
    GEN_ED10,
    MEM_ED10,
    OXY_ED30,
    POL_ED300,
    SXT_ED1_2,
    CTX_ED5,
    CTC_ED30
  ) %>%
  rename(
    Amoxicillin = AMX_ED10,
    Azithromycin = AZM_ED15,
    Ceftriaxone = CRO_ED30,
    Ciprofloxacin = CIP_ED5,
    Doxycycline = DOX_ED30,
    Florfenicol = FLR_ED30,
    Gentamicin = GEN_ED10,
    Meropenem = MEM_ED10,
    Oxytetracycline = OXY_ED30,
    `Polymyxin_B(PB)` = POL_ED300,
    `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
    Cefotaxime = CTX_ED5,
    `Cefotaxime/Clavulanic Acid` = CTC_ED30
  ) %>%
  mutate(
    across(
      all_of(abx_classes$Antimicrobial_substance),
      ~ as.numeric(trimws(.))
    )
  ) %>%
  filter(
    (
      (
        grepl("K", Isolate, ignore.case = TRUE) &
          grepl(
            "Klebsiella pneumoniae",
            VITEK_MS_Results,
            ignore.case = TRUE
          ) &
          grepl(
            "Metallic blue",
            `COLONY MORPHOLOGY ON C3GR`,
            ignore.case = TRUE
          )
      ) |
        (
          grepl("E", Isolate, ignore.case = TRUE) &
            grepl(
              "Klebsiella pneumoniae",
              VITEK_MS_Results,
              ignore.case = TRUE
            ) &
            grepl(
              "Pink|Red",
              `COLONY MORPHOLOGY ON C3GR`,
              ignore.case = TRUE
            )
        )
    ) &
      Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

message(
  paste(
    "ESCR K. pneumoniae isolates processed:",
    nrow(ESCR_Kpn_confirmed)
  )
)

##############################################################################
# Calculate resistance frequencies and zone diameter distributions
##############################################################################

ESCR_Kpn_long <- ESCR_Kpn_confirmed %>%
  pivot_longer(
    cols = all_of(abx_classes$Antimicrobial_substance),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm)) %>%
  left_join(
    Kpn_EPI_CUTOFF,
    by = "Antimicrobial_substance"
  ) %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    Has_Breakpoint = !is.na(R_limit),
    is_resistant = case_when(
      is.na(R_limit) ~ NA_real_,
      mm <= R_limit ~ 1,
      mm > R_limit ~ 0
    )
  )

##############################################################################
# Resistance Summary
##############################################################################

Summary_Resistance <- ESCR_Kpn_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    N_defined = sum(!is.na(is_resistant)),
    R_count = sum(is_resistant, na.rm = TRUE),
    Has_Breakpoint = any(Has_Breakpoint),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    `Resistance % [95% CI]` = {
      
      if (!Has_Breakpoint || N_defined == 0) {
        
        "ND"
        
      } else {
        
        ci <- binom.test(
          as.integer(R_count),
          as.integer(N_defined)
        )$conf.int
        
        sprintf(
          "%.1f%% [%.1f-%.1f]",
          100 * R_count / N_defined,
          ci[1] * 100,
          ci[2] * 100
        )
      }
    }
  ) %>%
  ungroup() %>%
  select(
    Antimicrobial_substance,
    N,
    `Resistance % [95% CI]`
  )

##############################################################################
# Frequency Distribution of Zone Diameters
##############################################################################

Percentage_Wide <- ESCR_Kpn_long %>%
  count(
    Antimicrobial_substance,
    mm
  ) %>%
  group_by(
    Antimicrobial_substance
  ) %>%
  mutate(
    Percentage = (n / sum(n)) * 100
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = Antimicrobial_substance,
    names_from = mm,
    values_from = Percentage,
    values_fill = 0
  )

##############################################################################
# Sort mm columns numerically
##############################################################################

mm_cols <- names(Percentage_Wide)[-1]
mm_cols <- mm_cols[
  order(as.numeric(mm_cols))
]

##############################################################################
# Combine Resistance Summary with mm Distribution
##############################################################################

Final_Kpn_Formatted <- Summary_Resistance %>%
  left_join(
    Percentage_Wide,
    by = "Antimicrobial_substance"
  ) %>%
  select(
    Antimicrobial_substance,
    N,
    `Resistance % [95% CI]`,
    all_of(mm_cols)
  ) %>%
  mutate(
    across(
      all_of(mm_cols),
      ~{
        res <- round(.x, 1)
        
        ifelse(
          res == 0,
          "-",
          format(
            res,
            nsmall = 1,
            trim = TRUE
          )
        )
      }
    )
  )

##############################################################################
# Create and Style Table
##############################################################################

final_ESCR_Kpn_table <- Final_Kpn_Formatted %>%
  flextable() %>%
  set_table_properties(
    layout = "autofit",
    width = 1
  ) %>%
  theme_booktabs() %>%
  set_caption(
    caption = paste0(
      "Table: Resistance Percentages and ZOI Frequency Distribution for ESCR K. pneumoniae isolates (N=",
      nrow(ESCR_Kpn_confirmed),
      ")."
    )
  ) %>%
  fontsize(
    size = 6.5,
    part = "all"
  ) %>%
  bold(part = "header") %>%
  padding(
    padding.top = 1,
    padding.bottom = 1,
    padding.left = 0,
    padding.right = 0,
    part = "all"
  ) %>%
  align(
    align = "center",
    j = 2:ncol(Final_Kpn_Formatted),
    part = "all"
  ) %>%
  bg(
    j = 3,
    bg = "#F2F2F2",
    part = "body"
  ) %>%
  valign(
    valign = "bottom",
    part = "header"
  ) %>%
  add_footer_lines(
    values = "ND = No ECOFF or interpretive breakpoint available."
  )

##############################################################################
# Export
##############################################################################

if (!dir.exists("Results")) {
  dir.create("Results")
}

doc <- read_docx() %>%
  body_add_flextable(
    value = final_ESCR_Kpn_table
  ) %>%
  body_end_section_landscape()

print(
  doc,
  target = "Results/Kpn_Resistance_Frequency_Table.docx"
)

save_as_docx(
  final_ESCR_Kpn_table,
  path = "Results/Kpn_Resistance_Frequency_Table.docx"
)

##########################################################################




##############################################################################
# Calculate resistance frequencies
##############################################################################

ESCR_long <- ESCR_confirmed_data %>%
  pivot_longer(
    cols = all_of(abx_classes$Antimicrobial_substance),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm)) %>%
  left_join(
    Ecoli_EPI_CUTOFF,
    by = "Antimicrobial_substance"
  ) %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    Has_Breakpoint = !is.na(R_limit),
    is_resistant = case_when(
      is.na(R_limit) ~ NA_real_,
      mm <= R_limit ~ 1,
      mm > R_limit ~ 0
    )
  )

##############################################################################
# Resistance summary
##############################################################################

Summary_Resistance <- ESCR_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N_defined = sum(!is.na(is_resistant)),
    R_count = sum(is_resistant, na.rm = TRUE),
    Has_Breakpoint = any(Has_Breakpoint),
    .groups = "drop"
  ) %>%
  filter(
    Has_Breakpoint,
    N_defined > 0
  ) %>%
  rowwise() %>%
  mutate(
    Prop = R_count / N_defined,
    Lower = binom.test(
      as.integer(R_count),
      as.integer(N_defined)
    )$conf.int[1] * 100,
    Upper = binom.test(
      as.integer(R_count),
      as.integer(N_defined)
    )$conf.int[2] * 100
  ) %>%
  ungroup() %>%
  filter(Prop > 0)

##############################################################################
# Generate refined plot
##############################################################################

AMR_CI_Plot <- ggplot(
  Summary_Resistance,
  aes(
    x = reorder(
      Antimicrobial_substance,
      Prop
    ),
    y = Prop * 100,
    fill = Antimicrobial_substance
  )
) +
  geom_col(
    color = "black",
    alpha = 0.9,
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.25,
    color = "gray20",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", Prop * 100)
    ),
    hjust = 1.2,
    colour = "white",
    size = 3.5,
    fontface = "bold"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Resistance Percentage (%)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

##############################################################################
# Save and export
##############################################################################

if (!dir.exists("Results")) {
  dir.create("Results")
}

ggsave(
  "Results/Resistance_Refined_Plot.png",
  plot = AMR_CI_Plot,
  width = 9,
  height = 6,
  dpi = 300
)

doc <- read_docx() %>%
  body_add_img(
    src = "Results/Resistance_Refined_Plot.png",
    width = 6.5,
    height = 4.3
  )

print(
  doc,
  target = "Results/ESCR_AMR_Refined_Histogram.docx"
)

################################################################################   
# Figure 2
##############################################################################
# Calculate resistance frequencies
##############################################################################

ESCR_Kpn_long <- ESCR_Kpn_confirmed %>%
  pivot_longer(
    cols = all_of(abx_classes$Antimicrobial_substance),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm)) %>%
  left_join(
    Kpn_EPI_CUTOFF,
    by = "Antimicrobial_substance"
  ) %>%
  mutate(
    mm = as.numeric(str_trim(mm)),
    R_limit = as.numeric(str_replace_all(R, "[^0-9.]", "")),
    Has_Breakpoint = !is.na(R_limit),
    is_resistant = case_when(
      is.na(R_limit) ~ NA_real_,
      mm <= R_limit ~ 1,
      mm > R_limit ~ 0
    )
  )

##############################################################################
# Resistance Summary
##############################################################################

Summary_Resistance <- ESCR_Kpn_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N_defined = sum(!is.na(is_resistant)),
    R_count = sum(is_resistant, na.rm = TRUE),
    Has_Breakpoint = any(Has_Breakpoint),
    .groups = "drop"
  ) %>%
  filter(
    Has_Breakpoint,
    N_defined > 0
  ) %>%
  rowwise() %>%
  mutate(
    Prop = R_count / N_defined,
    Lower = binom.test(
      as.integer(R_count),
      as.integer(N_defined)
    )$conf.int[1] * 100,
    Upper = binom.test(
      as.integer(R_count),
      as.integer(N_defined)
    )$conf.int[2] * 100
  ) %>%
  ungroup() %>%
  filter(Prop > 0)

##############################################################################
# Generate Refined Plot
##############################################################################

Kpn_AMR_CI_Plot <- ggplot(
  Summary_Resistance,
  aes(
    x = reorder(
      Antimicrobial_substance,
      Prop
    ),
    y = Prop * 100,
    fill = Antimicrobial_substance
  )
) +
  geom_col(
    color = "black",
    alpha = 0.9,
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.25,
    color = "gray20",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", Prop * 100)
    ),
    hjust = 1.2,
    colour = "white",
    size = 3.5,
    fontface = "bold"
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10)
  ) +
  labs(
    title = "Resistance Percentage (%)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

##############################################################################
# Save and Export
##############################################################################

if (!dir.exists("Results")) {
  dir.create("Results")
}

ggsave(
  "Results/Kpn_Resistance_Refined_Plot.png",
  plot = Kpn_AMR_CI_Plot,
  width = 9,
  height = 6,
  dpi = 300
)

doc <- read_docx() %>%
  body_add_img(
    src = "Results/Kpn_Resistance_Refined_Plot.png",
    width = 6.5,
    height = 4.3
  )

print(
  doc,
  target = "Results/ESCR_Kpn_AMR_Refined_Histogram.docx"
)
################################################################################
# ## Computing for ESBL count
# # 1. Create 'Results' folder if it does not exist
# if (!dir.exists("Results")) {
#   dir.create("Results")
# }
# 
# # Helper function: formats n (%) to 1 decimal place, stripping '.0' for whole numbers
# fmt_n_pct <- function(n, total) {
#   if (total == 0) return("0 (0%)")
#   pct <- (n / total) * 100
#   pct_str <- sub("\\.0$", "", sprintf("%.1f", pct)) # Replaces .0 with empty string
#   sprintf("%d (%s%%)", n, pct_str)
# }
# # General summarizer for Origin, Season, and Region
# summarize_variable <- function(df, var_col, var_label) {
#   var_sym <- sym(var_col)
#   
#   df %>%
#     filter(!is.na(!!var_sym)) %>%
#     group_by(Values = as.character(!!var_sym)) %>%
#     summarise(
#       n_ecoli = sum(Species_Group == "E. coli", na.rm = TRUE),
#       n_kpneumo = sum(Species_Group == "K. pneumoniae", na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(
#       `ESBL E. coli` = map2_chr(n_ecoli, total_ecoli, fmt_n_pct),
#       `ESBL K. pneumoniae` = map2_chr(n_kpneumo, total_kpneumo, fmt_n_pct),
#       Variable = var_label
#     ) %>%
#     select(Variable, Values, `ESBL E. coli`, `ESBL K. pneumoniae`)
# }
# # Dedicated District summarizer sorted by Region (Kilimanjaro -> Mwanza), then District alphabetically
# summarize_districts <- function(df, district_col, region_col) {
#   dist_sym <- sym(district_col)
#   reg_sym <- sym(region_col)
#   
#   df %>%
#     filter(!is.na(!!dist_sym), !is.na(!!reg_sym)) %>%
#     group_by(Region = as.character(!!reg_sym), District = as.character(!!dist_sym)) %>%
#     summarise(
#       n_ecoli = sum(Species_Group == "E. coli", na.rm = TRUE),
#       n_kpneumo = sum(Species_Group == "K. pneumoniae", na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     # Order Region so Kilimanjaro comes first, Mwanza second
#     mutate(Region = factor(Region, levels = c("Kilimanjaro", "Mwanza"))) %>%
#     # Sort by Region order, then District name alphabetically
#     arrange(Region, District) %>%
#     mutate(
#       `ESBL E. coli` = map2_chr(n_ecoli, total_ecoli, fmt_n_pct),
#       `ESBL K. pneumoniae` = map2_chr(n_kpneumo, total_kpneumo, fmt_n_pct),
#       Variable = "District",
#       Values = District
#     ) %>%
#     select(Variable, Values, `ESBL E. coli`, `ESBL K. pneumoniae`)
# }
# # 2. Filter dataset for presumptive ESBL positives & deduplicate by unique INIKA_ID
# # Define species filtering conditions and deduplicate
# esbl_data <- joined_data %>%
#   rename(
#     Amoxicillin = AMX_ED10,
#     Azithromycin = AZM_ED15,
#     Ceftriaxone = CRO_ED30,
#     Ciprofloxacin = CIP_ED5,
#     Doxycycline = DOX_ED30,
#     Florfenicol = FLR_ED30,
#     Gentamicin = GEN_ED10,
#     Meropenem = MEM_ED10,
#     Oxytetracycline = OXY_ED30,
#     `Polymyxin_B(PB)` = POL_ED300,
#     `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
#     Cefotaxime = CTX_ED5,
#     `Cefotaxime/Clavulanic Acid` = CTC_ED30) %>%
#   mutate(across(all_of(abx_cols), ~ as.numeric(trimws(.)))) %>%
#   # Basic phenotypic resistance criteria
#   filter( `Cefotaxime/Clavulanic Acid`- Cefotaxime >= 5,
#     #ESBL == 1,
#     Ceftriaxone < 23 
#   ) %>%
# 
#   # Apply strict microbiological identification & chromogenic colony rules
#   filter(
#     # Condition A: E. coli Rules 
#     (
#       (
#         grepl("E", Isolate, ignore.case = TRUE) &
#           grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
#           grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
#       ) |
#         (
#           grepl("K", Isolate, ignore.case = TRUE) &
#             grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) &
#             grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
#         )
#     ) |
#       # Condition B: K. pneumoniae Rule 
#       (
#         trimws(VITEK_MS_Results) == "Klebsiella pneumoniae"
#       )
#   ) %>%
#   # Assign definitive Species_Group for table aggregation
#   mutate(
#     Species_Group = case_when(
#       grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) ~ "E. coli",
#       trimws(VITEK_MS_Results) == "Klebsiella pneumoniae" ~ "K. pneumoniae",
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   # Sort deterministically and retain unique INIKA_ID
#   arrange(INIKA_ID, Isolate) %>%
#   distinct(INIKA_ID, .keep_all = TRUE)
# 
# # Print logging summary messages to console
# n_ecoli <- sum(esbl_data$Species_Group == "E. coli", na.rm = TRUE)
# n_kpneumo <- sum(esbl_data$Species_Group == "K. pneumoniae", na.rm = TRUE)
# n_total <- nrow(esbl_data)
# 
# message(sprintf("ESCR E. coli isolates retained: %d", n_ecoli))
# message(sprintf("ESBL K. pneumoniae isolates retained: %d", n_kpneumo))
# message(sprintf("Total unique isolates processed: %d", n_total))
# # Compute total counts per species across unique INIKA_ID
# total_ecoli <- sum(esbl_data$Species_Group == "E. coli")
# total_kpneumo <- sum(esbl_data$Species_Group == "K. pneumoniae")
# 
# # Helper function to summarize species distributions per variable category
# summarize_variable <- function(df, var_col, var_label) {
#   var_sym <- sym(var_col)
#   
#   df %>%
#     filter(!is.na(!!var_sym)) %>%
#     group_by(Values = as.character(!!var_sym)) %>%
#     summarise(
#       n_ecoli = sum(Species_Group == "E. coli"),
#       n_kpneumo = sum(Species_Group == "K. pneumoniae"),
#       .groups = "drop"
#     ) %>%
#     mutate(
#       `ESBL E. coli` = map2_chr(n_ecoli, total_ecoli, fmt_n_pct),
#       `ESBL K. pneumoniae` = map2_chr(n_kpneumo, total_kpneumo, fmt_n_pct),
#       Variable = var_label
#     ) %>%
#     select(Variable, Values, `ESBL E. coli`, `ESBL K. pneumoniae`)
# }
# 
# # 3. Aggregate across Origin of sample, Season, Region, and District
# table_df <- bind_rows(
#   tibble(
#     Variable = "Overall",
#     Values = "Total Unique INIKA_ID",
#     `ESBL E. coli` = fmt_n_pct(total_ecoli, total_ecoli),       
#     `ESBL K. pneumoniae` = fmt_n_pct(total_kpneumo, total_kpneumo) 
#   ),
#   summarize_variable(esbl_data, "ORIGIN_OF_SAMPLE", "Origin"),
#   summarize_variable(esbl_data, "SEASON.x", "Season"),
#   summarize_variable(esbl_data, "REGION.x", "Region"),
#   summarize_districts(esbl_data, district_col = "DISTRICT.x", region_col = "REGION.x")
# )
# 
# # 4. Construct Publication-Ready Flextable
# esbl_table <- table_df %>%
#   flextable() %>%
#   set_header_labels(
#     Variable = "Variable",
#     Values = "Category / Value",
#     `ESBL E. coli` = sprintf("ESBL E. coli\n(N = %d), n (%%)", total_ecoli),
#     `ESBL K. pneumoniae` = sprintf("ESBL K. pneumoniae\n(N = %d), n (%%)", total_kpneumo)
#   ) %>%
#   merge_v(j = "Variable") %>%
#   valign(j = "Variable", valign = "top") %>% 
#   theme_booktabs() %>%
#   bold(part = "header") %>%
#   bold(i = ~ Values == "Total Unique INIKA_ID") %>%
#   align(j = c("ESBL E. coli", "ESBL K. pneumoniae"), align = "center", part = "all") %>%
#   align(j = c("Variable", "Values"), align = "left", part = "all") %>%
#   add_header_lines(
#     values = "Table 1: Distribution of ESBL E. coli and ESBL K. pneumoniae Isolates by Sample Origin, Season, Region, and District"
#   ) %>%
#   bold(i = 1, part = "header") %>%
#   add_footer_lines(
#     sprintf(
#       "Note: Filtered by ESBL == 1, Ceftriaxone < 23 mm, and deduplicated by unique INIKA_ID. Percentages calculated out of total ESBL E. coli (N = %d) and ESBL K. pneumoniae (N = %d).",
#       total_ecoli, total_kpneumo
#     )
#   ) %>%
#   italic(j = c("ESBL E. coli", "ESBL K. pneumoniae"), part = "header") %>%
#   font(fontname = "Times New Roman", part = "all") %>%
#   fontsize(size = 10, part = "all") %>%
#   autofit()
# 
# # 5. Export table to Word document in Results/ directory
# doc <- read_docx() %>%
#   body_add_flextable(esbl_table)
# 
# print(doc, target = "Results/ESBL_Ecoli_Kpneumoniae_Summary_Table.docx")
#######################################################################
## Counting for the Carbapenem Resistant
#  Grouped count (if you want to see counts per Isolate or per Colony Morphology)
df_grouped_count <- joined_data %>%
  filter(
    str_detect(`COLONY MORPHOLOGY ON CARBA`, "(?i)Red|Pinkish|Metalic blue")
  ) %>%
  distinct(INIKA_ID, Isolate, `COLONY MORPHOLOGY ON CARBA`) %>%
  count(Isolate, name = "unique_INIKA_count") # E.coli 11
###############################################################################
#  Grouped breakdown by Isolate type 
Joined_data_results_by_isolate <- joined_data %>%
  filter(
    str_detect(`COLONY MORPHOLOGY ON CARBA`, "(?i)Red|Pinkish|Metalic blue")
  ) %>%
  group_by(Isolate) %>%
  summarise(
    unique_INIKA_total = n_distinct(INIKA_ID),
    with_VITEK_MS = n_distinct(INIKA_ID[!is.na(VITEK_MS_Results) & VITEK_MS_Results != ""]),
    with_Isolate_NVI = n_distinct(INIKA_ID[!is.na(Isolate_NVI) & Isolate_NVI != ""]),
    .groups = "drop"
  )
################################################################################
# ## Computing the CR Isolation rate,(n,%)
# # ------------------------------------------------------------------------------
# # Helper function for proportions and Wilson 95% CI calculation
# # ------------------------------------------------------------------------------
# calc_rate_ci <- function(n_pos, n_total) {
#   if (n_total == 0) {
#     return(data.frame(
#       cr_n_pct = "0 (0%)",
#       ci_95    = "0 - 0"
#     ))
#   }
# 
#   # Format helper: Outputs 1 decimal place, or integer if ending in .0
#   fmt_num <- function(val) {
#     val_rounded <- round(val, 1)
#     if (val_rounded %% 1 == 0) {
#       sprintf("%.0f", val_rounded)
#     } else {
#       sprintf("%.1f", val_rounded)
#     }
#   }
# 
#   ci_res <- binom::binom.confint(x = n_pos, n = n_total, methods = "wilson")
# 
#   pct   <- (n_pos / n_total) * 100
#   lower <- ci_res$lower * 100
#   upper <- ci_res$upper * 100
# 
#   data.frame(
#     cr_n_pct = sprintf("%d (%s%%)", n_pos, fmt_num(pct)),
#     ci_95    = sprintf("%s - %s", fmt_num(lower), fmt_num(upper))
#   )
# }
# 
# # ------------------------------------------------------------------------------
# # DATA CLEANING & RECODING
# # ------------------------------------------------------------------------------
# Joined_data_clean <- joined_data %>%
#   select(INIKA_ID,ORIGIN_OF_SAMPLE, SEASON.x,REGION.x,DISTRICT.x,`COLONY MORPHOLOGY ON CARBA`,
#          VITEK_MS_Results,Isolate,PROTOCOL,Isolate_NVI) %>%
#   rename('SEASON' = SEASON.x, 'REGION' = REGION.x,'DISTRICT' = DISTRICT.x)%>%
#   # Recode category values
#   rename(
#     `Sample origin` = ORIGIN_OF_SAMPLE,
#     `Season`        = SEASON,
#     `Region`        = REGION,
#     `District`      = DISTRICT
#   )
# 
# # Base unique dataset for total tested counts (Denominator)
# denom_df <- Joined_data_clean %>%
#   group_by(INIKA_ID) %>%
#   slice(1) %>%
#   ungroup()
# 
# # Target isolates count matching criteria (Numerator)
# cr_isolates_df <- Joined_data_clean %>%
#   filter(str_detect(`COLONY MORPHOLOGY ON CARBA`, "(?i)Red|Pinkish|Metalic blue")) %>%
#   group_by(INIKA_ID) %>%
#   slice(1) %>%
#   ungroup()
# 
# # Function to compute row metrics for a given variable and grouping column
# compute_strata_table <- function(data_denom, data_cr, var_label, group_col) {
#   categories <- data_denom %>%
#     pull({{ group_col }}) %>%
#     unique() %>%
#     na.omit()
# 
#   map_dfr(categories, function(cat) {
#     n_tot <- data_denom %>% filter({{ group_col }} == cat) %>% pull(INIKA_ID) %>% n_distinct()
#     n_pos <- data_cr %>% filter({{ group_col }} == cat) %>% pull(INIKA_ID) %>% n_distinct()
# 
#     ci_stats <- calc_rate_ci(n_pos, n_tot)
# 
#     tibble(
#       Variable = var_label,
#       Category = as.character(cat),
#       `Total test` = n_tot,
#       `CR E. coli (n, %)` = ci_stats$cr_n_pct,
#       `95% CI` = ci_stats$ci_95
#     )
#   })
# }
# 
# # ------------------------------------------------------------------------------
# # Overall Isolation Rate
# # ------------------------------------------------------------------------------
# total_tested_overall <- n_distinct(denom_df$INIKA_ID)
# total_cr_overall <- n_distinct(cr_isolates_df$INIKA_ID)
# overall_ci <- calc_rate_ci(total_cr_overall, total_tested_overall)
# 
# overall_row <- tibble(
#   Variable = "Overall",
#   Category = "Overall",
#   `Total test` = total_tested_overall,
#   `CR E. coli (n, %)` = overall_ci$cr_n_pct,
#   `95% CI` = overall_ci$ci_95
# )
# 
# # ------------------------------------------------------------------------------
# # Stratified Computations
# # ------------------------------------------------------------------------------
# sample_origin_rows <- compute_strata_table(denom_df, cr_isolates_df, "Sample origin", `Sample origin`)
# season_rows        <- compute_strata_table(denom_df, cr_isolates_df, "Season", Season)
# region_rows        <- compute_strata_table(denom_df, cr_isolates_df, "Region", Region)
# 
# # District sorting: Kilimanjaro districts first (alphabetical), then Mwanza districts (alphabetical)
# kilimanjaro_districts <- denom_df %>%
#   filter(Region == "Kilimanjaro") %>%
#   pull(District) %>%
#   unique() %>%
#   na.omit() %>%
#   sort()
# 
# mwanza_districts <- denom_df %>%
#   filter(Region == "Mwanza") %>%
#   pull(District) %>%
#   unique() %>%
#   na.omit() %>%
#   sort()
# 
# ordered_districts <- c(kilimanjaro_districts, mwanza_districts)
# 
# district_rows <- map_dfr(ordered_districts, function(dist) {
#   n_tot <- denom_df %>% filter(District == dist) %>% pull(INIKA_ID) %>% n_distinct()
#   n_pos <- cr_isolates_df %>% filter(District == dist) %>% pull(INIKA_ID) %>% n_distinct()
# 
#   ci_stats <- calc_rate_ci(n_pos, n_tot)
# 
#   tibble(
#     Variable = "District",
#     Category = as.character(dist),
#     `Total test` = n_tot,
#     `CR E. coli (n, %)` = ci_stats$cr_n_pct,
#     `95% CI` = ci_stats$ci_95
#   )
# })
# 
# # ------------------------------------------------------------------------------
# # Combine and Format Final Table
# # ------------------------------------------------------------------------------
# final_table_df <- bind_rows(
#   overall_row,
#   sample_origin_rows,
#   season_rows,
#   region_rows,
#   district_rows
# )
# 
# # Blank out repeated Variable labels for clean publication display
# final_table_formatted <- final_table_df %>%
#   mutate(
#     Variable = if_else(duplicated(Variable), "", Variable)
#   )
# 
# # ------------------------------------------------------------------------------
# # Build Flextable and Save to Word Document
# # ------------------------------------------------------------------------------
# if (!dir.exists("Results")) {
#   dir.create("Results")
# }
# 
# ft <- flextable(final_table_formatted) %>%
#   set_header_labels(
#     Variable = "Variable",
#     Category = "Category",
#     `Total test` = "Total test",
#     `CR E. coli (n, %)` = "CR E. coli (n, %)",
#     `95% CI` = "95% CI"
#   ) %>%
#   autofit() %>%
#   theme_booktabs() %>%
#   align(j = c("Total test", "CR E. coli (n, %)", "95% CI"), align = "center", part = "all") %>%
#   bold(i = ~ Variable != "", j = "Variable", bold = TRUE) %>%
#   italic(j = "CR E. coli (n, %)", part = "header") %>%  # Optional: Italics for species name
#   fontsize(size = 10, part = "all") %>%
#   font(fontname = "Times New Roman", part = "all")
# 
# # Write to Word file in 'Results' directory
# doc <- read_docx() %>%
#   body_add_par("Table: Carbapenem-Resistant (CR) E. coli Isolation Rates Stratified by Socio-Demographic and Geographic Characteristics", style = "heading 2") %>%
#   body_add_flextable(ft)
# 
# output_path <- file.path("Results", "CR_E_coli_Isolation_Rates_Table.docx")
# print(doc, target = output_path)
# 
# message("Table successfully compiled and saved to: ", output_path)
################################################################################
### Carbapenem-Resistant (CR) Enterobacterales Isolation and Species Confirmation 
# Helper function for proportions and Wilson 95% CI calculation
# ------------------------------------------------------------------------------
calc_rate_ci <- function(n_pos, n_total) {
  if (n_total == 0) {
    return(list(
      cr_n_pct = "0 (0%)",
      ci_95    = "0 - 0",
      combined = "0 (0%) [0 - 0]"
    ))
  }
  
  ci_res <- binom::binom.confint(x = n_pos, n = n_total, methods = "wilson")
  
  pct   <- (n_pos / n_total) * 100
  lower <- max(0, ci_res$lower * 100)
  upper <- min(100, ci_res$upper * 100)
  
  n_pct_str <- sprintf("%d (%s%%)", n_pos, fmt_num(pct))
  ci_str    <- sprintf("%s - %s", fmt_num(lower), fmt_num(upper))
  
  list(
    cr_n_pct = n_pct_str,
    ci_95    = ci_str,
    combined = sprintf("%s [%s]", n_pct_str, ci_str)
  )
}

# Rounds to 2 decimal places, omits .00 for whole numbers
fmt_num <- function(val) {
  if (is.na(val) || val == 0) return("0")
  
  val_rounded <- round(val, 2)
  
  # Format with up to 2 decimals, then trim unnecessary trailing zeros
  formatted <- sprintf("%.2f", val_rounded)
  formatted <- sub("\\.00$", "", formatted)          # e.g., 3.00 -> 3
  formatted <- sub("(\\.[1-9])0$", "\\1", formatted) # e.g., 3.50 -> 3.5
  
  return(formatted)
}

# ------------------------------------------------------------------------------
# DATA CLEANING & ISOLATE CONFIRMATION RECODING
# ------------------------------------------------------------------------------
Joined_data_clean <- joined_data %>%
  select(INIKA_ID, ORIGIN_OF_SAMPLE, SEASON.x, REGION.x, DISTRICT.x, `COLONY MORPHOLOGY ON CARBA`,
         VITEK_MS_Results, Isolate, PROTOCOL, Isolate_NVI) %>%
  rename('SEASON' = SEASON.x, 'REGION' = REGION.x, 'DISTRICT' = DISTRICT.x) %>%
  mutate(
    is_initial_ecoli = str_detect(coalesce(Isolate, ""), "(?i)E\\.?\\s*coli"),
    is_vitek_ecoli   = str_detect(coalesce(VITEK_MS_Results, ""), "(?i)Escherichia\\s+coli"),
    is_nvi_ecoli     = str_detect(coalesce(Isolate_NVI, ""), "(?i)E\\.?\\s*coli|Escherichia\\s+coli"),
    is_nvi_kpneuma   = str_detect(coalesce(Isolate_NVI, ""), "(?i)K\\.?\\s*pneumoniae|Klebsiella\\s+pneumoniae"),
    
    CONFIRMED_SPECIES = case_when(
      is_initial_ecoli & (is_vitek_ecoli | is_nvi_ecoli) ~ "CR Escherichia coli",
      is_initial_ecoli & is_nvi_kpneuma                  ~ "CR Klebsiella pneumoniae",
      TRUE ~ NA_character_
    )
  ) %>%
  rename(
    `Sample origin` = ORIGIN_OF_SAMPLE,  
    `Season`        = SEASON,
    `Region`        = REGION,
    `District`      = DISTRICT
  )

denom_df <- Joined_data_clean %>%
  group_by(INIKA_ID) %>%
  slice(1) %>%
  ungroup()

cr_isolates_df <- Joined_data_clean %>%
  filter(str_detect(coalesce(`COLONY MORPHOLOGY ON CARBA`, ""), "(?i)Red|Pinkish|Metalic blue")) %>%
  group_by(INIKA_ID) %>%
  slice(1) %>%
  ungroup()

# ------------------------------------------------------------------------------
# STRATA COMPUTATION FUNCTION
# ------------------------------------------------------------------------------
compute_strata_table <- function(data_denom, data_cr, var_label, group_col) {
  categories <- data_denom %>%
    pull({{ group_col }}) %>%
    unique() %>%
    na.omit()
  
  map_dfr(categories, function(cat) {
    n_tot <- data_denom %>% filter({{ group_col }} == cat) %>% pull(INIKA_ID) %>% n_distinct()
    
    data_cr_sub <- data_cr %>% filter({{ group_col }} == cat)
    n_pos <- data_cr_sub %>% pull(INIKA_ID) %>% n_distinct()
    
    ci_stats <- calc_rate_ci(n_pos, n_tot)
    
    n_ecoli <- data_cr_sub %>% filter(CONFIRMED_SPECIES == "CR Escherichia coli") %>% pull(INIKA_ID) %>% n_distinct()
    n_kpneu <- data_cr_sub %>% filter(CONFIRMED_SPECIES == "CR Klebsiella pneumoniae") %>% pull(INIKA_ID) %>% n_distinct()
    
    ci_ecoli <- calc_rate_ci(n_ecoli, n_tot)
    ci_kpneu <- calc_rate_ci(n_kpneu, n_tot)
    
    tibble(
      Variable = var_label,
      Category = as.character(cat),
      `Total test` = n_tot,
      `pCR Escherichia coli n (%) [95% CI]` = ci_stats$combined,
      `CR Escherichia coli n (%) [95% CI]` = ci_ecoli$combined,
      `CR Klebsiella pneumoniae n (%) [95% CI]` = ci_kpneu$combined
    )
  })
}

# ------------------------------------------------------------------------------
# COMPUTE TABLE SECTIONS
# ------------------------------------------------------------------------------
total_tested_overall <- n_distinct(denom_df$INIKA_ID)
total_cr_overall    <- n_distinct(cr_isolates_df$INIKA_ID)
overall_ci          <- calc_rate_ci(total_cr_overall, total_tested_overall)

overall_ecoli    <- cr_isolates_df %>% filter(CONFIRMED_SPECIES == "CR Escherichia coli") %>% pull(INIKA_ID) %>% n_distinct()
overall_kpneu    <- cr_isolates_df %>% filter(CONFIRMED_SPECIES == "CR Klebsiella pneumoniae") %>% pull(INIKA_ID) %>% n_distinct()

overall_ci_ecoli <- calc_rate_ci(overall_ecoli, total_tested_overall)
overall_ci_kpneu <- calc_rate_ci(overall_kpneu, total_tested_overall)

overall_row <- tibble(
  Variable = "Overall",
  Category = "Overall",
  `Total test` = total_tested_overall,
  `pCR Escherichia coli n (%) [95% CI]` = overall_ci$combined,
  `CR Escherichia coli n (%) [95% CI]` = overall_ci_ecoli$combined,
  `CR Klebsiella pneumoniae n (%) [95% CI]` = overall_ci_kpneu$combined
)

sample_origin_rows <- compute_strata_table(denom_df, cr_isolates_df, "Sample origin", `Sample origin`)
season_rows        <- compute_strata_table(denom_df, cr_isolates_df, "Season", Season)
region_rows        <- compute_strata_table(denom_df, cr_isolates_df, "Region", Region)

kilimanjaro_districts <- denom_df %>% filter(Region == "Kilimanjaro") %>% pull(District) %>% unique() %>% na.omit() %>% sort()
mwanza_districts      <- denom_df %>% filter(Region == "Mwanza") %>% pull(District) %>% unique() %>% na.omit() %>% sort()
ordered_districts     <- c(kilimanjaro_districts, mwanza_districts)

district_rows <- map_dfr(ordered_districts, function(dist) {
  n_tot <- denom_df %>% filter(District == dist) %>% pull(INIKA_ID) %>% n_distinct()
  data_cr_sub <- cr_isolates_df %>% filter(District == dist)
  n_pos <- data_cr_sub %>% pull(INIKA_ID) %>% n_distinct()
  
  ci_stats <- calc_rate_ci(n_pos, n_tot)
  
  n_ecoli <- data_cr_sub %>% filter(CONFIRMED_SPECIES == "CR Escherichia coli") %>% pull(INIKA_ID) %>% n_distinct()
  n_kpneu <- data_cr_sub %>% filter(CONFIRMED_SPECIES == "CR Klebsiella pneumoniae") %>% pull(INIKA_ID) %>% n_distinct()
  
  ci_ecoli <- calc_rate_ci(n_ecoli, n_tot)
  ci_kpneu <- calc_rate_ci(n_kpneu, n_tot)
  
  tibble(
    Variable = "District",
    Category = as.character(dist),
    `Total test` = n_tot,
    `pCR Escherichia coli n (%) [95% CI]` = ci_stats$combined,
    `CR Escherichia coli n (%) [95% CI]` = ci_ecoli$combined,
    `CR Klebsiella pneumoniae n (%) [95% CI]` = ci_kpneu$combined
  )
})

final_table_df <- bind_rows(
  overall_row,
  sample_origin_rows,
  season_rows,
  region_rows,
  district_rows
)

# Identify non-duplicated variable positions for formatting before blanking them out
bold_rows <- which(!duplicated(final_table_df$Variable))

final_table_formatted <- final_table_df %>%
  mutate(Variable = if_else(duplicated(Variable), "", Variable))

# ------------------------------------------------------------------------------
# FLEXTABLE GENERATION
# ------------------------------------------------------------------------------
if (!dir.exists("Results")) {
  dir.create("Results")
}

ft <- flextable(final_table_formatted) %>%
  set_header_labels(
    Variable = "Variable",
    Category = "Category",
    `Total test` = "Total test",
    `pCR Escherichia coli n (%) [95% CI]` = "pCR Escherichia coli\nn (%) [95% CI]",
    `CR Escherichia coli n (%) [95% CI]` = "CR Escherichia coli\nn (%) [95% CI]",
    `CR Klebsiella pneumoniae n (%) [95% CI]` = "CR Klebsiella pneumoniae\nn (%) [95% CI]"
  ) %>%
  add_header_row(
    top = TRUE,
    values = c(
      "Variable", 
      "Category", 
      "Total test", 
      "Presumptive CR Isolates", 
      "Confirmed isolate", 
      "Confirmed isolate"
    )
  ) %>%
  merge_h(part = "header") %>%
  merge_v(part = "header") %>%
  compose(
    i = 2, j = "pCR Escherichia coli n (%) [95% CI]",
    value = as_paragraph("pCR ", as_i("Escherichia coli"), "\nn (%) [95% CI]"),
    part = "header"
  ) %>%
  compose(
    i = 2, j = "CR Escherichia coli n (%) [95% CI]",
    value = as_paragraph("CR ", as_i("Escherichia coli"), "\nn (%) [95% CI]"),
    part = "header"
  ) %>%
  compose(
    i = 2, j = "CR Klebsiella pneumoniae n (%) [95% CI]",
    value = as_paragraph("CR ", as_i("Klebsiella pneumoniae"), "\nn (%) [95% CI]"),
    part = "header"
  ) %>%
  add_footer_lines("CR = Carbapenem-Resistant, pCR = Presumptive Carbapenem-Resistant; CI = Confidence Interval (calculated via Wilson score method).") %>%
  autofit() %>%
  theme_booktabs() %>%
  align(
    j = c("Total test", "pCR Escherichia coli n (%) [95% CI]", 
          "CR Escherichia coli n (%) [95% CI]", "CR Klebsiella pneumoniae n (%) [95% CI]"),
    align = "center", 
    part = "all"
  ) %>%
  bold(i = bold_rows, j = 1, bold = TRUE, part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  fontsize(size = 8, part = "footer") %>%
  font(fontname = "Times New Roman", part = "all")

# Save Word Document
doc <- read_docx() %>%
  body_add_par("Table: Carbapenem-Resistant (CR) Enterobacterales Isolation and Species Confirmation Stratified by Socio-Demographic and Geographic Characteristics", style = "heading 2") %>%
  body_add_flextable(ft)

output_path <- file.path("Results", "CR_E_coli_Isolation_Rates_Table.docx")
print(doc, target = output_path)

message("Table successfully updated with footer and revised species headers. Output saved to: ", output_path)
################################################################################
# Subsetting unique participant isolate records
cr_isolates_subset <- joined_data %>%
  filter(
    str_detect(coalesce(`COLONY MORPHOLOGY ON CARBA`, ""), "(?i)Red|Pinkish|Metalic blue")
  ) %>%
  distinct(INIKA_ID, Isolate, `COLONY MORPHOLOGY ON CARBA`, .keep_all = TRUE) %>%
  select(INIKA_ID,SEASON,REGION,DISTRICT, `COLONY MORPHOLOGY ON CARBA`,
         PROTOCOL,Isolate,VITEK_MS_Results, Isolate_NVI) # The isolate tested at TVLA was not turned to E.coli ~ Citrobacter
################################################################################
# ESBL Counting 2

# 1. Create 'Results' folder if it does not exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Helper function: formats n (%) to 1 decimal place, stripping '.0' for whole numbers
fmt_n_pct <- function(n, total) {
  if (is.na(total) || total == 0) return("0 (0%)")
  pct <- (n / total) * 100
  pct_str <- sub("\\.0$", "", sprintf("%.1f", pct))
  sprintf("%d (%s%%)", n, pct_str)
}

# 2. Rename & Convert Antibiotics Data
cleaned_data <- joined_data %>%
  rename(
    Amoxicillin = AMX_ED10,
    Azithromycin = AZM_ED15,
    Ceftriaxone = CRO_ED30,
    Ciprofloxacin = CIP_ED5,
    Doxycycline = DOX_ED30,
    Florfenicol = FLR_ED30,
    Gentamicin = GEN_ED10,
    Meropenem = MEM_ED10,
    Oxytetracycline = OXY_ED30,
    `Polymyxin_B(PB)` = POL_ED300,
    `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
    Cefotaxime = CTX_ED5,
    `Cefotaxime/Clavulanic Acid` = CTC_ED30
  ) %>%
  mutate(across(all_of(abx_cols), ~ as.numeric(trimws(.))))

# 3. Filter Confirmed 3GCR Isolates
c3gcr_data <- cleaned_data %>%
  filter(
    Ceftriaxone < 23,
    (
      (grepl("E|K", Isolate, ignore.case = TRUE) & grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE)) |
        (grepl("E|K", Isolate, ignore.case = TRUE) & trimws(VITEK_MS_Results) == "Klebsiella pneumoniae")
    )
  ) %>%
  mutate(
    C3GCR_Species = case_when(
      grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) ~ "Confirmed 3GCR Escherichia coli",
      trimws(VITEK_MS_Results) == "Klebsiella pneumoniae" ~ "Confirmed 3GCR Klebsiella pneumoniae",
      TRUE ~ NA_character_
    )
  )

# Calculate total Confirmed 3GCR per species (deduplicated by INIKA_ID)
total_c3gcr_ecoli <- c3gcr_data %>%
  filter(C3GCR_Species == "Confirmed 3GCR Escherichia coli") %>%
  distinct(INIKA_ID) %>%
  nrow()

total_c3gcr_kpneumo <- c3gcr_data %>%
  filter(C3GCR_Species == "Confirmed 3GCR Klebsiella pneumoniae") %>%
  distinct(INIKA_ID) %>%
  nrow()

# 4. Filter ESBL Positives from Confirmed 3GCR Data
esbl_data <- c3gcr_data %>%
  filter(`Cefotaxime/Clavulanic Acid` - Cefotaxime >= 5) %>%
  filter(
    (
      (grepl("E", Isolate, ignore.case = TRUE) & grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)) |
        (grepl("K", Isolate, ignore.case = TRUE) & grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) & grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE))
    ) |
      (trimws(VITEK_MS_Results) == "Klebsiella pneumoniae")
  ) %>%
  mutate(
    Species_Group = case_when(
      grepl("Escherichia coli", VITEK_MS_Results, ignore.case = TRUE) ~ "E. coli",
      trimws(VITEK_MS_Results) == "Klebsiella pneumoniae" ~ "K. pneumoniae",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(INIKA_ID, Isolate) %>%
  distinct(INIKA_ID, Species_Group, .keep_all = TRUE)

# Logging summary messages
n_esbl_ecoli <- sum(esbl_data$Species_Group == "E. coli", na.rm = TRUE)
n_esbl_kpneumo <- sum(esbl_data$Species_Group == "K. pneumoniae", na.rm = TRUE)

message(sprintf("Confirmed 3GCR E. coli total: %d", total_c3gcr_ecoli))
message(sprintf("Confirmed 3GCR K. pneumoniae total: %d", total_c3gcr_kpneumo))
message(sprintf("ESBL E. coli isolates retained: %d", n_esbl_ecoli))
message(sprintf("ESBL K. pneumoniae isolates retained: %d", n_esbl_kpneumo))

# 5. Summarizer Functions
summarize_variable <- function(c3g_df, esbl_df, var_col, var_label) {
  var_sym <- sym(var_col)
  
  c3g_counts <- c3g_df %>%
    filter(!is.na(!!var_sym)) %>%
    group_by(Values = as.character(!!var_sym)) %>%
    summarise(
      c3g_ecoli = n_distinct(INIKA_ID[C3GCR_Species == "Confirmed 3GCR Escherichia coli"]),
      c3g_kpneumo = n_distinct(INIKA_ID[C3GCR_Species == "Confirmed 3GCR Klebsiella pneumoniae"]),
      .groups = "drop"
    )
  
  esbl_counts <- esbl_df %>%
    filter(!is.na(!!var_sym)) %>%
    group_by(Values = as.character(!!var_sym)) %>%
    summarise(
      n_ecoli = n_distinct(INIKA_ID[Species_Group == "E. coli"]),
      n_kpneumo = n_distinct(INIKA_ID[Species_Group == "K. pneumoniae"]),
      .groups = "drop"
    )
  
  full_join(c3g_counts, esbl_counts, by = "Values") %>%
    mutate(
      across(c(c3g_ecoli, c3g_kpneumo, n_ecoli, n_kpneumo), ~ coalesce(., 0L)),
      `Total confirmed 3GCR Escherichia coli` = as.character(c3g_ecoli),
      `ESBL E. coli` = map2_chr(n_ecoli, c3g_ecoli, fmt_n_pct),
      `Confirmed 3GCR Klebsiella pneumoniae` = as.character(c3g_kpneumo),
      `ESBL K. pneumoniae` = map2_chr(n_kpneumo, c3g_kpneumo, fmt_n_pct),
      Variable = var_label
    ) %>%
    select(
      Variable, Values, 
      `Total confirmed 3GCR Escherichia coli`, `ESBL E. coli`, 
      `Confirmed 3GCR Klebsiella pneumoniae`, `ESBL K. pneumoniae`
    )
}

summarize_districts <- function(c3g_df, esbl_df, district_col, region_col) {
  dist_sym <- sym(district_col)
  reg_sym <- sym(region_col)
  
  c3g_counts <- c3g_df %>%
    filter(!is.na(!!dist_sym), !is.na(!!reg_sym)) %>%
    group_by(Region = as.character(!!reg_sym), District = as.character(!!dist_sym)) %>%
    summarise(
      c3g_ecoli = n_distinct(INIKA_ID[C3GCR_Species == "Confirmed 3GCR Escherichia coli"]),
      c3g_kpneumo = n_distinct(INIKA_ID[C3GCR_Species == "Confirmed 3GCR Klebsiella pneumoniae"]),
      .groups = "drop"
    )
  
  esbl_counts <- esbl_df %>%
    filter(!is.na(!!dist_sym), !is.na(!!reg_sym)) %>%
    group_by(Region = as.character(!!reg_sym), District = as.character(!!dist_sym)) %>%
    summarise(
      n_ecoli = n_distinct(INIKA_ID[Species_Group == "E. coli"]),
      n_kpneumo = n_distinct(INIKA_ID[Species_Group == "K. pneumoniae"]),
      .groups = "drop"
    )
  
  full_join(c3g_counts, esbl_counts, by = c("Region", "District")) %>%
    mutate(
      across(c(c3g_ecoli, c3g_kpneumo, n_ecoli, n_kpneumo), ~ coalesce(., 0L)),
      Region = factor(Region, levels = c("Kilimanjaro", "Mwanza"))
    ) %>%
    arrange(Region, District) %>%
    mutate(
      `Total confirmed 3GCR Escherichia coli` = as.character(c3g_ecoli),
      `ESBL E. coli` = map2_chr(n_ecoli, c3g_ecoli, fmt_n_pct),
      `Confirmed 3GCR Klebsiella pneumoniae` = as.character(c3g_kpneumo),
      `ESBL K. pneumoniae` = map2_chr(n_kpneumo, c3g_kpneumo, fmt_n_pct),
      Variable = "District",
      Values = District
    ) %>%
    select(
      Variable, Values, 
      `Total confirmed 3GCR Escherichia coli`, `ESBL E. coli`, 
      `Confirmed 3GCR Klebsiella pneumoniae`, `ESBL K. pneumoniae`
    )
}

# 6. Aggregate Table Data
table_df <- bind_rows(
  tibble(
    Variable = "Overall",
    Values = "Confirmed 3GCR",
    `Total confirmed 3GCR Escherichia coli` = as.character(total_c3gcr_ecoli),
    `ESBL E. coli` = fmt_n_pct(n_esbl_ecoli, total_c3gcr_ecoli),
    `Confirmed 3GCR Klebsiella pneumoniae` = as.character(total_c3gcr_kpneumo),
    `ESBL K. pneumoniae` = fmt_n_pct(n_esbl_kpneumo, total_c3gcr_kpneumo)
  ),
  summarize_variable(c3gcr_data, esbl_data, "ORIGIN_OF_SAMPLE", "Origin"),
  summarize_variable(c3gcr_data, esbl_data, "SEASON.x", "Season"),
  summarize_variable(c3gcr_data, esbl_data, "REGION.x", "Region"),
  summarize_districts(c3gcr_data, esbl_data, district_col = "DISTRICT.x", region_col = "REGION.x")
)

# 7. Construct Publication-Ready Flextable
border_style <- fp_border(color = "black", width = 1)
title_text <- "Table 1: Distribution of ESBL E. coli and ESBL K. pneumoniae Isolates by Sample Origin, Season, Region, and District"

esbl_table <- table_df %>%
  flextable() %>%
  # --- UPDATE: Set caption outside table grid ---
  set_caption(caption = title_text) %>%
  # ---------------------------------------------
set_header_labels(
  Variable = "Variable",
  Values = "Category / Value",
  `Total confirmed 3GCR Escherichia coli` = "Total confirmed 3GCR Escherichia coli",
  `ESBL E. coli` = sprintf("ESBL E. coli\n(N = %d), n (%%)", n_esbl_ecoli),
  `Confirmed 3GCR Klebsiella pneumoniae` = "Confirmed 3GCR Klebsiella pneumoniae",
  `ESBL K. pneumoniae` = sprintf("ESBL K. pneumoniae\n(N = %d), n (%%)", n_esbl_kpneumo)
) %>%
  merge_v(j = "Variable", part = "body") %>%
  valign(j = "Variable", valign = "top", part = "body") %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = ~ Values == "Confirmed 3GCR", part = "body") %>%
  align(j = c("Total confirmed 3GCR Escherichia coli", "ESBL E. coli", "Confirmed 3GCR Klebsiella pneumoniae", "ESBL K. pneumoniae"), align = "center", part = "all") %>%
  align(j = c("Variable", "Values"), align = "left", part = "all") %>%
  italic(j = c("Total confirmed 3GCR Escherichia coli", "ESBL E. coli", 
               "Confirmed 3GCR Klebsiella pneumoniae", "ESBL K. pneumoniae"), 
         part = "header") %>%
  add_footer_lines(
    sprintf(
      "Note: Percentages calculated out of total Confirmed 3GCR Escherichia coli (N = %d) and Confirmed 3GCR Klebsiella pneumoniae (N = %d).",
      total_c3gcr_ecoli, total_c3gcr_kpneumo
    )
  ) %>%
  italic(part = "footer") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()
# 8. Export to Word
doc <- read_docx() %>%
  body_add_flextable(esbl_table)

print(doc, target = "Results/ESBL_Ecoli_Kpneumoniae_Summary_Table.docx")
######################################################################
## Computing the AMR for the confirmed CR E.coli
# 1. Define Antibiotic Class Mapping
abx_classes <- tribble(
  ~Antimicrobial_substance,          ~Class,
  "Amoxicillin",                     "Penicillins",
  "Azithromycin",                    "Macrolides",
  "Ceftriaxone",                     "Cephalosporins (3rd Gen)",
  "Cefotaxime",                      "Cephalosporins (3rd Gen)",
  "Cefotaxime/Clavulanic Acid",      "β-lactam/Inhibitor",
  "Ciprofloxacin",                   "Fluoroquinolones",
  "Doxycycline",                     "Tetracyclines",
  "Oxytetracycline",                 "Tetracyclines",
  "Gentamicin",                      "Aminoglycosides",
  "Meropenem",                       "Carbapenems",
  "Florfenicol",                     "Phenicols",
  "Polymyxin_B(PB)",                 "Polymyxins",
  "Sulfamethoxazole/Trimethoprim",   "Sulfonamides"
)

abx_cols <- abx_classes$Antimicrobial_substance

# 2. Main Processing Function (Count ONLY for CR E. coli)
process_cr_ecoli_amr_results <- function(joined_data, Ecoli_EPI_CUTOFF) {
  
  # Step 1: Clean, confirm species, and filter strictly for CR E. coli
  Data_Clean <- joined_data %>%
    select(
      INIKA_ID, ORIGIN_OF_SAMPLE, SEASON.x, REGION.x, DISTRICT.x,
      `COLONY MORPHOLOGY ON CARBA`, VITEK_MS_Results, Isolate, Isolate_NVI,
      AMX_ED10, AZM_ED15, CRO_ED30, CIP_ED5, DOX_ED30, FLR_ED30,
      GEN_ED10, MEM_ED10, OXY_ED30, POL_ED300, SXT_ED1_2,
      CTX_ED5, CTC_ED30
    ) %>%
    rename(
      'SEASON' = SEASON.x, 
      'REGION' = REGION.x,
      'DISTRICT' = DISTRICT.x,
      Amoxicillin = AMX_ED10,
      Azithromycin = AZM_ED15,
      Ceftriaxone = CRO_ED30,
      Ciprofloxacin = CIP_ED5,
      Doxycycline = DOX_ED30,
      Florfenicol = FLR_ED30,
      Gentamicin = GEN_ED10,
      Meropenem = MEM_ED10,
      Oxytetracycline = OXY_ED30,
      `Polymyxin_B(PB)` = POL_ED300,
      `Sulfamethoxazole/Trimethoprim` = SXT_ED1_2,
      Cefotaxime = CTX_ED5,
      `Cefotaxime/Clavulanic Acid` = CTC_ED30
    ) %>%
    # Multi-level species confirmation logic
    mutate(
      # Check exclusively if Isolate_NVI indicates E. coli
      is_nvi_ecoli = str_detect(coalesce(Isolate_NVI, ""), "(?i)E\\.?\\s*coli|Escherichia\\s+coli|^E$"),
      
      CONFIRMED_SPECIES = case_when(
        is_nvi_ecoli ~ "CR Escherichia coli",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(across(all_of(abx_cols), ~ as.numeric(trimws(.)))) %>%

    filter(
      CONFIRMED_SPECIES == "CR Escherichia coli",
      str_detect(coalesce(`COLONY MORPHOLOGY ON CARBA`, ""), "(?i)Red|Pinkish|Pink|Metalic blue")
    ) %>%
    # Deduplicate strictly by INIKA_ID
    group_by(INIKA_ID) %>%
    slice(1) %>%
    ungroup()
  
  message(sprintf("Total confirmed CR E. coli unique isolates processed: %d", nrow(Data_Clean)))
  
  # Step 2: Pivot long and merge with ECOFF cutoffs
  Ecoli_Long <- Data_Clean %>%
    pivot_longer(
      cols = any_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    ) %>%
    left_join(Ecoli_EPI_CUTOFF, by = "Antimicrobial_substance") %>%
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
  
  # Step 3: Calculation helper logic
  calculate_summary <- function(df, group_var) {
    group_cols <- if (group_var == "Overall") "Antimicrobial_substance" else c(group_var, "Antimicrobial_substance")
    df %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        Total = n(),
        Resistant = sum(Categorical_Result == "R"),
        .groups = "drop"
      ) %>%
      mutate(
        Grouping_Variable = group_var,
        Grouping_Value = if (group_var == "Overall") "Overall" else as.character(.data[[group_var]])
      )
  }
  
  group_vars <- c("Overall", "ORIGIN_OF_SAMPLE", "REGION", "SEASON")
  Combined_Summaries <- map_dfr(group_vars, ~ calculate_summary(Ecoli_Long, .x))
  
  # Step 4: Formatting final output table
  Final_Table <- Combined_Summaries %>%
    rowwise() %>%
    mutate(
      Percentage = round((Resistant / Total) * 100, 1),
      CI_text = if (Grouping_Variable == "Overall") {
        test <- binom.test(Resistant, Total)
        paste0(" (", round(test$conf.int[1] * 100, 1), "–", round(test$conf.int[2] * 100, 1), ")")
      } else { "" },
      Display_Value = paste0(Percentage, "% (", Resistant, "/", Total, ")", CI_text)
    ) %>%
    ungroup() %>%
    select(Antimicrobial_substance, Grouping_Value, Display_Value) %>%
    pivot_wider(names_from = Grouping_Value, values_from = Display_Value) %>%
    left_join(abx_classes, by = "Antimicrobial_substance") %>%
    select(Class, Antimicrobial_substance, Overall, everything()) %>%
    arrange(Class, Antimicrobial_substance)
  
  return(Final_Table)
}

# 3. Run and Export Results
final_CR_ECO_AMR_table <- process_cr_ecoli_amr_results(joined_data, Ecoli_EPI_CUTOFF)

if (!dir.exists("Results")) dir.create("Results")

# Export to Excel
write_xlsx(final_CR_ECO_AMR_table, "Results/CR_Ecoli_AMR_Resistance_Table.xlsx")

# Export to Word Document
ft_compact <- flextable(final_CR_ECO_AMR_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding = 1, part = "all") %>%
  merge_v(j = ~ Class) %>%
  bold(part = "header") %>%
  set_table_properties(layout = "fixed") %>%
  width(width = 1.1) %>%
  width(j = 1:2, width = 1.5)

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5),
  type = "continuous"
)

save_as_docx(
  ft_compact, 
  path = "Results/CR_Ecoli_AMR_Results_Final_Fit.docx",
  pr_section = sect_properties
)
