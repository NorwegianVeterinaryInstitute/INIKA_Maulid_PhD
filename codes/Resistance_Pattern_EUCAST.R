

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

# Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
#   "data/ECOFF_E.coli_K.pneumoniae.xlsx",
#   sheet = 2
# )
# 
# Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT
# 
# 
# ECOFF_EUCAST_BREAK_POINT <- read_excel(
#   "data/ECOFF_E.coli_K.pneumoniae.xlsx"
# )
# 
#EPI_CUTOFF <- ECOFF_EUCAST_BREAK_POINT
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
group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE")

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


#Used later for: MDR calculation and Table grouping
#Step 3: Define ESCR E. coli
# Define the Main Processing Function 
process_ecoli_amr_results <- function(joined_data, Ecoli_EPI_CUTOFF) {
  
  # Cleaning and Filtering 
  
  abx_cols <- abx_classes$Antimicrobial_substance
  
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
  group_vars <- c("Overall", "REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE")
  
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
# Distribution of Zone Diameters and Resistance Frequency in ESCR E. coli
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