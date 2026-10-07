

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

########################

######################
if (!dir.exists("Results")) {
  dir.create("Results")
}

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


####################################
#Global functions
##############################
rename_antibiotics <- function(df) {
  
  df %>%
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
    )
}

convert_antibiotics <- function(df) {
  
  df %>%
    mutate(
      across(
        all_of(abx_cols),
        ~ as.numeric(trimws(.))
      )
    )
}
get_escr_ecoli <- function(df) {
  
  df %>%
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
    ) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
}

get_escr_kpn <- function(df) {
  
  df %>%
    filter(
      (
        (
          grepl("K", Isolate, ignore.case = TRUE) &
            grepl("Klebsiella pneumoniae", VITEK_MS_Results, ignore.case = TRUE) &
            grepl("Metallic blue", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
        ) |
          (
            grepl("E", Isolate, ignore.case = TRUE) &
              grepl("Klebsiella pneumoniae", VITEK_MS_Results, ignore.case = TRUE) &
              grepl("Pink|Red", `COLONY MORPHOLOGY ON C3GR`, ignore.case = TRUE)
          )
      ) &
        Ceftriaxone < 23
    ) %>%
    distinct(INIKA_ID, .keep_all = TRUE)
}

# Format helper: Outputs 1 decimal place, or integer if ending in .0
fmt_num <- function(val) {
  val_rounded <- round(val, 1)
  if (val_rounded %% 1 == 0) {
    sprintf("%.0f", val_rounded)
  } else {
    sprintf("%.1f", val_rounded)
  }
}



# Helper function: formats n (%) to 1 decimal place, stripping '.0' for whole numbers
fmt_n_pct <- function(n, total) {
  if (is.na(total) || total == 0) return("0 (0%)")
  pct <- (n / total) * 100
  pct_str <- sub("\\.0$", "", sprintf("%.1f", pct))
  sprintf("%d (%s%%)", n, pct_str)
}


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
# Supplementary Table 5
# AMR Frequencies for Confirmed ESCR E. coli
##############################################################################

# Antibiotic class mapping
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

##############################################################################
# Main Processing Function
##############################################################################

process_ecoli_amr_results <- function(joined_data, Ecoli_EPI_CUTOFF) {
  
  ###########################################################################
  # Cleaning and Filtering
  ###########################################################################
  
  Data_Clean <- joined_data %>%
    select(
      INIKA_ID,
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE,
      `COLONY MORPHOLOGY ON C3GR`,
      Isolate,
      VITEK_MS_Results,
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
    get_escr_ecoli()
  
  message(
    paste(
      "ESCR E. coli isolates retained:",
      nrow(Data_Clean)
    )
  )
  
  ###########################################################################
  # Pivot and Classify
  ###########################################################################
  
  Ecoli_Long <- Data_Clean %>%
    pivot_longer(
      cols = all_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "value"
    ) %>%
    left_join(
      Ecoli_EPI_CUTOFF,
      by = "Antimicrobial_substance"
    ) %>%
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
  
  ###########################################################################
  # Summary Function
  ###########################################################################
  
  calculate_summary <- function(df, group_var) {
    
    group_cols <- if (group_var == "Overall") {
      "Antimicrobial_substance"
    } else {
      c(group_var, "Antimicrobial_substance")
    }
    
    df %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        Total = n(),
        Resistant = sum(Categorical_Result == "R"),
        .groups = "drop"
      ) %>%
      mutate(
        Grouping_Variable = group_var,
        Grouping_Value = if (group_var == "Overall") {
          "Overall"
        } else {
          as.character(.data[[group_var]])
        }
      )
  }
  
  ###########################################################################
  # Generate Summaries
  ###########################################################################
  
  group_vars <- c(
    "Overall",
    "ORIGIN_OF_SAMPLE",
    "REGION.x",
    "SEASON.x"
  )
  
  Combined_Summaries <- map_dfr(
    group_vars,
    ~ calculate_summary(Ecoli_Long, .x)
  )
  
  ###########################################################################
  # Format Final Table
  ###########################################################################
  
  Final_Table <- Combined_Summaries %>%
    rowwise() %>%
    mutate(
      Percentage = round((Resistant / Total) * 100, 1),
      CI_text = if (Grouping_Variable == "Overall") {
        
        test <- binom.test(Resistant, Total)
        
        paste0(
          " (",
          round(test$conf.int[1] * 100, 1),
          "–",
          round(test$conf.int[2] * 100, 1),
          ")"
        )
        
      } else {
        ""
      },
      Display_Value = paste0(
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
# Run Analysis
##############################################################################

final_ESCR_ECO_AMR_table <- process_ecoli_amr_results(
  joined_data,
  Ecoli_EPI_CUTOFF
)

print(final_ESCR_ECO_AMR_table)

##############################################################################
# Export Excel
##############################################################################

write_xlsx(
  final_ESCR_ECO_AMR_table,
  "Results/Ecoli_ESCR_Resistance_Table.xlsx"
)

##############################################################################
# Export Word
##############################################################################

ft <- flextable(final_ESCR_ECO_AMR_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(size = 9, part = "all") %>%
  merge_v(j = ~ Class) %>%
  bold(part = "header")

sect_properties <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(
    bottom = 0.5,
    top = 0.5,
    right = 0.5,
    left = 0.5
  ),
  type = "continuous"
)

ft_compact <- ft %>%
  fontsize(size = 8, part = "all") %>%
  padding(padding = 1, part = "all") %>%
  set_table_properties(layout = "fixed") %>%
  width(width = 1.1) %>%
  width(j = 1:2, width = 1.5)

save_as_docx(
  ft_compact,
  path = "Results/Ecoli_AMR_Results_Final_Fit.docx",
  pr_section = sect_properties
)

##############################################################################
# Supplementary Table 6
# AMR Frequencies for Confirmed ESCR Klebsiella pneumoniae
##############################################################################

process_kpn_amr_analysis <- function(joined_data, Kpn_EPI_CUTOFF) {
  
  ###########################################################################
  # Data Cleaning & Filtering
  ###########################################################################
  
  colnames(joined_data) <- trimws(colnames(joined_data))
  
  cleaned_data <- joined_data %>%
    select(
      INIKA_ID,
      REGION.x,
      SEASON.x,
      ORIGIN_OF_SAMPLE,
      `COLONY MORPHOLOGY ON C3GR`,
      Isolate,
      VITEK_MS_Results,
      ESCR_KPN_presumptive,
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
    get_escr_kpn()
  
  message(
    paste(
      "ESCR K. pneumoniae isolates retained:",
      nrow(cleaned_data)
    )
  )
  
  ###########################################################################
  # Pivot & Susceptibility Classification
  ###########################################################################
  
  kpn_long <- cleaned_data %>%
    pivot_longer(
      cols = all_of(abx_classes$Antimicrobial_substance),
      names_to = "Antimicrobial_substance",
      values_to = "Measured_Zone"
    ) %>%
    left_join(
      Kpn_EPI_CUTOFF,
      by = "Antimicrobial_substance"
    ) %>%
    mutate(
      across(
        c(S, R),
        ~ as.numeric(str_remove_all(.x, "[^0-9.]")),
        .names = "{.col}_num"
      ),
      Measured_Zone = as.numeric(Measured_Zone),
      Result = case_when(
        Measured_Zone >= S_num ~ "S",
        Measured_Zone <= R_num ~ "R",
        Measured_Zone > R_num & Measured_Zone < S_num ~ "I",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Result))
  
  ###########################################################################
  # Summary Helper Function
  ###########################################################################
  
  get_summary <- function(df, var_name) {
    
    group_cols <-
      if (var_name == "Overall") {
        "Antimicrobial_substance"
      } else {
        c(var_name, "Antimicrobial_substance")
      }
    
    df %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        n_res = sum(Result == "R"),
        total = n(),
        .groups = "drop"
      ) %>%
      mutate(
        Grouping_Var = var_name,
        Grouping_Val = if (
          var_name == "Overall"
        ) {
          "Overall"
        } else {
          as.character(.data[[var_name]])
        }
      )
  }
  
  ###########################################################################
  # Generate Table
  ###########################################################################
  
  group_vars <- c(
    "Overall",
    "ORIGIN_OF_SAMPLE",
    "REGION.x",
    "SEASON.x"
  )
  
  final_tablekpn <- map_dfr(
    group_vars,
    ~ get_summary(kpn_long, .x)
  ) %>%
    rowwise() %>%
    mutate(
      perc = round(
        (n_res / total) * 100,
        1
      ),
      ci = if (Grouping_Var == "Overall") {
        
        bt <- binom.test(
          n_res,
          total
        )
        
        paste0(
          " (",
          round(bt$conf.int[1] * 100, 1),
          "–",
          round(bt$conf.int[2] * 100, 1),
          ")"
        )
        
      } else {
        ""
      },
      cell_text = paste0(
        perc,
        "% (",
        n_res,
        "/",
        total,
        ")",
        ci
      )
    ) %>%
    ungroup() %>%
    select(
      Antimicrobial_substance,
      Grouping_Val,
      cell_text
    ) %>%
    pivot_wider(
      names_from = Grouping_Val,
      values_from = cell_text
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
  
  ###########################################################################
  # Export Word
  ###########################################################################
  
  ft <- flextable(final_tablekpn) %>%
    theme_booktabs() %>%
    fontsize(size = 8, part = "all") %>%
    padding(padding = 1, part = "all") %>%
    merge_v(j = ~ Class) %>%
    bold(part = "header") %>%
    set_table_properties(layout = "autofit")
  
  sect_props <- prop_section(
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
    ft,
    path = "Results/Kpn_AMR_Table_Landscape.docx",
    pr_section = sect_props
  )
  
  ###########################################################################
  # Export Excel
  ###########################################################################
  
  write_xlsx(
    final_tablekpn,
    "Results/Kpn_AMR_Table.xlsx"
  )
  
  message(
    "Success: Table saved to Results folder as Word and Excel."
  )
  
  final_tablekpn
}

##############################################################################
# Execute Analysis
##############################################################################

kpn_amr_final_results <- process_kpn_amr_analysis(
  joined_data,
  Kpn_EPI_CUTOFF
)


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
    rename_antibiotics() %>%
    convert_antibiotics() %>%
    get_escr_ecoli()
  
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
  # MDR Calculation
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
  
  MDR_Data
}

##############################################################################
# Run Analysis
##############################################################################

ecoli_mdr_results <- process_ecoli_MDR_results(
  joined_data,
  Ecoli_EPI_CUTOFF,
  abx_classes
)

##############################################################################
# Summary Table
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
    `MDR n (%)` = sprintf(
      "%d (%.1f%%)",
      MDR_n,
      Percent
    ),
    `95% CI` = sprintf(
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
# Export
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

save_as_docx(
  doc_table,
  path = "Results/Ecoli_MDR_Final_Table.docx"
)
  
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
    rename_antibiotics() %>%
    convert_antibiotics() %>%
    get_escr_kpn()
  
  message(
    paste(
      "ESCR K. pneumoniae isolates for MDR analysis:",
      nrow(ESCR_confirmed_data)
    )
  )
  
  ###########################################################################
  # SIR Classification
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
  # MDR Calculation
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
  
  MDR_Data
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
    `MDR n (%)` = sprintf(
      "%d (%.1f%%)",
      MDR_n,
      Percent
    ),
    `95% CI` = sprintf(
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

save_as_docx(
  doc_table,
  path = "Results/Kpn_MDR_Final_Table.docx"
)

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
  rename_antibiotics() %>%
  convert_antibiotics() %>%
  get_escr_ecoli()

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


##############################################################################
# Supplementary Table 4
# Distribution of Zone Diameters and Resistance Frequency in ESCR K. pneumoniae
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
  rename_antibiotics() %>%
  convert_antibiotics() %>%
  get_escr_kpn()

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
##############################################################################
# Prepare plotting dataset
##############################################################################

Plot_Resistance <- ESCR_long %>%
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
  Plot_Resistance,
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

##############################################################################
# Computing ESBL Counts
##############################################################################

##############################################################################
# General summarizer for Origin, Season, and Region
##############################################################################

summarize_variable <- function(df, var_col, var_label) {
  
  var_sym <- sym(var_col)
  
  df %>%
    filter(!is.na(!!var_sym)) %>%
    group_by(Values = as.character(!!var_sym)) %>%
    summarise(
      n_ecoli = sum(Species_Group == "E. coli", na.rm = TRUE),
      n_kpneumo = sum(Species_Group == "K. pneumoniae", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      `ESBL E. coli` = map2_chr(
        n_ecoli,
        total_ecoli,
        fmt_n_pct
      ),
      `ESBL K. pneumoniae` = map2_chr(
        n_kpneumo,
        total_kpneumo,
        fmt_n_pct
      ),
      Variable = var_label
    ) %>%
    select(
      Variable,
      Values,
      `ESBL E. coli`,
      `ESBL K. pneumoniae`
    )
}

##############################################################################
# District summarizer
##############################################################################

summarize_districts <- function(
    df,
    district_col,
    region_col
) {
  
  dist_sym <- sym(district_col)
  reg_sym <- sym(region_col)
  
  df %>%
    filter(
      !is.na(!!dist_sym),
      !is.na(!!reg_sym)
    ) %>%
    group_by(
      Region = as.character(!!reg_sym),
      District = as.character(!!dist_sym)
    ) %>%
    summarise(
      n_ecoli = sum(
        Species_Group == "E. coli",
        na.rm = TRUE
      ),
      n_kpneumo = sum(
        Species_Group == "K. pneumoniae",
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Region = factor(
        Region,
        levels = c(
          "Kilimanjaro",
          "Mwanza"
        )
      )
    ) %>%
    arrange(
      Region,
      District
    ) %>%
    mutate(
      `ESBL E. coli` = map2_chr(
        n_ecoli,
        total_ecoli,
        fmt_n_pct
      ),
      `ESBL K. pneumoniae` = map2_chr(
        n_kpneumo,
        total_kpneumo,
        fmt_n_pct
      ),
      Variable = "District",
      Values = District
    ) %>%
    select(
      Variable,
      Values,
      `ESBL E. coli`,
      `ESBL K. pneumoniae`
    )
}

##############################################################################
# Create ESBL dataset
##############################################################################

esbl_data <- joined_data %>%
  rename_antibiotics() %>%
  convert_antibiotics() %>%
  filter(
    `Cefotaxime/Clavulanic Acid` - Cefotaxime >= 5,
    Ceftriaxone < 23
  ) %>%
  filter(
    (
      (
        grepl(
          "E",
          Isolate,
          ignore.case = TRUE
        ) &
          grepl(
            "Escherichia coli",
            VITEK_MS_Results,
            ignore.case = TRUE
          ) &
          grepl(
            "Pink|Red",
            `COLONY MORPHOLOGY ON C3GR`,
            ignore.case = TRUE
          )
      ) |
        (
          grepl(
            "K",
            Isolate,
            ignore.case = TRUE
          ) &
            grepl(
              "Escherichia coli",
              VITEK_MS_Results,
              ignore.case = TRUE
            ) &
            grepl(
              "Metallic blue",
              `COLONY MORPHOLOGY ON C3GR`,
              ignore.case = TRUE
            )
        )
    ) |
      (
        trimws(VITEK_MS_Results) ==
          "Klebsiella pneumoniae"
      )
  ) %>%
  mutate(
    Species_Group = case_when(
      grepl(
        "Escherichia coli",
        VITEK_MS_Results,
        ignore.case = TRUE
      ) ~ "E. coli",
      trimws(VITEK_MS_Results) ==
        "Klebsiella pneumoniae" ~
        "K. pneumoniae",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    INIKA_ID,
    Isolate
  ) %>%
  distinct(
    INIKA_ID,
    .keep_all = TRUE
  )

##############################################################################
# Logging
##############################################################################

n_ecoli <- sum(
  esbl_data$Species_Group == "E. coli",
  na.rm = TRUE
)

n_kpneumo <- sum(
  esbl_data$Species_Group == "K. pneumoniae",
  na.rm = TRUE
)

n_total <- nrow(esbl_data)

message(
  sprintf(
    "ESBL E. coli isolates retained: %d",
    n_ecoli
  )
)

message(
  sprintf(
    "ESBL K. pneumoniae isolates retained: %d",
    n_kpneumo
  )
)

message(
  sprintf(
    "Total unique isolates processed: %d",
    n_total
  )
)

##############################################################################
# Totals
##############################################################################

total_ecoli <- sum(
  esbl_data$Species_Group == "E. coli"
)

total_kpneumo <- sum(
  esbl_data$Species_Group == "K. pneumoniae"
)

##############################################################################
# Aggregate Table
##############################################################################

table_df <- bind_rows(
  
  tibble(
    Variable = "Overall",
    Values = "Total Unique INIKA_ID",
    `ESBL E. coli` =
      fmt_n_pct(
        total_ecoli,
        total_ecoli
      ),
    `ESBL K. pneumoniae` =
      fmt_n_pct(
        total_kpneumo,
        total_kpneumo
      )
  ),
  
  summarize_variable(
    esbl_data,
    "ORIGIN_OF_SAMPLE",
    "Origin"
  ),
  
  summarize_variable(
    esbl_data,
    "SEASON.x",
    "Season"
  ),
  
  summarize_variable(
    esbl_data,
    "REGION.x",
    "Region"
  ),
  
  summarize_districts(
    esbl_data,
    district_col = "DISTRICT.x",
    region_col = "REGION.x"
  )
)

##############################################################################
# Publication-ready table
##############################################################################

esbl_table <- table_df %>%
  flextable() %>%
  set_header_labels(
    Variable = "Variable",
    Values = "Category / Value",
    `ESBL E. coli` =
      sprintf(
        "ESBL E. coli\n(N = %d), n (%%)",
        total_ecoli
      ),
    `ESBL K. pneumoniae` =
      sprintf(
        "ESBL K. pneumoniae\n(N = %d), n (%%)",
        total_kpneumo
      )
  ) %>%
  merge_v(j = "Variable") %>%
  valign(
    j = "Variable",
    valign = "top"
  ) %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(
    i = ~ Values == "Total Unique INIKA_ID"
  ) %>%
  align(
    j = c(
      "ESBL E. coli",
      "ESBL K. pneumoniae"
    ),
    align = "center",
    part = "all"
  ) %>%
  align(
    j = c(
      "Variable",
      "Values"
    ),
    align = "left",
    part = "all"
  ) %>%
  add_header_lines(
    values =
      "Table 1: Distribution of ESBL E. coli and ESBL K. pneumoniae Isolates by Sample Origin, Season, Region, and District"
  ) %>%
  bold(
    i = 1,
    part = "header"
  ) %>%
  add_footer_lines(
    sprintf(
      paste0(
        "Note: Filtered using ",
        "Cefotaxime/Clavulanic Acid − Cefotaxime ≥ 5 mm, ",
        "Ceftriaxone < 23 mm, species confirmation, and ",
        "deduplicated by unique INIKA_ID. ",
        "Percentages calculated out of total ESBL E. coli (N = %d) ",
        "and ESBL K. pneumoniae (N = %d)."
      ),
      total_ecoli,
      total_kpneumo
    )
  ) %>%
  italic(
    j = c(
      "ESBL E. coli",
      "ESBL K. pneumoniae"
    ),
    part = "header"
  ) %>%
  font(
    fontname = "Times New Roman",
    part = "all"
  ) %>%
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  autofit()

##############################################################################
# Export
##############################################################################

doc <- read_docx() %>%
  body_add_flextable(
    esbl_table
  )

print(
  doc,
  target = "Results/ESBL_Ecoli_Kpneumoniae_Summary_Table.docx"
)
###############################################################################
# CARBAPENEM-RESISTANT (CR) ENTEROBACTERALES
###############################################################################

###############################################################################
# Initial descriptive summaries
###############################################################################

df_grouped_count <- joined_data %>%
  filter(
    str_detect(
      coalesce(`COLONY MORPHOLOGY ON CARBA`, ""),
      "(?i)Red|Pinkish|Metalic blue"
    )
  ) %>%
  distinct(
    INIKA_ID,
    Isolate,
    `COLONY MORPHOLOGY ON CARBA`
  ) %>%
  count(
    Isolate,
    name = "unique_INIKA_count"
  )

Joined_data_results_by_isolate <- joined_data %>%
  filter(
    str_detect(
      coalesce(`COLONY MORPHOLOGY ON CARBA`, ""),
      "(?i)Red|Pinkish|Metalic blue"
    )
  ) %>%
  group_by(Isolate) %>%
  summarise(
    unique_INIKA_total =
      n_distinct(INIKA_ID),
    
    with_VITEK_MS =
      n_distinct(
        INIKA_ID[
          !is.na(VITEK_MS_Results) &
            VITEK_MS_Results != ""
        ]
      ),
    
    with_Isolate_NVI =
      n_distinct(
        INIKA_ID[
          !is.na(Isolate_NVI) &
            Isolate_NVI != ""
        ]
      ),
    
    .groups = "drop"
  )

###############################################################################
# Wilson CI helper
###############################################################################

calc_rate_ci <- function(
    n_pos,
    n_total
) {
  
  if (n_total == 0) {
    
    return(
      list(
        cr_n_pct = "0 (0%)",
        ci_95 = "0 - 0",
        combined = "0 (0%) [0 - 0]"
      )
    )
  }
  
  ci_res <- binom::binom.confint(
    x = n_pos,
    n = n_total,
    methods = "wilson"
  )
  
  pct <- (n_pos / n_total) * 100
  
  lower <- max(
    0,
    ci_res$lower * 100
  )
  
  upper <- min(
    100,
    ci_res$upper * 100
  )
  
  n_pct_str <- sprintf(
    "%d (%s%%)",
    n_pos,
    fmt_num(pct)
  )
  
  ci_str <- sprintf(
    "%s - %s",
    fmt_num(lower),
    fmt_num(upper)
  )
  
  list(
    cr_n_pct = n_pct_str,
    ci_95 = ci_str,
    combined = sprintf(
      "%s [%s]",
      n_pct_str,
      ci_str
    )
  )
}

###############################################################################
# Data cleaning and isolate confirmation
###############################################################################

Joined_data_clean <- joined_data %>%
  select(
    INIKA_ID,
    ORIGIN_OF_SAMPLE,
    SEASON.x,
    REGION.x,
    DISTRICT.x,
    `COLONY MORPHOLOGY ON CARBA`,
    VITEK_MS_Results,
    Isolate,
    PROTOCOL,
    Isolate_NVI
  ) %>%
  rename(
    SEASON = SEASON.x,
    REGION = REGION.x,
    DISTRICT = DISTRICT.x
  ) %>%
  mutate(
    
    is_initial_ecoli =
      str_detect(
        coalesce(Isolate, ""),
        "(?i)E\\.?\\s*coli"
      ),
    
    is_vitek_ecoli =
      str_detect(
        coalesce(VITEK_MS_Results, ""),
        "(?i)Escherichia\\s+coli"
      ),
    
    is_nvi_ecoli =
      str_detect(
        coalesce(Isolate_NVI, ""),
        "(?i)E\\.?\\s*coli|Escherichia\\s+coli"
      ),
    
    is_nvi_kpneuma =
      str_detect(
        coalesce(Isolate_NVI, ""),
        "(?i)K\\.?\\s*pneumoniae|Klebsiella\\s+pneumoniae"
      ),
    
    CONFIRMED_SPECIES = case_when(
      is_initial_ecoli &
        (is_vitek_ecoli | is_nvi_ecoli) ~
        "CR Escherichia coli",
      
      is_initial_ecoli &
        is_nvi_kpneuma ~
        "CR Klebsiella pneumoniae",
      
      TRUE ~ NA_character_
    )
  ) %>%
  rename(
    `Sample origin` = ORIGIN_OF_SAMPLE,
    Season = SEASON,
    Region = REGION,
    District = DISTRICT
  )

###############################################################################
# Denominators and presumptive CR isolates
###############################################################################

denom_df <- Joined_data_clean %>%
  distinct(
    INIKA_ID,
    .keep_all = TRUE
  )

cr_isolates_df <- Joined_data_clean %>%
  filter(
    str_detect(
      coalesce(`COLONY MORPHOLOGY ON CARBA`, ""),
      "(?i)Red|Pinkish|Metalic blue"
    )
  ) %>%
  distinct(
    INIKA_ID,
    .keep_all = TRUE
  )

###############################################################################
# Strata computation
###############################################################################

compute_strata_table <- function(
    data_denom,
    data_cr,
    var_label,
    group_col
) {
  
  categories <- data_denom %>%
    pull({{group_col}}) %>%
    unique() %>%
    na.omit()
  
  map_dfr(categories, function(cat) {
    
    n_tot <- data_denom %>%
      filter({{group_col}} == cat) %>%
      pull(INIKA_ID) %>%
      n_distinct()
    
    data_cr_sub <- data_cr %>%
      filter({{group_col}} == cat)
    
    n_pos <- data_cr_sub %>%
      pull(INIKA_ID) %>%
      n_distinct()
    
    ci_stats <- calc_rate_ci(
      n_pos,
      n_tot
    )
    
    n_ecoli <- data_cr_sub %>%
      filter(
        CONFIRMED_SPECIES ==
          "CR Escherichia coli"
      ) %>%
      pull(INIKA_ID) %>%
      n_distinct()
    
    n_kpneu <- data_cr_sub %>%
      filter(
        CONFIRMED_SPECIES ==
          "CR Klebsiella pneumoniae"
      ) %>%
      pull(INIKA_ID) %>%
      n_distinct()
    
    ci_ecoli <- calc_rate_ci(
      n_ecoli,
      n_tot
    )
    
    ci_kpneu <- calc_rate_ci(
      n_kpneu,
      n_tot
    )
    
    tibble(
      Variable = var_label,
      Category = as.character(cat),
      
      `Total test` = n_tot,
      
      `pCR Isolates n (%) [95% CI]` =
        ci_stats$combined,
      
      `CR Escherichia coli n (%) [95% CI]` =
        ci_ecoli$combined,
      
      `CR Klebsiella pneumoniae n (%) [95% CI]` =
        ci_kpneu$combined
    )
  })
}

###############################################################################
# Overall row
###############################################################################

total_tested_overall <- n_distinct(
  denom_df$INIKA_ID
)

total_cr_overall <- n_distinct(
  cr_isolates_df$INIKA_ID
)

overall_ci <- calc_rate_ci(
  total_cr_overall,
  total_tested_overall
)

overall_ecoli <- cr_isolates_df %>%
  filter(
    CONFIRMED_SPECIES ==
      "CR Escherichia coli"
  ) %>%
  pull(INIKA_ID) %>%
  n_distinct()

overall_kpneu <- cr_isolates_df %>%
  filter(
    CONFIRMED_SPECIES ==
      "CR Klebsiella pneumoniae"
  ) %>%
  pull(INIKA_ID) %>%
  n_distinct()

overall_row <- tibble(
  Variable = "Overall",
  Category = "Overall",
  
  `Total test` =
    total_tested_overall,
  
  `pCR Isolates n (%) [95% CI]` =
    overall_ci$combined,
  
  `CR Escherichia coli n (%) [95% CI]` =
    calc_rate_ci(
      overall_ecoli,
      total_tested_overall
    )$combined,
  
  `CR Klebsiella pneumoniae n (%) [95% CI]` =
    calc_rate_ci(
      overall_kpneu,
      total_tested_overall
    )$combined
)

###############################################################################
# Stratified sections
###############################################################################

sample_origin_rows <- compute_strata_table(
  denom_df,
  cr_isolates_df,
  "Sample origin",
  `Sample origin`
)

season_rows <- compute_strata_table(
  denom_df,
  cr_isolates_df,
  "Season",
  Season
)

region_rows <- compute_strata_table(
  denom_df,
  cr_isolates_df,
  "Region",
  Region
)

kilimanjaro_districts <- denom_df %>%
  filter(
    Region == "Kilimanjaro"
  ) %>%
  pull(District) %>%
  unique() %>%
  na.omit() %>%
  sort()

mwanza_districts <- denom_df %>%
  filter(
    Region == "Mwanza"
  ) %>%
  pull(District) %>%
  unique() %>%
  na.omit() %>%
  sort()

ordered_districts <- c(
  kilimanjaro_districts,
  mwanza_districts
)

district_rows <- map_dfr(
  ordered_districts,
  function(dist) {
    
    compute_strata_table(
      denom_df %>% filter(District == dist),
      cr_isolates_df %>% filter(District == dist),
      "District",
      District
    )
  }
)

###############################################################################
# Final table
###############################################################################

final_table_df <- bind_rows(
  overall_row,
  sample_origin_rows,
  season_rows,
  region_rows,
  district_rows
) %>%
  mutate(
    Variable =
      if_else(
        duplicated(Variable),
        "",
        Variable
      )
  )

###############################################################################
# Flextable
###############################################################################

ft <- flextable(final_table_df) %>%
  theme_booktabs() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Table: Carbapenem-Resistant Enterobacterales Isolation and Species Confirmation",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft)

print(
  doc,
  target =
    "Results/CR_Enterobacterales_Confirmation_Table.docx"
)

###############################################################################
# Export isolate-level subset
###############################################################################

cr_isolates_subset <- joined_data %>%
  filter(
    str_detect(
      coalesce(`COLONY MORPHOLOGY ON CARBA`, ""),
      "(?i)Red|Pinkish|Metalic blue"
    )
  ) %>%
  distinct(
    INIKA_ID,
    Isolate,
    `COLONY MORPHOLOGY ON CARBA`,
    .keep_all = TRUE
  ) %>%
  select(
    INIKA_ID,
    SEASON,
    REGION,
    DISTRICT,
    `COLONY MORPHOLOGY ON CARBA`,
    PROTOCOL,
    Isolate,
    VITEK_MS_Results,
    Isolate_NVI
  )

  
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




##############################################################################
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
)