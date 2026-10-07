################################################################################
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(flextable)
library(officer)
################################################################################
# Importing the file
joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
spec(joined_data)
names(joined_data)

#################################################################################
# Computing for the demographic features of study participants

# Column Mapping Configuration
DATASET_COLUMN_NAMES <- list(
  id       = "INIKA_ID", 
  age      = "Age_yrs",
  region   = "REGION", 
  district = "DISTRICT", 
  origin   = "ORIGIN_OF_SAMPLE", 
  gender   = "GENDER",
  season   = "SEASON"
)

ID_COL       <- DATASET_COLUMN_NAMES$id
AGE_COL      <- DATASET_COLUMN_NAMES$age
REGION_COL   <- DATASET_COLUMN_NAMES$region
DISTRICT_COL <- DATASET_COLUMN_NAMES$district
ORIGIN_COL   <- DATASET_COLUMN_NAMES$origin
GENDER_COL   <- DATASET_COLUMN_NAMES$gender
SEASON_COL   <- DATASET_COLUMN_NAMES$season

# --- Pre-processing & Categorization ---

joined_data_processed <- joined_data %>%
  mutate(
    !!ID_COL       := as.character(!!sym(ID_COL)),
    !!REGION_COL   := as.character(!!sym(REGION_COL)),
    
    # Correcting spelling from 'Ilemala' to 'Ilemela' during column standardization
    !!DISTRICT_COL := if_else(as.character(!!sym(DISTRICT_COL)) == "Ilemala", "Ilemela", as.character(!!sym(DISTRICT_COL))),
    
    !!ORIGIN_COL   := as.character(!!sym(ORIGIN_COL)),
    !!GENDER_COL   := as.character(!!sym(GENDER_COL)),
    !!SEASON_COL   := as.character(!!sym(SEASON_COL)),
    
    Age_group = cut(
      !!sym(AGE_COL),
      breaks = c(10, 18, 28, 38, 48, 58, 68, Inf), 
      right = FALSE, 
      labels = c("10 - 17", "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 +"),
      exclude.lowest = TRUE 
    )
  )

# --- Summary Calculation Helper ---

total_unique_ids <- joined_data_processed %>%
  summarise(Total_Unique = n_distinct(!!sym(ID_COL), na.rm = TRUE)) %>% 
  pull(Total_Unique)

calculate_summary_by_id <- function(data, var_name, total_ids, id_col) {
  if (var_name == DISTRICT_COL) {
    res <- data %>%
      group_by(!!sym(REGION_COL), !!sym(var_name)) %>%
      summarise(Frequency = n_distinct(!!sym(id_col), na.rm = TRUE), .groups = 'drop') %>%
      rename(Value = !!sym(var_name), Region_Context = !!sym(REGION_COL)) %>%
      mutate(Variable = "District")
  } else {
    res <- data %>%
      group_by(!!sym(var_name)) %>%
      summarise(Frequency = n_distinct(!!sym(id_col), na.rm = TRUE), .groups = 'drop') %>%
      rename(Value = !!sym(var_name)) %>%
      mutate(Variable = var_name, Region_Context = NA_character_)
  }
  
  res %>%
    mutate(Raw_Percentage = (Frequency / total_ids) * 100) %>%
    select(Variable, Value, Frequency, Raw_Percentage, Region_Context)
}

# Run raw data aggregations
raw_summary <- bind_rows(
  calculate_summary_by_id(joined_data_processed, REGION_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, GENDER_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, ORIGIN_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, SEASON_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, DISTRICT_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, "Age_group", total_unique_ids, ID_COL)
)

# Smart decimal formatter (1 decimal place, whole number if trailing zero)
format_percentage <- function(val) {
  rounded <- round(val, 1)
  formatted <- sprintf("%.1f", rounded)
  sub("\\.0$", "", formatted)
}

# --- Ordering, Grouping & Formatting ---

Demographic_final_pivot_table <- raw_summary %>%
  mutate(
    Variable_Label = case_match(
      Variable,
      REGION_COL ~ "Region",
      GENDER_COL ~ "Gender",
      ORIGIN_COL ~ "Origin of sample",
      SEASON_COL ~ "Season",
      "District" ~ "District",
      "Age_group" ~ "Age group",
      .default = Variable
    ),
    # Hierarchy: Region -> Gender -> Origin -> Season -> District -> Age group
    Var_Order = factor(
      Variable_Label, 
      levels = c("Origin of sample", "Gender",  "Season", "Region","District", "Age group")
    ),
    # Region Sorting (Kilimanjaro first, Mwanza second)
    Region_Order = case_when(
      Variable_Label == "District" & grepl("Kilimanjaro", Region_Context, ignore.case = TRUE) ~ 1,
      Variable_Label == "District" & grepl("Mwanza", Region_Context, ignore.case = TRUE) ~ 2,
      Variable_Label == "Region" & grepl("Kilimanjaro", Value, ignore.case = TRUE) ~ 1,
      Variable_Label == "Region" & grepl("Mwanza", Value, ignore.case = TRUE) ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(Var_Order, Region_Order, Value) %>%
  mutate(
    Percentage = sapply(Raw_Percentage, format_percentage),
    Variable = ifelse(duplicated(Variable_Label), "", Variable_Label),
    Value = as.character(Value)
  ) %>%
  select(Variable, Category = Value, Frequency, Percentage)

# --- Multi-Format Export to "Results" Folder ---

if (!dir.exists("Results")) {
  dir.create("Results", recursive = TRUE)
}

# Plain Text Data Formats
write_csv(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.csv")
write_tsv(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.tsv")
saveRDS(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.rds")

#  Word Format (.docx) Table
border_line <- fp_border(color = "black", width = 1)
header_line <- fp_border(color = "black", width = 0.5)

ft <- flextable(Demographic_final_pivot_table) %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  bold(j = "Variable", part = "body") %>%
  bold(part = "header") %>%
  align(j = c("Variable", "Category"), align = "left", part = "all") %>%
  align(j = c("Frequency", "Percentage"), align = "right", part = "all") %>%
  set_header_labels(
    Variable = "Variable",
    Category = "Category / Level",
    Frequency = "Frequency (n)",
    Percentage = "Percentage (%)"
  ) %>%
  border_remove() %>%
  hline_top(border = border_line, part = "header") %>%
  hline_bottom(border = header_line, part = "header") %>%
  hline_bottom(border = border_line, part = "body") %>%
  autofit() %>%
  add_footer_lines(values = paste0("N = ", total_unique_ids, " unique participants. Percentages are formatted to one decimal place or whole numbers when trailing decimals are zero.")) %>%
  fontsize(size = 8, part = "footer")

# Build Word Document
doc <- read_docx()

caption_p <- fpar(
  ftext("Table 2 ", prop = fp_text(font.family = "Arial", font.size = 10, bold = TRUE)),
  ftext("Demographic and baseline characteristics of study participants", prop = fp_text(font.family = "Arial", font.size = 10))
)

doc <- doc %>%
  body_add_fpar(caption_p) %>%
  body_add_flextable(ft)

print(doc, target = "Results/Demographic_Summary_Table.docx")

cat("All files (.csv, .tsv, .rds, .docx) successfully formatted and saved to 'Results/'.\n")
################################################################################
## 30.09.2026 # MMJ
# Importing the file
joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")
spec(joined_data)
names(joined_data)

#################################################################################
# Computing for the demographic features of study participants

joined_data <- joined_data %>%
  select(
    INIKA_ID = INIKA_ID.x.x, SEASON = SEASON.x, REGION = REGION.x,
    DISTRICT = DISTRICT.x,AGE = Age_yrs,GENDER, ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    "COLONY MORPHOLOGY ON C3GR","COLONY MORPHOLOGY ON CARBA","COLONY MORPHOLOGY ON XLD",
    "COLONY MORPHOLOGY ON BGA","TSI Media_Slope","TSI Media_Butt",
    "TSI Media_Gas","TSIMedia_H2S","SIM Media_H2S","SIM Media Indole",
    "SIM Media Motility","CITRATE" ,"UREASE", "Isolate","Freezing_ID_No",
    "INIKA_prefix","Isolate_suffix",Isolate_ID_repeat_count = Isolate_ID_repeat_count.x,
    Isolate_ID_row_count = Isolate_ID_row_count.x,"TVLA_ID","VITEK_MS_Results",
    "ESC","PROTOCOL","ORGANISM","AMX_ED10","AZM_ED15","CRO_ED30","CIP_ED5","DOX_ED30",
    "FLR_ED30","GEN_ED10","MEM_ED10","POL_ED300","SXT_ED1_2","CTX_ED5","CTC_ED30",
    "CRO_ED5","ESBL_Selection","ESCR_ECO_presumptive","ESCR_KPN_presumptive",
    "ESBL_Presumptivefinal", "Løpenr","Isolate_NVI"
  )
  

# Column Mapping Configuration
DATASET_COLUMN_NAMES <- list(
  id       = "INIKA_ID", 
  age      = "AGE",
  region   = "REGION", 
  district = "DISTRICT", 
  origin   = "ORIGIN_OF_SAMPLE", 
  gender   = "GENDER",
  season   = "SEASON"
)

ID_COL       <- DATASET_COLUMN_NAMES$id
AGE_COL      <- DATASET_COLUMN_NAMES$age
REGION_COL   <- DATASET_COLUMN_NAMES$region
DISTRICT_COL <- DATASET_COLUMN_NAMES$district
ORIGIN_COL   <- DATASET_COLUMN_NAMES$origin
GENDER_COL   <- DATASET_COLUMN_NAMES$gender
SEASON_COL   <- DATASET_COLUMN_NAMES$season

# --- Pre-processing & Categorization ---

joined_data_processed <- joined_data %>%
  mutate(
    !!ID_COL       := as.character(!!sym(ID_COL)),
    !!REGION_COL   := as.character(!!sym(REGION_COL)),
    
    # Correcting spelling from 'Ilemala' to 'Ilemela' during column standardization
    !!DISTRICT_COL := if_else(as.character(!!sym(DISTRICT_COL)) == "Ilemala", "Ilemela", as.character(!!sym(DISTRICT_COL))),
    
    !!ORIGIN_COL   := as.character(!!sym(ORIGIN_COL)),
    !!GENDER_COL   := as.character(!!sym(GENDER_COL)),
    !!SEASON_COL   := as.character(!!sym(SEASON_COL)),
    
    Age_group = cut(
      !!sym(AGE_COL),
      breaks = c(10, 18, 28, 38, 48, 58, 68, Inf), 
      right = FALSE, 
      labels = c("10 - 17", "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 +"),
      exclude.lowest = TRUE 
    )
  )

# --- Summary Calculation Helper ---

total_unique_ids <- joined_data_processed %>%
  summarise(Total_Unique = n_distinct(!!sym(ID_COL), na.rm = TRUE)) %>% 
  pull(Total_Unique)

calculate_summary_by_id <- function(data, var_name, total_ids, id_col) {
  if (var_name == DISTRICT_COL) {
    res <- data %>%
      group_by(!!sym(REGION_COL), !!sym(var_name)) %>%
      summarise(Frequency = n_distinct(!!sym(id_col), na.rm = TRUE), .groups = 'drop') %>%
      rename(Value = !!sym(var_name), Region_Context = !!sym(REGION_COL)) %>%
      mutate(Variable = "District")
  } else {
    res <- data %>%
      group_by(!!sym(var_name)) %>%
      summarise(Frequency = n_distinct(!!sym(id_col), na.rm = TRUE), .groups = 'drop') %>%
      rename(Value = !!sym(var_name)) %>%
      mutate(Variable = var_name, Region_Context = NA_character_)
  }
  
  res %>%
    mutate(Raw_Percentage = (Frequency / total_ids) * 100) %>%
    select(Variable, Value, Frequency, Raw_Percentage, Region_Context)
}

# Run raw data aggregations
raw_summary <- bind_rows(
  calculate_summary_by_id(joined_data_processed, REGION_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, GENDER_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, ORIGIN_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, SEASON_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, DISTRICT_COL, total_unique_ids, ID_COL),
  calculate_summary_by_id(joined_data_processed, "Age_group", total_unique_ids, ID_COL)
)

# Smart decimal formatter (1 decimal place, whole number if trailing zero)
format_percentage <- function(val) {
  rounded <- round(val, 1)
  formatted <- sprintf("%.1f", rounded)
  sub("\\.0$", "", formatted)
}

# --- Ordering, Grouping & Formatting ---

Demographic_final_pivot_table <- raw_summary %>%
  mutate(
    Variable_Label = case_match(
      Variable,
      REGION_COL ~ "Region",
      GENDER_COL ~ "Gender",
      ORIGIN_COL ~ "Origin of sample",
      SEASON_COL ~ "Season",
      "District" ~ "District",
      "Age_group" ~ "Age group",
      .default = Variable
    ),
    # Hierarchy: Region -> Gender -> Origin -> Season -> District -> Age group
    Var_Order = factor(
      Variable_Label, 
      levels = c("Origin of sample", "Gender",  "Season","Region", "District", "Age group")
    ),
    # Region Sorting (Kilimanjaro first, Mwanza second)
    Region_Order = case_when(
      Variable_Label == "District" & grepl("Kilimanjaro", Region_Context, ignore.case = TRUE) ~ 1,
      Variable_Label == "District" & grepl("Mwanza", Region_Context, ignore.case = TRUE) ~ 2,
      Variable_Label == "Region" & grepl("Kilimanjaro", Value, ignore.case = TRUE) ~ 1,
      Variable_Label == "Region" & grepl("Mwanza", Value, ignore.case = TRUE) ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(Var_Order, Region_Order, Value) %>%
  mutate(
    Percentage = sapply(Raw_Percentage, format_percentage),
    Variable = ifelse(duplicated(Variable_Label), "", Variable_Label),
    Value = as.character(Value)
  ) %>%
  select(Variable, Category = Value, Frequency, Percentage)

# --- Multi-Format Export to "Results" Folder ---

if (!dir.exists("Results")) {
  dir.create("Results", recursive = TRUE)
}

# Plain Text Data Formats
write_csv(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.csv")
write_tsv(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.tsv")
saveRDS(Demographic_final_pivot_table, file = "Results/Demographic_Summary_Table.rds")

# Define border styles following publication standards
border_thick <- fp_border(color = "black", width = 1.2)
border_thin  <- fp_border(color = "black", width = 0.5)
border_grey  <- fp_border(color = "#D3D3D3", width = 0.5) # Soft divider between variables

ft <- flextable(Demographic_final_pivot_table) %>%
  # 1. Typography & Font Styling
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 9.5, part = "header") %>%
  bold(j = "Variable", part = "body") %>%
  bold(part = "header") %>%
  
  # 2. Variable Merging & Deduplication (Appears once per block)
  merge_v(j = "Variable") %>%
  valign(j = "Variable", val = "top") %>%
  
  # 3. Column Alignment
  align(j = c("Variable", "Category"), align = "left", part = "all") %>%
  align(j = c("Frequency", "Percentage"), align = "right", part = "all") %>%
  
  # 4. Standard Column Headers
  set_header_labels(
    Variable = "Demographic Variable",
    Category = "Category / Level",
    Frequency = "Frequency (n)",
    Percentage = "Percentage (%)"
  ) %>%
  
  # 5. Cell Padding & Spacing
  padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
  
  # 6. Academic Booktabs Borders
  border_remove() %>%
  hline_top(border = border_thick, part = "header") %>%
  hline_bottom(border = border_thin, part = "header") %>%
  hline_bottom(border = border_thick, part = "body") %>%
  # Adds light horizontal lines between different demographic variables
  hline(i = ~ lead(Variable) != Variable, border = border_grey, part = "body") %>%
  
  # 7. Layout & Footer
  autofit() %>%
  add_footer_lines(
    values = paste0(
      "N = ", total_unique_ids, " unique study participants. ",
      "Percentages are presented relative to valid responses. Values rounded to 1 decimal place."
    )
  ) %>%
  fontsize(size = 8, part = "footer") %>%
  italic(part = "footer")

# Ensure target directory exists
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Build Word Document
doc <- read_docx()

# Format Caption according to APA / Publication Style
caption_p <- fpar(
  ftext("Table 2. ", prop = fp_text(font.family = "Arial", font.size = 10, bold = TRUE)),
  ftext("Demographic and baseline characteristics of study participants.", prop = fp_text(font.family = "Arial", font.size = 10, italic = TRUE))
)

doc <- doc %>%
  body_add_fpar(caption_p) %>%
  body_add_par("", style = "Normal") %>% # Spacer paragraph
  body_add_flextable(ft)

print(doc, target = "Results/Demographic_Summary_Table.docx")

cat("Publication-ready demographic table successfully saved to 'Results/Demographic_Summary_Table.docx'\n")
#############################################################