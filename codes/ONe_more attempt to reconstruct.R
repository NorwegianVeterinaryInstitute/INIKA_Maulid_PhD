
##############################################################################
# SETUP
##############################################################################

library(tidyverse)
library(readxl)
library(writexl)
library(flextable)
library(officer)
library(ggplot2)

##############################################################################
# MASTER DATA
##############################################################################

joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")

spec(joined_data)
names(joined_data)[grep("INIKA", names(joined_data))]

joined_data <- joined_data %>%
  rename(
    INIKA_ID = INIKA_ID.x.x,
    ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x,
    ESBL = ESBL_Presumptivefinal
  )%>%
  select(INIKA_ID, REGION, SEASON, ORIGIN_OF_SAMPLE, Isolate,PROTOCOL,
         VITEK_MS_Results, ESBL, AMX_ED10, AZM_ED15, CRO_ED30,
         CIP_ED5, DOX_ED30, FLR_ED30, GEN_ED10, MEM_ED10,
         OXY_ED30, POL_ED300, SXT_ED1_2, CTX_ED5, CTC_ED30) %>%
  rename(
    Amoxicillin = "AMX_ED10", Azithromycin = "AZM_ED15",
    Ceftriaxone = "CRO_ED30", Ciprofloxacin = "CIP_ED5",
    Doxycycline = "DOX_ED30", Florfenicol = "FLR_ED30",
    Gentamicin = "GEN_ED10", Meropenem = "MEM_ED10",
    Oxytetracycline = "OXY_ED30",
    "Popymyxin_B(PB)" = "POL_ED300",
    "Sulfamethoxazole/Trimethoprim" = "SXT_ED1_2",
    Cefotaxime = "CTX_ED5",
    "Cefotaxime/ClavulanicAcid" = "CTC_ED30"
  ) 

names(joined_data)



##############################################################################
# ECOFF TABLES
##############################################################################

ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx"
)

EPI_CUTOFF <- ECOFF_EUCAST_BREAK_POINT

Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx",
  sheet = 2
)

Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT


ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx"
)

EPI_CUTOFF <- ECOFF_EUCAST_BREAK_POINT

#Import ECOFF table with the sheet for Klebisella
Kpn_ECOFF_EUCAST_BREAK_POINT <- read_excel(
  "data/ECOFF_E.coli_K.pneumoniae.xlsx",
  sheet = "K.pneumoniae_ECOFF_EUCAST"
)

Kpn_EPI_CUTOFF <- Kpn_ECOFF_EUCAST_BREAK_POINT


##############################################################################
# ESCR E. coli Resistance Analysis
##############################################################################

##############################################################################
# 1. DEFINE ESCR E. coli POPULATION
##############################################################################

ESCR_ECO_confirmed <- joined_data %>%
  filter(
    PROTOCOL == "CGR3",
    VITEK_MS_Results == "Escherichia coli",
    Ceftriaxone < 23
  ) %>%
  distinct(INIKA_ID, .keep_all = TRUE)

message(
  paste(
    "Number of ESCR E. coli isolates:",
    nrow(ESCR_ECO_confirmed)
  )
)

##############################################################################
# 2. DEFINE AST VARIABLES
##############################################################################

tested_antibiotics <- c(
  "Amoxicillin",
  "Azithromycin",
  "Ceftriaxone",
  "Ciprofloxacin",
  "Doxycycline",
  "Florfenicol",
  "Gentamicin",
  "Meropenem",
  "Oxytetracycline",
  "Popymyxin_B(PB)",
  "Sulfamethoxazole/Trimethoprim",
  "Cefotaxime",
  "Cefotaxime/ClavulanicAcid"
)

##############################################################################
# 3. CONVERT WIDE TO LONG FORMAT
##############################################################################

ESCR_long <- ESCR_ECO_confirmed %>%
  pivot_longer(
    cols = all_of(tested_antibiotics),
    names_to = "Antimicrobial_substance",
    values_to = "mm"
  ) %>%
  filter(!is.na(mm))

##############################################################################
# 4. HARMONIZE NAMES WITH ECOFF TABLE
##############################################################################

ESCR_long <- ESCR_long %>%
  mutate(
    Antimicrobial_substance = recode(
      Antimicrobial_substance,
      "Popymyxin_B(PB)" = "Polymyxin B",
      "Sulfamethoxazole/Trimethoprim" =
        "Sulfamethoxazole_Trimethoprim",
      "Cefotaxime/ClavulanicAcid" =
        "Cefotaxime _Clavulanic Acid"
    )
  )

##############################################################################
# 5. CHECK NAME MATCHES
##############################################################################

cat("\nNames in AST data missing from ECOFF:\n")

setdiff(
  unique(ESCR_long$Antimicrobial_substance),
  unique(EPI_CUTOFF$Antimicrobial_substance)
)

##############################################################################
# 6. JOIN ECOFF TABLE
##############################################################################

ESCR_long <- ESCR_long %>%
  left_join(
    EPI_CUTOFF,
    by = "Antimicrobial_substance"
  )

##############################################################################
# 7. VERIFY JOIN SUCCESS
##############################################################################

ESCR_long %>%
  summarise(
    Missing_S = sum(is.na(S)),
    Missing_R = sum(is.na(R))
  ) %>%
  print()

##############################################################################
# 8. CALCULATE RESISTANCE
##############################################################################

ESCR_long <- ESCR_long %>%
  mutate(
    mm = as.numeric(mm),
    R_limit = as.numeric(
      str_replace_all(
        as.character(R),
        "[^0-9.]",
        ""
      )
    ),
    is_resistant =
      if_else(
        mm <= R_limit,
        1,
        0,
        missing = 0
      )
  )

##############################################################################
# 9. RESISTANCE SUMMARY
##############################################################################

Summary_Resistance <- ESCR_long %>%
  group_by(Antimicrobial_substance) %>%
  summarise(
    N = n(),
    R_count = sum(is_resistant),
    Prop = R_count / N,
    .groups = "drop"
  ) %>%
  mutate(
    CI = map2(
      R_count,
      N,
      ~ binom.test(.x, .y)$conf.int * 100
    ),
    Lower = map_dbl(CI, 1),
    Upper = map_dbl(CI, 2),
    `Resistance % [95% CI]` =
      sprintf(
        "%.1f%% [%0.1f-%0.1f]",
        Prop * 100,
        Lower,
        Upper
      )
  ) %>%
  select(
    Antimicrobial_substance,
    N,
    R_count,
    `Resistance % [95% CI]`
  )

##############################################################################
# 10. PRINT RESULTS
##############################################################################

print(Summary_Resistance)

##############################################################################
# 11. OPTIONAL ANTIBIOTIC CLASS MAPPING
##############################################################################

drug_class_map <- tibble::tribble(
  ~Antimicrobial_substance,       ~Antibiotic_Class,
  "Amoxicillin",                  "Penicillins",
  "Ceftriaxone",                  "3rd Gen Cephalosporins",
  "Cefotaxime",                   "3rd Gen Cephalosporins",
  "Cefotaxime _Clavulanic Acid",  "Beta-lactam/Inhibitor combinations",
  "Meropenem",                    "