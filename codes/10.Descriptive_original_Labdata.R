
library(dplyr)
library(rlang)
library(binom)
library(purrr)
library(tibble)
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(tidyr)

# This is the function for the calculation of prevalence that makes it automatically, and gives you a nice table in the end
# Only need to run this once, so I place it above all other codes!

calculate_prevalence <- function(data, target_var, value, group_vars = NULL, conf_level = 0.95) {
  target_var <- ensym(target_var)
  
  if (!is.null(group_vars)) {
    group_syms <- syms(group_vars)
    
    result <- data %>%
      group_by(!!!group_syms) %>%
      summarise(
        total = n(),
        count = sum(!!target_var == value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        ci = pmap(list(count, total), ~ binom::binom.confint(..1, ..2, conf.level = conf_level, methods = "wilson")),
        prevalence_percent = round((count / total) * 100, 2),
        ci_lower = map_dbl(ci, ~ round(.x$lower[1] * 100, 2)),
        ci_upper = map_dbl(ci, ~ round(.x$upper[1] * 100, 2))
      ) %>%
      select(-ci)
  } else {
    total <- nrow(data)
    count <- sum(data[[as_string(target_var)]] == value, na.rm = TRUE)
    ci <- binom::binom.confint(count, total, conf.level = conf_level, methods = "wilson")
    
    result <- tibble(
      total = total,
      count = count,
      prevalence_percent = round((count / total) * 100, 2),
      ci_lower = round(ci$lower[1] * 100, 2),
      ci_upper = round(ci$upper[1] * 100, 2)
    )
  }
  
  return(result)
}


#We will work on the  total dataset, but need to make sure we don't have duplications of results
#for any of the INIKA_IDs (Samples) in one column, so we need to make a new column for Klebsiella_isolates and for Eco_isolates as well as for Salmonellla_isolates,
#the same will we do for the VITEK_MS_Results, 
#All Isolates not E.coli will in that column be negative,
# the same will apply for Klebsiella, we can also do one column containing ESC_Isolate; regardless of Klebsiella or E.coli,
#but not for the isolate you will not have the AST results at all, so we will miss that criteria for a lot of them
# Here you have to select or not select the columns needed for the analyses

## Importing the joined_data file
library(readr)
joined_data <- read_csv("data/CLEANED_DATA/joined_data.csv")
View(joined_data)

Descriptive_Lab <- joined_data%>%
  mutate(Results_ECO = case_when(
    Isolate.x == "E.coli" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_KLEBP = case_when(
    Isolate.x == "K.pneumoniae" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_SAM = case_when(
    Isolate.x == "S.typhimurium" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_Styphi = case_when(
    Isolate.x == "S.typhi" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_K.aerogenes = case_when(
    Isolate.x == "K.aerogenes" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_K.oxytoca = case_when(
    Isolate.x == "K.oxytoca" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_C.freundii = case_when(
    Isolate.x == "C.freundii" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_S.paratyphiA = case_when(
    Isolate.x == "S.paratyphi A" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_P.aeruginosa = case_when(
    Isolate.x == "P.aeruginosa" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_Y.enterocolitica = case_when(
    Isolate.x == "Yersinia enterocolitica" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_UnidentifiedGNR = case_when(
    Isolate.x == "Unidentified GNR" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_Acinetobacterspp = case_when(
    Isolate.x == "Acinetobacter spp" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_A.hydrophilia = case_when(
    Isolate.x == "Aeromonas hydrophilia" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_C.koseri = case_when(
    Isolate.x == "C.koseri" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_E.cloacae = case_when(
    Isolate.x == "E.cloacae" ~ "1",
    TRUE ~ "0"))%>%
  mutate(Results_Shigellaspp = case_when(
    Isolate.x == "Shigella spp" ~ "1",
    TRUE ~ "0"))%>%
  select(INIKA_ID.x, Results_ECO, Results_KLEBP , Results_K.oxytoca, Results_K.aerogenes, 
         Results_Styphi, Results_SAM,Results_Shigellaspp,
         Results_E.cloacae,Results_C.koseri,Results_A.hydrophilia,Results_Acinetobacterspp,
         Results_UnidentifiedGNR,Results_Y.enterocolitica,Results_S.paratyphiA,
         Results_C.freundii,Results_P.aeruginosa,
         REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x,
  ) 

UNIKDESC<-unique(Descriptive_Lab )
  # need to select for each Isolate separately
 DescriptECO<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_ECO, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 # Note we are getting duplicates, need to remove negative results
 

 
 # Step 1: Count how many times each ID appears
 
 DescriptECO_flagged <- DescriptECO %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptECO_final <- DescriptECO_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_ECO)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
# Step 3. Renaming the columns
 
 unique_DATA <- DescriptECO_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 #Make sure to select them first in the dataset joined_data ( which was done in the 9.Joining...that you will work on!
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_ECO,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_ECO <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_ECO, "Results/prevalence_table_ECO.rds")
 write_tsv(prevalence_table_ECO, "Results/prevalence_table_ECO.tsv")
 
 # This can now be repeated for each of the bacteria using the datasets provided above:
 #Results_KLEBP (Klebsiella pneumoniae , Salmonella Typhimurium, Klebseilla oxytoca and so on!) 
# From Step 1, Replace with the actualk datasets needed
            
#############################################################################          
 
 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptKLEBP<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_KLEBP, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 

 # Step 1: Count how many times each ID appears
 
 DescriptKLEBP_flagged <- DescriptKLEBP %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptKLEBP_final <- DescriptKLEBP_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_KLEBP)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptKLEBP_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_KLEBP,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_KLEBP <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_KLEBP, "Results/prevalence_table_KLEBP.rds")
 write_tsv(prevalence_table_KLEBP, "Results/prevalence_table_KLEBP.tsv")
 
######################################################################################### 

 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptK.oxytoca<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_K.oxytoca, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 

 # Step 1: Count how many times each ID appears
 
 DescriptK.oxytoca_flagged <- DescriptK.oxytoca %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptK.oxytoca_final <- DescriptK.oxytoca_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_K.oxytoca)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptK.oxytoca_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_K.oxytoca,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_K.oxytoca <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.oxytoca, "Results/prevalence_table_K.oxytoca.rds")
 write_tsv(prevalence_table_K.oxytoca, "Results/prevalence_table_K.oxytoca.tsv")
 
######################################################################################### 
  
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptK.aerogenes<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_K.aerogenes, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptK.aerogenes_flagged <- DescriptK.aerogenes %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptK.aerogenes_final <- DescriptK.aerogenes_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_K.aerogenes)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptK.aerogenes_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_K.aerogenes,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_K.aerogenes <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.aerogenes, "Results/prevalence_table_K.aerogenes.rds")
 write_tsv(prevalence_table_K.aerogenes, "Results/prevalence_table_K.aerogenes.tsv")
#################################################################################### 

 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptSAM<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_SAM, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptSAM_flagged <- DescriptSAM %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptSAM_final <- DescriptSAM_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_SAM)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptSAM_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_SAM,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_SAM <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_SAM, "Results/prevalence_table_SAM.rds")
 write_tsv(prevalence_table_SAM, "Results/prevalence_table_SAM.tsv") 
################################################################################# 
 
 UNIKDESC<-unique(Descriptive_Lab)
 # need to select for each Isolate separately
 DescriptStyphi<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Styphi, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 

 
 # Step 1: Count how many times each ID appears
 
 DescriptStyphi_flagged <- DescriptStyphi %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptStyphi_final <- DescriptStyphi_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Styphi)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptStyphi_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_Styphi,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_Styphi <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Styphi, "Results/prevalence_table_Styphi.rds")
 write_tsv(prevalence_table_Styphi, "Results/prevalence_table_Styphi.tsv")
################################################################################# 

 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptStyphi<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Styphi, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptStyphi_flagged <- DescriptStyphi %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptStyphi_final <- DescriptStyphi_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Styphi)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptStyphi_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_Styphi,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_Styphi <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Styphi, "Results/prevalence_table_Styphi.rds")
 write_tsv(prevalence_table_Styphi, "Results/prevalence_table_Styphi.tsv")
################################################################################# 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptS.paratyphiA<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_S.paratyphiA, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)

 
 # Step 1: Count how many times each ID appears
 
 DescriptS.paratyphiA_flagged <- DescriptS.paratyphiA %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptS.paratyphiA_final <- DescriptS.paratyphiA_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_S.paratyphiA)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptS.paratyphiA_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_S.paratyphiA,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_S.paratyphiA <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_S.paratyphiA, "Results/prevalence_table_S.paratyphiA.rds")
 write_tsv(prevalence_table_S.paratyphiA, "Results/prevalence_table_S.paratyphiA.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptC.freundii<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_C.freundii, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptC.freundii_flagged <- DescriptC.freundii %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptC.freundii_final <- DescriptC.freundii_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_C.freundii)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptC.freundii_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_C.freundii,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_C.freundii <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_C.freundii, "Results/prevalence_table_C.freundii.rds")
 write_tsv(prevalence_table_C.freundii, "Results/prevalence_table_C.freundii.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptC.koseri<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_C.koseri, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptC.koseri_flagged <- DescriptC.koseri %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptC.koseri_final <- DescriptC.koseri_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_C.koseri)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptC.koseri_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_C.koseri,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_C.koseri <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_C.koseri, "Results/prevalence_table_C.koseri.rds")
 write_tsv(prevalence_table_C.koseri, "Results/prevalence_table_C.koseri.tsv") 
################################################################################
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptP.aeruginosa<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_P.aeruginosa, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptP.aeruginosa_flagged <- DescriptP.aeruginosa %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptP.aeruginosa_final <- DescriptP.aeruginosa_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_P.aeruginosa)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptP.aeruginosa_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_P.aeruginosa,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_P.aeruginosa <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_P.aeruginosa, "Results/prevalence_table_P.aeruginosa.rds")
 write_tsv(prevalence_table_P.aeruginosa, "Results/prevalence_table_P.aeruginosa.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptShigellaspp<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Shigellaspp, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptShigellaspp_flagged <- DescriptShigellaspp %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptShigellaspp_final <- DescriptShigellaspp_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Shigellaspp)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptShigellaspp_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_Shigellaspp,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_Shigellaspp <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Shigellaspp, "Results/prevalence_table_Shigellaspp.rds")
 write_tsv(prevalence_table_Shigellaspp, "Results/prevalence_table_Shigellaspp.tsv")
 ################################################################################
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptY.enterocolitica<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Y.enterocolitica, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)

 
 # Step 1: Count how many times each ID appears
 
 DescriptY.enterocolitica_flagged <- DescriptY.enterocolitica %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptY.enterocolitica_final <- DescriptY.enterocolitica_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Y.enterocolitica)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptY.enterocolitica_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_Y.enterocolitica,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_Y.enterocolitica <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Y.enterocolitica, "Results/prevalence_table_Y.enterocoliticap.rds")
 write_tsv(prevalence_table_Y.enterocolitica, "Results/prevalence_table_Y.enterocolitica.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptA.hydrophilia<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_A.hydrophilia, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptA.hydrophilia_flagged <- DescriptA.hydrophilia %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptA.hydrophilia_final <- DescriptA.hydrophilia_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_A.hydrophilia)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptA.hydrophilia_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_A.hydrophilia,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_A.hydrophilia <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_A.hydrophilia, "Results/prevalence_table_A.hydrophilia.rds")
 write_tsv(prevalence_table_A.hydrophilia, "Results/prevalence_table_A.hydrophilia.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptAcinetobacterspp<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Acinetobacterspp, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptAcinetobacterspp_flagged <- DescriptAcinetobacterspp %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptAcinetobacterspp_final <- DescriptAcinetobacterspp_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Acinetobacterspp)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptAcinetobacterspp_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_Acinetobacterspp,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_Acinetobacterspp <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Acinetobacterspp, "Results/prevalence_table_Acinetobacterspp.rds")
 write_tsv(prevalence_table_Acinetobacterspp, "Results/prevalence_table_Acinetobacterspp.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptE.cloacae<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_E.cloacae, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptE.cloacae_flagged <- DescriptE.cloacae %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptE.cloacae_final <- DescriptE.cloacae_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_E.cloacae)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptE.cloacae_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_E.cloacae,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_E.cloacae <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_E.cloacae, "Results/prevalence_table_E.cloacae.rds")
 write_tsv(prevalence_table_E.cloacae, "Results/prevalence_table_E.cloacae.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptUnidentifiedGNR<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_UnidentifiedGNR, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 

 # Step 1: Count how many times each ID appears
 
 DescriptUnidentifiedGNR_flagged <- DescriptUnidentifiedGNR %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptUnidentifiedGNR_final <- DescriptUnidentifiedGNR_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_UnidentifiedGNR)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA <- DescriptUnidentifiedGNR_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     unique_DATA,
     target_var = Results_UnidentifiedGNR,
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_UnidentifiedGNR <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_UnidentifiedGNR, "Results/prevalence_table_UnidentifiedGNR.rds")
 write_tsv(prevalence_table_UnidentifiedGNR, "Results/prevalence_table_UnidentifiedGNR.tsv")
################################################################################ 
 # Occurrences of ESCR E.coli and K.pneumoniae
 
 ESCR_Descriptive_Lab <- joined_data%>%
   mutate(Results_E.coli = case_when(
     VITEK_MS_Results == "Escherichia coli" ~ "1",
     TRUE ~ "0"))%>%
   mutate(Results_K.pneumoniae = case_when(
     VITEK_MS_Results == "Klebsiella pneumoniae" ~ "1",
     TRUE ~ "0"))%>%
   
   select(TVLA_ID,ESBL_Selection, ESBL, Results_E.coli, Results_K.pneumoniae,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x,
   ) 
 
 UNIKESCR<-unique(ESCR_Descriptive_Lab )
 # need to select for each Isolate separately
 ESCR_E.coli<-ESCR_Descriptive_Lab%>% 
   select(TVLA_ID, Results_E.coli, TVLA_ID,ESBL_Selection, ESBL,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 # Note we are getting duplicates, need to remove negative results
 
 
 
 # Step 1: Count how many times each ID appears
 
 ESCR_DescriptiveE.coli_flagged <- ESCR_E.coli %>%
   group_by(TVLA_ID) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 ESCR_DescriptiveE.coli_final <- ESCR_DescriptiveE.coli_flagged %>%
   group_by(TVLA_ID) %>%
   arrange(desc(Results_E.coli)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 ESCR_unique_DATA <- ESCR_DescriptiveE.coli_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 #Make sure to select them first in the dataset joined_data ( which was done in the 9.Joining...that you will work on!
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     ESCR_unique_DATA, 
     target_var = Results_E.coli, 
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_E.coli <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_E.coli, "Results/prevalence_table_E.coli.rds")
 write_tsv(prevalence_table_E.coli, "Results/prevalence_table_E.coli.tsv")
 
################################################################################
 # Occurrences of ESCR K.pneumoniae
 
 ESCR_Descriptive_Lab <- joined_data%>%
   mutate(Results_E.coli = case_when(
     VITEK_MS_Results == "Escherichia coli" ~ "1",
     TRUE ~ "0"))%>%
   mutate(Results_K.pneumoniae = case_when(
     VITEK_MS_Results == "Klebsiella pneumoniae" ~ "1",
     TRUE ~ "0"))%>%
   
   select(TVLA_ID,ESBL_Selection, ESBL, Results_E.coli, Results_K.pneumoniae,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x,
   ) 
 
 UNIKESCR<-unique(ESCR_Descriptive_Lab )
 # need to select for each Isolate separately
 ESCR_K.pneumoniae<-ESCR_Descriptive_Lab%>% 
   select(TVLA_ID, Results_K.pneumoniae, TVLA_ID,ESBL_Selection, ESBL,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 # Note we are getting duplicates, need to remove negative results
 
 
 
 # Step 1: Count how many times each ID appears
 
 ESCR_DescriptiveK.pneumoniae_flagged <- ESCR_K.pneumoniae %>%
   group_by(TVLA_ID) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 ESCR_DescriptiveK.pneumoniae_final <- ESCR_DescriptiveK.pneumoniae_flagged %>%
   group_by(TVLA_ID) %>%
   arrange(desc(Results_K.pneumoniae)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 ESCR_unique_DATA <- ESCR_DescriptiveK.pneumoniae_final %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE= ORIGIN_OF_SAMPLE.x
   )
 
 
 # Step 4. Loop over grouping variables_ You can add more variables once you have control!,
 #Make sure to select them first in the dataset joined_data ( which was done in the 9.Joining...that you will work on!
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 for (var in group_vars_list) {
   result <- calculate_prevalence(
     ESCR_unique_DATA, 
     target_var = Results_K.pneumoniae, 
     value = "1",
     group_vars = c(var)
   )
   prevalence_results[[var]] <- result
 }
 # Combine all results into one table
 prevalence_table <- bind_rows(
   lapply(names(prevalence_results), function(var) {
     prevalence_results[[var]] %>%
       mutate(grouping_variable = var) %>%
       rename(Group_Value = !!sym(var))
   }),
   .id = "group_id"
 )
 
 # Final table one for each bacteria.
 prevalence_table_K.pneumoniae <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.pneumoniae, "Results/prevalence_table_K.pneumoniae.rds")
 write_tsv(prevalence_table_K.pneumoniae, "Results/prevalence_table_K.pneumoniae.tsv")
######################################################################## 
 ## ESBL prevalence of E.coli
 # Define the list of grouping variables
 group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x", "DISTRICT.x")
 
 # --- Step 1: Prepare the Base Data (Filter E. coli and ensure unique isolates) ---
 Ecoli_data <- joined_data %>%
   # Filter only Escherichia coli isolates
   filter(VITEK_MS_Results == "Escherichia coli") %>%
   # Use TVLA_ID to ensure each unique isolate is counted only once
   distinct(TVLA_ID, .keep_all = TRUE)
 
 # --- Step 2: Define the Prevalence Calculation Function ---
 calculate_prevalence_ci <- function(data, group_var) {
   
   # Calculate success (ESBL positive) and total (all E.coli) counts
   counts <- data %>%
     group_by(!!sym(group_var)) %>%
     summarise(
       N_ESBL_Positive = sum(ESBL == "1", na.rm = TRUE), 
       N_Total_Ecoli = n(),                          
       .groups = 'drop'
     ) %>%
     filter(N_Total_Ecoli > 0)
   
   # Apply prop.test to each row (group) to get CI
   ci_results <- counts %>%
     rowwise() %>%
     mutate(
       Prop_Test = list(broom::tidy(prop.test(x = N_ESBL_Positive, n = N_Total_Ecoli, conf.level = 0.95))),
       # Add the Group_Type column here, before unnesting
       Group_Type = group_var 
     ) %>%
     ungroup() %>%
     tidyr::unnest(Prop_Test) %>% 
     
     # Calculate Prevalence and format to Percentage (Rounded to 1 decimal place)
     mutate(
       Prevalence_Percent = round((N_ESBL_Positive / N_Total_Ecoli) * 100, 1),
       Lower_CI_95_Percent = round(conf.low * 100, 1),
       Upper_CI_95_Percent = round(conf.high * 100, 1)
     ) %>%
     
     # Select and rename final columns, ensuring Group_Type is included
     select(
       Group_Type, # <--- Correctly select the column created above
       Grouping_Value = !!sym(group_var),
       N_ESBL_Positive,
       N_Total_Ecoli,
       Prevalence_Percent,
       Lower_CI_95_Percent,
       Upper_CI_95_Percent
     ) %>%
     # Relocate Group_Type to the first position (optional, but good practice)
     relocate(Group_Type, .before = 1)
   
   return(ci_results)
 }
 # --- Step 3: Apply the function across all grouping variables and combine results ---
 prevalence_results_list <- purrr::map(group_vars_list, ~calculate_prevalence_ci(Ecoli_data, .x))
 
 ESBL_Prevalence_Table_Ecoli_Pct <- bind_rows(prevalence_results_list)
 
 # View the final table
 print(head(ESBL_Prevalence_Table_Ecoli_Pct))
  write.csv(ESBL_Prevalence_Table_Ecoli_Pct, "Results/ESBL_Ecoli_Prevalence_Percent_by_Group.csv")
  write_tsv(ESBL_Prevalence_Table_Ecoli_Pct, "Results/ESBL_Ecoli_Prevalence_Percent_by_Group.tsv")
###############################################################################
  ## ESBL prevalence of K.pneumoniae
  
  # Define the list of grouping variables
  group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x", "DISTRICT.x")
  
  # --- Step 1: Prepare the Base Data (Filter E. coli and ensure unique isolates) ---
  Kpn_data <- joined_data %>%
    # Filter only Escherichia coli isolates
    filter(VITEK_MS_Results == "Klebsiella pneumoniae") %>%
    # Use TVLA_ID to ensure each unique isolate is counted only once
    distinct(TVLA_ID, .keep_all = TRUE)
  
  # --- Step 2: Define the Prevalence Calculation Function ---
  calculate_prevalence_ci <- function(data, group_var) {
    
    # Calculate success (ESBL positive) and total (all Kpn) counts
    counts <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(
        N_ESBL_Positive = sum(ESBL == "1", na.rm = TRUE), 
        N_Total_Kpn = n(),                          
        .groups = 'drop'
      ) %>%
      # Use the column name created in the summarise step:
      filter(N_Total_Kpn > 0)                       
    
    # Apply prop.test to each row (group) to get CI
    ci_results <- counts %>%
      rowwise() %>%
      mutate(
        Prop_Test = list(broom::tidy(prop.test(x = N_ESBL_Positive, n = N_Total_Kpn, conf.level = 0.95))),
        # Add the Group_Type column here, before unnesting
        Group_Type = group_var 
      ) %>%
      ungroup() %>%
      tidyr::unnest(Prop_Test) %>% 
      
      # Calculate Prevalence and format to Percentage (Rounded to 1 decimal place)
      mutate(
        Prevalence_Percent = round((N_ESBL_Positive / N_Total_Kpn) * 100, 1),
        Lower_CI_95_Percent = round(conf.low * 100, 1),
        Upper_CI_95_Percent = round(conf.high * 100, 1)
      ) %>%
      
      # Select and rename final columns, ensuring Group_Type is included
      select(
        Group_Type, # <--- Correctly select the column created above
        Grouping_Value = !!sym(group_var),
        N_ESBL_Positive,
        N_Total_Kpn,
        Prevalence_Percent,
        Lower_CI_95_Percent,
        Upper_CI_95_Percent
      ) %>%
      # Relocate Group_Type to the first position (optional, but good practice)
      relocate(Group_Type, .before = 1)
    
    return(ci_results)
  }
  # --- Step 3: Apply the function across all grouping variables and combine results ---
  prevalence_results_list <- purrr::map(group_vars_list, ~calculate_prevalence_ci(Kpn_data, .x))
  
  ESBL_Prevalence_Table_Kpn_Pct <- bind_rows(prevalence_results_list)
  
  # View the final table
  print(head(ESBL_Prevalence_Table_Kpn_Pct))
  write.csv(ESBL_Prevalence_Table_Kpn_Pct, "Results/ESBL_Kpn_Prevalence_Percent_by_Group.csv")
  write_tsv(ESBL_Prevalence_Table_Kpn_Pct, "Results/ESBL_Kpn_Prevalence_Percent_by_Group.tsv")
######################################################################################
  ## Concordant results between Biochemical and VITEK MS 
  
  library(dplyr)
  library(stringr)
  library(dplyr)
  library(stringr)
  
  # --- Step 1: Calculate the total number of unique isolates (Denominator) ---
  Total_Unique_Isolates <- joined_data %>%
    distinct(TVLA_ID) %>%
    nrow()
  
  # --- Step 2: Define the Concordance Calculation Function ---
  
  calculate_organism_concordance <- function(data, isolate_name, vitek_name) {
    
    # Count the unique TVLA_IDs where both columns match for the specific organism
    N_Concordant <- data %>%
      filter(
        Isolate.x == isolate_name & VITEK_MS_Results == vitek_name
      ) %>%
      distinct(TVLA_ID) %>%
      nrow()
    
    # Calculate Percentage
    Percentage <- round((N_Concordant / Total_Unique_Isolates) * 100, 2)
    
    return(
      data.frame(
        Organism = vitek_name,
        N_Concordant = N_Concordant,
        Total_Unique_Isolates = Total_Unique_Isolates,
        Concordance_Percentage = Percentage
      )
    )
  }
  
  # --- Step 3: Calculate results for E. coli and K. pneumoniae separately ---
  
  # E. coli Concordance
  Ecoli_Concordance <- calculate_organism_concordance(
    joined_data, 
    isolate_name = "E.coli", 
    vitek_name = "Escherichia coli"
  )
  
  # K. pneumoniae Concordance
  Kpneumoniae_Concordance <- calculate_organism_concordance(
    joined_data, 
    isolate_name = "K.pneumoniae", 
    vitek_name = "Klebsiella pneumoniae"
  )
  
  # --- Step 4: Combine and Print the Final Results ---
  
  Separate_Concordance_Summary <- bind_rows(Ecoli_Concordance, Kpneumoniae_Concordance)
  
  print(Separate_Concordance_Summary)
  
  # Save file 
  write.csv(Separate_Concordance_Summary, "Results/Separate_Concordance_Summary.csv")
  write_tsv(Separate_Concordance_Summary, "Results/Separate_Concordance_Summary.tsv")
#####################################################################################
## Frequencies
  
  get_frequencies <- function(df) {
    freq_list <- lapply(d, function(col) {
      freq_table <- as.data.frame(table(col, useNA = "ifany"))
      colnames(freq_table) <- c("Value", "Frequency")
      return(freq_table)
    })
    names(freq_list) <- colnames(df)
    return(freq_list)
  }

  
 # Below a  quick way of getting the frequencies of the variables
  #Apply the function to the dataset
  
  frequencies <- get_frequencies(Descriptive_Lab)
  
  # Combine the frequencies into a single dataframe
  combined_frequencies <- bind_rows(
    lapply(names(frequencies), function(colname) {
      freq_table <- frequencies[[colname]]
      freq_table <- mutate(freq_table, Column = colname)
    }),
    .id = "id"
  )
  
  
  
 # writing out the table to have as a helping file- opened in Excel, where all the antibiotic names were identified and corrected below
  
  
  
  write_tsv(combined_frequencies, "Results/combined_frequencies.tsv")
  
###################################################################################
  library(dplyr)
  library(purrr) # Needed for the map functions, often used with dplyr
  
  # 1. Define the function (as provided)
  get_frequencies <- function(df) {
    freq_list <- lapply(df, function(col) {
      freq_table <- as.data.frame(table(col, useNA = "ifany"))
      colnames(freq_table) <- c("Value", "Frequency")
      return(freq_table)
    })
    names(freq_list) <- colnames(df)
    return(freq_list)
  }
  
  # 2. Apply the function to the 'joined_data' dataset
  frequencies <- get_frequencies(joined_data)
  
  # 3. Combine the frequencies into a single dataframe
  combined_frequencies <- bind_rows(
    lapply(names(frequencies), function(colname) {
      freq_table <- frequencies[[colname]]
      freq_table <- mutate(freq_table, Column = colname)
    }),
    .id = "id"
  ) %>%
    # Clean up the final columns
    select(Column, Value, Frequency)
  
  # Save the file
  write_csv(combined_frequencies, "Results/combined_frequencies.tsv")
  write_csv(frequencies, "Results/combined_frequencies.tsv")
  
#################################################################