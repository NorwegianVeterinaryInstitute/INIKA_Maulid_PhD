
library(tidyverse)  # for data manipulation and visualization
library(readxl)     # for reading Excel files)
library(rlang)
library(binom)




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

joined_data <- read_csv("data/CLEANED_DATA/UniqueData.csv")
spec(joined_data)

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
 DescriptECO1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_ECO, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptECO1_flagged <- DescriptECO1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptECO1_final <- DescriptECO1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_ECO)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
# Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptECO1_final %>%
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
     unique_DATA1,
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
 prevalence_table_ECO1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_ECO1, "Results/prevalence_table_ECO1.rds")
 write_tsv(prevalence_table_ECO1, "Results/prevalence_table_ECO1.tsv")
#############################################################################          
 
 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptKLEBP1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_KLEBP, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 

 # Step 1: Count how many times each ID appears
 
 DescriptKLEBP1_flagged <- DescriptKLEBP1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptKLEBP1_final <- DescriptKLEBP1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_KLEBP)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptKLEBP1_final %>%
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
     unique_DATA1,
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
 prevalence_table_KLEBP1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_KLEBP1, "Results/prevalence_table_KLEBP1.rds")
 write_tsv(prevalence_table_KLEBP1, "Results/prevalence_table_KLEBP1.tsv")
 
######################################################################################### 

 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptK.oxytoca1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_K.oxytoca, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 

 # Step 1: Count how many times each ID appears
 
 DescriptK.oxytoca1_flagged <- DescriptK.oxytoca1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptK.oxytoca1_final <- DescriptK.oxytoca1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_K.oxytoca)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptK.oxytoca1_final %>%
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
     unique_DATA1,
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
 prevalence_table_K.oxytoca1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.oxytoca1, "Results/prevalence_table_K.oxytoca1.rds")
 write_tsv(prevalence_table_K.oxytoca1, "Results/prevalence_table_K.oxytoca1.tsv")
 
######################################################################################### 
  
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptK.aerogenes1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_K.aerogenes, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptK.aerogenes1_flagged <- DescriptK.aerogenes1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptK.aerogenes1_final <- DescriptK.aerogenes1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_K.aerogenes)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptK.aerogenes1_final %>%
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
     unique_DATA1,
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
 prevalence_table_K.aerogenes1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.aerogenes1, "Results/prevalence_table_K.aerogenes1.rds")
 write_tsv(prevalence_table_K.aerogenes1, "Results/prevalence_table_K.aerogenes1.tsv")
#################################################################################### 

 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptSAM1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_SAM, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptSAM1_flagged <- DescriptSAM1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptSAM1_final <- DescriptSAM1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_SAM)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptSAM1_final %>%
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
     unique_DATA1,
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
 prevalence_table_SAM1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_SAM1, "Results/prevalence_table_SAM1.rds")
 write_tsv(prevalence_table_SAM1, "Results/prevalence_table_SAM1.tsv") 
################################################################################# 
 
 UNIKDESC<-unique(Descriptive_Lab)
 # need to select for each Isolate separately
 DescriptStyphi1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Styphi, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 

 
 # Step 1: Count how many times each ID appears
 
 DescriptStyphi1_flagged <- DescriptStyphi1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptStyphi1_final <- DescriptStyphi1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Styphi)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptStyphi1_final %>%
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
     unique_DATA1,
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
 prevalence_table_Styphi1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Styphi1, "Results/prevalence_table_Styphi1.rds")
 write_tsv(prevalence_table_Styphi1, "Results/prevalence_table_Styphi1.tsv")
################################################################################# 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptS.paratyphiA1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_S.paratyphiA, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)

 
 # Step 1: Count how many times each ID appears
 
 DescriptS.paratyphiA1_flagged <- DescriptS.paratyphiA1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptS.paratyphiA1_final <- DescriptS.paratyphiA1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_S.paratyphiA)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptS.paratyphiA1_final %>%
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
     unique_DATA1,
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
 prevalence_table_S.paratyphiA1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_S.paratyphiA1, "Results/prevalence_table_S.paratyphiA1.rds")
 write_tsv(prevalence_table_S.paratyphiA1, "Results/prevalence_table_S.paratyphiA1.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptC.freundii1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_C.freundii, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptC.freundii1_flagged <- DescriptC.freundii1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptC.freundii1_final <- DescriptC.freundii1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_C.freundii)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptC.freundii1_final %>%
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
     unique_DATA1,
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
 prevalence_table_C.freundii1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_C.freundii1, "Results/prevalence_table_C.freundii1.rds")
 write_tsv(prevalence_table_C.freundii1, "Results/prevalence_table_C.freundii1.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptC.koseri1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_C.koseri, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptC.koseri1_flagged <- DescriptC.koseri1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptC.koseri1_final <- DescriptC.koseri1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_C.koseri)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptC.koseri1_final %>%
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
     unique_DATA1,
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
 prevalence_table_C.koseri1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_C.koseri1, "Results/prevalence_table_C.koseri1.rds")
 write_tsv(prevalence_table_C.koseri1, "Results/prevalence_table_C.koseri1.tsv") 
################################################################################
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptP.aeruginosa1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_P.aeruginosa, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptP.aeruginosa1_flagged <- DescriptP.aeruginosa1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptP.aeruginosa1_final <- DescriptP.aeruginosa1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_P.aeruginosa)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptP.aeruginosa1_final %>%
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
     unique_DATA1,
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
 prevalence_table_P.aeruginosa1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_P.aeruginosa1, "Results/prevalence_table_P.aeruginosa1.rds")
 write_tsv(prevalence_table_P.aeruginosa1, "Results/prevalence_table_P.aeruginosa1.tsv") 
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptShigellaspp1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Shigellaspp, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptShigellaspp1_flagged <- DescriptShigellaspp1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptShigellaspp1_final <- DescriptShigellaspp1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Shigellaspp)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptShigellaspp1_final %>%
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
     unique_DATA1,
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
 prevalence_table_Shigellaspp1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Shigellaspp1, "Results/prevalence_table_Shigellaspp1.rds")
 write_tsv(prevalence_table_Shigellaspp1, "Results/prevalence_table_Shigellaspp1.tsv")
 ################################################################################
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptY.enterocolitica1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Y.enterocolitica, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)

 
 # Step 1: Count how many times each ID appears
 
 DescriptY.enterocolitica1_flagged <- DescriptY.enterocolitica1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptY.enterocolitica1_final <- DescriptY.enterocolitica1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Y.enterocolitica)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptY.enterocolitica1_final %>%
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
     unique_DATA1,
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
 prevalence_table_Y.enterocolitica1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Y.enterocolitica1, "Results/prevalence_table_Y.enterocolitica1.rds")
 write_tsv(prevalence_table_Y.enterocolitica1, "Results/prevalence_table_Y.enterocolitica1.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptA.hydrophilia1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_A.hydrophilia, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptA.hydrophilia1_flagged <- DescriptA.hydrophilia1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptA.hydrophilia1_final <- DescriptA.hydrophilia1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_A.hydrophilia)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptA.hydrophilia1_final %>%
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
     unique_DATA1,
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
 prevalence_table_A.hydrophilia1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_A.hydrophilia1, "Results/prevalence_table_A.hydrophilia1.rds")
 write_tsv(prevalence_table_A.hydrophilia1, "Results/prevalence_table_A.hydrophilia1.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptAcinetobacterspp1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_Acinetobacterspp, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptAcinetobacterspp1_flagged <- DescriptAcinetobacterspp1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptAcinetobacterspp1_final <- DescriptAcinetobacterspp1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_Acinetobacterspp)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptAcinetobacterspp1_final %>%
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
     unique_DATA1,
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
 prevalence_table_Acinetobacterspp1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_Acinetobacterspp1, "Results/prevalence_table_Acinetobacterspp1.rds")
 write_tsv(prevalence_table_Acinetobacterspp1, "Results/prevalence_table_Acinetobacterspp1.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptE.cloacae1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_E.cloacae, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 
 
 
 # Step 1: Count how many times each ID appears
 
 DescriptE.cloacae1_flagged <- DescriptE.cloacae1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptE.cloacae1_final <- DescriptE.cloacae1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_E.cloacae)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptE.cloacae1_final %>%
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
     unique_DATA1,
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
 prevalence_table_E.cloacae1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_E.cloacae1, "Results/prevalence_table_E.cloacae1.rds")
 write_tsv(prevalence_table_E.cloacae1, "Results/prevalence_table_E.cloacae1.tsv")
################################################################################ 
 UNIKDESC<-unique(Descriptive_Lab )
 # need to select for each Isolate separately
 DescriptUnidentifiedGNR1<-Descriptive_Lab%>% 
   select(INIKA_ID.x, Results_UnidentifiedGNR, 
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 

 # Step 1: Count how many times each ID appears
 
 DescriptUnidentifiedGNR1_flagged <- DescriptUnidentifiedGNR1 %>%
   group_by(INIKA_ID.x) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 DescriptUnidentifiedGNR1_final <- DescriptUnidentifiedGNR1_flagged %>%
   group_by(INIKA_ID.x) %>%
   arrange(desc(Results_UnidentifiedGNR)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 unique_DATA1 <- DescriptUnidentifiedGNR1_final %>%
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
     unique_DATA1,
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
 prevalence_table_UnidentifiedGNR1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_UnidentifiedGNR1, "Results/prevalence_table_UnidentifiedGNR1.rds")
 write_tsv(prevalence_table_UnidentifiedGNR1, "Results/prevalence_table_UnidentifiedGNR1.tsv")
################################################################################ 
 # ESCR for E.coli Trial
 
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(stats)
 
 # ====================================================================
 # !!! IMPORTANT: REPLACE THIS SYNTHETIC DATA SETUP WITH YOUR DATA !!!
 # ====================================================================
 # Assuming 'joined_data' is your source data frame containing:
 # TVLA_ID, Isolate.x, VITEK_MS_Results, ESBL, REGION.x, DISTRICT.x, 
 # ORIGIN_OF_SAMPLE.x, and SEASON.x.
 
 # --- 1. Custom Function to Calculate Prevalence and 95% CI ---
 
 # Computes prevalence and CI for given numerator (x) and denominator (n).
 calculate_prevalence_CI <- function(x, n) {
   if (n == 0) {
     return(data.frame(
       N = 0, ESBL_N = 0, Prevalence = NA, Lower_CI = NA, Upper_CI = NA
     ))
   }
   
   # Use binom.test for the exact 95% CI (Clopper-Pearson)
   ci_result <- binom.test(x = x, n = n, conf.level = 0.95)
   
   return(data.frame(
     N = n,
     ESBL_N = x,
     Prevalence = round(ci_result$estimate * 100, 2), # As percentage
     Lower_CI = round(ci_result$conf.int[1] * 100, 2), # As percentage
     Upper_CI = round(ci_result$conf.int[2] * 100, 2)  # As percentage
   ))
 }
 
 # --------------------------------------------------------------------
 # --- 2. DATA PREPARATION: DEFINING NUMERATOR AND DENOMINATOR ---
 # --------------------------------------------------------------------
 
 # 2.1. Define the Denominator Population
 # Denominator: All unique isolates confirmed as E. coli by VITEK-MS.
 Denominator_Data <- joined_data %>%
   filter(VITEK_MS_Results == "Escherichia coli") %>%
   # Select grouping variables and TVLA_ID
   select(TVLA_ID, REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x) %>%
   # Deduplicate to count only unique isolates
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>% 
   ungroup() %>%
   rename(
     REGION = REGION.x,
     DISTRICT = DISTRICT.x,
     SEASON = SEASON.x,
     ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x
   )
 
 # 2.2. Define the Numerator Population
 # Numerator: Isolate.x in (E.coli, K.pneumoniae) AND VITEK_MS_Results is E. coli AND ESBL == 1.
 numerator_acceptable_isolates <- c("E.coli", "K.pneumoniae") 
 
 Numerator_Data <- joined_data %>%
   filter(
     Isolate.x %in% numerator_acceptable_isolates & 
       VITEK_MS_Results == "Escherichia coli" &       
       as.numeric(ESBL) == 1                          
   ) %>%
   select(TVLA_ID) %>%
   # Deduplicate to count only unique isolates in the numerator
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>%
   ungroup() %>%
   mutate(is_numerator = 1) # Flag for counting ESBL positive
 
 # 2.3. Merge the two datasets for aggregation.
 Prevalence_Data_Merged <- Denominator_Data %>%
   left_join(Numerator_Data, by = "TVLA_ID") %>%
   # If an isolate is in the denominator but not the numerator, its flag is 0
   mutate(is_numerator = replace_na(is_numerator, 0))
 
 
 # --------------------------------------------------------------------
 # --- 3. LOOP AND GROUPED CALCULATION ---
 # --------------------------------------------------------------------
 
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 prevalence_results <- list()
 
 for (var in group_vars_list) {
   # Aggregate the numerator and denominator counts within the current group
   count_summary <- Prevalence_Data_Merged %>%
     group_by(!!sym(var)) %>%
     summarise(
       x_successes = sum(is_numerator), # Numerator count (ESBL positive from specific initial isolates)
       n_total = n()                   # Denominator count (Total E. coli in group)
     ) %>%
     ungroup()
   
   # Apply the custom CI function to the aggregated counts
   result <- count_summary %>%
     rowwise() %>%
     mutate(prevalence_data = list(calculate_prevalence_CI(x = x_successes, n = n_total))) %>%
     unnest(prevalence_data) %>%
     ungroup() %>%
     # Format for final table
     mutate(grouping_variable = var) %>%
     rename(Group_Value = !!sym(var)) %>%
     select(-x_successes, -n_total) 
   
   prevalence_results[[var]] <- result
 }
 
 # --------------------------------------------------------------------
 # --- 4. COMBINE AND CLEAN THE FINAL TABLE ---
 # --------------------------------------------------------------------
 
 final_prevalence_table <- bind_rows(prevalence_results) %>%
   relocate(grouping_variable, Group_Value, .before = 1) %>%
   select(
     `Grouping Variable` = grouping_variable, 
     `Group Value` = Group_Value, 
     `Total E.coli (N)` = N, 
     `ESBL Positive (N)` = ESBL_N, 
     `Prevalence (%)` = Prevalence, 
     `95% CI Lower (%)` = Lower_CI, 
     `95% CI Upper (%)` = Upper_CI
   ) %>%
   arrange(`Grouping Variable`, `Group Value`)
 
 # Print the final result table
 print(final_prevalence_table)
 # Save the file
 write_tsv(final_prevalence_table, "Results/162_presumptiveESBL_Ecoli_final_prevalence_table.tsv")
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
 ESCR_E.coli1<-ESCR_Descriptive_Lab%>%
   select(TVLA_ID, Results_E.coli, TVLA_ID,ESBL_Selection, ESBL,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 # Note we are getting duplicates, need to remove negative results



 # Step 1: Count how many times each ID appears

 ESCR_DescriptiveE.coli1_flagged <- ESCR_E.coli1 %>%
   group_by(TVLA_ID) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 ESCR_DescriptiveE.coli1_final <- ESCR_DescriptiveE.coli1_flagged %>%
   group_by(TVLA_ID) %>%
   arrange(desc(Results_E.coli)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()

 # Step 3. Renaming the columns

 ESCR_unique_DATA1 <- ESCR_DescriptiveE.coli1_final %>%
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
     ESCR_unique_DATA1,
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
 prevalence_table_E.coli1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_E.coli1, "Results/prevalence_table_E.coli1.rds")
 write_tsv(prevalence_table_E.coli1, "Results/prevalence_table_E.coli1.tsv")
 
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
 ESCR_K.pneumoniae1<-ESCR_Descriptive_Lab%>% 
   select(TVLA_ID, Results_K.pneumoniae, TVLA_ID,ESBL_Selection, ESBL,
          REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x)
 # Note we are getting duplicates, need to remove negative results
 
 
 
 # Step 1: Count how many times each ID appears
 
 ESCR_DescriptiveK.pneumoniae1_flagged <- ESCR_K.pneumoniae1 %>%
   group_by(TVLA_ID) %>%
   mutate(ID_repeat_count = n()) %>%
   ungroup()
 # Step 2: Create a filtered version that keeps only one row per ID based on your logic
 ESCR_DescriptiveK.pneumoniae1_final <- ESCR_DescriptiveK.pneumoniae1_flagged %>%
   group_by(TVLA_ID) %>%
   arrange(desc(Results_K.pneumoniae)) %>%  # Prioritize "1" over "0"
   slice_head(n = 1) %>%           # Keep only the first row per ID
   ungroup()
 
 # Step 3. Renaming the columns
 
 ESCR_unique_DATA1 <- ESCR_DescriptiveK.pneumoniae1_final %>%
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
     ESCR_unique_DATA1, 
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
 prevalence_table_K.pneumoniae1 <- prevalence_table %>%
   relocate(grouping_variable, Group_Value, .before = 1)
 ## Save the file
 saveRDS(prevalence_table_K.pneumoniae1, "Results/prevalence_table_K.pneumoniae1.rds")
 write_tsv(prevalence_table_K.pneumoniae1, "Results/prevalence_table_K.pneumoniae1.tsv")
######################################################################## 
 # Prevalence of ESBL(ESCR) E.coli
 group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x", "DISTRICT.x")
 
 # --- Step 1: Prepare the Base Data (Filter E. coli AND K. pneumoniae isolates) ---
 Ecoli_Kpneumoniae_data <- joined_data %>%
   # Filter based on the OR condition in Isolate.x AND the VITEK condition
   filter(
     Isolate.x %in% c("E.coli", "K.pneumoniae") &
       VITEK_MS_Results == "Escherichia coli"
   ) %>%
   # Use TVLA_ID to ensure each unique isolate is counted only once
   distinct(TVLA_ID, .keep_all = TRUE)
 
 # --- Step 2: Define the Prevalence Calculation Function (Prevalence ONLY) ---
 calculate_prevalence <- function(data, group_var) {
   
   # Calculate success (ESBL positive) and total (all E.coli) counts
   prevalence_table <- data %>%
     group_by(!!sym(group_var)) %>%
     summarise(
       N_ESBL_Positive = sum(ESBL == "1", na.rm = TRUE),
       N_Total_Ecoli = n(),
       .groups = 'drop'
     ) %>%
     filter(N_Total_Ecoli > 0) %>%
     
     # Calculate Prevalence and format to Percentage (Rounded to 1 decimal place)
     mutate(
       Group_Type = group_var,
       Prevalence_Percent = round((N_ESBL_Positive / N_Total_Ecoli) * 100, 1)
     ) %>%
     
     # Select and rename final columns
     select(
       Group_Type,
       Grouping_Value = !!sym(group_var),
       N_ESBL_Positive,
       N_Total_Ecoli,
       Prevalence_Percent
     ) %>%
     relocate(Group_Type, .before = 1)
   
   return(prevalence_table)
 }
 
 # --- Step 3: Apply the function across all grouping variables and combine results ---
 prevalence_results_list <- purrr::map(group_vars_list, ~calculate_prevalence(Ecoli_Kpneumoniae_data, .x))
 
 ESBL_Prevalence_Table_Ecoli_Pct <- bind_rows(prevalence_results_list)
 
 # Display the final table (optional)
  print(ESBL_Prevalence_Table_Ecoli_Pct)
 
 # Save the fie
  write_tsv(ESBL_Prevalence_Table_Ecoli_Pct,"Results/131_Obs_ECO_ESBL_Prevalence_Table_Ecoli_Pct.tsv")
################################################################################ 
 
 ## ESBL prevalence of E.coli with 95%CI
 # Define the list of grouping variables
 group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x", "DISTRICT.x")
 
 # --- Step 1: Prepare the Base Data (Filter E. coli AND K. pneumoniae isolates) ---
 Ecoli_Kpneumoniae_data <- joined_data %>%
   # Filter based on the OR condition in Isolate.x AND the VITEK condition
   filter(
     Isolate.x %in% c("E.coli", "K.pneumoniae") &
       VITEK_MS_Results == "Escherichia coli"
   ) %>%
   # Use TVLA_ID to ensure each unique isolate is counted only once
   distinct(TVLA_ID, .keep_all = TRUE)
 
 # Step 2: Define the Prevalence Calculation Function (Function remains UNCHANGED) ---
 calculate_prevalence_ci <- function(data, group_var) {
   
   # Calculate success (ESBL positive) and total (all E.coli/K.pneumoniae) counts
   counts <- data %>%
     group_by(!!sym(group_var)) %>%
     summarise(
       # The ESBL column is already filtered to include 1s
       N_ESBL_Positive = sum(ESBL == "1", na.rm = TRUE),
       N_Total_Ecoli = n(), # Renamed to N_Total_Isolates for clarity, but kept N_Total_Ecoli for consistency with original code
       .groups = 'drop'
     ) %>%
     filter(N_Total_Ecoli > 0)
   
   # Apply prop.test to each row (group) to get CI
   ci_results <- counts %>%
     rowwise() %>%
     mutate(
       Prop_Test = list(broom::tidy(prop.test(x = N_ESBL_Positive, n = N_Total_Ecoli, conf.level = 0.95))),
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
       Group_Type,
       Grouping_Value = !!sym(group_var),
       N_ESBL_Positive,
       N_Total_Ecoli,
       Prevalence_Percent,
       Lower_CI_95_Percent,
       Upper_CI_95_Percent
     ) %>%
     relocate(Group_Type, .before = 1)
   
   return(ci_results)
 }
 
 # --- Step 3: Apply the function across all grouping variables and combine results ---
 # Note the change to use the new data frame name
 prevalence_results_list <- purrr::map(group_vars_list, ~calculate_prevalence_ci(Ecoli_Kpneumoniae_data, .x))
 
 ESBL_Prevalence_Table_Ecoli_Pct <- bind_rows(prevalence_results_list)
 
 print(head(ESBL_Prevalence_Table_Ecoli_Pct))
 write.csv(ESBL_Prevalence_Table_Ecoli_Pct, "Results/131_Obs_ESBL_Ecoli_Prevalence_Percent_by_Group.csv")
 write_tsv(ESBL_Prevalence_Table_Ecoli_Pct, "Results/131_Obs_ESBL_Ecoli_Prevalence_Percent_by_Group.tsv")
###############################################################################
 # Chi-square calculation for statistical Isolation comparison for E.coli
 
 # FUNCTION FOR PREVALENCE AND 95% CI ---
 # This function calculates CI but returns a list of values to be expanded row-wise.
 calculate_prevalence_CI <- function(x, n) {
   if (n == 0) {
     return(list(
       Prevalence = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_
     ))
   }
   
   # Use binom.test for the exact 95% CI (Clopper-Pearson)
   ci_result <- binom.test(x = x, n = n, conf.level = 0.95)
   
   return(list(
     Prevalence = round(ci_result$estimate * 100, 2), 
     Lower_CI = round(ci_result$conf.int[1] * 100, 2), 
     Upper_CI = round(ci_result$conf.int[2] * 100, 2)
   ))
 }
 
 # DATA PREPARATION: DEFINING NUMERATOR AND DENOMINATOR ---

 # Define the list of grouping variables
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 
 # Define the Denominator Population 
 Denominator_Data <- joined_data %>%
   filter(VITEK_MS_Results == "Escherichia coli") %>%
   select(TVLA_ID, REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x) %>%
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>%
   ungroup() %>%
   rename(
     REGION = REGION.x, DISTRICT = DISTRICT.x, SEASON = SEASON.x, ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x
   )
 
 # Define the Numerator Population (ESBL Positive)
 numerator_acceptable_isolates <- c("E.coli", "K.pneumoniae")
 
 Numerator_Data <- joined_data %>%
   filter(
     Isolate.x %in% numerator_acceptable_isolates &
       VITEK_MS_Results == "Escherichia coli" &
       as.numeric(ESBL) == 1
   ) %>%
   select(TVLA_ID) %>%
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>%
   ungroup() %>%
   mutate(is_numerator = 1)
 
 #  Merge the two datasets for aggregation.
 Prevalence_Data_Merged <- Denominator_Data %>%
   left_join(Numerator_Data, by = "TVLA_ID") %>%
   mutate(is_numerator = replace_na(is_numerator, 0))
 
 # AGGREGATION AND STATISTICAL CALCULATION ---
 
 all_results <- list()
 i <- 1
 
 ## Overall Total Calculation 
 
 # Aggregate for the Overall Total
 overall_summary <- Prevalence_Data_Merged %>%
   summarise(
     x_successes = sum(is_numerator),
     n_total = n()
   )
 
 # Apply CI function and calculate working columns
 overall_result_ci <- overall_summary %>%
   rowwise() %>%
   mutate(ci_data = list(calculate_prevalence_CI(x = x_successes, n = n_total))) %>%
   unnest_wider(ci_data) %>%
   ungroup() %>%
   mutate(
     `Grouping Variable` = "OVERALL_TOTAL",
     `Group Value` = "Total",
     ChiSq_P_Value = NA_real_,
     ChiSq_Test_Result = "Overall Prevalence (No comparison possible)"
   ) %>%
   #  Standardize Final Column Selection and Naming ---
   select(
     `Grouping Variable`, 
     `Group Value`, 
     `Total E.coli (N)` = n_total, 
     `ESBL Positive (N)` = x_successes,
     `Prevalence (%)` = Prevalence, 
     `95% CI Lower (%)` = Lower_CI, 
     `95% CI Upper (%)` = Upper_CI, 
     ChiSq_P_Value, 
     ChiSq_Test_Result
   )
 
 all_results[[i]] <- overall_result_ci
 i <- i + 1
 
 ## Segregated Calculations with Chi-Square Comparison 
 
 for (var in group_vars_list) {
   # Aggregate the numerator and denominator counts within the current group
   count_summary <- Prevalence_Data_Merged %>%
     group_by(!!sym(var)) %>%
     summarise(
       x_successes = sum(is_numerator),
       n_total = n()
     ) %>%
     ungroup()
   
   # Apply the custom CI function and calculate working columns
   result_ci <- count_summary %>%
     rowwise() %>%
     mutate(ci_data = list(calculate_prevalence_CI(x = x_successes, n = n_total))) %>%
     unnest_wider(ci_data) %>%
     ungroup() %>%
     # Calculate the column needed for the Chi-Square test:
     mutate(ESBL_Negative_N = n_total - x_successes) 
   
   # Perform Chi-Square Test (Only if >= 2 categories exist)
   if (nrow(result_ci) >= 2) {
     # Create the 2xN contingency matrix
     contingency_matrix <- result_ci %>%
       select(x_successes, ESBL_Negative_N) %>%
       as.matrix()
     
     # Run prop.test (equivalent to Chi-square for 2xN table)
     chisq_result <- tryCatch({
       prop.test(contingency_matrix)
     }, error = function(e) {
       warning(paste("Chi-square failed for", var, ":", e$message))
       return(NULL)
     })
     
     if (!is.null(chisq_result)) {
       p_value <- chisq_result$p.value
       test_result <- paste0(
         "X2(", chisq_result$parameter, ", N=", sum(result_ci$n_total), ") = ",
         round(chisq_result$statistic, 2), ", p=", round(chisq_result$p.value, 4)
       )
     } else {
       p_value <- NA_real_
       test_result <- "Chi-Square Test Failed/Warning"
     }
   } else {
     p_value <- NA_real_
     test_result <- "Insufficient categories (N < 2 for ChiSq)"
   }
   
   #  Format and Combine Results
   result_formatted <- result_ci %>%
     mutate(
       `Grouping Variable` = var,
       `Group Value` = !!sym(var),
       ChiSq_P_Value = p_value,
       ChiSq_Test_Result = test_result
     ) %>%
     # Standardize Final Column Selection and Naming ---
     select(
       `Grouping Variable`,
       `Group Value`,
       `Total E.coli (N)` = n_total,
       `ESBL Positive (N)` = x_successes,
       `Prevalence (%)` = Prevalence,
       `95% CI Lower (%)` = Lower_CI,
       `95% CI Upper (%)` = Upper_CI,
       ChiSq_P_Value,
       ChiSq_Test_Result
     ) %>%
     arrange(`Group Value`)
   
   all_results[[i]] <- result_formatted
   i <- i + 1
 }
 
 
 #  COMBINE AND CLEAN THE FINAL TABLE ---
 
 ESCR_ECO_prevalence_ChiSq_P_Value <- bind_rows(all_results)
 
 
 # Print the final result table
 print(ESCR_ECO_prevalence_ChiSq_P_Value)
 
 # Save the file
 write_tsv(ESCR_ECO_prevalence_ChiSq_P_Value, "Results/ESCR_ECO_prevalence_ChiSq_P_Value.tsv")
############################################################################### 
 # Chi-square calculation for statistical Isolation comparison for Klebsiella pneumoniae
 
 # FUNCTION FOR PREVALENCE AND 95% CI ---
 # This function calculates CI but returns a list of values to be expanded row-wise.
 calculate_prevalence_CI <- function(x, n) {
   if (n == 0) {
     return(list(
       Prevalence = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_
     ))
   }
   
   # Use binom.test for the exact 95% CI (Clopper-Pearson)
   ci_result <- binom.test(x = x, n = n, conf.level = 0.95)
   
   return(list(
     Prevalence = round(ci_result$estimate * 100, 2), 
     Lower_CI = round(ci_result$conf.int[1] * 100, 2), 
     Upper_CI = round(ci_result$conf.int[2] * 100, 2)
   ))
 }
 
 # DATA PREPARATION: DEFINING NUMERATOR AND DENOMINATOR ---
 
 # Define the list of grouping variables
 group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE")
 
 # Define the Denominator Population 
 Denominator_Data <- joined_data %>%
   filter(VITEK_MS_Results == "Klebsiella pneumoniae") %>%
   select(TVLA_ID, REGION.x, DISTRICT.x, ORIGIN_OF_SAMPLE.x, SEASON.x) %>%
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>%
   ungroup() %>%
   rename(
     REGION = REGION.x, DISTRICT = DISTRICT.x, SEASON = SEASON.x, ORIGIN_OF_SAMPLE = ORIGIN_OF_SAMPLE.x
   )
 
 # Define the Numerator Population (ESBL Positive)
 numerator_acceptable_isolates <- c("E.coli", "K.pneumoniae")
 
 Numerator_Data <- joined_data %>%
   filter(
     Isolate.x %in% numerator_acceptable_isolates &
       VITEK_MS_Results == "Klebsiella pneumoniae" &
       as.numeric(ESBL) == 1
   ) %>%
   select(TVLA_ID) %>%
   group_by(TVLA_ID) %>%
   slice_head(n = 1) %>%
   ungroup() %>%
   mutate(is_numerator = 1)
 
 #  Merge the two datasets for aggregation.
 Prevalence_Data_Merged <- Denominator_Data %>%
   left_join(Numerator_Data, by = "TVLA_ID") %>%
   mutate(is_numerator = replace_na(is_numerator, 0))
 
 # AGGREGATION AND STATISTICAL CALCULATION ---
 
 all_results <- list()
 i <- 1
 
 ## Overall Total Calculation 
 
 # Aggregate for the Overall Total
 overall_summary <- Prevalence_Data_Merged %>%
   summarise(
     x_successes = sum(is_numerator),
     n_total = n()
   )
 
 # Apply CI function and calculate working columns
 overall_result_ci <- overall_summary %>%
   rowwise() %>%
   mutate(ci_data = list(calculate_prevalence_CI(x = x_successes, n = n_total))) %>%
   unnest_wider(ci_data) %>%
   ungroup() %>%
   mutate(
     `Grouping Variable` = "OVERALL_TOTAL",
     `Group Value` = "Total",
     ChiSq_P_Value = NA_real_,
     ChiSq_Test_Result = "Overall Prevalence (No comparison possible)"
   ) %>%
   #  Standardize Final Column Selection and Naming ---
   select(
     `Grouping Variable`, 
     `Group Value`, 
     `Total E.coli (N)` = n_total, 
     `ESBL Positive (N)` = x_successes,
     `Prevalence (%)` = Prevalence, 
     `95% CI Lower (%)` = Lower_CI, 
     `95% CI Upper (%)` = Upper_CI, 
     ChiSq_P_Value, 
     ChiSq_Test_Result
   )
 
 all_results[[i]] <- overall_result_ci
 i <- i + 1
 
 ## Segregated Calculations with Chi-Square Comparison 
 
 for (var in group_vars_list) {
   # Aggregate the numerator and denominator counts within the current group
   count_summary <- Prevalence_Data_Merged %>%
     group_by(!!sym(var)) %>%
     summarise(
       x_successes = sum(is_numerator),
       n_total = n()
     ) %>%
     ungroup()
   
   # Apply the custom CI function and calculate working columns
   result_ci <- count_summary %>%
     rowwise() %>%
     mutate(ci_data = list(calculate_prevalence_CI(x = x_successes, n = n_total))) %>%
     unnest_wider(ci_data) %>%
     ungroup() %>%
     # Calculate the column needed for the Chi-Square test:
     mutate(ESBL_Negative_N = n_total - x_successes) 
   
   # Perform Chi-Square Test (Only if >= 2 categories exist)
   if (nrow(result_ci) >= 2) {
     # Create the 2xN contingency matrix
     contingency_matrix <- result_ci %>%
       select(x_successes, ESBL_Negative_N) %>%
       as.matrix()
     
     # Run prop.test (equivalent to Chi-square for 2xN table)
     chisq_result <- tryCatch({
       prop.test(contingency_matrix)
     }, error = function(e) {
       warning(paste("Chi-square failed for", var, ":", e$message))
       return(NULL)
     })
     
     if (!is.null(chisq_result)) {
       p_value <- chisq_result$p.value
       test_result <- paste0(
         "X2(", chisq_result$parameter, ", N=", sum(result_ci$n_total), ") = ",
         round(chisq_result$statistic, 2), ", p=", round(chisq_result$p.value, 4)
       )
     } else {
       p_value <- NA_real_
       test_result <- "Chi-Square Test Failed/Warning"
     }
   } else {
     p_value <- NA_real_
     test_result <- "Insufficient categories (N < 2 for ChiSq)"
   }
   
   #  Format and Combine Results
   result_formatted <- result_ci %>%
     mutate(
       `Grouping Variable` = var,
       `Group Value` = !!sym(var),
       ChiSq_P_Value = p_value,
       ChiSq_Test_Result = test_result
     ) %>%
     # Standardize Final Column Selection and Naming ---
     select(
       `Grouping Variable`,
       `Group Value`,
       `Total E.coli (N)` = n_total,
       `ESBL Positive (N)` = x_successes,
       `Prevalence (%)` = Prevalence,
       `95% CI Lower (%)` = Lower_CI,
       `95% CI Upper (%)` = Upper_CI,
       ChiSq_P_Value,
       ChiSq_Test_Result
     ) %>%
     arrange(`Group Value`)
   
   all_results[[i]] <- result_formatted
   i <- i + 1
 }
 
 
 #  COMBINE AND CLEAN THE FINAL TABLE ---
 
 ESCR_Kpn_prevalence_ChiSq_P_Value <- bind_rows(all_results)
 
 
 # Print the final result table
 print(ESCR_Kpn_prevalence_ChiSq_P_Value)
 
 # Save the file
 write_tsv(ESCR_Kpn_prevalence_ChiSq_P_Value, "Results/ESCR_Kpn_prevalence_ChiSq_P_Value.tsv")
################################################################################  
  ## ESBL prevalence of K.pneumoniae
  
 group_vars_list <- c("REGION.x", "SEASON.x", "ORIGIN_OF_SAMPLE.x", "DISTRICT.x")
 
 # --- Step 1: Prepare the Base Data (Filter E. coli AND K. pneumoniae isolates) ---
 Ecoli_Kpneumoniae_data <- joined_data %>%
   # Filter based on the OR condition in Isolate.x AND the VITEK condition
   filter(
     Isolate.x %in% c("E.coli", "K.pneumoniae") &
       VITEK_MS_Results == "Klebsiella pneumoniae"
   ) %>%
   # Use TVLA_ID to ensure each unique isolate is counted only once
   distinct(TVLA_ID, .keep_all = TRUE)
 
 # Step 2: Define the Prevalence Calculation Function (Function remains UNCHANGED) ---
 calculate_prevalence_ci <- function(data, group_var) {
   
   # Calculate success (ESBLKpn positive) and total (all E.coli/K.pneumoniae) counts
   counts <- data %>%
     group_by(!!sym(group_var)) %>%
     summarise(
       # The ESBL column is already filtered to include 1s
       N_ESBL_Positive = sum(ESBL == "1", na.rm = TRUE),
       N_Total_Kpn = n(), # Renamed to N_Total_Isolates for clarity, but kept N_Total_Ecoli for consistency with original code
       .groups = 'drop'
     ) %>%
     filter(N_Total_Kpn > 0)
   
   # Apply prop.test to each row (group) to get CI
   ci_results <- counts %>%
     rowwise() %>%
     mutate(
       Prop_Test = list(broom::tidy(prop.test(x = N_ESBL_Positive, n = N_Total_Kpn, conf.level = 0.95))),
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
       Group_Type,
       Grouping_Value = !!sym(group_var),
       N_ESBL_Positive,
       N_Total_Kpn,
       Prevalence_Percent,
       Lower_CI_95_Percent,
       Upper_CI_95_Percent
     ) %>%
     relocate(Group_Type, .before = 1)
   
   return(ci_results)
 }
 
 # --- Step 3: Apply the function across all grouping variables and combine results ---
 # Note the change to use the new data frame name
 prevalence_results_list <- purrr::map(group_vars_list, ~calculate_prevalence_ci(Ecoli_Kpneumoniae_data, .x))
 
 ESBL_Prevalence_Table_Kpn_Pct <- bind_rows(prevalence_results_list)
 
 
 write.csv(ESBL_Prevalence_Table_Kpn_Pct, "Results/54_Obs_ESBL_Kpn_Prevalence_Percent_by_Group.csv")
 write_tsv(ESBL_Prevalence_Table_Kpn_Pct, "Results/54_Obs_ESBL_Kpn_Prevalence_Percent_by_Group.tsv")
 
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
