

#Dataset from Kobotoolbox One duplicate is removed as it is impossible to know which ID to keep!
Kobo<-Demographic_cleaned_data%>%
filter(!INIKA_ID =="21267")

#Dataset from labdata ( This originates from the script 1.Cleaned....) also here the duplication detected in Kobotoolbox are removed!
Labdata<-Cleaned_Labdata_Original%>%
  filter(!INIKA_ID == "21267")

#join Kobo Labdata

joined_Kobo_Labdata<-full_join(Labdata,Kobo, by="INIKA_ID")

Checknotjoined<-anti_join(Labdata,Kobo, by="INIKA_ID")



# we also have the Isolate_ID in Labdata ,
# Madelaine_explain: Here we join by left keeping all data in Labdata and then assigning targeted bacteria to one and those with no growth to 0, all others we assign to EX so we can filter these away for further analyses.

joined_Kobo_Labdata_left<-left_join(Labdata,Kobo, by="INIKA_ID")%>%
  mutate(Results_Bio=case_when(
    Isolate %in% c("E.coli","K.pneumoniae", "S.typhimurium") ~ "1",
    Isolate=="No growth" ~ "0",
                               TRUE ~ "EX"))%>%
  #select(INIKA_ID,Isolate_ID, Isolate_ID_Corrected, Isolate, Results_Bio)%>%
  filter(!Results_Bio == "EX")


NOTJOINED<-anti_join(Labdata,Kobo, by="INIKA_ID")
# None, all joined above!

# try to use the Malditof_all instead of the VITEKcombinedBio
# We need to use the VITEKresults where all confirmed and not confirmed isolates are included
# We first make the negative confirmed that have been assigned as NA to be 0.
# 1.12.25 Note we assign all E.coli and Klebsiella as ESC, this is not quite true as these might not be ESC, they could be CPE or they could be normal E.coli
# depending if they are resistant or not to Cefotaxome or Ceftriaxone.

VITEK<-MALDITOF_RESULTS_cleaned_data%>%
  mutate(
    ESC = if_else(
      VITEK_MS_Results %in% c("Escherichia coli", "Klebsiella pneumoniae"),
      1,
      0
    )
  )

# 1.12.25 Check each step below Madelaine
VITECECO<-VITEK%>%
  filter(VITEK_MS_Results=="Escherichia coli")

VITEKKLEB<-VITEK%>%
  filter(VITEK_MS_Results=="Klebsiella pneumoniae")

#VITEK<-Malditof_all%>%
#  select(Isolate_ID,TVLA_ID, VITEK_MS_Results)

NewselectedKOBO<-joined_Kobo_Labdata_left
  #select(INIKA_ID, Isolate_ID,Isolate_ID_Corrected, Results_Bio)

# When join the one below with the NewselectedKOBO you will miss the ones that are true negative

# Note 1.12.25 the isolate_ID below is not join the selected isolates for conformation at all,
# I have # before to not do anything not needed However I detected that this isolate_ID belongs to the INIKA_ID=2156 and therefore needs to be changed in the 1.2

joinedVITEK_KOBO_LAB_LEFT <- left_join(NewselectedKOBO, VITEK, by = "Isolate_ID") %>%
#  mutate(Isolate_ID = case_when(
#    Isolate_ID == "2146_1_D" ~ "22146_1_D",
#   TRUE ~ Isolate_ID
#  )) %>%
  mutate(ESC = as.character(ESC)) %>%  # Ensure ESC is character
  mutate(ESC = case_when(
    is.na(ESC) & Results_Bio == "0" ~ "0",
    TRUE ~ ESC
  )) %>%
  filter(ESC %in% c("1", "0"))


Negatives<-joinedVITEK_KOBO_LAB_LEFT%>%
  filter(ESC=="0")

# from this we can do the calculation for the occurrence of ESC (Klebsiella  and E.coli , separately eventually together)


# check dublicates 
joinedVITEK_KOBO_LAB_LEFT_duplicates <- joinedVITEK_KOBO_LAB_LEFT$Isolate_ID[duplicated(joinedVITEK_KOBO_LAB_LEFT$Isolate_ID)]
view(joinedVITEK_KOBO_LAB_LEFT_duplicates)

# None left at 1.12.25!
# There were some duplicates where the one isolate ID is a bacteria not of interest or result not identified, which we need to exclude before continuing with sumamrising the results.
# 21228_1_R, 23224_1_R, 23248_2_R, 2377_2_R, 23370_1_R, 23372_1_R


  # Define what counts as "not interesting"
  not_interesting <- c("Enterobacter asburiae  and Enterobacter cloacae", "Enterobacter cloacae and Enterobacter asburiae", "No Identification","Citrobacter Werkmanii" )

# Filter logic
filtered_df <- joinedVITEK_KOBO_LAB_LEFT %>%
  group_by(Isolate_ID) %>%
  filter(
    # Keep all rows if there's only one per isolate_ID
    n() == 1 |
      # Otherwise, keep only rows that are NOT "not interesting"
      !(VITEK_MS_Results %in% not_interesting)
  ) %>%
  ungroup()

filtered_df_duplicates <- filtered_df$Isolate_ID[duplicated(filtered_df$Isolate_ID)]
view(filtered_df_duplicates)
#None left at 1.12.25!
# Now there is only one real duplicate left, 23277_2_R Klebsiella pneumoniae, will be removed when doing unique!

unik_DATA<-unique(filtered_df)
# This data set has only 809 observations left on 1.12.25 ( Note have changed!), But note one sample might have several isolates so we need to remove those that are duplicates if there are more than one isolate per INIKA_ID, but first we need to check how many that can be as if there are some that needs to stay like Klebisella pneumoniae or E.coli then we need to keep them on that condition
 
# Making a dataset only containing the INIKA_IDs
INIKA_ID <-unik_DATA%>%
   select(INIKA_ID.x)

UNIQUE_INIKA_ID <- unique(INIKA_ID )
# These are 797 samples! ( observations ("persons")) This is the real denominator! 

not_interesting <- c("0")

  # Filter logic
  filtered_unique <- unik_DATA%>%
  group_by(INIKA_ID) %>%
  filter(
    # Keep all rows if there's only one per isolate_ID
    n() == 1 |
      # Otherwise, keep only rows that are NOT "not interesting"
      !(ESC %in% not_interesting)
  ) %>%
  ungroup()


  filtered_unique_counts<-filtered_unique%>%
    add_count(INIKA_ID, name = "Isolate_Count")

  # Total observations 784: From this it will only be 781 unique persons ( samples)
  #Note there are 2 isolates both Klebsiella and E.coli from the same person, 
  #but there are also one person with two E.coli 23272-1-R and 23272-2-R?- The last one can only be counted one time!
  



ESCall_summary<-unik_DATA%>%
filter(ESC %in% c("1", "0")) %>%
  
  summarise(
    total = n(),
    count_1 = sum(ESC == "1"),
    percentage_1 = (count_1 / total) * 100
  )

print(ESCall_summary)

ESCnegatives<-joinedVITEK_KOBO_LAB_LEFT%>%
  filter(ESC=="0")

ESCECO_summary<-joinedVITEK_KOBO_LAB_LEFT%>%
  filter(VITEK_MS_Results == "Escherichia coli")

ESCKleb_summary<-joinedVITEK_KOBO_LAB_LEFT%>%
  filter(VITEK_MS_Results == "Klebsiella pneumoniae" )

# The denominator should be ? - It does not work to filter out before.




library(dplyr)
library(rlang)
library(binom)
library(purrr)
library(tibble)


unique_DATA <- filtered_unique %>%
  rename(
    REGION = REGION.x,
    DISTRICT = DISTRICT.x,
    SEASON = SEASON.x
  )

# 1.12.25 Here things needs to be done and reported with care, even if the function here is called prevalence, the Occurrence of the true prevalence will not be possible with the dataset we have as we haven't confirmed all ESC resistant E.coli.
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


group_vars_list <- c("REGION", "SEASON", "DISTRICT", "ORIGIN_OF_SAMPLE", "GENDER" )
prevalence_results <- list()

for (var in group_vars_list) {
  result <- calculate_prevalence(
    unique_DATA,
    target_var = VITEK_MS_Results,
    value = "Klebsiella pneumoniae",
    group_vars = c(var)
  )
  prevalence_results[[var]] <- result
}

# Combine all results into one table with a new column indicating the grouping variable
prevalence_table <- bind_rows(
  lapply(names(prevalence_results), function(var) {
    prevalence_results[[var]] %>%
      mutate(grouping_variable = var)
  }),
  .id = "group_id"
)
prevalence_table_Kleb <- prevalence_table %>%
  relocate(grouping_variable, .before = 1)

write.csv(prevalence_table_Kleb, "Kleb_prevalence_summary.csv", row.names = FALSE)

# View the table
print(prevalence_table_Kleb )


prevalence_results <- list()

for (var in group_vars_list) {
  result <- calculate_prevalence(
    unique_DATA,
    target_var = VITEK_MS_Results,
    value = "Escherichia coli",
    group_vars = c(var)
  )
  prevalence_results[[var]] <- result
}

# Combine all results into one table with a new column indicating the grouping variable
prevalence_table <- bind_rows(
  lapply(names(prevalence_results), function(var) {
    prevalence_results[[var]] %>%
      mutate(grouping_variable = var)
  }),
  .id = "group_id"
)
prevalence_table_ECO <- prevalence_table %>%
  relocate(grouping_variable, .before = 1)

write.csv(prevalence_table_ECO, "ECO_prevalence_summary.csv", row.names = FALSE)

# View the table
print(ECO_prevalence_table)
  