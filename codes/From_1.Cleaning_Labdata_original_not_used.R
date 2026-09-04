





# we probably don't need thos below
## Filter No growth from Isolate column 
Cleaned_Labdata_NoGrowth <- 
  Cleaned_Labdata_Original %>%
  filter(Isolate == "No growth")

## Save the file N=413
# Saving the file as tsv to be opened in excel format
write_tsv(Cleaned_Labdata_NoGrowth , "data/CLEANED_DATA/Cleaned_Labdata_NoGrowth-2025-10-09.tsv")
# Save as rds to be used further work in R
saveRDS(Cleaned_Labdata_NoGrowth , "data/CLEANED_DATA/Cleaned_Labdata_NoGrowth-2025-10-09.rds")

