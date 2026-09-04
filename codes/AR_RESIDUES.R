## Cleaning data sets for Antibiotics residues
# Load libraries
library(tidyverse)
library(readxl)
library(dplyr)
library(broom)
library(pROC)
library(ResourceSelection)
library(flextable) 
library(officer)

# Importing the AR results for cleaning
AR_Results <- read_excel("data/Updated_stool_muscles_ARs_data_edited.xlsx", 
                                                    sheet = "Conc in 2g (stool), 15g(muscle)")
View(AR_Results)
names(AR_Results)
# Selecting relevant column and lines

AR_Results <- AR_Results %>% 
  select(-c(1, 3)) %>%          
  slice(-c(1:7,9:33))

# Set the column names to be the content of the first row
colnames(AR_Results) <- AR_Results[1, ]

# Remove that first row so it isn't duplicated in the data
AR_Results <- AR_Results %>% 
  slice(-1)

# Rename the first column and Create the new column by extracting only the numbers
AR_Results <- AR_Results %>% 
  rename(SUA_SAMPLE_ID = `Sample codes`) %>% 
  mutate(
    AR_RESIDUES_SAMPLE = as.numeric(str_extract(SUA_SAMPLE_ID, "\\d+")),
    .after = SUA_SAMPLE_ID
  )

# Import another dataset for joining
joined_data <- read_csv("data/CLEANED_DATA/joined_data_16.3.26.csv")

# Join the joined_data with AR_Results by AR_RESIDUES_SAMPLE
combined_data <- joined_data %>% 
  left_join(AR_Results, by = c("AR_RESIDUES_SAMPLE" = "AR_RESIDUES_SAMPLE"))
names(combined_data)
combined_data <- combined_data %>% 
  select(INIKA_ID.x.x,Age_yrs,GENDER,SEASON.x,REGION.x,DISTRICT.x,ORIGIN_OF_SAMPLE.x,
         AR_RESIDUES_SAMPLE,"COLONY MORPHOLOGY ON C3GR","COLONY MORPHOLOGY ON CARBA",
         "COLONY MORPHOLOGY ON XLD","COLONY MORPHOLOGY ON BGA","TSI Media_Slope" ,
         "TSI Media_Butt","TSI Media_Gas",TSIMedia_H2S, "SIM Media_H2S" ,"SIM Media Indole" ,
         "SIM Media Motility",CITRATE,UREASE,Isolate,Isolate_ID,TVLA_ID,
         VITEK_MS_Results,ESC,PROTOCOL,AMX_ED10,AZM_ED15,CRO_ED30,CIP_ED5,DOX_ED30,
         FLR_ED30,GEN_ED10,MEM_ED10,OXY_ED30,POL_ED300,SXT_ED1_2,CTX_ED5,
         CTC_ED30,CRO_ED5,ESBL_Selection,ESCR_ECO_presumptive,ESBL_Presumptivefinal,
         Isolate_NVI,SUA_SAMPLE_ID,"Sufapyridine","Sulfamethoxazole","Trimethoprim","Ciprofloxacin",
         "Florfenicol","Amoxicilin","Meropenen","Tetracycline","Azithromycin","Tylosin","Doxycycline",
         "Gentamycin","Ceftriaxone","Polymxin B")

#  Define variables and prepare the data 
antibiotic_list <- c(
  "Sufapyridine", "Sulfamethoxazole", "Trimethoprim", "Ciprofloxacin",
  "Florfenicol", "Amoxicilin", "Meropenen", "Tetracycline", "Azithromycin",
  "Tylosin", "Doxycycline", "Gentamycin", "Ceftriaxone", "Polymxin B"
)

cleaned_data <- combined_data %>%
  mutate(across(all_of(antibiotic_list), ~ as.numeric(as.character(.x))))

unique_samples <- cleaned_data %>%
  group_by(AR_RESIDUES_SAMPLE, REGION.x, ORIGIN_OF_SAMPLE.x, SEASON.x) %>%
  summarise(across(all_of(antibiotic_list), ~ max(.x, na.rm = TRUE)), .groups = "drop")

long_samples <- unique_samples %>%
  pivot_longer(cols = all_of(antibiotic_list), names_to = "Antibiotic", values_to = "Concentration") %>%
  filter(!is.na(Concentration), Concentration > 0)

# Build the publication-grade histogram 
histogram_plot <- ggplot(long_samples, aes(x = Concentration, fill = SEASON.x)) +
  geom_histogram(bins = 15, position = "identity", alpha = 0.6, color = "#ffffff", linewidth = 0.2) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_grid(REGION.x ~ ORIGIN_OF_SAMPLE.x, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = expression(paste("Antibiotic Concentration (", mu, "g/g - Log Scale)")),
    y = "Count of Detections",
    fill = "Season"
  ) +
  theme_bw(base_size = 11) + 
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

#  Save the plot locally as a temporary high-res image 
dir.create("Results", showWarnings = FALSE)
temp_img_path <- "Results/temp_histogram.png"

ggsave(
  filename = temp_img_path,
  plot = histogram_plot,
  width = 6.5,  # Fits perfectly inside standard 1-inch Word margins
  height = 4.5,
  dpi = 300
)

#  Use 'officer' to inject the chart into a Word Document 
doc <- read_docx() %>%
  # Add a title to the document
  body_add_par("Antibiotic Residues Analysis Results", style = "heading 1") %>%
  body_add_par("The faceted histogram below shows the distribution of antibiotic concentrations across different regions, sample origins, and seasons.", style = "Normal") %>%
  body_add_break() %>%
  # Insert the generated high-resolution plot image
  body_add_img(src = temp_img_path, width = 6.5, height = 4.5) %>%
  body_add_par("Figure 1: Stratified histogram of antibiotic concentrations (Log Scale).", style = "Image Caption")

#  Write the Word Document out to file 
print(doc, target = "Results/Antibiotic_Residues_Histogram_Report.docx")

# Clean up: delete the temporary image file if you only want the Word Doc
file.remove(temp_img_path)



