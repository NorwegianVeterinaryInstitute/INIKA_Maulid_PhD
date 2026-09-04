# Load required packages
library(dplyr)
library(ggplot2)
library(binom)



# Calculate prevalence and 95% CI using Wilson method
summary_df <- AMR_Data_Classified %>%
  group_by(antibiotic_name) %>%
  summarise(
    n_total = n(),
    n_R = sum(phenotype == "R"),
    prevalence = n_R / n_total
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(binom.confint(n_R, n_total, method = "wilson")),
    CI_lower = ci[[1]]$lower,
    CI_upper = ci[[1]]$upper
  ) %>%
  select(-ci)

print(summary_df)

# Plot: Bar chart with error bars
ggplot(summary_df, aes(x = antibiotic_name, y = prevalence)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  labs(title = "Prevalence of Resistance by Antibiotic",
       x = "Antibiotic",
       y = "Prevalence (R)") +
  theme_minimal()