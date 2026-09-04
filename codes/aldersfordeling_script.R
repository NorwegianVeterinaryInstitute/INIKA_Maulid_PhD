# Aldersfordeling i R
# Bytt ut df og age med dine egne datasett-/variabelnavn

# Grunnleggende beskrivelse/ Basic description
summary(Original_joined$Age_yrs)

# Antall manglende verdier/ Number of missing values
sum(is.na(Original_joined$Age_yrs))

# Histogram
hist(Original_joined$Age_yrs,
     breaks = 20,
     main = "Aldersfordeling",
     xlab = "Alder",
     col = "lightblue",
     border = "white")

# Tetthetskurve/ Density curve
plot(density(Original_joined $Age_yrs, na.rm = TRUE),
     main = "Density of Age",
     xlab = "Age")

# Kvantiler (nyttig for gruppering)/ Quantiles_ Useful for grouping
quantile(Original_joined$Age_yrs,
         probs = c(0, 0.25, 0.50, 0.75, 1),
         na.rm = TRUE)

# Lage kvartilgrupper/ Making Quartiles
Original_joined $age_grp <- cut(
  Original_joined$Age_yrs,
  breaks = quantile(Original_joined$Age_yrs,
                    probs = c(0, 0.25, 0.50, 0.75, 1),
                    na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

# Sjekk gruppestørrelser/ Check group sizes
table(Original_joined$age_grp)

# Hvis du bruker ggplot2/ Using ggplot
library(ggplot2)

ggplot(Original_joined, aes(x = Age_yrs)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  labs(title = "Agedistribution",
       x = "Age",
       y = "Number") +
  theme_minimal()
