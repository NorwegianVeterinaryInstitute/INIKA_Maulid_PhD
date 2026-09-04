# Model comparison helper for Model7 vs Model8

# 1. Compare nested models using likelihood ratio test
anova(model_Mv7, model_Mv6, test = "Chisq")

# 2. Calculate odds ratios and 95% confidence intervals
exp(cbind(
  OR = coef(model_Mv6),
  confint(model_Mv6)
))

# 3. Akaike weights (relative support for each model)
AIC7 <- 713.11
AIC6 <- 711.47

Delta7 <- AIC7 - min(AIC7, AIC6)
Delta6 <- AIC6 - min(AIC7, AIC6)

weight7 <- exp(-0.5 * Delta7) / (exp(-0.5 * Delta7) + exp(-0.5 * Delta7))
weight6 <- exp(-0.5 * Delta6) / (exp(-0.5 * Delta7) + exp(-0.5 * Delta6))

cat("Model 7 weight:", round(weight7, 3), "
")
cat("Model 6 weight:", round(weight6, 3), "
")

# Interpretation notes:
# - Delta AIC = 1.64, therefore both models have similar support.
# - Model 6 has the lower AIC and is therefore slightly preferred.
# - ORIGIN_OF_SAMPLE has p = 0.0566, which is borderline rather than clearly non-significant.
# - Removing ORIGIN_OF_SAMPLE changes the significance of Heard_AMR,
#   suggesting that ORIGIN_OF_SAMPLE contributes information to the model.
# - Unless there is a strong substantive reason to remove it,
#   Model 6 would usually be retained.
