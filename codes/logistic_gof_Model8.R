# Goodness-of-fit checks for logistic regression (Model_6)

# -------------------------------------------------
# 1. Hosmer-Lemeshow goodness-of-fit test
# -------------------------------------------------
# Install once if needed:
# install.packages("ResourceSelection")

library(ResourceSelection)

# Predicted probabilities model_Mv6
pred <- fitted(model_Mv6)

# Hosmer-Lemeshow test (10 groups)
hoslem.test(data_model$Case, pred, g = 10)

# Interpretation:
# p > 0.05 : no evidence of poor fit
# p < 0.05 : model may not fit well


# -------------------------------------------------
# 2. Pseudo R-squared
# -------------------------------------------------
# install.packages("pscl")

library(pscl)

pR2(model_Mv6)


# -------------------------------------------------
# 3. Classification table
# -------------------------------------------------

pred_class <- ifelse(pred > 0.5, 1, 0)

table(
  Observed = data_model$Case,
  Predicted = pred_class
)


# -------------------------------------------------
# 4. ROC curve and AUC
# -------------------------------------------------
# install.packages("pROC")

library(pROC)

roc_obj <- roc(data_model$Case, pred)

plot(roc_obj,
     main = "ROC curve for Model_6")

auc(roc_obj)

# Interpretation:
# 0.5  = no discrimination
# 0.7+ = acceptable
# 0.8+ = good
# 0.9+ = excellent


# -------------------------------------------------
# 5. Calibration plot
# -------------------------------------------------

decile <- cut(pred,
              breaks = quantile(pred,
                                probs = seq(0,1,0.1),
                                na.rm = TRUE),
              include.lowest = TRUE)

obs <- tapply(data_model$Case, decile, mean)

exp <- tapply(pred, decile, mean)

plot(exp, obs,
     xlab = "Predicted probability",
     ylab = "Observed proportion",
     main = "Calibration plot")

abline(0,1,col="red",lty=2)


# -------------------------------------------------
# 6. Influence diagnostics
# -------------------------------------------------

plot(cooks.distance(model_Mv6),
     type = "h",
     main = "Cook's distance")


# -------------------------------------------------
# 7. Odds ratios with 95% CI
# -------------------------------------------------

exp(cbind(
  OR = coef(model_Mv6),
  confint(model_Mv6)
))
