## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github("BarkleyBG/causalgeex")

## ------------------------------------------------------------------------
library(causalgeex)

n <- 200
data <- data.frame(
  Covar1 = rnorm(n),
  Covar2 = rnorm(n)
)
trtprobs <- plogis(0.2  + 0.5*data$Covar1 + 0.1*data$Covar2 + 0.2*data$Covar1 * data$Covar2)
data$BinaryTrt <- rbinom(n, 1, trtprobs)

outprobs <- plogis(0.5  + 0.2*data$Covar1 - 0.2*data$Covar2 + 0.2*data$Covar1 * data$Covar2 -0.2*data$BinaryTrt)
data$BinaryOutcome <- rbinom(n, 1, outprobs)

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
outcome_regression_formula <- BinaryOutcome ~ BinaryTrt + Covar1*Covar2
treatment_model_formula <- BinaryTrt ~ Covar1 * Covar2

DRIPTW <- estimateDRIPTW(
  data = data,
  outcome_formula = outcome_regression_formula,
  treatment_formula = treatment_model_formula,
  outcome_model_method = "logistic",
  treatment_model_method = "logistic",
  deriv_control = geex::setup_deriv_control(method="simple")
)

(ests <- DRIPTW@estimates)
(vcov <- DRIPTW@vcov)

## ------------------------------------------------------------------------
ests[10]
ests[10] + stats::qnorm(c(0.025,0.975))*vcov[10,10]

