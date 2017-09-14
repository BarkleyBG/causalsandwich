## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github("BarkleyBG/causalsandwich")

## ------------------------------------------------------------------------
library(causalsandwich)

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
formula <- BinaryOutcome | BinaryTrt ~ Covar1 * Covar2

IPTW <- estimateIPTW(
  data = data, 
  formula = formula 
) 

## ------------------------------------------------------------------------
(ests <- IPTW@estimates)
(vcov <- IPTW@vcov) 

## ------------------------------------------------------------------------
param_num <- 5
ests[param_num]
ests[param_num] + stats::qnorm(c(0.025,0.975))*vcov[param_num,param_num]

## ------------------------------------------------------------------------
outcome_regression_formula <- BinaryOutcome ~ BinaryTrt + Covar1 * Covar2

GF <- estimateGF(
  data = data,
  treatment_var_name = "BinaryTrt",
  formula = outcome_regression_formula, 
  model_method = "logistic" 
) 
 
(ests <- GF@estimates)
(vcov <- GF@vcov) 

## ------------------------------------------------------------------------
param_num <- 6
ests[param_num]
ests[param_num] + stats::qnorm(c(0.025,0.975))*vcov[param_num,param_num]

## ------------------------------------------------------------------------
outcome_regression_formula <- BinaryOutcome ~ BinaryTrt + Covar1 * Covar2
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
param_num <- 10
ests[param_num]
ests[param_num] + stats::qnorm(c(0.025,0.975))*vcov[param_num,param_num]

