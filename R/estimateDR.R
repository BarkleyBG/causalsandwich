
#' Estimate ATE with IPTW
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param outcome_formula Two part formula: Outcome ~ Covars + Trt
#' @param outcome_model_method currently only supported "logistic" for logit-link binomial GLM.
#' @param treatment_formula Two part formula: Trt ~ Covars
#' @param treatment_model_method currently only supported "logistic" for logit-link binomial GLM.
#' @param data the dataframe. Will be coerced from "tbl_df" to data.frame.
#' @param weight_type Currently only supports "unstabilized"
#' @param ... additional args
#'
#' @export
estimateDRIPTW <- function(
  data,
  outcome_formula,
  outcome_model_method="logistic",
  treatment_formula,
  treatment_model_method="logistic",
  weight_type="unstabilized",
  ...
){

  ## tibbles not allowed
  if ( "tbl_df" %in% class(data) ) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  outcome_formula <- Formula::as.Formula(outcome_formula)
  outcome_formula_terms <- stats::terms(outcome_formula)
  outcome_var_name <- row.names(attr(outcome_formula_terms, "factors"))[1]

  treatment_formula <- Formula::as.Formula(treatment_formula)
  treatment_formula_terms <- stats::terms(treatment_formula)
  treatment_var_name <- row.names(attr(treatment_formula_terms, "factors"))[1]


  outcome_vector <- data[[outcome_var_name]]
  treatment_vector <- data[[treatment_var_name]]


  outcome_model_fit <- fitOutcomeModel(
    data = data,
    formula = outcome_formula,
    model_method = outcome_model_method,
    treatment_var_name = treatment_var_name
  )


  treatment_model_fit <- fitTreatmentModel(
    data = data,
    formula = treatment_formula,
    model_method = treatment_model_method
  )

  args_list <- append(
    outcome_model_fit,
    treatment_model_fit
    )

  num_outcome_params <- length(args_list$outcome_param_ests)
  idx_outcome_params <- 1:num_outcome_params
  num_treatment_params <- length(args_list$treatment_param_ests)
  idx_treatment_params <- num_outcome_params+(1:num_treatment_params)

  args_list <- append(
    args_list,
    list(
      outcome = outcome_vector,
      outcome_var_name = outcome_var_name,
      num_outcome_params = num_outcome_params,
      idx_outcome_params = idx_outcome_params,

      treatment = treatment_vector,
      treatment_var_name = treatment_var_name,
      num_treatment_params = num_treatment_params,
      idx_treatment_params = idx_treatment_params
    )
  )
  if (weight_type == "unstabilized") {
    # calcFunDRIPTW <- match.fun('unstabilizedDRIPTW')
    calcFunDRIPTW <- unstabilizedDRIPTW
  } else {
    # if (weight_type!= "unstabilized") {
    stop("only unstabilized weights implemented")}

  args_list$calcFunDRIPTW <- calcFunDRIPTW


  IPTWs <- do.call(calcFunDRIPTW, args = args_list)
  estimated_delta <- mean(IPTWs)
  theta_hat <- c(
    args_list$outcome_param_ests,
    args_list$treatment_param_ests,
    estimated_delta
  )

  outer_args_list <- args_list[names(formals( eeFunDRIPTW ))[-1]]

  estimates <- geex::m_estimate(
    estFUN = eeFunDRIPTW,
    data = data,
    compute_roots = FALSE,
    roots = theta_hat,
    outer_args = outer_args_list,
    ...
  )
}


#' Function to Compute Horwitz-Thompson
#'
#' Computes estimates for each individual by IPTW (unstabilized) methods.
#'
#' @param outcome vector of outcome values
#' @param treatment vector of treatment values
#' @param prob_treated vector of propensity scores
#' @param prob_outcome1 vector of E(Y=1 | Trt=1, covars)
#' @param prob_outcome0 vector of E(Y=1 | Trt=0, covars)
#' @param ... dots
#'
#' @export
unstabilizedDRIPTW <- function(
  outcome, treatment,
  prob_treated,
  prob_outcome1,prob_outcome0,
  ...
){

  expected_y1 <- (
    (treatment*outcome) - (treatment - prob_treated)*prob_outcome1
  ) /  (prob_treated)

  expected_y0 <- (
    ((1-treatment)*outcome) + (treatment - prob_treated)*prob_outcome0
  ) /  (1-prob_treated)

  expected_y1-expected_y0
}


#' Estimating Function for DRIPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @param outcome_model_obj The fitted outcome model object (usually a glm).
#' @param outcome_var_name The name of the column in the dataframe indicating outcome of interest
#' @param num_outcome_params n_trt
#' @param idx_outcome_params vector of indices
#' @param trt_model_obj The fitted treatment model object (usually a glm).
#' @param treatment_var_name The name of the column in the dataframe indicating treatment
#' @param num_treatment_params n_trt
#' @param idx_treatment_params vector of indices
#' @inheritParams estimateDRIPTW
#' @param calcFunDRIPTW this is a function object specified in the weight_type argument
#' @param predictOutcome the function to predict outcome probs
#' @param predictTreatment the function to predict propensity score
#'
#' @export
eeFunDRIPTW <- function(
  data,
  outcome_model_obj, trt_model_obj,
  num_outcome_params, num_treatment_params,
  idx_outcome_params, idx_treatment_params,
  outcome_var_name,treatment_var_name,
  calcFunDRIPTW,
  predictOutcome,predictTreatment
){

  closureOutcomeModel <- grab_psiFUN_glm(data=data,object = outcome_model_obj)
  closureTreatmentModel <- grab_psiFUN_glm(data=data,object = trt_model_obj)

  outcome_model_matrices <- makeOutcomeModelMats(
    data = data,
    outcome_model_obj = outcome_model_obj,
    treatment_var_name = treatment_var_name
  )


  treatment_model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)
  treatment <- data[[treatment_var_name]]
  outcome <- data[[outcome_var_name]]


  ## Estimating function for IPTW (unstabilized)
  closureDRIPTW <- function(theta){

    num_params  <- length(theta)
    stopifnot((num_outcome_params + num_treatment_params + 1) == num_params)
    prob_outcomes <- lapply(
      outcome_model_matrices,
      predictOutcome,
      theta = theta[idx_outcome_params]
    )
    prob_treated <- predictTreatment(
      model_matrix = treatment_model_matrix,
      theta = theta[idx_treatment_params]
    )


    DRIPTWs <-
      calcFunDRIPTW(
        outcome = outcome,
        treatment = treatment,
        prob_treated = prob_treated,
        prob_outcome1 = prob_outcomes[[1]],
        prob_outcome0 = prob_outcomes[[2]]
      )

    DRIPTWs - theta[num_params]
  }

  closureStacked <- function(theta){

    c(
      closureOutcomeModel(theta[idx_outcome_params]),
      closureTreatmentModel(theta[idx_treatment_params]),
      closureDRIPTW(theta) ## Returns a singleton
    ) ## Returns a vector of len=length(theta)
  }
}
