
#' Estimate ATE with IPTW
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Three-part formula: Outcome | Treatment ~ Covars.
#' @param model_method Currently only supports "logistic" for logit-link binomial GLM.
#' @param weight_type Currently only supports "unstabilized"
#' @inheritParams estimateDRIPTW
#'
#' @export
estimateIPTW <- function(
  data, formula, model_method="logistic", weight_type="unstabilized", ...
){

  ## tibbles not allowed
  if ( "tbl_df" %in% class(data) ) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  formula <- Formula::as.Formula(formula)
  len_lhs_formula <- length(formula)[1]
  formula_terms <- stats::terms(formula)
  formula_factors <- attr(formula_terms, "factors")
  outcome_var_name <- row.names(formula_factors)[1]
  treatment_var_name <- row.names(formula_factors)[2]
  modeling_formula <- formula(
    stats::terms(formula, lhs = len_lhs_formula, rhs = -2)
  )

  treatment_vector <- data[[treatment_var_name]]
  outcome_vector <- data[[outcome_var_name]]


  treatment_model_fit <- fitTreatmentModel(
    data = data,
    formula = modeling_formula,
    model_method = model_method
  )
  args_list <- treatment_model_fit

  args_list <- append(
    args_list,
    list(
      outcome = outcome_vector,
      outcome_var_name = outcome_var_name,
      treatment = treatment_vector,
      treatment_var_name = treatment_var_name
    )
  )
  if (weight_type == "unstabilized") {
    # calcFunIPTW <- match.fun('unstabilizedIPTW')
    calcFunIPTW <- unstabilizedIPTW
  } else {
    # if (weight_type!= "unstabilized") {
    stop("only unstabilized weights implemented")}

  args_list$calcFunIPTW <- calcFunIPTW


  IPTWs <- do.call(calcFunIPTW, args = args_list)
  estimated_delta <- mean(IPTWs)
  theta_hat <- c(args_list$treatment_param_ests,estimated_delta )

  outer_args_list <- args_list[names(formals( eeFunIPTW ))[-1]]

  estimates <- geex::m_estimate(
    estFUN = eeFunIPTW,
    data = data,
    compute_roots = FALSE,
    roots = theta_hat,
    outer_args = outer_args_list,
    ...
  )
}


## #' Function to Compute Horwitz-Thompson
## #'
## #' Computes estimates for each individual by IPTW (unstabilized) methods.
## #'
## #' @param outcome vector of outcome values
## #' @param treatment vector of treatment values
## #' @param prob_treated vectof of propensity scores
## #' @inheritParams unstabilizedDRIPTW
## #'
## #' @export
unstabilizedIPTW <- function(
  outcome, treatment, prob_treated,...
){
  signed_trt_div_ps <- ifelse(treatment, 1/prob_treated, -1/(1-prob_treated))
  # ){ ##this takes care of PS=1 or PS = 0
  outcome*signed_trt_div_ps
}


### #' Estimating Function for IPTW
### #'
### #' This function is to be passed into geex::m_estimate
### #'
### #' @param trt_model_obj The fitted model object (usually a glm).
### #' @param outcome_var_name The name of the column in the dataframe indicating outcome of interest
### #' @inheritParams estimateIPTW
### #' @param calcFunIPTW this is a function object specified in the weight_type argument
### #' @inheritParams eeFunDRIPTW
### #'
### #' @export
eeFunIPTW <- function(
 data, trt_model_obj,
 outcome_var_name,
 treatment_var_name,
  calcFunIPTW,
  predictTreatment
){

  # warning("need to use new geex grab_psifun")
  closureModel <- grab_psiFUN_glm(data=data,object = trt_model_obj)

  model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)
  outcome <- data[[outcome_var_name]]
  treatment <- data[[treatment_var_name]]


  ## Estimating function for IPTW (unstabilized)
  closureIPTW <- function(theta){
    num_params  <- length(theta)

    prob_treated <- predictTreatment(
      model_matrix = model_matrix,
      theta = theta[-num_params]
    )

    IPTWs <-
      calcFunIPTW(
        outcome = outcome,
        treatment = treatment,
        prob_treated = prob_treated
      )
    IPTWs - theta[num_params]
  }

  closureStacked <- function(theta){

    c(
      closureModel(theta[-length(theta)]),## Returns vector of len=length(theta)-1
      closureIPTW(theta) ## Returns a single
    ) ## Returns a vector of len=length(theta)
  }
}


