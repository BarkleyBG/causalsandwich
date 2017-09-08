

#' Estimate ATE with IPTW
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Three-part formula: Outcome | Treatment ~ model_predictors. Will be coerced to object of type Formula.
#' @param data the dataframe. Will be coerced from "tbl_df" to data.frame.
#' @param model_method currently only supported "logistic" for logit-link binomial GLM.
#' @param weight_type Currently only supports "unstabilized"
#' @param ... additional args
#'
#' @export
estimateIPTW <- function(
  data, formula, model_method="logistic", weight_type="unstabilized", ...
){

  ## tibbles not allowed
  if ( "tbl_df" %in% class(data) ) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  # if (weight_type!= "unstabilized") {stop("only unstabilized weights implemented")}

  formula <- Formula::as.Formula(formula)
  len_lhs_formula <- length(formula)[1]
  formula_terms <- stats::terms(formula)
  outcome_var_name <- attr(formula_terms, "term.labels")[1]
  modeling_formula <- formula(
    stats::terms(formula, lhs = len_lhs_formula, rhs = -2)
  )

  ## Fit GLM model
  if (model_method == "logistic") {
    trt_model_obj  <- stats::glm(
      formula = modeling_formula,
      data = data,
      family = stats::binomial()
    )

     # browser()
    ##doing the following for estimation.
    # model_matrix <- stats::model.matrix(trt_model_obj, data = data)
    treatment <- stats::model.response(stats::model.frame(
      trt_model_obj$formula, data = data))
    prob_treated <- stats::predict(trt_model_obj, type = "response")
  } else { stop("only model_method='logistic' implemented") }



  if (weight_type == "unstabilized") {
    calcFunIPTW <- match.fun('unstabilizedIPTW')
  } else {
    # if (weight_type!= "unstabilized") {
    stop("only unstabilized weights implemented")}

  IPTWs <- calcFunIPTW(
    outcome = data[[outcome_var_name]],
    treatment = treatment,
    prob_treated = prob_treated
  )
  estimated_delta <- mean(IPTWs)
  theta_hat <- c(stats::coef(trt_model_obj), estimated_delta)

  warning("need to use new geex grab_psifun in eeFunIPTW")

  estimates <- geex::m_estimate(
    estFUN = eeFunIPTW,
    data = data,
    compute_roots = FALSE,
    roots = theta_hat,
    # root_control = setup_root_control(start = c(coef(trt_model_obj), 0)),
    # inner_args = list(
    #   weight_type = weight_type
    # ),
    outer_args = list(
      trt_model_obj = trt_model_obj,
      outcome_var_name = outcome_var_name,
      calcFunIPTW = calcFunIPTW

    )
  )
}

# print('got LHS formula')
# len_rhs_formula <- length(formula)[2]

# grouping_var_name <- attr( stats::terms(formula),
#  lhs = 0, rhs = len_rhs_formula), 'term.labels')
# model.matrix(modeling_formula,data=data)
# mf1 <- model.frame(formula=modeling_formula, data=data)
#
# # formula <- as.formula()
#  model.response(mf1)
#

#' Function to Compute Horwitz-Thompson
#'
#' Computes estimates for each individual by IPTW (unstabilized) methods.
#'
#' @param outcome vector of outcome values
#' @param treatment vector of treatment values
#' @param prob_treated vectof of propensity scores
#'
#' @export
unstabilizedIPTW <- function(
  outcome, treatment, prob_treated
){
  signed_trt_div_ps <- ifelse(treatment, 1/prob_treated, -1/(1-prob_treated))
  # ){ ##this takes care of PS=1 or PS = 0
  outcome*signed_trt_div_ps
  #   return( outcome * (-(!treatment)) )
  # }
  # outcome * (
  #   (treatment/prob_treated) - ((1-treatment)/(1-prob_treated))
  # )
}

# options(error='recover')

#' Estimating Function for IPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @param trt_model_obj The fitted model object (usually a glm).
#' @param outcome_var_name The name of the column in the dataframe indicating outcome of interest
#' @inheritParams estimateIPTW
#' @param calcFunIPTW this is a function object specified in the weight_type argument
#'
#' @export
eeFunIPTW <- function(
  data, trt_model_obj, outcome_var_name,calcFunIPTW
){

  # warning("need to use new geex grab_psifun")
  closureModel <- grab_psiFUN_glm(data=data,object = trt_model_obj)

  model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)
  # Y         <- as.numeric(stats::model.frame(geex::grab_response_formula(object), data = data)[[1]])
  treatment <- stats::model.response(stats::model.frame(trt_model_obj$formula, data = data))
  outcome <- data[[outcome_var_name]]


  ## Estimating function for IPTW (unstabilized)
  closureIPTW <- function(theta){
    num_params  <- length(theta)
    # num_model_params <- length(coef(trt_model_obj)) ##num_params-1
    linear_predictor  <- model_matrix %*% theta[-num_params]
    prob_treated <- stats::plogis(linear_predictor)

    # if (weight_type == "unstabilized") {
    ## Y*A/PS - Y*(1-A)/(1-PS)
    IPTWs <-
      calcFunIPTW(
        # unstabilizedIPTW(
        outcome = outcome,
        treatment = treatment,
        prob_treated = prob_treated
      )
    # } else {
    #   stop("weight_type == 'unstabilized' only implemented. Not stabilized yet")
    # }

    IPTWs - theta[num_params]
  }

  closureStacked <- function(theta){

    c(
      closureModel(theta[-length(theta)]),## Returns vector of len=length(theta)-1
      closureIPTW(theta) ## Returns a single
    ) ## Returns a vector of len=length(theta)
  }
}


