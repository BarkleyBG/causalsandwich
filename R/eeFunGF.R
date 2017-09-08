
# options(error='recover')

#' Estimating Function for IPTW
#'
#' This function is to be passed into geex::m_estimate
#'
#' @param outcome_model_obj The fitted model object (usually a glm).
#' @inheritParams estimateGF
#'
#' @export
eeFunGF <- function(
  data,
  outcome_model_obj,
  # outcome_var_name,
  treatment_var_name
){
  # ##' @param outcome_var_name The name of the column in the dataframe indicating outcome of interest

  # warning("need to use new geex grab_psifun")
  closureModel <- grab_psiFUN_glm(data=data,object = outcome_model_obj)
  data1 <- data0 <- data
  data1[[treatment_var_name]] <- 1
  data0[[treatment_var_name]] <- 0
  model_matrix1 <- stats::model.matrix(outcome_model_obj$formula,data = data1)
  model_matrix0 <- stats::model.matrix(outcome_model_obj$formula,data = data0)
  # Y         <- as.numeric(stats::model.frame(geex::grab_response_formula(object), data = data)[[1]])
  # treatment <- stats::model.response(stats::model.frame(trt_model_obj, data = data))
  # outcome <- data[[outcome_var_name]]


  ## Estimating function for IPTW (unstabilized)
  closureGF <- function(theta){
    num_params  <- length(theta)
    # num_model_params <- length(coef(trt_model_obj)) ##num_params-1
    linear_predictor1  <- model_matrix1 %*% theta[-num_params]
    prob_outcome1 <- stats::plogis(linear_predictor1)
    linear_predictor0  <- model_matrix0 %*% theta[-num_params]
    prob_outcome0 <- stats::plogis(linear_predictor0)

    (prob_outcome1-prob_outcome0) - theta[num_params]
  }

  closureStacked <- function(theta){

    c(
      closureModel(theta[-length(theta)]),## Returns vector of len=length(theta)-1
      closureGF(theta) ## Returns a single
    ) ## Returns a vector of len=length(theta)
  }
}


#' Calcualte G-formula-based Causal Effect Estimates
#'
#' @inheritParams eeFunGF
#'
calcGF <- function(
  data, outcome_model_obj,treatment_var_name
){
  data1 <- data0 <- data
  data1[[treatment_var_name]] <- 1
  data0[[treatment_var_name]] <- 0

  pred1 <- stats::predict(outcome_model_obj, newdata=data1, type="response")
  pred0 <- stats::predict(outcome_model_obj, newdata=data0, type="response")
  # mean(pred1)-mean(pred0)
  # mean(pred1-pred0)
  pred_delta <- pred1-pred0

}
#' Estimate ATE with IPTW
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Three-part formula: Outcome | Treatment ~ model_predictors. Will be coerced to object of type Formula.
#' @param data the dataframe. Will be coerced from "tbl_df" to data.frame.
#' @param model_method currently only supported "logistic" for logit-link binomial GLM.
#' @param treatment_var_name the name of the treatment variable. Must line up with formula.
#'
#' @export
estimateGF <- function(
  data, formula, model_method="logistic", treatment_var_name
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
    outcome_model_obj  <- stats::glm(
      formula = modeling_formula,
      data = data,
      family = stats::binomial()
    )

  } else { stop("only model_method='logistic' implemented") }

  est_GF <- calcGF(
    data= data,
    outcome_model_obj = outcome_model_obj,
    treatment_var_name = treatment_var_name
  )
  estimated_delta <- mean(est_GF)

  theta_hat <- c(stats::coef(outcome_model_obj), estimated_delta)

  estimates <- geex::m_estimate(
    estFUN = eeFunGF,
    data = data,
    compute_roots = FALSE,
    roots = theta_hat,
    outer_args = list(
      outcome_model_obj = outcome_model_obj,
      treatment_var_name = treatment_var_name
    )
  )
}



