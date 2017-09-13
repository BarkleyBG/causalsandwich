

#' Estimate ATE with G-Formula
#'
#' Fits a parametric model and estimates ATE via IPTW with Wald-type confidence
#' intervals from the empirical sandwich standard error estimates.
#'
#' @param formula Two-part formula: Outcome ~ Treatment + Covariates.
#' @param model_method Currently only supported "logistic" for logit-link binomial GLM.
#' @param treatment_var_name The name of the treatment variable's column in `data`
#' @inheritParams estimateDRIPTW
#'
#' @export
estimateGF <- function(
  data, formula, model_method="logistic", treatment_var_name,...
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
  modeling_formula <- formula(
    stats::terms(formula, lhs = len_lhs_formula, rhs = -2)
  )



  outcome_vector <- data[[outcome_var_name]]
  treatment_vector <- data[[treatment_var_name]]

  outcome_model_fit <- fitOutcomeModel(
    data = data,
    formula =  modeling_formula,
    model_method =  model_method,
    treatment_var_name = treatment_var_name
  )
  args_list <- outcome_model_fit

  args_list <- append(
    args_list,
    list(
      outcome = outcome_vector,
      outcome_var_name = outcome_var_name,
      treatment = treatment_vector,
      treatment_var_name = treatment_var_name
    )
  )

  estimated_delta <- mean(args_list$prob_outcome1-args_list$prob_outcome0)
  theta_hat <- c(args_list$outcome_param_ests,estimated_delta)

  outer_args_list <- args_list[names(formals( eeFunGF ))[-1]]

  estimates <- geex::m_estimate(
    estFUN = eeFunGF,
    data = data,
    compute_roots = FALSE,
    roots = theta_hat,
    outer_args = outer_args_list,
    ...
  )
}



## #' Estimating Function for IPTW
## #'
## #' This function is to be passed into geex::m_estimate
## #'
## #' @param outcome_model_obj The fitted model object (usually a glm).
## #' @inheritParams estimateGF
## #' @inheritParams eeFunDRIPTW
## #'
## #' @export
eeFunGF <- function(
  data,
  outcome_model_obj,
  predictOutcome,
  treatment_var_name
){

  closureModel <- grab_psiFUN_glm(data=data,object = outcome_model_obj)

  outcome_model_matrices <- makeOutcomeModelMats(
    data = data,
    outcome_model_obj = outcome_model_obj,
    treatment_var_name = treatment_var_name
  )
  model_matrix1 <- outcome_model_matrices[[1]]
  model_matrix0 <- outcome_model_matrices[[2]]


  closureGF <- function(theta){
    num_params  <- length(theta)

    prob_outcomes <- lapply(
      outcome_model_matrices,
      predictOutcome,
      theta = theta[-num_params]
    )

    (prob_outcomes[[1]] - prob_outcomes[[2]]) - theta[num_params]
  }

  closureStacked <- function(theta){

    c(
      closureModel(theta[-length(theta)]),## Returns vector of len=length(theta)-1
      closureGF(theta) ## Returns a single
    ) ## Returns a vector of len=length(theta)
  }
}


