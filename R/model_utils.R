
# #' Prep Outcome Model Matrices for Prediction
# #'
# #' @inheritParams eeFunDRIPTW
# #'
# #' @export
makeOutcomeModelMats <- function(
  data,
  outcome_model_obj,
  treatment_var_name
){
  lapply(1:0, function(trt_val){
    data_trt_val <- data
    data_trt_val[[treatment_var_name]] <- trt_val
    model_matrix_trt_val <-
      stats::model.matrix(outcome_model_obj$formula,data = data_trt_val)
  })
}


## #' Prep Outcome Model Matrices for Prediction
## #'
## #' @inheritParams eeFunDRIPTW
## #' @inheritParams estimateGF
## #'
## #' @export
fitOutcomeModel <- function(
  data,
  formula,
  model_method,
  treatment_var_name
  # outcome_var_name,

){
  ## Fit GLM model
  if (model_method == "logistic") {
    outcome_model_obj  <- stats::glm(
      formula = formula,
      data = data,
      family = stats::binomial()
    )

    # outcome_vector <- stats::model.response(stats::model.frame(
    #   outcome_model_obj$formula, data = data))
    # stopifnot(all(outcome_vector == data[[outcome_var_name]]))

    outcome_model_matrices <- makeOutcomeModelMats(
      data = data,
      outcome_model_obj = outcome_model_obj,
      treatment_var_name = treatment_var_name
    )
    outcome_param_ests <- outcome_model_obj$coefficients
    predictOutcome <- function(model_matrix,theta){
      stats::plogis(model_matrix %*% theta)
    }
    # prob_outcome1 <- stats::plogis(outcome_model_matrices[[1]] %*% outcome_param_ests)
    # prob_outcome0 <- stats::plogis(outcome_model_matrices[[2]] %*% outcome_param_ests)
    prob_outcomes <- lapply(
      outcome_model_matrices,
      predictOutcome,
      theta = outcome_param_ests
    )

    out <- list(
      outcome_model_obj = outcome_model_obj,
      # outcome_vector = outcome_vector,

      outcome_param_ests = outcome_param_ests,
      prob_outcome1 = prob_outcomes[[1]],
      prob_outcome0 = prob_outcomes[[2]],
      predictOutcome = predictOutcome

    )
  } else { stop("only model_method='logistic' implemented") }

  out
}

## #' Prep Outcome Model Matrices for Prediction
## #'
## #' @inheritParams fitOutcomeModel
## #'
## #' @export
fitTreatmentModel <- function(
  data,
  formula,
  model_method
){
  if (model_method == "logistic") {
    trt_model_obj  <- stats::glm(
      formula = formula,
      data = data,
      family = stats::binomial()
    )

    # treatment_vector <- stats::model.response(stats::model.frame(
    #   trt_model_obj$formula, data = data))
    # stopifnot(all(treatment_vector == data[[treatment_var_name]]))


    treatment_model_matrix <- stats::model.matrix(trt_model_obj$formula,data = data)
    treatment_param_ests <-  trt_model_obj$coefficients
    # prob_treated <- stats::plogis(treatment_model_matrix %*% treatment_param_ests)
    predictTreatment <- function(model_matrix,theta){
      stats::plogis(model_matrix %*% theta)
    }
    prob_treated <- predictTreatment(
      model_matrix = treatment_model_matrix,
      theta = treatment_param_ests
    )
    out <- list(
      trt_model_obj = trt_model_obj,
      treatment_param_ests  = treatment_param_ests,
      prob_treated = prob_treated,
      predictTreatment = predictTreatment
    )
  } else { stop("only model_method='logistic' implemented") }
  out
}
