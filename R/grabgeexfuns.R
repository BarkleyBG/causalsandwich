
#' Temporary Function Until geex>1.0.3
#'
#' @param object stats::glm model object
#' @param data data
#' @param weights 1
#' @param ... more stuffs
#'
#' @export
grab_psiFUN_glm <- function(object, data, weights = 1, ...){

  X  <- stats::model.matrix(object$formula, data = data)
  Y  <- as.numeric(stats::model.frame(geex::grab_response_formula(object), data = data)[[1]])
  n  <- length(Y)
  p  <- length(stats::coef(object))
  phi    <- as.numeric(summary(object)$dispersion[1])
  W      <- weights
  family <- object$family$family
  link   <- object$family$link
  invlnk <- object$family$linkinv
  family_link <- paste(family, link, sep = '_')

  stopifnot(length(W) == 1 | length(W) == n)
  if(length(W) == 1){
    W <- rep(W, n)
  }

  function(theta){
    lp <- X %*% theta # linear predictor
    f  <- as.numeric(invlnk(lp))  # fitted values
    r  <- Y - f       # residuals

    ### TODO: this is cludgy and needs to be reworked to be more general
    if(family_link == 'gaussian_identity'){
      D <- X
      V <- phi * diag(1, nrow = n, ncol = n)
    } else if(family_link == 'binomial_logit'){
      D <- apply(X, 2, function(x) x * exp(lp)/((1+exp(lp))^2) )
      if (n==1){D <- t(D)} ##takes care of base case
      V <- phi * diag(f * (1 - f), ncol = length(f) )/length(f)
    }

    t(D) %*% solve(V) %*% diag(W, nrow = n, ncol = n) %*% (r)
  }
}
