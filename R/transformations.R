#' @rdname log_transform
#' @title Logarithmic transformation 
#' @description Computes logarithmic transformation for parameters with positivity constraints.
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters 
#'@return a vector of transformed parameters 
#'@export
log_transform <- function(model, theta){
  log_theta <- theta
  log_theta[model$theta_positivity] <- log(theta[model$theta_positivity])
  return(log_theta)
}

#' @rdname exp_transform
#' @title Exponential transformation 
#' @description Computes exponential transformation for parameters with positivity constraints.
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param log_theta a vector of parameters 
#'@return a vector of transformed parameters 
#'@export
exp_transform <- function(model, log_theta){
  theta <- log_theta
  theta[model$theta_positivity] <- exp(log_theta[model$theta_positivity])
  return(theta)
}
