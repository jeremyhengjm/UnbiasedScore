#' @rdname unbiased_gradient
#' @title Unbiased estimator of the gradient of the log-likelihood at a discretization level
#' @description Estimates the expectation of a functional with respect to the smoothing distribution at a discretization level
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param nof_niterations number of iterations
#' @return the MCMC estimate
#' @export
mcmc_estimator <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1, nof_iterations){
  
  # initialize chains
  chain_state <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
  
  # initialize estimators computation
  mcmcestimator <- model$functional(theta, discretization, chain_state, observations)
  
  for (iter in 2:nof_iterations){
    chain_state <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, chain_state)$new_trajectory
    mcmcestimator <- mcmcestimator + model$functional(theta, discretization, chain_state, observations)
  }
  
  # compute mcmc gradient estimator
  mcmcestimator <- mcmcestimator / nof_iterations
  
  return(mcmcestimator)
}
