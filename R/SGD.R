#' @rdname SGD
#' @title Stochastic gradient descent
#' @description Run a stochastic gradient descent using unbiased estimator of the gradient of the log-likelihood
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta_initial an initial vector of parameters 
#' @param observations a matrix of observations of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled2_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_coupled_residuals}}
#' @param coupled4_resampling a 4-way coupled resampling scheme, such as \code{\link{coupled4_maximal_coupled_residuals}}
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param minimum_level coarsest discretization level 
#' @param maximum_level finest discretization level 
#' @param level_distribution list containing mass_function and tail_function that specify the distribution of levels, 
#' e.g. by calling \code{\link{compute_level_distribution}} 
#' @param learning_rate stepsize of the SGD algorithm
#' @param stopping_threshold criterion to terminate iterations
#' @param max_iterations maximum number of SGD iterations
#' @return a list with objects such as 
#' theta parameters at the last SGD iteration
#' trajectory parameters across the SGD iterations
#' @export
SGD <- function(model, theta_initial, observations, nparticles, resampling_threshold = 1, coupled2_resampling, coupled4_resampling, 
                k = 0, m = 1, minimum_level, maximum_level, level_distribution, 
                learning_rate = 1e-3, stopping_threshold = 1e-4, max_iterations = 1e6){
  
  # initialize 
  iter <- 1
  continue <- TRUE
  theta <- theta_initial
  
  # Trace of each parameter through iterations
  trajectory <- matrix(0, model$theta_dimension, max_iterations)
  trajectory[, 1] <- theta
  
  # iterate until convergence or maximum number of iterations
  while (iter < max_iterations && continue){
    iter <- iter + 1
    
    # compute unbiased estimator of gradient 
    gradient <- coupled_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                            k = k, m = m, minimum_level, maximum_level, level_distribution)$unbiasedestimator
    
    # parameter update
    theta_new <- theta + learning_rate * gradient
    trajectory[, iter] <- theta_new
    
    # check for convergence
    difference <- sqrt(sum((theta_new - theta)^2))
    continue <- (difference > stopping_threshold)
    
    # update theta values
    theta <- theta_new
    
    # print output
    cat("Iteration:", iter, "\n", "Parameters:", theta, "\n")
  }
  
  # truncate trace
  trajectory <- trajectory[, 1:iter]
  
  return(list(theta = theta, trajectory = trajectory))
  
}
  