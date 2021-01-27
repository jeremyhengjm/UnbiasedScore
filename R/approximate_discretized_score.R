#' @rdname approximate_discretized_score
#' @title MCMC estimator of score function at a discretization level
#' @description Computes MCMC estimator of the time-discretized score function
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param initialization choice of distribution to initialize chains, such as \code{dynamics} or the default \code{particlefilter} 
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param max_iterations number of MCMC iterations
#' @return a list with objects such as 
#' mcmcestimator is the MCMC estimator of the discretized score 
#' cost is the cost of the algorithm 
#' elapsedtime is the time taken by the algorithm
#' @export
approximate_discretized_score <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1, 
                                          initialization = "particlefilter", algorithm = "CPF", max_iterations){
  
  # initialize chains
  if (initialization == "dynamics"){
    chain_state <- simulate_SDE(model, theta, discretization)$new_trajectory
  }
  if (initialization == "particlefilter"){
    chain_state <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
  }
  
  # pre-allocate
  theta_dimension <- model$theta_dimension
  mcmcestimator <- matrix(0, nrow = max_iterations + 1, ncol = theta_dimension)
  cost <- rep(0, max_iterations + 1)
  elapsedtime <- rep(0, max_iterations + 1)
  
  # initialize MCMC estimator
  functional_cumsum <- model$functional(theta, discretization, chain_state, observations)
  mcmcestimator[1, ] <- functional_cumsum
  
  # MCMC iterations
  for (i in 1:max_iterations){
    # start timer
    tic() 
    
    # sample from MCMC kernel
    chain_state <- kernel(model, theta, discretization, observations, nparticles, resampling_threshold, chain_state, algorithm)$new_trajectory
    
    # evaluate functional
    functional_cumsum <- functional_cumsum + model$functional(theta, discretization, chain_state, observations)
    
    # update MCMC estimator
    mcmcestimator[i+1, ] <- functional_cumsum / (i+1)
    
    # update cost
    cost[i+1] <- cost[i] + nparticles * discretization$nsteps 
    
    # end timer and update elapsed time
    timer <- toc(quiet = TRUE)
    elapsedtime[i+1] <- elapsedtime[i] + (timer$toc - timer$tic)
    
  }
  
  return(list(mcmcestimator = mcmcestimator, cost = cost, elapsedtime = elapsedtime))
}