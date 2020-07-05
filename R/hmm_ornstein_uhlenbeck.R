#' @rdname hmm_ornstein_uhlenbeck
#' @title Construct hidden markov model obtained by discretizing Ornstein-Uhlenbeck process
#' @param nobservations number of observations or time interval of interest
#' @description This function returns a list with objects such as
#' * xdimension is the dimension of the latent process
#' * ydimension is the dimension of the observation process
#' * theta_dimension is the dimension of the parameter space
#' * construct_discretization outputs a list containing stepsize, nsteps, statelength and obstimes
#' * construct_successive_discretization outputs lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' * sigma is the diffusion coefficient of the process
#' * rinit to sample from the initial distribution
#' * rtransition to sample from the Markov transition
#' * dtransition to evaluate the transition density
#' * dmeasurement to evaluate the measurement density
#' * functional is the smoothing function to compute gradient of log-likelihood 
#' @return A list 
#' @export
hmm_ornstein_uhlenbeck <- function(times){
  # model dimensions
  xdimension <- 1
  ydimension <- 1
  
  # parameters of the model
  theta_dimension <- 3
  
  # time intervals
  nobservations <- length(times) # nobservations at unit times
  nintervals <- nobservations - 1 
  
  # construct discretization
  construct_discretization <- function(level){
    stepsize <- 2^(-level)
    nsteps <- nintervals * 2^level
    all_stepsizes <- rep(stepsize, nsteps)
    statelength <- nsteps + 1
    obstimes <- rep(FALSE, statelength)
    obs_index <- seq(1, statelength, by = 2^level)
    obstimes[obs_index] <- TRUE
    
    return(list(stepsize = all_stepsizes, nsteps = nsteps, 
                statelength = statelength, obstimes = obstimes))
  }
  
  construct_successive_discretization <- function(level_fine){
    # construct fine discretization
    fine <- construct_discretization(level_fine)
    
    # construct coarse discretization
    level_coarse <- level_fine - 1
    coarse <- construct_discretization(level_coarse)
    
    # index coarse times
    coarsetimes <- rep(FALSE, fine$statelength)
    coarse_index <- seq(1, fine$statelength, by = 2)
    coarsetimes[coarse_index] <- TRUE
    
    return(list(fine = fine, coarse = coarse, coarsetimes = coarsetimes))
  }
  
  # drift 
  drift <- function(theta, x) theta[1] * (theta[2] - x) # drift
  jacobian_drift <- function(theta, x) c(theta[2] - x, theta[1], 0)
  
  # diffusivity 
  sigma <- 1
  Sigma <- sigma^2
  Omega <- sigma^(-2)
  
  # sample from initial distribution
  rinit <- function(nparticles){
    return(matrix(rnorm(nparticles, mean = 0, sd = 1), nrow = 1)) # returns 1 x N
  }
  
  # sample from Markov transition kernel
  rtransition <- function(theta, stepsize, xparticles, rand){ 
    # theta is a vector of size 3
    # xparticles is a vector of size N
    # should take xparticles of size 1 x N and output new particles of size 1 x N
    # rand is a vector of size N following a standard normal distribution
    # output new particles in a vector of size N
    
    return(xparticles + stepsize * drift(theta, xparticles) + sigma * sqrt(stepsize) *  rand) 
  }
  
  # evaluate Markov transition density 
  dtransition <- function(theta, stepsize, xparticles, next_xparticles){
    # theta is a vector of size 3
    # xparticles and next_xparticles are vectors of size N
    
    particle_mean <- xparticles + stepsize * drift(theta, xparticles) 
    return(dnorm(next_xparticles, mean = particle_mean, sd = sigma * sqrt(stepsize), log = TRUE))
  }
  
  # evaluate observation density
  dmeasurement <- function(theta, xparticles, observation){
    # xparticles is a vector of size N
    # theta is a vector of size 3
    # observation is a number
    return(dnorm(observation, mean = xparticles, sd = sqrt(theta[3]), log = TRUE))
  }
  
  # evaluate gradient of log-observation density
  gradient_dmeasurement <- function(theta, xstate, observation){
    # theta is a vector of size 3
    # xstate is a number 
    # observation is a number
    partial_derivative_variance <- - 1 / (2 * theta[3]) + (observation - xstate)^2 / (2 * theta[3]^2)
    return(c(0, 0, partial_derivative_variance))
  }
  
  functional <- function(theta, discretization, xtrajectory, observations){
    # theta is a vector of size 3
    # discretization is a list created by construct_discretization
    # xtrajectory is a vector of length statelength 
    # observations is a matrix of size nobservations x 1
    
    # discretization objects
    stepsize <- discretization$stepsize # vector of length nsteps
    nsteps <- discretization$nsteps 
    statelength <- discretization$statelength  
    obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1
    
    # initialize 
    output <- rep(0, theta_dimension)
    index_obs <- 1
    xstate <- xtrajectory[1]
    output <- output + gradient_dmeasurement(theta, xstate, observations[index_obs, ])
    
    # loop over time 
    for (k in 1:nsteps){
      jacobian <- jacobian_drift(theta, xtrajectory[k])
      output <- output - stepsize[k] * Omega * drift(theta, xtrajectory[k]) * jacobian  # JHo: coeff 0.5 removed
      output <- output + Omega * (xtrajectory[k+1] - xtrajectory[k]) * jacobian
      
      # observation time
      if (obstimes[k+1]){
        index_obs <- index_obs + 1
        xstate <- xtrajectory[k+1]
        output <- output + gradient_dmeasurement(theta, xstate, observations[index_obs, ])
      }
    }
    
    return(output)
  }
  
  model <- list(xdimension = xdimension,
                ydimension = ydimension,
                theta_dimension = theta_dimension,
                construct_discretization = construct_discretization,
                construct_successive_discretization = construct_successive_discretization,
                sigma = sigma,
                rinit = rinit, 
                rtransition = rtransition, 
                dtransition = dtransition, 
                dmeasurement = dmeasurement,
                functional = functional)
  return(model)
  
}