#' @rdname hmm_logistic_diffusion_full
#' @title Construct hidden Markov model obtained by discretizing Lamperti transformed logistic diffusion process
#' @description Define a hidden Markov model obtained by discretizing Lamperti transformed logistic diffusion process. 
#' @param times vector specifying observation times
#' @return a list with objects such as: 
#' \code{xdimension} is the dimension of the latent process; 
#' \code{ydimension} is the dimension of the observation process; 
#' \code{theta_dimension} is the dimension of the parameter space; 
#' \code{theta_names} is a vector of characters to enumerate the parameters;
#' \code{theta_positivity} is a vector logicals to index parameters with positivity constraints;
#' \code{construct_discretization} outputs a list containing stepsize, nsteps, statelength and obstimes; 
#' \code{construct_successive_discretization} outputs lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level;
#' \code{sigma} is the diffusion coefficient of the process; 
#' \code{rinit} to sample from the initial distribution; 
#' \code{rtransition} to sample from the Markov transition; 
#' \code{dtransition} to evaluate the transition density; 
#' \code{dmeasurement} to evaluate the measurement density; 
#' \code{functional} is the smoothing function to compute gradient of log-likelihood. 
#' @export
hmm_logistic_diffusion_full <- function(times){
  # model dimensions
  xdimension <- 1
  ydimension <- 2
  
  # parameters of the model
  theta_dimension <- 4 # inferring diffusivity parameter
  theta_names <- NULL
  for (j in 1:theta_dimension){
    theta_names <- c(theta_names, paste("theta", j, sep = ""))
  }
  theta_positivity <- c(FALSE, TRUE, TRUE, TRUE)
  
  is_discrete_observation <- TRUE
  
  # time intervals
  time_intervals <- diff(times)
  nintervals <- length(time_intervals)
  min_interval <- min(time_intervals)
  
  # construct discretization
  construct_discretization <- function(level){
    stepsize <- min_interval * 2^(-level) 
    all_stepsizes <- NULL
    obstimes <- TRUE # resample at initial time
    nsteps <- 0 
    for (p in 1:nintervals){
      multiple <- floor(time_intervals[p] / stepsize)
      remainder <- time_intervals[p] %% stepsize
      if (remainder == 0){
        add_stepsizes <- rep(stepsize, multiple)
        nsteps <- nsteps + multiple
        add_obstimes <- c(rep(FALSE, multiple-1), TRUE)
      } else {
        add_stepsizes <- c(rep(stepsize, multiple), remainder)
        nsteps <- nsteps + multiple + 1
        add_obstimes <- c(rep(FALSE, multiple), TRUE)
      }
      all_stepsizes <- c(all_stepsizes, add_stepsizes)
      obstimes <- c(obstimes, add_obstimes)
    }
    statelength <- nsteps + 1
    
    return(list(stepsize = all_stepsizes, nsteps = nsteps, 
                statelength = statelength, obstimes = obstimes))
  }
  
  construct_successive_discretization <- function(level_fine){
    # define step sizes
    stepsize_fine <- min_interval * 2^(-level_fine) 
    stepsize_coarse <- 2 * stepsize_fine
    
    # construct fine discretization
    all_stepsizes_fine <- NULL
    obstimes_fine <- TRUE # resample at initial time
    nsteps_fine <- 0 
    for (p in 1:nintervals){
      multiple <- floor(time_intervals[p] / stepsize_fine)
      remainder <- time_intervals[p] %% stepsize_fine
      if (remainder == 0){
        add_stepsizes <- rep(stepsize_fine, multiple)
        nsteps_fine <- nsteps_fine + multiple
        add_obstimes <- c(rep(FALSE, multiple-1), TRUE)
      } else {
        add_stepsizes <- c(rep(stepsize_fine, multiple), remainder)
        nsteps_fine <- nsteps_fine + multiple + 1
        add_obstimes <- c(rep(FALSE, multiple), TRUE)
      }
      all_stepsizes_fine <- c(all_stepsizes_fine, add_stepsizes)
      obstimes_fine <- c(obstimes_fine, add_obstimes)
    }
    statelength_fine <- nsteps_fine + 1
    fine <- list(stepsize = all_stepsizes_fine, nsteps = nsteps_fine, 
                 statelength = statelength_fine, obstimes = obstimes_fine)
    
    # construct coarse discretization
    all_stepsizes_coarse <- NULL
    obstimes_coarse <- TRUE # resample at initial time
    nsteps_coarse <- 0 
    coarsetimes <- TRUE # index coarse times
    for (p in 1:nintervals){
      multiple <- floor(time_intervals[p] / stepsize_coarse)
      remainder <- time_intervals[p] %% stepsize_coarse
      if (remainder == 0){
        add_stepsizes <- rep(stepsize_coarse, multiple)
        nsteps_coarse <- nsteps_coarse + multiple
        add_obstimes <- c(rep(FALSE, multiple-1), TRUE)
        add_coarsetimes <- rep(c(FALSE, TRUE), multiple)
      } else {
        add_stepsizes <- c(rep(stepsize_coarse, multiple), remainder)
        nsteps_coarse <- nsteps_coarse + multiple + 1
        add_obstimes <- c(rep(FALSE, multiple), TRUE)
        if (floor(remainder / stepsize_fine) == 0){
          add_coarsetimes <- c(rep(c(FALSE, TRUE), multiple), TRUE)
        } else {
          add_coarsetimes <- c(rep(c(FALSE, TRUE), multiple), FALSE, TRUE)
        }
      }
      all_stepsizes_coarse <- c(all_stepsizes_coarse, add_stepsizes)
      obstimes_coarse <- c(obstimes_coarse, add_obstimes)
      coarsetimes <- c(coarsetimes, add_coarsetimes)
    }
    statelength_coarse <- nsteps_coarse + 1
    coarse <- list(stepsize = all_stepsizes_coarse, nsteps = nsteps_coarse, 
                   statelength = statelength_coarse, obstimes = obstimes_coarse)
    
    return(list(fine = fine, coarse = coarse, coarsetimes = coarsetimes))
  }
  
  # drift 
  drift <- function(theta, x) theta[1] / theta[3] - theta[2] * exp(theta[3] * x) / theta[3]
  jacobian_drift <- function(theta, x) c(1 / theta[3], 
                                         - exp(theta[3] * x) / theta[3], 
                                         - theta[1] / theta[3]^2 - theta[2] * exp(theta[3] * x) * (theta[3] * x - 1) / theta[3]^2, 
                                         0)
  
  # diffusivity 
  constant_sigma <- TRUE
  sigma <- 1
  Sigma <- 1
  Omega <- 1
  
  # sample from initial distribution
  init_mean <- 5
  init_sd <- 10
  rinit <- function(nparticles){
    return(matrix(rnorm(nparticles, mean = init_mean, sd = init_sd) / theta[3], nrow = 1)) # returns 1 x N
  }
  
  gradient_init <- function(theta, xstate){
    # theta is a vector of size 3
    # xstate is a number 
    partial_derivative_theta3 <- 1 / theta[3] - theta[3] * (xstate - init_mean / theta[3])^2 / init_sd^2 - 
      init_mean * (xstate - init_mean / theta[3]) / init_sd^2
    
    return(c(0, 0, partial_derivative_theta3, 0))
  }
  
  # sample from Markov transition kernel
  rtransition <- function(theta, stepsize, xparticles, rand){ 
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
    # xparticles is a vector of size N
    # should take xparticles of size 1 x N and output new particles of size 1 x N
    # rand is a vector of size N following a standard normal distribution
    # output new particles in a vector of size N
    
    return(xparticles + stepsize * drift(theta, xparticles) + sigma * sqrt(stepsize) *  rand) 
  }
  
  # evaluate Markov transition density 
  dtransition <- function(theta, stepsize, xparticles, next_xparticles){
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
    # xparticles and next_xparticles are vectors of size N
    
    particle_mean <- xparticles + stepsize * drift(theta, xparticles) 
    return(dnorm(next_xparticles, mean = particle_mean, sd = sigma * sqrt(stepsize), log = TRUE))
  }
  
  dmeasurement <- function(theta, stepsize, xparticles, observation){
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
    # xparticles is a vector of size N
    # observation is a vector of size 2 
    
    obsdensity <- dnbinom(observation[1], size = theta[4], mu = exp(theta[3] * xparticles), log = TRUE) + 
      dnbinom(observation[2], size = theta[4], mu = exp(theta[3] * xparticles), log = TRUE)
    return(obsdensity)
  }
  
  # evaluate gradient of log-observation density
  gradient_dmeasurement <- function(theta, stepsize, xstate, observation){
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
    # xstate is a number 
    # observation is a vector of size 2
    
    exp_theta3_xstate <- exp(theta[3] * xstate)
    
    partial_derivative_theta3 <- - 2 * theta[4] * xstate * exp_theta3_xstate / (theta[4] + exp_theta3_xstate) + 
      (observation[1] + observation[2]) * xstate * (1 - exp_theta3_xstate /  (theta[4] + exp_theta3_xstate))
    
    partial_derivative_theta4 <- digamma(observation[1] + theta[4]) + digamma(observation[2] + theta[4]) - 
      2 * digamma(theta[4]) + 2 * (log(theta[4]) - log(theta[4] + exp_theta3_xstate)) + 
      2 * (1 - theta[4] / (theta[4] + exp_theta3_xstate)) - (observation[1] + observation[2]) / (theta[4] + exp_theta3_xstate)
    
    return(c(0, 0, partial_derivative_theta3, partial_derivative_theta4))
  }
  
  functional <- function(theta, discretization, xtrajectory, observations){
    # theta is a vector of size 3
    # discretization is a list created by construct_discretization
    # xtrajectory is a vector of length statelength 
    # observations is a matrix of size nobservations x 2
    
    # discretization objects
    stepsize <- discretization$stepsize # vector of length nsteps
    nsteps <- discretization$nsteps 
    statelength <- discretization$statelength  
    obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1
    
    # initialize
    output <- rep(0, theta_dimension)
    index_obs <- 1
    xstate <- xtrajectory[1]
    output <- output + gradient_init(theta, xstate)
    output <- output + gradient_dmeasurement(theta, stepsize, xstate, observations[index_obs, ])
    
    # loop over time 
    for (k in 1:nsteps){
      jacobian <- jacobian_drift(theta, xtrajectory[k])
      output <- output - stepsize[k] * Omega * drift(theta, xtrajectory[k]) * jacobian  
      output <- output + Omega * (xtrajectory[k+1] - xtrajectory[k]) * jacobian
      
      # observation time
      if (obstimes[k+1]){
        index_obs <- index_obs + 1
        xstate <- xtrajectory[k+1]
        output <- output + gradient_dmeasurement(theta, stepsize, xstate, observations[index_obs, ])
      }
    }
    return(output)
  }
  
  model <- list(xdimension = xdimension,
                ydimension = ydimension,
                theta_dimension = theta_dimension,
                theta_names = theta_names,
                theta_positivity = theta_positivity,
                is_discrete_observation = is_discrete_observation,
                construct_discretization = construct_discretization,
                construct_successive_discretization = construct_successive_discretization,
                constant_sigma = constant_sigma,
                sigma = sigma,
                rinit = rinit, 
                rtransition = rtransition, 
                dtransition = dtransition, 
                dmeasurement = dmeasurement,
                functional = functional)
  return(model)
  
}