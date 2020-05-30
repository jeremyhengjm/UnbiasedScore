#' @rdname hmm_ornstein_uhlenbeck
#' @title Construct hidden markov model obtained by discretizing Ornstein-Uhlenbeck process
#' @param terminal_time number specifying time interval of interest
#' @description This function returns a list with objects such as
#' * xdimension is the dimension of the latent process
#' * ydimension is the dimension of the observation process
#' * theta_dimension is the dimension of the parameter space
#' * statelength is the length of the discretized process at each level
#' * obstimes is a vector of logicals indexing the observation times
#' * sigma is the diffusion coefficient of the process
#' * rinit to sample from the initial distribution
#' * rtransition to sample from the Markov transition
#' * dtransition to evaluate the transition density
#' * dmeasurement to evaluate the measurement density
#' * functional is the smoothing function to compute gradient of log-likelihood 
#' @return A list 
#' @export
hmm_ornstein_uhlenbeck <- function(terminal_time){
  # model dimensions
  xdimension <- 1
  ydimension <- 1
  
  # parameters of the model
  theta_dimension <- 3
  
  # drift 
  drift <- function(theta, x) theta[1] * (theta[2] - x) # drift
  jacobian_drift <- function(theta, x) c(theta[2] - x, theta[1], 0)
  
  # diffusivity 
  sigma <- 1
  Sigma <- sigma^2
  Omega <- sigma^(-2)
  
  # compute length of latent states
  statelength <- function(level){
    return(terminal_time * 2^level + 1)
  }
  
  # sample from initial distribution
  x_star <- 0
  rinit <- function(nparticles){
    return(matrix(rep(x_star, nparticles), nrow = 1)) # returns 1 x N
  }
  
  # sample from Markov transition kernel
  rtransition <- function(theta, level, xparticles, rand){ 
    # theta is a vector of size 3
    # xparticles is a vector of size N
    # should take xparticles of size 1 x N and output new particles of size 1 x N
    # rand is a vector of size N following a standard normal distribution
    # output new particles in a vector of size N
    stepsize <- 2^(-level) # stepsize of discretization
    return(xparticles + stepsize * drift(theta, xparticles) + sigma * sqrt(stepsize) *  rand) 
  }
  
  # evaluate Markov transition density 
  dtransition <- function(theta, level, xparticles, next_xparticles){
    # theta is a vector of size 3
    # xparticles and next_xparticles are vectors of size N
    stepsize <- 2^(-level) # stepsize of discretization
    particle_mean <- xparticles + stepsize * drift(theta, xparticles) 
    return(dnorm(next_xparticles, mean = particle_mean, sd = sigma * sqrt(stepsize), log = TRUE))
  }
  
  # observations are at unit times
  obstimes <- function(level){
    obs_times <- rep(FALSE, terminal_time * 2^level + 1)
    obs_index <- seq(2^level + 1, terminal_time * 2^level + 1, by = 2^level)
    obs_times[obs_index] <- TRUE
    return(obs_times)
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
  
  functional <- function(theta, level, xtrajectory, observations){
    # theta is a vector of size 3
    # xtrajectory is a vector of length statelength = terminal_time * 2^level + 1
    # observations is a vector of length terminal_time
    stepsize <- 2^(-level) # stepsize of discretization
    nsteps <- terminal_time * 2^level
    output <- rep(0, theta_dimension)
    for (k in 1:nsteps){
      jacobian <- jacobian_drift(theta, xtrajectory[k])
      output <- output - stepsize * Omega * drift(theta, xtrajectory[k]) * jacobian  # JHo: coeff 0.5 removed
      output <- output + Omega * (xtrajectory[k+1] - xtrajectory[k]) * jacobian
    }
    for (p in 1:terminal_time){
      xstate <- xtrajectory[p * 2^level + 1]
      output <- output + gradient_dmeasurement(theta, xstate, observations[p])
    }
    return(output)
  }
  
  model <- list(xdimension = xdimension,
                ydimension = ydimension,
                theta_dimension = theta_dimension,
                statelength = statelength,
                obstimes = obstimes,
                sigma = sigma,
                rinit = rinit, 
                rtransition = rtransition, 
                dtransition = dtransition, 
                dmeasurement = dmeasurement,
                functional = functional)
  return(model)
}