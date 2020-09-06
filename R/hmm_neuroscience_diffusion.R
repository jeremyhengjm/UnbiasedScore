#' @rdname hmm_neuroscience_diffusion
#' @title Construct hidden markov model obtained by discretizing neuroscience diffusion
#' @param spiketimes1 vector specifying cell spike times for the first grid cell
#' @param spiketimes2 vector specifying cell spike times for the second grid cell
#' @description This function returns a list with objects such as
#' * xdimension is the dimension of the latent process, 
#' * ydimension is the dimension of the observation process, 
#' * theta_dimension is the dimension of the parameter space, 
#' * compute_observations compute observation counts given spike times and discretization level, 
#' * construct_discretization outputs a list containing stepsize, nsteps, statelength and obstimes, 
#' * construct_successive_discretization outputs lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' * sigma is the diffusion coefficient of the process, 
#' * rinit to sample from the initial distribution, 
#' * rtransition to sample from the Markov transition, 
#' * dmeasurement to evaluate the measurement density, 
#' * functional is the smoothing function to compute gradient of log-likelihood. 
#' @return A list 
#' @export
hmm_neuroscience_diffusion <- function(spiketimes1, spiketimes2){
  
  # model dimensions
  xdimension <- 2
  ydimension <- 2
  
  # number of parameters to be inferred
  theta_dimension <- 12
  
  # compute observation counts given spike times
  start_time <- floor(min(c(min(spiketimes1), min(spiketimes2))))
  end_time <- ceiling(max(c(max(spiketimes1), max(spiketimes2))))
  
  compute_observations <- function(level){
    nbins <- 2^level # same as time discretization
    times <- seq(start_time, end_time, length.out = nbins + 1)
    observations <- matrix(0, nrow = nbins, ncol = ydimension)
  
    for (i in 1:nbins){
      bin_start <- times[i]
      bin_end <- times[i+1]
      observations[i, 1] <- sum((spiketimes1 > bin_start) & (spiketimes1 <= bin_end))
      observations[i, 2] <- sum((spiketimes2 > bin_start) & (spiketimes2 <= bin_end))
    }
    
    return(list(observations = observations, nbins = nbins, times = times))
  }
  
  # construct discretization
  construct_discretization <- function(level){
    stepsize <- 2^(-level)
    nsteps <- 2^level
    all_stepsizes <- rep(stepsize, nsteps)
    statelength <- nsteps + 1
    obstimes <- rep(TRUE, statelength) # same discretization as the latent process
    obstimes[1] <- FALSE # deterministic initialization 
    
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
  drift <- function(theta, x){
    # theta is a vector of length theta_dimension
    # x is a matrix of size xdimension x nparticles
    # output is a matrix of size xdimension x nparticles
    
    # parameters for cell 1
    alpha1 <- theta[1]
    beta1 <- theta[2]
    gamma1 <- theta[3]
    delta1 <- theta[4]
    sigma1 <- theta[5]
    kappa1 <- theta[6]
    
    # parameters for cell 2
    alpha2 <- theta[7]
    beta2 <- theta[8]
    gamma2 <- theta[9]
    delta2 <- theta[10]
    sigma2 <- theta[11]
    kappa2 <- theta[12]
    
    # compute drift
    nparticles <- ncol(x)
    output <- matrix(0, nrow = xdimension, ncol = nparticles)
    output[1, ] <- alpha1 * tanh(beta1 * sigma2 * x[2, ] + gamma1) / sigma1 - delta1 * x[1, ]
    output[2, ] <- alpha2 * tanh(beta2 * sigma1 * x[1, ] + gamma2) / sigma2 - delta2 * x[2, ]
    
    return(output)
  }
  
  # jacobian of drift
  jacobian_drift <- function(theta, x){
    # theta is a vector of length theta_dimension
    # x is a vector of length xdimension 
    # output is a matrix of size xdimension x theta_dimension
    
    # parameters for cell 1
    alpha1 <- theta[1]
    beta1 <- theta[2]
    gamma1 <- theta[3]
    delta1 <- theta[4]
    sigma1 <- theta[5]
    kappa1 <- theta[6]
    
    # parameters for cell 2
    alpha2 <- theta[7]
    beta2 <- theta[8]
    gamma2 <- theta[9]
    delta2 <- theta[10]
    sigma2 <- theta[11]
    kappa2 <- theta[12]
    
    # compute jacobian matrix
    output <- matrix(0, nrow = xdimension, ncol = theta_dimension)
    
    # non-zero partial derivatives for cell 1
    argument <- beta1 * sigma2 * x[2] + gamma1
    output[1, 1] <- tanh(argument) / sigma1 # w.r.t. alpha1
    output[1, 2] <- alpha1 * sigma2 * x[2] * (1 - tanh(argument)^2) / sigma1 # w.r.t. beta1
    output[1, 3] <- alpha1 * (1 - tanh(argument)^2) / sigma1 # w.r.t. gamma1
    output[1, 4] <- - x[1] # w.r.t. delta1
    output[1, 5] <- - alpha1 * tanh(argument) / sigma1^2 # w.r.t. sigma1
    output[1, 11] <- alpha1 * beta1 * x[2] * (1 - tanh(argument)^2) / sigma1 # w.r.t. sigma2
    
    # non-zero partial derivatives for cell 2
    argument <- beta2 * sigma1 * x[1] + gamma2
    output[2, 5] <- alpha2 * beta2 * x[1] * (1 - tanh(argument)^2) / sigma2 # w.r.t. sigma1
    output[2, 7] <- tanh(argument) / sigma2 # w.r.t. alpha2
    output[2, 8] <- alpha2  * sigma1 * x[1] * (1 - tanh(argument)^2) / sigma2 # w.r.t. beta2
    output[2, 9] <- alpha2 * (1 - tanh(argument)^2) / sigma2 # w.r.t. gamma2
    output[2, 10] <- - x[2] # w.r.t. delta2
    output[2, 11] <- - alpha2 * tanh(argument) / sigma2^2 # w.r.t sigma2
    
    return(output)
  } 
  
  # # diffusivity 
  sigma <- 1
  # Sigma <- 1
  # Omega <- 1
  
  # sample from initial distribution
  rinit <- function(nparticles){
    # output is a matrix of size xdimension x nparticles
    return(matrix(0, nrow = xdimension, ncol = nparticles)) 
  }
  
  # sample from Markov transition kernel
  rtransition <- function(theta, stepsize, xparticles, rand){ 
    # theta is a vector of length theta_dimension
    # stepsize is the time discretization step size 
    # xparticles is a matrix of size xdimension x nparticles
    # rand is a matrix of size xdimension x nparticles following the standard normal distribution
    # outputs new particles in a matrix of size xdimension x nparticles
    
    return(xparticles + stepsize * drift(theta, xparticles) + sqrt(stepsize) *  rand) 
  }
  
  dmeasurement <- function(theta, stepsize, xparticles, observation){
    # theta is a vector of size theta_dimension
    # stepsize is the time discretization step size 
    # xparticles is a matrix of size xdimension x nparticles
    # observation is a vector of size ydimension  
    # output is a vector of size nparticles
    
    # parameters in observation model
    kappa1 <- theta[6]
    kappa2 <- theta[12]
    
    # evaluate observation density
    obsdensity <- observation[1] * (log(stepsize) + kappa1 + xparticles[1, ]) - stepsize * exp(kappa1 + xparticles[1, ]) - lfactorial(observation[1]) + 
                  observation[2] * (log(stepsize) + kappa2 + xparticles[2, ]) - stepsize * exp(kappa2 + xparticles[2, ]) - lfactorial(observation[2])
    
    return(obsdensity)
  }
  
  # evaluate gradient of log-observation density
  gradient_dmeasurement <- function(theta, stepsize, xstate, observation){
    # theta is a vector of size theta_dimension
    # stepsize is the time discretization step size 
    # xstate is a vector of size xdimension
    # observation is a vector of size ydimension  
    # output is a vector of size theta_dimension
    
    # parameters in observation model
    kappa1 <- theta[6]
    kappa2 <- theta[12]
    
    # compute gradient of log-observation density 
    output <- rep(0, theta_dimension)
    output[6] <- observation[1] - stepsize * exp(kappa1 + xstate[1]) # w.r.t. kappa1
    output[12] <- observation[2] - stepsize * exp(kappa2 + xstate[2]) # w.r.t. kappa2 
    
    return(output)
  }
  
  functional <- function(theta, discretization, xtrajectory, observations){
    # theta is a vector of size theta_dimension
    # discretization is a list created by construct_discretization
    # xtrajectory is a matrix of size xdimension x statelength
    # observations is a matrix of size nobservations x ydimension
    
    # discretization objects
    stepsize <- discretization$stepsize # vector of length nsteps
    nsteps <- discretization$nsteps 
    statelength <- discretization$statelength  
    obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1
    
    # initialize
    output <- rep(0, theta_dimension)
    index_obs <- 0 # deterministic initialization
  
    # loop over time 
    for (k in 1:nsteps){
      xstate <- xtrajectory[, k]
      jacobian <- jacobian_drift(theta, xstate) # matrix of size xdimension x theta_dimension
      current_drift <- drift(theta, matrix(xstate, nrow = xdimension, ncol = 1)) # matrix of size xdimension x 1
      output <- output - stepsize[k] * as.numeric(t(jacobian) %*% current_drift)
      
      next_state <- xtrajectory[, k+1]
      xincrement <- matrix(next_state - xstate, nrow = xdimension, ncol = 1) # matrix of size xdimension x 1
      output <- output + as.numeric(t(jacobian) %*% xincrement)
      
      # observation time
      if (obstimes[k+1]){
        index_obs <- index_obs + 1
        xstate <- xtrajectory[, k+1]
        output <- output + gradient_dmeasurement(theta, stepsize[k], xstate, observations[index_obs, ])
      }
    }
    return(output)
  }
  
  model <- list(xdimension = xdimension,
                ydimension = ydimension,
                theta_dimension = theta_dimension,
                compute_observations = compute_observations, 
                construct_discretization = construct_discretization,
                construct_successive_discretization = construct_successive_discretization,
                sigma = sigma, 
                rinit = rinit, 
                rtransition = rtransition, 
                dmeasurement = dmeasurement,
                functional = functional)
  return(model)
  
}