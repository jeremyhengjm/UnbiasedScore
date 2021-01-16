#' @rdname hmm_ornstein_uhlenbeck
#' @title Construct hidden markov model obtained by discretizing Ornstein-Uhlenbeck process
#' @param times observation times
#' @description This function returns a list with objects such as
#' * xdimension is the dimension of the latent process, 
#' * ydimension is the dimension of the observation process, 
#' * theta_dimension is the dimension of the parameter space, 
#' * construct_discretization outputs a list containing stepsize, nsteps, statelength and obstimes,
#' * construct_successive_discretization outputs lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level, 
#' * sigma is the diffusion coefficient of the process, 
#' * rinit to sample from the initial distribution, 
#' * rtransition to sample from the Markov transition, 
#' * dtransition to evaluate the transition density, 
#' * dmeasurement to evaluate the measurement density,
#' * functional is the smoothing function to compute gradient of log-likelihood.
#' @return A list 
#' @export
hmm_ornstein_uhlenbeck <- function(times){
  # model dimensions
  xdimension <- 1
  ydimension <- 1
  
  # parameters of the model
  theta_dimension <- 3
  
  # type of observation model
  is_discrete_observation <- TRUE
  
  # time intervals
  nobservations <- length(times) # nobservations at unit times
  
  # construct discretization
  construct_discretization <- function(level){
    stepsize <- 2^(-level)
    nsteps <- nobservations * 2^level
    all_stepsizes <- rep(stepsize, nsteps)
    statelength <- nsteps + 1
    obstimes <- rep(FALSE, statelength)
    # No observation at first time step for deterministic initialization
    obs_index <- seq(2^level+1, statelength, by = 2^level)
    # obs_index <- seq(1, statelength, by = 2^level)
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
  constant_sigma <- TRUE
  sigma <- 1
  Sigma <- sigma^2
  Omega <- sigma^(-2)
  
  # sample from initial distribution
  x_star <- 0 # deterministic init
  rinit <- function(nparticles){
    #return(matrix(rnorm(nparticles, mean = 0, sd = 1), nrow = 1)) # returns 1 x N
    return(matrix(x_star, nrow = xdimension, ncol = nparticles))
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
  
  # evaluate observation density
  dmeasurement <- function(theta, stepsize, xparticles, observation){
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
    # xparticles is a vector of size N
    # observation is a number
    return(dnorm(observation, mean = xparticles, sd = sqrt(theta[3]), log = TRUE))
  }
  
  # evaluate gradient of log-observation density
  gradient_dmeasurement <- function(theta, stepsize, xstate, observation){
    # theta is a vector of size 3
    # stepsize is the time discretization step size 
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
    index_obs <- 0 
    
    # loop over time 
    for (k in 1:nsteps){
      jacobian <- jacobian_drift(theta, xtrajectory[k])
      output <- output - stepsize[k] * Omega * drift(theta, xtrajectory[k]) * jacobian  # JHo: coeff 0.5 removed
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
  
  # compute moments using Kalman smoothing
  smoothing_moments <- function(theta, observations){
    # theta is a vector of size 3
    # observations is a matrix of size nobservations x 1
    
    # unit step size
    delta <- 1
    
    ## Compute smoothing distribution via the Kalman filter ##
    
    # theta-dependent quantities
    e1 <- exp(- delta * theta[1])
    var_dyn <- sigma ** 2 * (1 - exp(- 2 * delta * theta[1])) / (2 * theta[1])
    var_obs <- theta[3]
    
    # Build prior mean and variance jointly for all time steps
    prior_mean <- matrix(0, nrow = nobservations, ncol=1)
    prior_var <- matrix(0, nrow = nobservations, ncol=nobservations)
    
    for (k in 1:nobservations){
      if (k == 1){
        m_prev <- x_star
        P_prev <- 0
      } else {
        m_prev <- prior_mean[k-1]
        P_prev <- prior_var[k-1, k-1]
      }
      
      prior_mean[k] <- theta[2] + (m_prev - theta[2]) * e1
      prior_var[k, k] <- P_prev * e1 ** 2 + var_dyn
      
      # Correlations between time steps
      for (l in (k+1):nobservations){
        if (l > nobservations) break
        prior_var[k, l] <- e1 ** (l-k) * prior_var[k, k]
        prior_var[l, k] <- prior_var[k, l]
      }
    }
    
    # Observation matrix
    H <- diag(nobservations)
    
    # Apply one Kalman filter step with all observations at once
    innov <- observations - H %*% prior_mean
    cov_innov <- H %*% prior_var %*% t(H) + var_obs * diag(nobservations)
    gain <- prior_var %*% t(H) %*% solve(cov_innov)
    post_mean <- prior_mean + gain %*% innov
    post_var <- (diag(nobservations) - gain %*% H) %*% prior_var
    
    return(list(post_mean = post_mean, post_var = post_var, prior_mean = prior_mean))
  }
  
  compute_gradients <- function(theta, observations){
    # theta is a vector of size 3
    # observations is a matrix of size nobservations x 1
    
    # unit step size
    delta <- 1
    
    e1 <- exp(- delta * theta[1])
    var_dyn <- sigma ** 2 * (1 - exp(- 2 * delta * theta[1])) / (2 * theta[1])
    
    # \nabla_{\theta} \Sigma_{\theta} / \Sigma_{\theta}
    dlogS = ((2 * delta * theta[1] + 1) * exp(- 2 * delta * theta[1]) - 1) / ( theta[1] * (1 - exp(- 2 * delta * theta[1])) )
    c1 = dlogS / (2 * var_dyn)
    
    # Constants of different polynomial terms in the 1st component of \nabla_{\theta} \log p_{\theta}
    c1x2 = c1
    c1xxp = -(delta / var_dyn + 2 * c1) * e1
    c1xp2 = (delta / var_dyn + c1) * e1**2
    c1x = (-2 * c1 * (1 - e1) + delta * e1 / var_dyn) * theta[2]
    c1xp =  (delta * (1 - e1) / var_dyn - delta * e1 / var_dyn + 2 * c1 * (1 - e1)) * e1 * theta[2]
    c1c = c1 * theta[2]**2 * (1 - e1)**2 - (1 - e1) * delta * e1 * theta[2]**2 / var_dyn
    
    moments <- smoothing_moments(theta, observations)
    
    grad <- matrix(0, nrow = nobservations, ncol=3)
    for (k in 1:nobservations){
      if (k == 1){
        Exp <- x_star # Expectation x'
        Vxp <- 0      # Variance x'
        Vxxp <- 0     # Correlation x and x'
      } else {
        Exp <- moments$post_mean[k-1]
        Vxp <- moments$post_var[k-1, k-1]
        Vxxp <- moments$post_var[k-1, k]
      }
      
      # Second moment
      Exp2 <- Vxp + Exp**2
      Ex <- moments$post_mean[k]
      Ex2 <- moments$post_var[k, k] + Ex**2
      Exxp <- Vxxp + Ex * Exp
      
      grad[k, 1] <- c1x2 * Ex2 + c1xxp * Exxp + c1xp2 * Exp2 + c1x * Ex + c1xp * Exp + c1c - dlogS / 2
      grad[k, 2] <- (1 - e1) * (Ex - e1 * Exp + (e1 - 1) * theta[2]) / var_dyn
      grad[k, 3] <- (Ex2 - 2 * observations[k] * Ex + observations[k]**2) / (2 * theta[3]**2) - 1 / (2 * theta[3])
    }
    
    gradient <- colSums(grad)
    return(gradient)
  }
  
  model <- list(xdimension = xdimension,
                ydimension = ydimension,
                theta_dimension = theta_dimension,
                is_discrete_observation = is_discrete_observation,
                construct_discretization = construct_discretization,
                construct_successive_discretization = construct_successive_discretization,
                constant_sigma = constant_sigma,
                sigma = sigma,
                rinit = rinit, 
                rtransition = rtransition, 
                dtransition = dtransition, 
                dmeasurement = dmeasurement,
                functional = functional,
                smoothing_moments = smoothing_moments,
                compute_gradients = compute_gradients)
  return(model)
  
}