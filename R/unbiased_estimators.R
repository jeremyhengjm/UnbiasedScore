#' @rdname unbiased_gradient
#' @title Unbiased estimator of the gradient of the log-likelihood at a discretization level
#' @description Estimates the expectation of a functional with respect to the smoothing distribution at a discretization level
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param level time discretization level
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param coupled_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_coupled_residuals}}
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param max_niterations iteration at which to stop the while loop (default to infinity)
#' @return a list with objects such as 
#' mcmcestimator is the MCMC estimator of the gradient at level
#' unbiasedestimator is an unbiased estimator of the gradient at level
#' meetingtime is the meeting time of the two chains at level
#' iteration is the number of iterations taken
#' finished indicates if the algorithm has completed successfully
#' @export
unbiased_gradient <- function(model, theta, level, observations, nparticles, coupled_resampling, 
                              k = 0, m = 1, max_iterations = Inf){
  
  # initialize chains
  chain_state1 <- CPF(model, theta, level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  chain_state2 <- CPF(model, theta, level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  
  # initialize estimators computation
  mcmcestimator <- model$functional(theta, level, chain_state1, observations)
  theta_dimension <- model$theta_dimension
  if (k > 0){
    mcmcestimator <- rep(0, theta_dimension)
  }
  
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, theta_dimension)
  chain_state1 <- CPF(model, theta, level, observations, nparticles, chain_state1)$new_trajectory
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * 
      (model$functional(theta, level, chain_state1, observations) - 
         model$functional(theta, level, chain_state2, observations))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + model$functional(theta, level, chain_state1, observations)
  }
  
  # initialize
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    # increment counter
    iter <- iter + 1
    
    # sample from 2-way coupled CPF kernel
    coupled2CPF <- coupled2_CPF(model, theta, level, observations, nparticles, coupled_resampling, 
                                chain_state1, chain_state2)
    chain_state1 <- coupled2CPF$new_trajectory1
    chain_state2 <- coupled2CPF$new_trajectory2
    
    # update gradient estimators
    if (meet){
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + model$functional(theta, level, chain_state1, observations)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + model$functional(theta, level, chain_state1, observations)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * 
          (model$functional(theta, level, chain_state1, observations) - 
             model$functional(theta, level, chain_state2, observations))
      }
    }
    
    # check if meeting occurs 
    if (all(chain_state1 == chain_state2) && !meet){
      meet <- TRUE # recording meeting time tau
      meetingtime <- iter
    }
    
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
    
  }
  
  # compute mcmc gradient estimator
  mcmcestimator <- mcmcestimator / (m - k + 1)
  
  # compute unbiased gradient estimator
  unbiasedestimator <- mcmcestimator + correction
  
  return(list(mcmcestimator = mcmcestimator, unbiasedestimator = unbiasedestimator, 
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

#' @rdname unbiased_gradient_increment
#' @title Unbiased estimator of the difference of the gradient of log-likelihood at two discretization levels
#' @description Estimates the difference of the expectation of a functional with respect to the smoothing distribution at two discretization levels
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param level time discretization level
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param coupled_resampling a 4-way coupled resampling scheme, such as \code{\link{coupled4_maximal_coupled_residuals}}
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param max_niterations iteration at which to stop the while loop (default to infinity)
#' @return a list with objects such as 
#' mcmcestimator_coarse is the MCMC estimator of the gradient at level-1
#' mcmcestimator_fine is the MCMC estimator of the gradient at level
#' unbiasedestimator_coarse is an unbiased estimator of the gradient at level-1
#' unbiasedestimator_fine is an unbiased estimator of the gradient at level
#' mcmcestimator is the MCMC estimator of the gradient increment between the two discretization levels
#' unbiasedestimator is an unbiased estimator of the gradient increment between the two discretization levels
#' meetingtime_coarse is the meeting time of the two chains at level-1
#' meetingtime_fine is the meeting time of the two chains at level
#' iteration is the number of iterations taken
#' finished indicates if the algorithm has completed successfully
#' @export
unbiased_gradient_increment <- function(model, theta, level, observations, nparticles, coupled_resampling, 
                                        k = 0, m = 1, max_iterations = Inf){
  
  # initialize chains
  chain_state_coarse1 <- CPF(model, theta, level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  chain_state_coarse2 <- CPF(model, theta, level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  chain_state_fine1 <- CPF(model, theta, level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  chain_state_fine2 <- CPF(model, theta, level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
  
  
  # initialize coarse and fine estimators computation
  mcmcestimator_coarse <- model$functional(theta, level-1, chain_state_coarse1, observations)
  mcmcestimator_fine <- model$functional(theta, level, chain_state_fine1, observations)
  theta_dimension <- model$theta_dimension
  if (k > 0){
    mcmcestimator_coarse <- rep(0, theta_dimension)
    mcmcestimator_fine <- rep(0, theta_dimension)
  }
  
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction_coarse <- rep(0, theta_dimension)
  correction_fine <- rep(0, theta_dimension)
  chain_state_coarse1 <- CPF(model, theta, level-1, observations, nparticles, chain_state_coarse1)$new_trajectory
  chain_state_fine1 <- CPF(model, theta, level, observations, nparticles, chain_state_fine1)$new_trajectory
  if (k == 0){
    correction_coarse <- correction_coarse + (min(1, (0 - k + 1)/(m - k + 1))) *
      (model$functional(theta, level-1, chain_state_coarse1, observations) - 
         model$functional(theta, level-1, chain_state_coarse2, observations))
    correction_fine <- correction_fine + (min(1, (0 - k + 1)/(m - k + 1))) * 
      (model$functional(theta, level, chain_state_fine1, observations) - 
         model$functional(theta, level, chain_state_fine2, observations))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator_coarse <- mcmcestimator_coarse + model$functional(theta, level-1, chain_state_coarse1, observations)
    mcmcestimator_fine <- mcmcestimator_fine + model$functional(theta, level, chain_state_fine1, observations)
  }
  
  # initialize
  iter <- 1
  meet_coarse <- FALSE
  meet_fine <- FALSE
  finished <- FALSE
  meetingtime_coarse <- Inf
  meetingtime_fine <- Inf
  
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    # increment counter
    iter <- iter + 1
    
    # sample from 4-way coupled CPF kernel
    coupled4CPF <- coupled4_CPF(model, theta, level, observations, nparticles, coupled_resampling,
                                chain_state_coarse1, chain_state_coarse2,
                                chain_state_fine1, chain_state_fine2)
    chain_state_coarse1 <- coupled4CPF$new_trajectory_coarse1
    chain_state_coarse2 <- coupled4CPF$new_trajectory_coarse2
    chain_state_fine1 <- coupled4CPF$new_trajectory_fine1
    chain_state_fine2 <- coupled4CPF$new_trajectory_fine2
    
    # update gradient estimators for coarse discretization level
    if (meet_coarse){
      if (k <= iter && iter <= m){
        mcmcestimator_coarse <- mcmcestimator_coarse + model$functional(theta, level-1, chain_state_coarse1, observations)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_coarse <- mcmcestimator_coarse + model$functional(theta, level-1, chain_state_coarse1, observations)
        }
        correction_coarse <- correction_coarse + (min(1, (iter-1 - k + 1)/(m - k + 1))) * 
          (model$functional(theta, level-1, chain_state_coarse1, observations) - 
             model$functional(theta, level-1, chain_state_coarse2, observations))
      }
    }
    
    # update gradient estimators for fine discretization level
    if (meet_fine){
      if (k <= iter && iter <= m){
        mcmcestimator_fine <- mcmcestimator_fine + model$functional(theta, level, chain_state_fine1, observations)
      }
    } else {
      if (k <= iter){
        if (iter <= m){
          mcmcestimator_fine <- mcmcestimator_fine + model$functional(theta, level, chain_state_fine1, observations)
        }
        correction_fine <- correction_fine + (min(1, (iter-1 - k + 1)/(m - k + 1))) * 
          (model$functional(theta, level, chain_state_fine1, observations) - 
             model$functional(theta, level, chain_state_fine2, observations))
      }
    }
    
    # check if meeting occurs for coarse discretization level
    if (all(chain_state_coarse1 == chain_state_coarse2) && !meet_coarse){
      meet_coarse <- TRUE # recording meeting time tau_coarse
      meetingtime_coarse <- iter
    }
    
    # check if meeting occurs for fine discretization level
    if (all(chain_state_fine1 == chain_state_fine2) && !meet_fine){
      meet_fine <- TRUE # recording meeting time tau_fine
      meetingtime_fine <- iter
    }
    
    # stop after max(m, tau_coarse, tau_fine) steps
    if (iter >= max(meetingtime_coarse, meetingtime_fine, m)){
      finished <- TRUE
    }
    
  }
  
  # compute mcmc gradient estimators and their difference 
  mcmcestimator_coarse <- mcmcestimator_coarse / (m - k + 1)
  mcmcestimator_fine <- mcmcestimator_fine / (m - k + 1)
  mcmcestimator <- mcmcestimator_fine - mcmcestimator_coarse
  
  # compute unbiased gradient estimators and their difference
  unbiasedestimator_coarse <- mcmcestimator_coarse + correction_coarse
  unbiasedestimator_fine <- mcmcestimator_fine + correction_fine
  unbiasedestimator <- unbiasedestimator_fine - unbiasedestimator_coarse 
  
  return(list(mcmcestimator_coarse = mcmcestimator_coarse, mcmcestimator_fine = mcmcestimator_fine, 
              unbiasedestimator_coarse = unbiasedestimator_coarse, unbiasedestimator_fine = unbiasedestimator_fine, 
              mcmcestimator = mcmcestimator, unbiasedestimator = unbiasedestimator, 
              meetingtime_coarse = meetingtime_coarse, meetingtime_fine = meetingtime_fine, 
              iteration = iter, finished = finished))
}


