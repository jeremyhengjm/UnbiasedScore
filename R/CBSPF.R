#' @rdname CBSPF
#' @title Conditional Backward Sampling Particle Filter
#' @description Runs a conditional particle filter with backward sampling (Whiteley, 2010)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param ref_trajectory a matrix of reference trajectory, of size xdimension x statelength
#'@return a matrix containing a new trajectory of size xdimension x statelength
#'@export
CBSPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1,
                  ref_trajectory){
  
  # get model/problem settings 
  nobservations <- nrow(observations)
  xdimension <- model$xdimension
  ydimension <- model$ydimension
  
  # discretization
  stepsize <- discretization$stepsize # vector of length nsteps
  statelength <- discretization$statelength
  nsteps <- discretization$nsteps
  obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1
  
  # create tree representation of the trajectories or store all states and ancestors
  xtrajectory <- array(0, dim = c(statelength, xdimension, nparticles))
  ancestries <- matrix(0, nrow = statelength, ncol = nparticles)
  store_logweights <- matrix(0, nrow = nobservations, ncol = nparticles)
  
  # initialization
  xparticles <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles[, nparticles] <- ref_trajectory[, 1]
  xtrajectory[1, , ] <- xparticles
  logweights <- rep(0, nparticles)
  ess <- rep(nparticles, statelength)
  log_ratio_normconst <- 0
  log_normconst_previous <- 0 # log normalizing constant at previous resampling time
  log_normconst <- rep(0, nobservations)
  index_obs <- 0 
  ancestors <- 1:nparticles
  last_obs <- 1 # index last observation
  
  # random initialization (only for discrete observations)
  if (obstimes[1]){
    # compute weights
    index_obs <- index_obs + 1 
    observation <- observations[index_obs, ] # 1 x ydimension 
    logweights <- model$dmeasurement(theta, stepsize[1], xparticles, observation)
    store_logweights[index_obs, ] <- logweights
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[1] <- 1 / sum(normweights^2)
    
    # compute normalizing constant
    log_ratio_normconst <- log(mean(weights)) + maxlogweights  
    log_normconst_previous <- log_ratio_normconst
    log_normconst[1] <- log_ratio_normconst
    
    # resampling
    if (ess[1] < resampling_threshold * nparticles){
      rand <- runif(nparticles)
      ancestors <- multinomial_resampling(normweights, nparticles, rand) 
      xparticles <- xparticles[, ancestors]
      logweights <- rep(0, nparticles)
    }
  }
  
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    # randn <- matrix(Rfast::Rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    xparticles <- model$rtransition(theta, stepsize[k], xparticles, randn) # size: xdimension x nparticles
    xparticles[, nparticles] <- ref_trajectory[, k+1] 
    ancestors[nparticles] <- nparticles
    
    # update tree storage
    xtrajectory[k+1, , ] <- xparticles
    ancestries[k, ] <- ancestors
    ancestors <- 1:nparticles
    
    if (obstimes[k+1]){
      # compute weights
      index_obs <- index_obs + 1 
      observation <- observations[index_obs, ] # 1 x ydimension 
      # dealing with continuous observational model [new code]
      if (model$is_discrete_observation){
        logweights <- logweights + model$dmeasurement(theta, stepsize[k], xparticles, observation)
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory <- array(0, dim = c(k-last_obs+1, xdimension, nparticles))
        x_sub_trajectory[ , , ] <- xtrajectory[last_obs:k, , ]
        logweights <- logweights + model$dmeasurement(theta, stepsize[k], x_sub_trajectory, observation)
      }
      last_obs <- k+1 # index last observation
      # [end of new code]
      store_logweights[index_obs, ] <- logweights
      maxlogweights <- max(logweights)
      weights <- exp(logweights - maxlogweights)
      normweights <- weights / sum(weights)
      ess[k+1] <- 1 / sum(normweights^2)
      
      # resampling
      if (k < nsteps && ess[k+1] < resampling_threshold * nparticles){
        rand <- runif(nparticles)
        ancestors <- multinomial_resampling(normweights, nparticles, rand) 
        xparticles <- xparticles[, ancestors]
        logweights <- rep(0, nparticles)
        
        # compute normalizing constant and update log_normconst_previous
        log_ratio_normconst <- log_normconst_previous + log(mean(weights)) + maxlogweights  
        log_normconst_previous <- log_ratio_normconst # log normalizing constant at previous resampling time
      } else {
        # compute normalizing constant without updating log_normconst_previous
        log_ratio_normconst <- log_normconst_previous + log(mean(weights)) + maxlogweights  
      }
      
      # store normalizing constant
      log_normconst[index_obs] <- log_ratio_normconst
    }
  }
  
  # draw a trajectory
  # ancestor <- multinomial_resampling(normweights, 1, runif(1))
  # new_trajectory <- get_path_bs(model, discretization, xtrajectory, ancestries, ancestor, store_logweights)
  
  # trace ancestry to get a path
  new_trajectory <- matrix(0, nrow = xdimension, ncol = statelength)
  particle_index <- multinomial_resampling(normweights, 1, runif(1))
  new_trajectory[, statelength] <- xtrajectory[statelength, , particle_index]
  index_obs <- sum(obstimes)
  
  # last observation has already been taken into account if it is at the last time step
  if (obstimes[statelength]){
    index_obs <- index_obs - 1
  }
  
  for (k in nsteps:1){
    if (obstimes[k]){
      logweights <- store_logweights[index_obs, ]
      bs_logweights <- logweights + model$dtransition(theta, stepsize[k], xtrajectory[k, , ], new_trajectory[, k+1])
      bs_maxlogweights <- max(bs_logweights)
      bs_weights <- exp(bs_logweights - bs_maxlogweights)
      bs_normweights <- bs_weights / sum(bs_weights)
      index_obs <- index_obs - 1
      
      particle_index <- multinomial_resampling(bs_normweights, 1, runif(1))
    } else {
      particle_index <- ancestries[k, particle_index]
    }
    new_trajectory[, k] <- xtrajectory[k, , particle_index]
  }
  
  return(list(new_trajectory = new_trajectory, ess = ess, log_normconst = log_normconst))
  
}
