#' @rdname multilevel_CBSPF
#' @title Multilevel Conditional Backward Sampling Particle Filter
#' @description Runs two coupled conditional particle filters (one at each discretization level) with backward sampling (Whiteley, 2010)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param ref_trajectory_coarse a matrix of reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_fine a matrix of reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @return two new trajectories stored as matrices of size xdimension x statelength_coarse/fine
#' @export
multilevel_CBSPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                             ref_trajectory_coarse = NULL, ref_trajectory_fine = NULL){
  
  if (is.null(ref_trajectory_coarse) || is.null(ref_trajectory_fine)){
    cat('multilevel_CBSPF must be called with non-null trajectories\n')
  }
  
  # get model/problem settings 
  nobservations <- nrow(observations)
  xdimension <- model$xdimension
  ydimension <- model$ydimension
  
  # discretization
  nsteps <- discretization$fine$nsteps
  stepsize_fine <- discretization$fine$stepsize
  stepsize_coarse <- discretization$coarse$stepsize
  statelength_fine <- discretization$fine$statelength
  statelength_coarse <- discretization$coarse$statelength
  obstimes <- discretization$fine$obstimes # vector of length nsteps+1 indexing observation times
  coarsetimes <- discretization$coarsetimes # vector of length nsteps+1 indexing coarse times
  
  # create tree representation of the trajectories or store all states and ancestors
  xtrajectory_coarse <- array(0, dim = c(statelength_coarse, xdimension, nparticles))
  ancestries_coarse <- matrix(0, nrow = statelength_coarse, ncol = nparticles)
  store_logweights_coarse <- matrix(0, nrow = nobservations, ncol = nparticles)
  xtrajectory_fine <- array(0, dim = c(statelength_fine, xdimension, nparticles))
  ancestries_fine <- matrix(0, nrow = statelength_fine, ncol = nparticles)
  store_logweights_fine <- matrix(0, nrow = nobservations, ncol = nparticles)
  
  # initialization
  xparticles_coarse <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles_fine <- xparticles_coarse
  if (!is.null(ref_trajectory_coarse)){
    xparticles_coarse[, nparticles] <- ref_trajectory_coarse[, 1]
  }
  if (!is.null(ref_trajectory_fine)){
    xparticles_fine[, nparticles] <- ref_trajectory_fine[, 1]
  }
  
  # initialize tree to storage trajectories
  xtrajectory_coarse[1, , ] <- xparticles_coarse
  xtrajectory_fine[1, , ] <- xparticles_fine
  
  logweights_coarse <- rep(0, nparticles)
  logweights_fine <- rep(0, nparticles)
  index_obs <- 0 
  ancestors_coarse <- 1:nparticles
  ancestors_fine <- 1:nparticles
  
  # index last observation
  last_obs_coarse <- 1
  last_obs_fine <- 1
  
  # random initialization
  if (obstimes[1]){
    # compute weights
    index_obs <- index_obs + 1
    observation <- observations[index_obs, ] # 1 x ydimension 
    logweights_fine <- model$dmeasurement(theta, stepsize_fine[1], xparticles_fine, observation)
    store_logweights_fine[index_obs, ] <- logweights_fine
    maxlogweights_fine <- max(logweights_fine)
    weights_fine <- exp(logweights_fine - maxlogweights_fine)
    normweights_fine <- weights_fine / sum(weights_fine)
    ess_fine <- 1 / sum(normweights_fine^2)
    
    logweights_coarse <- model$dmeasurement(theta, stepsize_coarse[1], xparticles_coarse, observation)
    store_logweights_coarse[index_obs, ] <- logweights_coarse
    maxlogweights_coarse <- max(logweights_coarse)
    weights_coarse <- exp(logweights_coarse - maxlogweights_coarse)
    normweights_coarse <- weights_coarse / sum(weights_coarse)
    ess_coarse <- 1 / sum(normweights_coarse^2)
    
    # resampling
    min_ess <- min(ess_fine, ess_coarse)
    if (min_ess < resampling_threshold * nparticles){
      rand <- runif(nparticles)
      ancestors <- coupled_resampling(normweights_coarse, normweights_fine, nparticles, rand)
      ancestors_coarse <- ancestors[, 1]
      ancestors_fine <- ancestors[, 2]
      xparticles_coarse <- xparticles_coarse[, ancestors_coarse]
      xparticles_fine <- xparticles_fine[, ancestors_fine]
      
      # reset weights
      logweights_coarse <- rep(0, nparticles)
      logweights_fine <- rep(0, nparticles)
    }
  }
  
  index_coarse <- 0 # index coarse times
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    xparticles_fine <- model$rtransition(theta, stepsize_fine[k], xparticles_fine, randn) # size: xdimension x nparticles
    if (!is.null(ref_trajectory_fine)){
      xparticles_fine[, nparticles] <- ref_trajectory_fine[, k+1]
      ancestors_fine[nparticles] <- nparticles
    }
    
    if (coarsetimes[k+1]){
      # increment number of coarse steps
      index_coarse <- index_coarse + 1 
      
      # combine Brownian increments
      if (coarsetimes[k]){
        combined_randn <- randn # use same Brownian increment if we cannot fit a fine stepsize in the remainder of coarse level
      } else {
        combined_randn <- (previous_randn + randn) / sqrt(2)
      }
      
      # propagate under latent dynamics
      xparticles_coarse <- model$rtransition(theta, stepsize_coarse[index_coarse], xparticles_coarse, combined_randn) # size: xdimension x nparticles
      if (!is.null(ref_trajectory_coarse)){
        xparticles_coarse[, nparticles] <- ref_trajectory_coarse[, index_coarse+1]
        ancestors_coarse[nparticles] <- nparticles
      }
    }
    previous_randn <- randn
    
    # update tree storage
    xtrajectory_fine[k+1, , ] <- xparticles_fine
    ancestries_fine[k, ] <- ancestors_fine
    ancestors_fine <- 1:nparticles
    
    if (coarsetimes[k+1]){
      xtrajectory_coarse[index_coarse+1, , ] <- xparticles_coarse
      ancestries_coarse[index_coarse, ] <- ancestors_coarse
      ancestors_coarse <- 1:nparticles
    }
    
    if (obstimes[k+1]){
      # compute weights
      index_obs <- index_obs + 1
      observation <- observations[index_obs, ] # 1 x ydimension 
      
      if (model$is_discrete_observation){
        # observation only depends on current particles in a discrete model
        logweights_fine <- logweights_fine + model$dmeasurement(theta, stepsize_fine[k], xparticles_fine, observation)
        
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory <- array(0, dim = c(k-last_obs_fine+1, xdimension, nparticles))
        x_sub_trajectory[ , , ] <- xtrajectory_fine[last_obs_fine:k, , ]
        logweights_fine <- logweights_fine + model$dmeasurement(theta, stepsize_fine[k], x_sub_trajectory, observation)
      }
      
      store_logweights_fine[index_obs, ] <- logweights_fine
      
      # index last observation
      last_obs_fine <- k + 1
      maxlogweights_fine <- max(logweights_fine)
      weights_fine <- exp(logweights_fine - maxlogweights_fine)
      normweights_fine <- weights_fine / sum(weights_fine)
      ess_fine <- 1 / sum(normweights_fine^2)
      
      if (model$is_discrete_observation){
        # observation only depends on current particles in a discrete model
        logweights_coarse <- logweights_coarse + model$dmeasurement(theta, stepsize_coarse[index_coarse], xparticles_coarse, observation)
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory <- array(0, dim = c(index_coarse-last_obs_coarse+1, xdimension, nparticles))
        x_sub_trajectory[ , , ] <- xtrajectory_coarse[last_obs_coarse:index_coarse, , ]
        logweights_coarse <- logweights_coarse + model$dmeasurement(theta, stepsize_coarse[index_coarse], x_sub_trajectory, observation)
      }
      
      store_logweights_coarse[index_obs, ] <- logweights_coarse
      
      # index last observation
      last_obs_coarse <- index_coarse + 1
      maxlogweights_coarse <- max(logweights_coarse)
      weights_coarse <- exp(logweights_coarse - maxlogweights_coarse)
      normweights_coarse <- weights_coarse / sum(weights_coarse)
      ess_coarse <- 1 / sum(normweights_coarse^2)
      
      # resampling
      min_ess <- min(ess_fine, ess_coarse)
      if (k < nsteps && min_ess < resampling_threshold * nparticles){
        rand <- runif(nparticles)
        ancestors <- coupled_resampling(normweights_coarse, normweights_fine, nparticles, rand)
        ancestors_coarse <- ancestors[, 1]
        ancestors_fine <- ancestors[, 2]
        xparticles_coarse <- xparticles_coarse[, ancestors_coarse]
        xparticles_fine <- xparticles_fine[, ancestors_fine]
        
        # reset weights
        logweights_coarse <- rep(0, nparticles)
        logweights_fine <- rep(0, nparticles)
      }
    }
  }
  
  # draw a pair of trajectories 
  # rand <- runif(1)
  # ancestor <- coupled_resampling(normweights_coarse, normweights_fine, 1, rand)
  # ancestor_coarse <- ancestor[, 1]
  # ancestor_fine <- ancestor[, 2]
  
  # new_trajectory_coarse <- get_path(model, discretization$coarse, xtrajectory_coarse, ancestries_coarse, ancestor_coarse)
  # new_trajectory_fine <- get_path(model, discretization$fine, xtrajectory_fine, ancestries_fine, ancestor_fine)
  
  
  # trace ancestry to get a path
  new_trajectory_coarse <- matrix(0, nrow = xdimension, ncol = statelength_coarse)
  new_trajectory_fine <- matrix(0, nrow = xdimension, ncol = statelength_fine)
  particle_index <- coupled_resampling(normweights_coarse, normweights_fine, 1, runif(1))
  particle_index_coarse <- particle_index[, 1]
  particle_index_fine <- particle_index[, 2]
  new_trajectory_coarse[, statelength_coarse] <- xtrajectory_coarse[statelength_coarse, , particle_index_coarse]
  new_trajectory_fine[, statelength_fine] <- xtrajectory_fine[statelength_fine, , particle_index_fine]
  index_obs <- sum(obstimes)
  index_coarse <- sum(coarsetimes) - 1
  
  # last observation has already been taken into account if it is at the last time step
  if (obstimes[statelength_fine]){
    index_obs <- index_obs - 1
  }
  
  for (k in nsteps:1){
    if (obstimes[k]){
      logweights_coarse <- store_logweights_coarse[index_obs, ]
      bs_logweights_coarse <- logweights_coarse + model$dtransition(theta, stepsize_coarse[index_coarse],
                                                                    xtrajectory_coarse[index_coarse, , ], new_trajectory_coarse[, index_coarse+1])
      bs_maxlogweights_coarse <- max(bs_logweights_coarse)
      bs_weights_coarse <- exp(bs_logweights_coarse - bs_maxlogweights_coarse)
      bs_normweights_coarse <- bs_weights_coarse / sum(bs_weights_coarse)
      
      logweights_fine <- store_logweights_fine[index_obs, ]
      bs_logweights_fine <- logweights_fine + model$dtransition(theta, stepsize_fine[k], xtrajectory_fine[k, , ], new_trajectory_fine[, k+1])
      bs_maxlogweights_fine <- max(bs_logweights_fine)
      bs_weights_fine <- exp(bs_logweights_fine - bs_maxlogweights_fine)
      bs_normweights_fine <- bs_weights_fine / sum(bs_weights_fine)
      
      index_obs <- index_obs - 1
      
      particle_index <- coupled_resampling(bs_normweights_coarse, bs_normweights_fine, 1, runif(1))
      particle_index_coarse <- particle_index[, 1]
      particle_index_fine <- particle_index[, 2]
    } else {
      if (coarsetimes[k]){
        particle_index_coarse <- ancestries_coarse[index_coarse, particle_index_coarse]
      }
      particle_index_fine <- ancestries_fine[k, particle_index_fine]
    }
    if (coarsetimes[k]){
      new_trajectory_coarse[, index_coarse] <- xtrajectory_coarse[index_coarse, , particle_index_coarse]
      # decrement index of coarse steps
      index_coarse <- index_coarse - 1 
    }
    new_trajectory_fine[, k] <- xtrajectory_fine[k, , particle_index_fine]
  }
  
  return(list(new_trajectory_coarse = new_trajectory_coarse, new_trajectory_fine = new_trajectory_fine))
  
}