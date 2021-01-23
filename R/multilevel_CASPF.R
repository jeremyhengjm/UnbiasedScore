#' @rdname multilevel_CASPF
#' @title Multilevel Conditional Ancestor Sampling Particle Filter 
#' @description Runs two coupled conditional particle filters (one at each discretization level) with ancestor sampling (Lindsten, Jordan and Schon, 2014)
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
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
#' if missing, this function store all states and ancestors
#' @return two new trajectories stored as matrices of size xdimension x statelength_coarse/fine
#' @export
multilevel_CASPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                           ref_trajectory_coarse = NULL, ref_trajectory_fine = NULL, treestorage = FALSE){
  
  if (is.null(ref_trajectory_coarse) || is.null(ref_trajectory_fine)){
    cat('multilevel_CASPF must be called with non-null trajectories\n')
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
  if (treestorage){
    Tree_coarse <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    Tree_fine <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) 
  } else {
    xtrajectory_coarse <- array(0, dim = c(statelength_coarse, xdimension, nparticles))
    ancestries_coarse <- matrix(0, nrow = statelength_coarse, ncol = nparticles)
    xtrajectory_fine <- array(0, dim = c(statelength_fine, xdimension, nparticles))
    ancestries_fine <- matrix(0, nrow = statelength_fine, ncol = nparticles)
  }
  
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
  if (treestorage){
    Tree_coarse$init(xparticles_coarse) 
    Tree_fine$init(xparticles_fine) 
  } else {
    xtrajectory_coarse[1, , ] <- xparticles_coarse
    xtrajectory_fine[1, , ] <- xparticles_fine
  }
  
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
    maxlogweights_fine <- max(logweights_fine)
    weights_fine <- exp(logweights_fine - maxlogweights_fine)
    normweights_fine <- weights_fine / sum(weights_fine)
    ess_fine <- 1 / sum(normweights_fine^2)
    
    logweights_coarse <- model$dmeasurement(theta, stepsize_coarse[1], xparticles_coarse, observation)
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
      
      # ancestor sampling
      as_logweights_coarse <- logweights_coarse + model$dtransition(theta, stepsize_coarse[1], xparticles_coarse, ref_trajectory_coarse[, 2])
      as_logweights_fine <- logweights_fine + model$dtransition(theta, stepsize_fine[1], xparticles_fine, ref_trajectory_fine[, 2])
      as_maxlogweights_coarse <- max(as_logweights_coarse)
      as_maxlogweights_fine <- max(as_logweights_fine)
      as_weights_coarse <- exp(as_logweights_coarse - as_maxlogweights_coarse)
      as_weights_fine <- exp(as_logweights_fine - as_maxlogweights_fine)
      as_normweights_coarse <- as_weights_coarse / sum(as_weights_coarse)
      as_normweights_fine <- as_weights_fine / sum(as_weights_fine)
      as_ancestors <- coupled_resampling(as_normweights_coarse, as_normweights_fine, 1, rand[nparticles])
      ancestors_coarse[nparticles] <- as_ancestors[, 1]
      ancestors_fine[nparticles] <- as_ancestors[, 2]
      
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
    if (treestorage){
      Tree_fine$update(xparticles_fine, ancestors_fine - 1) 
    } else {
      xtrajectory_fine[k+1, , ] <- xparticles_fine
      ancestries_fine[k, ] <- ancestors_fine
    }
    ancestors_fine <- 1:nparticles
    
    if (coarsetimes[k+1]){
      if (treestorage){
        Tree_coarse$update(xparticles_coarse, ancestors_coarse - 1) 
      } else {
        xtrajectory_coarse[index_coarse+1, , ] <- xparticles_coarse
        ancestries_coarse[index_coarse, ] <- ancestors_coarse
      }
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
        
        # ancestor sampling
        as_logweights_coarse <- logweights_coarse + model$dtransition(theta, stepsize_coarse[index_coarse+1], xparticles_coarse, ref_trajectory_coarse[, index_coarse+2])
        as_logweights_fine <- logweights_fine + model$dtransition(theta, stepsize_fine[k+1], xparticles_fine, ref_trajectory_fine[, k+2])
        as_maxlogweights_coarse <- max(as_logweights_coarse)
        as_maxlogweights_fine <- max(as_logweights_fine)
        as_weights_coarse <- exp(as_logweights_coarse - as_maxlogweights_coarse)
        as_weights_fine <- exp(as_logweights_fine - as_maxlogweights_fine)
        as_normweights_coarse <- as_weights_coarse / sum(as_weights_coarse)
        as_normweights_fine <- as_weights_fine / sum(as_weights_fine)
        as_ancestors <- coupled_resampling(as_normweights_coarse, as_normweights_fine, 1, rand[nparticles])
        ancestors_coarse[nparticles] <- as_ancestors[, 1]
        ancestors_fine[nparticles] <- as_ancestors[, 2]
        
        xparticles_coarse <- xparticles_coarse[, ancestors_coarse]
        xparticles_fine <- xparticles_fine[, ancestors_fine]
        
        # reset weights
        logweights_coarse <- rep(0, nparticles)
        logweights_fine <- rep(0, nparticles)
      }
    }
  }
  
  # draw a pair of trajectories 
  rand <- runif(1)
  ancestor <- coupled_resampling(normweights_coarse, normweights_fine, 1, rand)
  ancestor_coarse <- ancestor[, 1]
  ancestor_fine <- ancestor[, 2]
  
  if (treestorage){
    new_trajectory_coarse <- Tree_coarse$get_path(ancestor_coarse - 1)
    new_trajectory_fine <- Tree_fine$get_path(ancestor_fine - 1)
  } else {
    new_trajectory_coarse <- get_path(model, discretization$coarse, xtrajectory_coarse, ancestries_coarse, ancestor_coarse)
    new_trajectory_fine <- get_path(model, discretization$fine, xtrajectory_fine, ancestries_fine, ancestor_fine)
  }
  
  return(list(new_trajectory_coarse = new_trajectory_coarse, new_trajectory_fine = new_trajectory_fine))
  
}