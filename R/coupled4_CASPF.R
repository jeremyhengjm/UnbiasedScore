#' @rdname coupled4_CASPF
#' @title 4-way Coupled Ancestor Sampling Conditional Particle Filter
#' @description Runs four coupled conditional particle filters (two at each discretization level)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 4-way coupled resampling scheme, such as \code{\link{coupled4_maximalchains_maximallevels_independent_residuals}}
#' @param ref_trajectory_coarse1 a matrix of first reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_coarse2 a matrix of second reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_fine1 a matrix of first reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param ref_trajectory_fine2 a matrix of second reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
#' if missing, this function store all states and ancestors
#' @return four new trajectories stored as matrices of size xdimension x statelength_coarse/fine
#' @export
coupled4_CASPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                         ref_trajectory_coarse1, ref_trajectory_coarse2,
                         ref_trajectory_fine1, ref_trajectory_fine2, treestorage = FALSE){
  
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
  
  # check if trajectories are equal
  meet_coarse <- all(ref_trajectory_coarse1 == ref_trajectory_coarse2)
  meet_fine <- all(ref_trajectory_fine1 == ref_trajectory_fine2)
  
  # create tree representation of the trajectories or store all states and ancestors
  if (meet_coarse){
    if (treestorage){
      Tree_coarse1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) # only store one tree in this case
    } else {
      xtrajectory_coarse1 <- array(0, dim = c(statelength_coarse, xdimension, nparticles))
      ancestries_coarse1 <- matrix(0, nrow = statelength_coarse, ncol = nparticles)
    }
  } else {
    if (treestorage){
      Tree_coarse1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
      Tree_coarse2 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    } else {
      xtrajectory_coarse1 <- array(0, dim = c(statelength_coarse, xdimension, nparticles))
      xtrajectory_coarse2 <- array(0, dim = c(statelength_coarse, xdimension, nparticles))
      ancestries_coarse1 <- matrix(0, nrow = statelength_coarse, ncol = nparticles)
      ancestries_coarse2 <- matrix(0, nrow = statelength_coarse, ncol = nparticles)
    }
  }
  
  if (meet_fine){
    if (treestorage){
      Tree_fine1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) # only store one tree in this case
    } else {
      xtrajectory_fine1 <- array(0, dim = c(statelength_fine, xdimension, nparticles))
      ancestries_fine1 <- matrix(0, nrow = statelength_fine, ncol = nparticles)
    }
  } else {
    if (treestorage){
      Tree_fine1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
      Tree_fine2 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    } else {
      xtrajectory_fine1 <- array(0, dim = c(statelength_fine, xdimension, nparticles))
      xtrajectory_fine2 <- array(0, dim = c(statelength_fine, xdimension, nparticles))
      ancestries_fine1 <- matrix(0, nrow = statelength_fine, ncol = nparticles)
      ancestries_fine2 <- matrix(0, nrow = statelength_fine, ncol = nparticles)
    }
  }
  
  # initialization
  xparticles_coarse1 <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles_coarse2 <- xparticles_coarse1
  xparticles_fine1 <- xparticles_coarse1
  xparticles_fine2 <- xparticles_coarse1 
  xparticles_coarse1[, nparticles] <- ref_trajectory_coarse1[, 1]
  xparticles_coarse2[, nparticles] <- ref_trajectory_coarse2[, 1]
  xparticles_fine1[, nparticles] <- ref_trajectory_fine1[, 1]
  xparticles_fine2[, nparticles] <- ref_trajectory_fine2[, 1]
  
  # initialize tree to storage trajectories
  if (meet_coarse){
    if (treestorage){
      Tree_coarse1$init(xparticles_coarse1) 
    } else {
      xtrajectory_coarse1[1, , ] <- xparticles_coarse1
    }
  } else {
    if (treestorage){
      Tree_coarse1$init(xparticles_coarse1) 
      Tree_coarse2$init(xparticles_coarse2) 
    } else {
      xtrajectory_coarse1[1, , ] <- xparticles_coarse1
      xtrajectory_coarse2[1, , ] <- xparticles_coarse2
    }
  }
  
  if (meet_fine){
    if (treestorage){
      Tree_fine1$init(xparticles_fine1) 
    } else {
      xtrajectory_fine1[1, , ] <- xparticles_fine1
    }
  } else {
    if (treestorage){
      Tree_fine1$init(xparticles_fine1) 
      Tree_fine2$init(xparticles_fine2) 
    } else {
      xtrajectory_fine1[1, , ] <- xparticles_fine1
      xtrajectory_fine2[1, , ] <- xparticles_fine2
    }
  }
  logweights_coarse1 <- rep(0, nparticles)
  logweights_coarse2 <- rep(0, nparticles)
  logweights_fine1 <- rep(0, nparticles)
  logweights_fine2 <- rep(0, nparticles)
  index_obs <- 0 
  ancestors_coarse1 <- 1:nparticles
  ancestors_coarse2 <- 1:nparticles
  ancestors_fine1 <- 1:nparticles
  ancestors_fine2 <- 1:nparticles
  
  # index last observation
  last_obs_fine <- 1
  last_obs_coarse <- 1
  
  # random initialization
  if (obstimes[1]){
    # compute weights
    index_obs <- index_obs + 1
    observation <- observations[index_obs, ] # 1 x ydimension 
    if (meet_fine){
      logweights_fine1 <- model$dmeasurement(theta, stepsize_fine[1], xparticles_fine1, observation)
      logweights_fine2 <- logweights_fine1
    } else {
      logweights_fine1 <- model$dmeasurement(theta, stepsize_fine[1], xparticles_fine1, observation)
      logweights_fine2 <- model$dmeasurement(theta, stepsize_fine[1], xparticles_fine2, observation)
    }
    maxlogweights_fine1 <- max(logweights_fine1)
    maxlogweights_fine2 <- max(logweights_fine2)
    weights_fine1 <- exp(logweights_fine1 - maxlogweights_fine1)
    weights_fine2 <- exp(logweights_fine2 - maxlogweights_fine2)
    normweights_fine1 <- weights_fine1 / sum(weights_fine1)
    normweights_fine2 <- weights_fine2 / sum(weights_fine2)
    ess_fine1 <- 1 / sum(normweights_fine1^2)
    ess_fine2 <- 1 / sum(normweights_fine2^2)
    
    if (meet_coarse){
      logweights_coarse1 <- model$dmeasurement(theta, stepsize_coarse[1], xparticles_coarse1, observation)
      logweights_coarse2 <- logweights_coarse1
    } else {
      logweights_coarse1 <- model$dmeasurement(theta, stepsize_coarse[1], xparticles_coarse1, observation)
      logweights_coarse2 <- model$dmeasurement(theta, stepsize_coarse[1], xparticles_coarse2, observation)
    }
    maxlogweights_coarse1 <- max(logweights_coarse1)
    maxlogweights_coarse2 <- max(logweights_coarse2)
    weights_coarse1 <- exp(logweights_coarse1 - maxlogweights_coarse1)
    weights_coarse2 <- exp(logweights_coarse2 - maxlogweights_coarse2)
    normweights_coarse1 <- weights_coarse1 / sum(weights_coarse1)
    normweights_coarse2 <- weights_coarse2 / sum(weights_coarse2)
    ess_coarse1 <- 1 / sum(normweights_coarse1^2)
    ess_coarse2 <- 1 / sum(normweights_coarse2^2)
    
    # resampling
    min_ess <- min(ess_fine1, ess_fine2, ess_coarse1, ess_coarse2)
    if (min_ess < resampling_threshold * nparticles){
      rand <- runif(nparticles)
      ancestors <- coupled_resampling(normweights_coarse1, normweights_coarse2,
                                      normweights_fine1, normweights_fine2,
                                      nparticles, rand)
      ancestors_coarse1 <- ancestors[, 1]
      ancestors_coarse2 <- ancestors[, 2]
      ancestors_fine1 <- ancestors[, 3]
      ancestors_fine2 <- ancestors[, 4]
      
      # ancestor sampling
      as_logweights_fine1 <- logweights_fine1 + model$dtransition(theta, stepsize_fine[1], xparticles_fine1, ref_trajectory_fine1[, 2])
      as_logweights_fine2 <- logweights_fine2 + model$dtransition(theta, stepsize_fine[1], xparticles_fine2, ref_trajectory_fine2[, 2])
      as_maxlogweights_fine1 <- max(as_logweights_fine1)
      as_maxlogweights_fine2 <- max(as_logweights_fine2)
      as_weights_fine1 <- exp(as_logweights_fine1 - as_maxlogweights_fine1)
      as_weights_fine2 <- exp(as_logweights_fine2 - as_maxlogweights_fine2)
      as_normweights_fine1 <- as_weights_fine1 / sum(as_weights_fine1)
      as_normweights_fine2 <- as_weights_fine2 / sum(as_weights_fine2)
      
      as_logweights_coarse1 <- logweights_coarse1 + model$dtransition(theta, stepsize_coarse[1], xparticles_coarse1, ref_trajectory_coarse1[, 2])
      as_logweights_coarse2 <- logweights_coarse2 + model$dtransition(theta, stepsize_coarse[1], xparticles_coarse2, ref_trajectory_coarse2[, 2])
      as_maxlogweights_coarse1 <- max(as_logweights_coarse1)
      as_maxlogweights_coarse2 <- max(as_logweights_coarse2)
      as_weights_coarse1 <- exp(as_logweights_coarse1 - as_maxlogweights_coarse1)
      as_weights_coarse2 <- exp(as_logweights_coarse2 - as_maxlogweights_coarse2)
      as_normweights_coarse1 <- as_weights_coarse1 / sum(as_weights_coarse1)
      as_normweights_coarse2 <- as_weights_coarse2 / sum(as_weights_coarse2)
      
      as_ancestors <- coupled_resampling(as_normweights_coarse1, as_normweights_coarse2,
                                         as_normweights_fine1, as_normweights_fine2,
                                         1, rand[nparticles])
      
      ancestors_coarse1[nparticles] <- as_ancestors[, 1]
      ancestors_coarse2[nparticles] <- as_ancestors[, 2]
      ancestors_fine1[nparticles] <- as_ancestors[, 3]
      ancestors_fine2[nparticles] <- as_ancestors[, 4]
      
      xparticles_coarse1 <- xparticles_coarse1[, ancestors_coarse1]
      xparticles_coarse2 <- xparticles_coarse2[, ancestors_coarse2]
      xparticles_fine1 <- xparticles_fine1[, ancestors_fine1]
      xparticles_fine2 <- xparticles_fine2[, ancestors_fine2]
      
      # reset weights
      logweights_coarse1 <- rep(0, nparticles)
      logweights_coarse2 <- rep(0, nparticles)
      logweights_fine1 <- rep(0, nparticles)
      logweights_fine2 <- rep(0, nparticles)
    }
  }
  
  index_coarse <- 0 # index coarse times
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    if (meet_fine){
      xparticles_fine1 <- model$rtransition(theta, stepsize_fine[k], xparticles_fine1, randn) # size: xdimension x nparticles
      xparticles_fine2 <- xparticles_fine1
    } else {
      xparticles_fine1 <- model$rtransition(theta, stepsize_fine[k], xparticles_fine1, randn) # size: xdimension x nparticles
      xparticles_fine2 <- model$rtransition(theta, stepsize_fine[k], xparticles_fine2, randn) # size: xdimension x nparticles
    }
    xparticles_fine1[, nparticles] <- ref_trajectory_fine1[, k+1]
    xparticles_fine2[, nparticles] <- ref_trajectory_fine2[, k+1]
    
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
      if (meet_coarse){
        xparticles_coarse1 <- model$rtransition(theta, stepsize_coarse[index_coarse], xparticles_coarse1, combined_randn) # size: xdimension x nparticles
        xparticles_coarse2 <- xparticles_coarse1
      } else {
        xparticles_coarse1 <- model$rtransition(theta, stepsize_coarse[index_coarse], xparticles_coarse1, combined_randn) # size: xdimension x nparticles
        xparticles_coarse2 <- model$rtransition(theta, stepsize_coarse[index_coarse], xparticles_coarse2, combined_randn) 
      }
      xparticles_coarse1[, nparticles] <- ref_trajectory_coarse1[, index_coarse+1]
      xparticles_coarse2[, nparticles] <- ref_trajectory_coarse2[, index_coarse+1]
    }
    previous_randn <- randn
    
    # update tree storage
    if (meet_fine){
      if (treestorage){
        Tree_fine1$update(xparticles_fine1, ancestors_fine1 - 1) 
      } else {
        xtrajectory_fine1[k+1, , ] <- xparticles_fine1
        ancestries_fine1[k, ] <- ancestors_fine1
      }
    } else {
      if (treestorage){
        Tree_fine1$update(xparticles_fine1, ancestors_fine1 - 1) 
        Tree_fine2$update(xparticles_fine2, ancestors_fine2 - 1) 
      } else {
        xtrajectory_fine1[k+1, , ] <- xparticles_fine1
        xtrajectory_fine2[k+1, , ] <- xparticles_fine2
        ancestries_fine1[k, ] <- ancestors_fine1
        ancestries_fine2[k, ] <- ancestors_fine2
      }
    }
    ancestors_fine1 <- 1:nparticles
    ancestors_fine2 <- 1:nparticles
    
    if (coarsetimes[k+1]){
      if (meet_coarse){
        if (treestorage){
          Tree_coarse1$update(xparticles_coarse1, ancestors_coarse1 - 1) 
        } else {
          xtrajectory_coarse1[index_coarse+1, , ] <- xparticles_coarse1
          ancestries_coarse1[index_coarse, ] <- ancestors_coarse1
        }
      } else {
        if (treestorage){
          Tree_coarse1$update(xparticles_coarse1, ancestors_coarse1 - 1) 
          Tree_coarse2$update(xparticles_coarse2, ancestors_coarse2 - 1) 
        } else {
          xtrajectory_coarse1[index_coarse+1, , ] <- xparticles_coarse1
          xtrajectory_coarse2[index_coarse+1, , ] <- xparticles_coarse2
          ancestries_coarse1[index_coarse, ] <- ancestors_coarse1
          ancestries_coarse2[index_coarse, ] <- ancestors_coarse2
        }
      }
      ancestors_coarse1 <- 1:nparticles
      ancestors_coarse2 <- 1:nparticles
    }
    
    if (obstimes[k+1]){
      # compute weights
      index_obs <- index_obs + 1
      observation <- observations[index_obs, ] # 1 x ydimension 
      
      if (model$is_discrete_observation){
        # observation only depends on current particles in a discrete model
        if (meet_fine){
          logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, stepsize_fine[k], xparticles_fine1, observation)
          logweights_fine2 <- logweights_fine1
        } else {
          logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, stepsize_fine[k], xparticles_fine1, observation)
          logweights_fine2 <- logweights_fine2 + model$dmeasurement(theta, stepsize_fine[k], xparticles_fine2, observation)
        }
        
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory1 <- array(0, dim = c(k-last_obs_fine+1, xdimension, nparticles))
        x_sub_trajectory1[ , , ] <- xtrajectory_fine1[last_obs_fine:k, , ]
        
        if (meet_fine){
          logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, stepsize_fine[k], x_sub_trajectory1, observation)
          logweights_fine2 <- logweights_fine1
        } else {
          x_sub_trajectory2 <- array(0, dim = c(k-last_obs_fine+1, xdimension, nparticles))
          x_sub_trajectory2[ , , ] <- xtrajectory_fine2[last_obs_fine:k, , ]
          logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, stepsize_fine[k], x_sub_trajectory1, observation)
          logweights_fine2 <- logweights_fine2 + model$dmeasurement(theta, stepsize_fine[k], x_sub_trajectory2, observation)
        }
        
      }
      
      # index last observation
      last_obs_fine <- k + 1
      
      maxlogweights_fine1 <- max(logweights_fine1)
      maxlogweights_fine2 <- max(logweights_fine2)
      weights_fine1 <- exp(logweights_fine1 - maxlogweights_fine1)
      weights_fine2 <- exp(logweights_fine2 - maxlogweights_fine2)
      normweights_fine1 <- weights_fine1 / sum(weights_fine1)
      normweights_fine2 <- weights_fine2 / sum(weights_fine2)
      ess_fine1 <- 1 / sum(normweights_fine1^2)
      ess_fine2 <- 1 / sum(normweights_fine2^2)
      
      if (model$is_discrete_observation){
        # observation only depends on current particles in a discrete model
        if (meet_coarse){
          logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, stepsize_coarse[index_coarse], xparticles_coarse1, observation)
          logweights_coarse2 <- logweights_coarse1
        } else {
          logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, stepsize_coarse[index_coarse], xparticles_coarse1, observation)
          logweights_coarse2 <- logweights_coarse2 + model$dmeasurement(theta, stepsize_coarse[index_coarse], xparticles_coarse2, observation)
        }
        
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory1 <- array(0, dim = c(index_coarse-last_obs_coarse+1, xdimension, nparticles))
        x_sub_trajectory1[ , , ] <- xtrajectory_coarse1[last_obs_coarse:index_coarse, , ]
        
        if (meet_coarse){
          logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, stepsize_coarse[index_coarse], x_sub_trajectory1, observation)
          logweights_coarse2 <- logweights_coarse1
        } else {
          x_sub_trajectory2 <- array(0, dim = c(index_coarse-last_obs_coarse+1, xdimension, nparticles))
          x_sub_trajectory2[ , , ] <- xtrajectory_coarse2[last_obs_coarse:index_coarse, , ]
          logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, stepsize_coarse[index_coarse], x_sub_trajectory1, observation)
          logweights_coarse2 <- logweights_coarse2 + model$dmeasurement(theta, stepsize_coarse[index_coarse], x_sub_trajectory2, observation)
        }
        
      }
      
      # index last observation
      last_obs_coarse <- index_coarse + 1
      
      maxlogweights_coarse1 <- max(logweights_coarse1)
      maxlogweights_coarse2 <- max(logweights_coarse2)
      weights_coarse1 <- exp(logweights_coarse1 - maxlogweights_coarse1)
      weights_coarse2 <- exp(logweights_coarse2 - maxlogweights_coarse2)
      normweights_coarse1 <- weights_coarse1 / sum(weights_coarse1)
      normweights_coarse2 <- weights_coarse2 / sum(weights_coarse2)
      ess_coarse1 <- 1 / sum(normweights_coarse1^2)
      ess_coarse2 <- 1 / sum(normweights_coarse2^2)
      
      # resampling
      min_ess <- min(ess_fine1, ess_fine2, ess_coarse1, ess_coarse2)
      if (k < nsteps && min_ess < resampling_threshold * nparticles){
        rand <- runif(nparticles)
        ancestors <- coupled_resampling(normweights_coarse1, normweights_coarse2,
                                        normweights_fine1, normweights_fine2,
                                        nparticles, rand)
        ancestors_coarse1 <- ancestors[, 1]
        ancestors_coarse2 <- ancestors[, 2]
        ancestors_fine1 <- ancestors[, 3]
        ancestors_fine2 <- ancestors[, 4]
        
        # ancestor sampling
        as_logweights_fine1 <- logweights_fine1 + model$dtransition(theta, stepsize_fine[k+1], xparticles_fine1, ref_trajectory_fine1[, k+2])
        as_logweights_fine2 <- logweights_fine2 + model$dtransition(theta, stepsize_fine[k+1], xparticles_fine2, ref_trajectory_fine2[, k+2])
        as_maxlogweights_fine1 <- max(as_logweights_fine1)
        as_maxlogweights_fine2 <- max(as_logweights_fine2)
        as_weights_fine1 <- exp(as_logweights_fine1 - as_maxlogweights_fine1)
        as_weights_fine2 <- exp(as_logweights_fine2 - as_maxlogweights_fine2)
        as_normweights_fine1 <- as_weights_fine1 / sum(as_weights_fine1)
        as_normweights_fine2 <- as_weights_fine2 / sum(as_weights_fine2)
        
        as_logweights_coarse1 <- logweights_coarse1 + model$dtransition(theta, stepsize_coarse[index_coarse+1], xparticles_coarse1, ref_trajectory_coarse1[, index_coarse+2])
        as_logweights_coarse2 <- logweights_coarse2 + model$dtransition(theta, stepsize_coarse[index_coarse+1], xparticles_coarse2, ref_trajectory_coarse2[, index_coarse+2])
        as_maxlogweights_coarse1 <- max(as_logweights_coarse1)
        as_maxlogweights_coarse2 <- max(as_logweights_coarse2)
        as_weights_coarse1 <- exp(as_logweights_coarse1 - as_maxlogweights_coarse1)
        as_weights_coarse2 <- exp(as_logweights_coarse2 - as_maxlogweights_coarse2)
        as_normweights_coarse1 <- as_weights_coarse1 / sum(as_weights_coarse1)
        as_normweights_coarse2 <- as_weights_coarse2 / sum(as_weights_coarse2)
        
        as_ancestors <- coupled_resampling(as_normweights_coarse1, as_normweights_coarse2,
                                           as_normweights_fine1, as_normweights_fine2,
                                           1, rand[nparticles])
        
        ancestors_coarse1[nparticles] <- as_ancestors[, 1]
        ancestors_coarse2[nparticles] <- as_ancestors[, 2]
        ancestors_fine1[nparticles] <- as_ancestors[, 3]
        ancestors_fine2[nparticles] <- as_ancestors[, 4]
        
        xparticles_coarse1 <- xparticles_coarse1[, ancestors_coarse1]
        xparticles_coarse2 <- xparticles_coarse2[, ancestors_coarse2]
        xparticles_fine1 <- xparticles_fine1[, ancestors_fine1]
        xparticles_fine2 <- xparticles_fine2[, ancestors_fine2]
        
        # reset weights
        logweights_coarse1 <- rep(0, nparticles)
        logweights_coarse2 <- rep(0, nparticles)
        logweights_fine1 <- rep(0, nparticles)
        logweights_fine2 <- rep(0, nparticles)
      }
    }
  }
  
  # draw a pair of trajectories 
  rand <- runif(1)
  ancestor <- coupled_resampling(normweights_coarse1, normweights_coarse2, 
                                 normweights_fine1, normweights_fine2, 
                                 1, rand)
  ancestor_coarse1 <- ancestor[, 1]
  ancestor_coarse2 <- ancestor[, 2]
  ancestor_fine1 <- ancestor[, 3]
  ancestor_fine2 <- ancestor[, 4]

  if (meet_coarse){
    if (treestorage){
      new_trajectory_coarse1 <- Tree_coarse1$get_path(ancestor_coarse1 - 1)
    } else {
      new_trajectory_coarse1 <- get_path(model, discretization$coarse, xtrajectory_coarse1, ancestries_coarse1, ancestor_coarse1)
    }
    new_trajectory_coarse2 <- new_trajectory_coarse1
  } else {
    if (treestorage){
      new_trajectory_coarse1 <- Tree_coarse1$get_path(ancestor_coarse1 - 1)
      new_trajectory_coarse2 <- Tree_coarse2$get_path(ancestor_coarse2 - 1)
    } else {
      new_trajectory_coarse1 <- get_path(model, discretization$coarse, xtrajectory_coarse1, ancestries_coarse1, ancestor_coarse1)
      new_trajectory_coarse2 <- get_path(model, discretization$coarse, xtrajectory_coarse2, ancestries_coarse2, ancestor_coarse2)
    }
  }
  
  if (meet_fine){
    if (treestorage){
      new_trajectory_fine1 <- Tree_fine1$get_path(ancestor_fine1 - 1)
    } else {
      new_trajectory_fine1 <- get_path(model, discretization$fine, xtrajectory_fine1, ancestries_fine1, ancestor_fine1)
    }
    new_trajectory_fine2 <- new_trajectory_fine1
  } else {
    if (treestorage){
      new_trajectory_fine1 <- Tree_fine1$get_path(ancestor_fine1 - 1)
      new_trajectory_fine2 <- Tree_fine2$get_path(ancestor_fine2 - 1)
    } else {
      new_trajectory_fine1 <- get_path(model, discretization$fine, xtrajectory_fine1, ancestries_fine1, ancestor_fine1)
      new_trajectory_fine2 <- get_path(model, discretization$fine, xtrajectory_fine2, ancestries_fine2, ancestor_fine2)
    }
  }
  
  return(list(new_trajectory_coarse1 = new_trajectory_coarse1, new_trajectory_coarse2 = new_trajectory_coarse2, 
              new_trajectory_fine1 = new_trajectory_fine1, new_trajectory_fine2 = new_trajectory_fine2))
  
}