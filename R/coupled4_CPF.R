#' @rdname coupled4_CPF
#' @title 4-wayCoupled Coupled Conditional Particle Filter
#' @description Runs four coupled conditional particle filters (two at each discretization level)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param level_fine level of the finer time discretization
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param coupled_resampling a 4-way coupled resampling scheme, such as \code{\link{coupled4_maximal_coupled_residuals}}
#' @param ref_trajectory_coarse1 a matrix of first reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_coarse2 a matrix of second reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_fine1 a matrix of first reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param ref_trajectory_fine2 a matrix of second reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @return four new trajectories stored as matrices of size xdimension x statelength_coarse/fine
#' @export
coupled4_CPF <- function(model, theta, level_fine, observations, nparticles, coupled_resampling,
                         ref_trajectory_coarse1, ref_trajectory_coarse2,
                         ref_trajectory_fine1, ref_trajectory_fine2){
  
  # get model/problem settings 
  statelength_fine <- model$statelength(level_fine)
  nsteps <- statelength_fine - 1
  nsteps_interval <- 2^level_fine
  xdimension <- model$xdimension
  ydimension <- model$ydimension

  # check if trajectories are equal
  meet_coarse <- all(ref_trajectory_coarse1 == ref_trajectory_coarse2)
  meet_fine <- all(ref_trajectory_fine1 == ref_trajectory_fine2)
  
  # create tree representation of the trajectories
  if (meet_coarse){
    Tree_coarse1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) # only store one tree in this case
  } else {
    Tree_coarse1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    Tree_coarse2 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
  }
  if (meet_fine){
    Tree_fine1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) # only store one tree in this case
  } else {
    Tree_fine1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    Tree_fine2 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
  }
  
  # initialization
  xparticles_coarse1 <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles_coarse2 <- model$rinit(nparticles) 
  xparticles_fine1 <- model$rinit(nparticles) 
  xparticles_fine2 <- model$rinit(nparticles) 
  xparticles_coarse1[, nparticles] <- ref_trajectory_coarse1[, 1]
  xparticles_coarse2[, nparticles] <- ref_trajectory_coarse2[, 1]
  xparticles_fine1[, nparticles] <- ref_trajectory_fine1[, 1]
  xparticles_fine2[, nparticles] <- ref_trajectory_fine2[, 1]
  
  # initialize tree to storage trajectories
  if (meet_coarse){
    Tree_coarse1$init(xparticles_coarse1) 
  } else {
    Tree_coarse1$init(xparticles_coarse1) 
    Tree_coarse2$init(xparticles_coarse2) 
  }
  if (meet_fine){
    Tree_fine1$init(xparticles_fine1) 
  } else {
    Tree_fine1$init(xparticles_fine1) 
    Tree_fine2$init(xparticles_fine2) 
  }
  
  # pre-allocate
  normweights_coarse1 <- rep(1 / nparticles, nparticles) # no need to resample at first step
  normweights_coarse2 <- rep(1 / nparticles, nparticles) 
  normweights_fine1 <- rep(1 / nparticles, nparticles) 
  normweights_fine2 <- rep(1 / nparticles, nparticles) 
  logweights_coarse1 <- rep(0, nparticles)
  logweights_coarse2 <- rep(0, nparticles)
  logweights_fine1 <- rep(0, nparticles)
  logweights_fine2 <- rep(0, nparticles)
  index_obs <- 0
  ancestors_coarse1 <- 1:nparticles
  ancestors_coarse2 <- 1:nparticles
  ancestors_fine1 <- 1:nparticles
  ancestors_fine2 <- 1:nparticles
  
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    if (meet_fine){
      xparticles_fine1 <- model$rtransition(theta, level_fine, xparticles_fine1, randn) # size: xdimension x nparticles
      xparticles_fine2 <- xparticles_fine1
    } else {
      xparticles_fine1 <- model$rtransition(theta, level_fine, xparticles_fine1, randn) # size: xdimension x nparticles
      xparticles_fine2 <- model$rtransition(theta, level_fine, xparticles_fine2, randn) # size: xdimension x nparticles
    }
    xparticles_fine1[, nparticles] <- ref_trajectory_fine1[, k+1]
    xparticles_fine2[, nparticles] <- ref_trajectory_fine2[, k+1]
    ancestors_fine1[nparticles] <- nparticles
    ancestors_fine2[nparticles] <- nparticles
    
    if (k %% 2 == 0){
      combined_randn <- (previous_randn + randn) / sqrt(2) 
      if (meet_coarse){
        xparticles_coarse1 <- model$rtransition(theta, level_fine-1, xparticles_coarse1, combined_randn) # size: xdimension x nparticles
        xparticles_coarse2 <- xparticles_coarse1
      } else {
        xparticles_coarse1 <- model$rtransition(theta, level_fine-1, xparticles_coarse1, combined_randn) # size: xdimension x nparticles
        xparticles_coarse2 <- model$rtransition(theta, level_fine-1, xparticles_coarse2, combined_randn) 
      }
      xparticles_coarse1[, nparticles] <- ref_trajectory_coarse1[, k/2+1]
      xparticles_coarse2[, nparticles] <- ref_trajectory_coarse2[, k/2+1]
      ancestors_coarse1[nparticles] <- nparticles
      ancestors_coarse2[nparticles] <- nparticles
    }
    previous_randn <- randn
    
    # update tree storage
    if (meet_fine){
      Tree_fine1$update(xparticles_fine1, ancestors_fine1 - 1) 
    } else {
      Tree_fine1$update(xparticles_fine1, ancestors_fine1 - 1) 
      Tree_fine2$update(xparticles_fine2, ancestors_fine2 - 1) 
    }
    ancestors_fine1 <- 1:nparticles
    ancestors_fine2 <- 1:nparticles
    
    if (k %% 2 == 0){
      if (meet_coarse){
        Tree_coarse1$update(xparticles_coarse1, ancestors_coarse1 - 1) 
      } else {
        Tree_coarse1$update(xparticles_coarse1, ancestors_coarse1 - 1) 
        Tree_coarse2$update(xparticles_coarse2, ancestors_coarse2 - 1) 
      }
      ancestors_coarse1 <- 1:nparticles
      ancestors_coarse2 <- 1:nparticles
    }
    
    if (k %% nsteps_interval == 0){
      # compute weights
      index_obs <- index_obs + 1
      observation <- observations[index_obs, ] # 1 x ydimension 
      if (meet_fine){
        logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, xparticles_fine1, observation)
        logweights_fine2 <- logweights_fine1
      } else {
        logweights_fine1 <- logweights_fine1 + model$dmeasurement(theta, xparticles_fine1, observation)
        logweights_fine2 <- logweights_fine2 + model$dmeasurement(theta, xparticles_fine2, observation)
      }
      maxlogweights_fine1 <- max(logweights_fine1)
      maxlogweights_fine2 <- max(logweights_fine2)
      weights_fine1 <- exp(logweights_fine1 - maxlogweights_fine1)
      weights_fine2 <- exp(logweights_fine2 - maxlogweights_fine2)
      normweights_fine1 <- weights_fine1 / sum(weights_fine1)
      normweights_fine2 <- weights_fine2 / sum(weights_fine2)
      
      if (meet_coarse){
        logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, xparticles_coarse1, observation)
        logweights_coarse2 <- logweights_coarse1
      } else {
        logweights_coarse1 <- logweights_coarse1 + model$dmeasurement(theta, xparticles_coarse1, observation)
        logweights_coarse2 <- logweights_coarse2 + model$dmeasurement(theta, xparticles_coarse2, observation)
      }
      maxlogweights_coarse1 <- max(logweights_coarse1)
      maxlogweights_coarse2 <- max(logweights_coarse2)
      weights_coarse1 <- exp(logweights_coarse1 - maxlogweights_coarse1)
      weights_coarse2 <- exp(logweights_coarse2 - maxlogweights_coarse2)
      normweights_coarse1 <- weights_coarse1 / sum(weights_coarse1)
      normweights_coarse2 <- weights_coarse2 / sum(weights_coarse2)
      
      # resampling
      if (k < nsteps){
        rand <- runif(nparticles)
        ancestors <- coupled_resampling(normweights_coarse1, normweights_coarse2,
                                        normweights_fine1, normweights_fine2,
                                        nparticles, rand)
        ancestors_coarse1 <- ancestors[, 1]
        ancestors_coarse2 <- ancestors[, 2]
        ancestors_fine1 <- ancestors[, 3]
        ancestors_fine2 <- ancestors[, 4]
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
  
  # draw a pair of trajectories using multinomial sampling with common random variables 
  rand <- runif(1)
  ancestor <- coupled_resampling(normweights_coarse1, normweights_coarse2, 
                                 normweights_fine1, normweights_fine2, 
                                 1, rand)
  ancestor_coarse1 <- ancestor[, 1]
  ancestor_coarse2 <- ancestor[, 2]
  ancestor_fine1 <- ancestor[, 3]
  ancestor_fine2 <- ancestor[, 4]

  if (meet_coarse){
    new_trajectory_coarse1 <- Tree_coarse1$get_path(ancestor_coarse1 - 1)
    new_trajectory_coarse2 <- new_trajectory_coarse1
  } else {
    new_trajectory_coarse1 <- Tree_coarse1$get_path(ancestor_coarse1 - 1)
    new_trajectory_coarse2 <- Tree_coarse2$get_path(ancestor_coarse2 - 1)
  }
  if (meet_fine){
    new_trajectory_fine1 <- Tree_fine1$get_path(ancestor_fine1 - 1)
    new_trajectory_fine2 <- new_trajectory_fine1
  } else {
    new_trajectory_fine1 <- Tree_fine1$get_path(ancestor_fine1 - 1)
    new_trajectory_fine2 <- Tree_fine2$get_path(ancestor_fine2 - 1)
  }
  
  return(list(new_trajectory_coarse1 = new_trajectory_coarse1, new_trajectory_coarse2 = new_trajectory_coarse2, 
              new_trajectory_fine1 = new_trajectory_fine1, new_trajectory_fine2 = new_trajectory_fine2))
  
}