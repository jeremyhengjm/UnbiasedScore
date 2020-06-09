#' @rdname coupled2_CPF
#' @title 2-way Coupled Conditional Particle Filter
#' @description Runs two coupled conditional particle filters (at each discretization level)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param level level of the time discretization
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param coupled_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_coupled_residuals}}
#' @param ref_trajectory1 a matrix of first reference trajectory, of size xdimension x statelength
#' @param ref_trajectory2 a matrix of second reference trajectory, of size xdimension x statelength
#' @return a pair of new trajectories stored as matrices of size xdimension x statelength
#' @export
coupled2_CPF <- function(model, theta, level, observations, nparticles, coupled_resampling, 
                         ref_trajectory1, ref_trajectory2){
  
  # get model/problem settings
  statelength <- model$statelength(level)
  nsteps <- statelength - 1 
  nsteps_interval <- 2^level
  xdimension <- model$xdimension
  ydimension <- model$ydimension

  # check if chains have met
  meet <- all(ref_trajectory1 == ref_trajectory2)
  
  # create tree representation of the trajectories
  if (meet){
    Tree1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension) # only store one tree in this case
  } else {
    Tree1 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
    Tree2 <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
  }
  
  # initialization
  xparticles1 <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles2 <- model$rinit(nparticles) 
  xparticles1[, nparticles] <- ref_trajectory1[, 1]
  xparticles2[, nparticles] <- ref_trajectory2[, 1]
  
  # initialize tree to storage trajectories
  if (meet){
    Tree1$init(xparticles1) 
  } else {
    Tree1$init(xparticles1)
    Tree2$init(xparticles2) 
  }
  
  # pre-allocate
  normweights1 <- rep(1 / nparticles, nparticles) # no need to resample at first step
  normweights2 <- rep(1 / nparticles, nparticles)
  logweights1 <- rep(0, nparticles)
  logweights2 <- rep(0, nparticles)
  index_obs <- 0
  ancestors1 <- 1:nparticles
  ancestors2 <- 1:nparticles
  
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    if (meet){
      xparticles1 <- model$rtransition(theta, level, xparticles1, randn) # size: xdimension x nparticles
      xparticles2 <- xparticles1
    } else {
      xparticles1 <- model$rtransition(theta, level, xparticles1, randn) 
      xparticles2 <- model$rtransition(theta, level, xparticles2, randn) # size: xdimension x nparticles
    }
    xparticles1[, nparticles] <- ref_trajectory1[, k+1]
    xparticles2[, nparticles] <- ref_trajectory2[, k+1]
    ancestors1[nparticles] <- nparticles
    ancestors2[nparticles] <- nparticles
    
    # update tree storage
    if (meet){
      Tree1$update(xparticles1, ancestors1 - 1)    
    } else {
      Tree1$update(xparticles1, ancestors1 - 1)    
      Tree2$update(xparticles2, ancestors2 - 1)    
    }
    ancestors1 <- 1:nparticles
    ancestors2 <- 1:nparticles
    
    if (k %% nsteps_interval == 0){
      # compute weights
      index_obs <- index_obs + 1 
      observation <- observations[index_obs, ] # 1 x ydimension 
      if (meet){
        logweights1 <- logweights1 + model$dmeasurement(theta, xparticles1, observation)
        logweights2 <- logweights1
      } else {
        logweights1 <- logweights1 + model$dmeasurement(theta, xparticles1, observation)
        logweights2 <- logweights2 + model$dmeasurement(theta, xparticles2, observation)
      }
      maxlogweights1 <- max(logweights1)
      maxlogweights2 <- max(logweights2)
      weights1 <- exp(logweights1 - maxlogweights1)
      weights2 <- exp(logweights2 - maxlogweights2)
      normweights1 <- weights1 / sum(weights1)
      normweights2 <- weights2 / sum(weights2)
      
      # resampling
      if (k < nsteps){
        rand <- runif(nparticles)
        ancestors <- coupled_resampling(normweights1, normweights2, nparticles, rand)
        ancestors1 <- ancestors[, 1]
        ancestors2 <- ancestors[, 2]
        xparticles1 <- xparticles1[, ancestors1]
        xparticles2 <- xparticles2[, ancestors2]
        logweights1 <- rep(0, nparticles)
        logweights2 <- rep(0, nparticles)
      }
    }
  }
  
  # draw a pair of trajectories 
  rand <- runif(1)
  ancestor <- coupled_resampling(normweights1, normweights2, 1, rand)
  ancestor1 <- ancestor[, 1]
  ancestor2 <- ancestor[, 2]
  
  if (meet){
    new_trajectory1 <- Tree1$get_path(ancestor1 - 1)
    new_trajectory2 <- new_trajectory1
  } else {
    new_trajectory1 <- Tree1$get_path(ancestor1 - 1)
    new_trajectory2 <- Tree2$get_path(ancestor2 - 1)
  }

  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2))
  
}