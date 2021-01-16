#' @rdname coupled2_CBSPF
#' @title 2-way Coupled Conditional Particle Filter with Backward Sampling
#' @description Runs two coupled conditional particle filters (at each discretization level)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param ref_trajectory1 a matrix of first reference trajectory, of size xdimension x statelength
#' @param ref_trajectory2 a matrix of second reference trajectory, of size xdimension x statelength
#' @return a pair of new trajectories stored as matrices of size xdimension x statelength
#' @export
coupled2_CBSPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1, coupled_resampling, 
                         ref_trajectory1, ref_trajectory2){
  
  # get model/problem settings 
  nobservations <- nrow(observations)
  xdimension <- model$xdimension
  ydimension <- model$ydimension
  
  # discretization
  stepsize <- discretization$stepsize # vector of length nsteps
  statelength <- discretization$statelength
  nsteps <- discretization$nsteps
  obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1  

  # check if chains have met
  meet <- all(ref_trajectory1 == ref_trajectory2)
  
  # create tree representation of the trajectories or store all states and ancestors
  xtrajectory1 <- array(0, dim = c(statelength, xdimension, nparticles))
  xtrajectory2 <- array(0, dim = c(statelength, xdimension, nparticles))
  ancestries1 <- matrix(0, nrow = statelength, ncol = nparticles)
  ancestries2 <- matrix(0, nrow = statelength, ncol = nparticles)
  store_logweights1 <- matrix(0, nrow = nobservations, ncol = nparticles)
  store_logweights2 <- matrix(0, nrow = nobservations, ncol = nparticles)
  
  # initialization
  xparticles1 <- model$rinit(nparticles) # size: xdimension x nparticles
  xparticles2 <- xparticles1
  xparticles1[, nparticles] <- ref_trajectory1[, 1]
  xparticles2[, nparticles] <- ref_trajectory2[, 1]
  
  # initialize tree to storage trajectories
  xtrajectory1[1, , ] <- xparticles1
  xtrajectory2[1, , ] <- xparticles2
  logweights1 <- rep(0, nparticles)
  logweights2 <- rep(0, nparticles)
  index_obs <- 0
  ancestors1 <- 1:nparticles
  ancestors2 <- 1:nparticles
  
  # index last observation
  last_obs <- 1
  
  # random initialization
  if (obstimes[1]){
    # compute weights
    index_obs <- index_obs + 1
    observation <- observations[index_obs, ] # 1 x ydimension 
    if (meet){
      logweights1 <- model$dmeasurement(theta, stepsize[1], xparticles1, observation)
      logweights2 <- logweights1
      store_logweights1[index_obs, ] <- logweights1
      store_logweights2[index_obs, ] <- logweights2
    } else {
      logweights1 <- model$dmeasurement(theta, stepsize[1], xparticles1, observation)
      logweights2 <- model$dmeasurement(theta, stepsize[1], xparticles2, observation)
      store_logweights1[index_obs, ] <- logweights1
      store_logweights2[index_obs, ] <- logweights2
    }
    maxlogweights1 <- max(logweights1)
    maxlogweights2 <- max(logweights2)
    weights1 <- exp(logweights1 - maxlogweights1)
    weights2 <- exp(logweights2 - maxlogweights2)
    normweights1 <- weights1 / sum(weights1)
    normweights2 <- weights2 / sum(weights2)
    ess1 <- 1 / sum(normweights1^2)
    ess2 <- 1 / sum(normweights2^2)
    
    # resampling
    if (min(ess1, ess2) < resampling_threshold * nparticles){
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
  
  
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    if (meet){
      xparticles1 <- model$rtransition(theta, stepsize[k], xparticles1, randn) # size: xdimension x nparticles
      xparticles2 <- xparticles1
    } else {
      xparticles1 <- model$rtransition(theta, stepsize[k], xparticles1, randn) 
      xparticles2 <- model$rtransition(theta, stepsize[k], xparticles2, randn) # size: xdimension x nparticles
    }
    xparticles1[, nparticles] <- ref_trajectory1[, k+1]
    xparticles2[, nparticles] <- ref_trajectory2[, k+1]
    ancestors1[nparticles] <- nparticles
    ancestors2[nparticles] <- nparticles
    
    # update tree storage
    xtrajectory1[k+1, , ] <- xparticles1
    xtrajectory2[k+1, , ] <- xparticles2
    ancestries1[k, ] <- ancestors1
    ancestries2[k, ] <- ancestors2
    ancestors1 <- 1:nparticles
    ancestors2 <- 1:nparticles
    
    if (obstimes[k+1]){
      # compute weights
      index_obs <- index_obs + 1 
      observation <- observations[index_obs, ] # 1 x ydimension 
      
      if (model$is_discrete_observation){
        # observation only depends on current particles in a discrete model
        if (meet){
          logweights1 <- logweights1 + model$dmeasurement(theta, stepsize[k], xparticles1, observation)
          logweights2 <- logweights1
        } else {
          logweights1 <- logweights1 + model$dmeasurement(theta, stepsize[k], xparticles1, observation)
          logweights2 <- logweights2 + model$dmeasurement(theta, stepsize[k], xparticles2, observation)
        }
        
      } else {
        # observation depends on inter-observation states
        x_sub_trajectory1 <- array(0, dim = c(k-last_obs+1, xdimension, nparticles))
        x_sub_trajectory1[ , , ] <- xtrajectory1[last_obs:k, , ]
        
        if (meet){
          logweights1 <- logweights1 + model$dmeasurement(theta, stepsize[k], x_sub_trajectory1, observation)
          logweights2 <- logweights1
        } else {
          x_sub_trajectory2 <- array(0, dim = c(k-last_obs+1, xdimension, nparticles))
          x_sub_trajectory2[ , , ] <- xtrajectory2[last_obs:k, , ]
          logweights1 <- logweights1 + model$dmeasurement(theta, stepsize[k], x_sub_trajectory1, observation)
          logweights2 <- logweights2 + model$dmeasurement(theta, stepsize[k], x_sub_trajectory2, observation)
        }
      }
      
      # index last observation
      last_obs <- k+1
      
      store_logweights1[index_obs, ] <- logweights1
      store_logweights2[index_obs, ] <- logweights2
      
      maxlogweights1 <- max(logweights1)
      maxlogweights2 <- max(logweights2)
      weights1 <- exp(logweights1 - maxlogweights1)
      weights2 <- exp(logweights2 - maxlogweights2)
      normweights1 <- weights1 / sum(weights1)
      normweights2 <- weights2 / sum(weights2)
      ess1 <- 1 / sum(normweights1^2)
      ess2 <- 1 / sum(normweights2^2)
      
      # resampling
      if (k < nsteps && min(ess1, ess2) < resampling_threshold * nparticles){
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
  # ancestor <- coupled_resampling(normweights1, normweights2, 1, runif(1))
  # ancestor1 <- ancestor[, 1]
  # ancestor2 <- ancestor[, 2]
  
  # if (meet){
  #   new_trajectory1 <- get_path_bs(model, discretization, xtrajectory1, ancestries1, ancestor1, store_logweights1)
  #   new_trajectory2 <- new_trajectory1
  # } else {
  #   new_trajectory1 <- get_path_bs(model, discretization, xtrajectory1, ancestries1, ancestor1, store_logweights1)
  #   new_trajectory2 <- get_path_bs(model, discretization, xtrajectory2, ancestries2, ancestor2, store_logweights2)
  # }
  
  
  # trace ancestry to get a path
  new_trajectory1 <- matrix(0, nrow = xdimension, ncol = statelength)
  new_trajectory2 <- matrix(0, nrow = xdimension, ncol = statelength)
  particle_index <- coupled_resampling(normweights1, normweights2, 1, runif(1))
  particle_index1 <- particle_index[, 1]
  particle_index2 <- particle_index[, 2]
  new_trajectory1[, statelength] <- xtrajectory1[statelength, , particle_index1]
  new_trajectory2[, statelength] <- xtrajectory2[statelength, , particle_index2]
  index_obs <- sum(obstimes)

  # last observation has already been taken into account if it is at the last time step
  if (obstimes[statelength]){
    index_obs <- index_obs - 1
  }

  for (k in nsteps:1){
    if (obstimes[k]){
      logweights1 <- store_logweights1[index_obs, ]
      bs_logweights1 <- logweights1 + model$dtransition(theta, stepsize[k], xtrajectory1[k, , ], ref_trajectory1[, k+1])
      bs_maxlogweights1 <- max(bs_logweights1)
      bs_weights1 <- exp(bs_logweights1 - bs_maxlogweights1)
      bs_normweights1 <- bs_weights1 / sum(bs_weights1)

      logweights2 <- store_logweights2[index_obs, ]
      bs_logweights2 <- logweights2 + model$dtransition(theta, stepsize[k], xtrajectory2[k, , ], ref_trajectory2[, k+1])
      bs_maxlogweights2 <- max(bs_logweights2)
      bs_weights2 <- exp(bs_logweights2 - bs_maxlogweights2)
      bs_normweights2 <- bs_weights2 / sum(bs_weights2)

      index_obs <- index_obs - 1

      particle_index <- coupled_resampling(bs_normweights1, bs_normweights2, 1, runif(1))
      particle_index1 <- particle_index[, 1]
      particle_index2 <- particle_index[, 2]
    } else {
      particle_index1 <- ancestries1[k, particle_index1]
      particle_index2 <- ancestries2[k, particle_index2]
    }
    new_trajectory1[, k] <- xtrajectory1[k, , particle_index1]
    new_trajectory2[, k] <- xtrajectory2[k, , particle_index2]
  }
  

  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2))
  
}