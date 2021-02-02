#' @rdname simulate_multilevel_SDE
#' @title Simulate two time-discretized process following a stochastic differential equation
#' @description Simulate two trajectories following a stochastic differential equation using 
#' Euler-Maruyama at two successive discretization levels.
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' @return two new trajectories stored as matrices of size xdimension x statelength_coarse/fine.
#' @export
simulate_multilevel_SDE <- function(model, theta, discretization){
  
  # get model/problem settings 
  xdimension <- model$xdimension
  
  # discretization
  nsteps <- discretization$fine$nsteps
  stepsize_fine <- discretization$fine$stepsize
  stepsize_coarse <- discretization$coarse$stepsize
  statelength_fine <- discretization$fine$statelength
  statelength_coarse <- discretization$coarse$statelength
  coarsetimes <- discretization$coarsetimes # vector of length nsteps+1 indexing coarse times
  
  # pre-allocate
  xtrajectory_coarse <- matrix(0, nrow = xdimension, ncol = statelength_coarse)
  xtrajectory_fine <- matrix(0, nrow = xdimension, ncol = statelength_fine)
  
  # initialization
  xparticles_coarse <- model$rinit(1) # size: xdimension x 1
  xparticles_fine <- xparticles_coarse
  xtrajectory_coarse[, 1] <- xparticles_coarse
  xtrajectory_fine[, 1] <- xparticles_fine
  
  # Euler-Maruyama time-discretization
  index_coarse <- 0 # index coarse times
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * 1), nrow = xdimension, ncol = 1) # size: xdimension x 1
    xparticles_fine <- model$rtransition(theta, stepsize_fine[k], xparticles_fine, randn) # size: xdimension x 1
    xtrajectory_fine[, k+1] <- xparticles_fine
    
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
      xparticles_coarse <- model$rtransition(theta, stepsize_coarse[index_coarse], xparticles_coarse, combined_randn) # size: xdimension x 1
      xtrajectory_coarse[, index_coarse+1] <- xparticles_coarse
    }
    previous_randn <- randn
    
  }
  
  return(list(new_trajectory_coarse = xtrajectory_coarse, new_trajectory_fine = xtrajectory_fine))
  
}