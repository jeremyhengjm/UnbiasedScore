#' @rdname simulate_SDE
#' @title Simulate time-discretized process following a stochastic differential equation
#' @description Simulate a trajectory following a stochastic differential equation using Euler-Maruyama discretization.
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps and statelength 
#'@return a matrix containing a new trajectory of size xdimension x statelength.
#'@export
simulate_SDE <- function(model, theta, discretization){
  
  # get model/problem settings 
  xdimension <- model$xdimension
  
  # discretization
  stepsize <- discretization$stepsize # vector of length nsteps
  statelength <- discretization$statelength
  nsteps <- discretization$nsteps
  
  # pre-allocate
  xtrajectory <- matrix(0, nrow = xdimension, ncol = statelength)

  # initialization
  xparticles <- model$rinit(1) # size: xdimension x 1
  xtrajectory[, 1] <- xparticles
  
  # Euler-Maruyama time-discretization
  for (k in 1:nsteps){
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension), nrow = xdimension, ncol = 1) # size: xdimension x 1
    xparticles <- model$rtransition(theta, stepsize[k], xparticles, randn) # size: xdimension x 1
    xtrajectory[, k+1] <- xparticles
  }
  
  return(list(new_trajectory = xtrajectory))
  
}