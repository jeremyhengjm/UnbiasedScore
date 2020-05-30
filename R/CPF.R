#' @rdname CPF
#' @title Conditional Particle Filter
#' @description Runs a conditional particle filter
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param level level of the time discretization
#' @param observations a matrix of observations of size terminal_time x ydimension
#' @param nparticles number of particles
#' @param ref_trajectory a matrix of reference trajectory, of size xdimension x statelength; 
#'if missing, this function runs a standard particle filter
#'@return a matrix containing a new trajectory of size xdimension x statelength
#'@export
CPF <- function(model, theta, level, observations, nparticles, ref_trajectory = NULL){
  
  # get model/problem settings 
  statelength <- model$statelength(level)
  nsteps <- statelength - 1 
  xdimension <- model$xdimension
  ydimension <- model$ydimension
  obs_times <- model$obstimes(level)
  
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
  
  # initialization
  xparticles <- model$rinit(nparticles) # size: xdimension x nparticles
  if (!is.null(ref_trajectory)){
    xparticles[, nparticles] <- ref_trajectory[, 1]
  }
  Tree$init(xparticles) # tree to storage trajectories
  normweights <- rep(1 / nparticles, nparticles) # no need to resample at first step
  logweights <- rep(0, nparticles)
  ess <- rep(nparticles, nsteps)
  index_obs <- 0
  ancestors <- 1:nparticles
  log_normconst <- rep(0, nsteps)
  log_ratio_normconst <- 0
  
  for (k in 1:nsteps){
    
    # propagate under latent dynamics
    randn <- matrix(rnorm(xdimension * nparticles), nrow = xdimension, ncol = nparticles) # size: xdimension x nparticles
    xparticles <- model$rtransition(theta, level, xparticles, randn) # size: xdimension x nparticles
    if (!is.null(ref_trajectory)){
      xparticles[, nparticles] <- ref_trajectory[, k+1]
      ancestors[nparticles] <- nparticles
    }
    
    # update tree storage
    Tree$update(xparticles, ancestors - 1)    
    ancestors <- 1:nparticles
    
    if (obs_times[k+1]){
      # compute weights
      index_obs <- index_obs + 1 
      observation <- observations[index_obs, ] # 1 x ydimension 
      logweights <- logweights + model$dmeasurement(theta, xparticles, observation)
      maxlogweights <- max(logweights)
      weights <- exp(logweights - maxlogweights)
      normweights <- weights / sum(weights)
      ess[k] <- 1 / sum(normweights^2)
      
      # compute normalizing constant
      log_ratio_normconst <- log_ratio_normconst + log(mean(weights)) + maxlogweights  
      log_normconst[k] <- log_ratio_normconst
      
      # resampling
      if (k < nsteps){
        rand <- runif(nparticles)
        ancestors <- multinomial_resampling(normweights, nparticles, rand) 
        xparticles <- xparticles[, ancestors]
        logweights <- rep(0, nparticles)
      }
    }
  }
  
  # draw a trajectory
  ancestor <- multinomial_resampling(normweights, 1, runif(1))
  new_trajectory <- Tree$get_path(ancestor - 1)
  return(list(new_trajectory = new_trajectory, ess = ess, log_normconst = log_normconst))
  
}