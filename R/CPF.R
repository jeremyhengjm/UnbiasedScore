#' @rdname CPF
#' @title Conditional Particle Filter
#' @description Runs a conditional particle filter
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param ref_trajectory a matrix of reference trajectory, of size xdimension x statelength; 
#' if missing, this function runs a standard particle filter
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013); 
#' if missing, this function store all states and ancestors
#'@return a matrix containing a new trajectory of size xdimension x statelength
#'@export
CPF <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1,
                ref_trajectory = NULL, treestorage = FALSE){
  
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
  if (treestorage){
    Tree <- new(TreeClass, nparticles, 10*nparticles*xdimension, xdimension)
  } else {
    xtrajectory <- array(0, dim = c(statelength, xdimension, nparticles))
    ancestries <- matrix(0, nrow = statelength, ncol = nparticles)
  }
  
  # initialization
  xparticles <- model$rinit(nparticles) # size: xdimension x nparticles
  if (!is.null(ref_trajectory)){
    xparticles[, nparticles] <- ref_trajectory[, 1]
  }
  if (treestorage){
    Tree$init(xparticles) # tree to storage trajectories
  } else {
    xtrajectory[1, , ] <- xparticles
  }
  logweights <- rep(0, nparticles)
  ess <- rep(nparticles, statelength)
  log_ratio_normconst <- 0
  log_normconst_previous <- 0 # log normalizing constant at previous resampling time
  log_normconst <- rep(0, nobservations)
  index_obs <- 0 
  ancestors <- 1:nparticles
  
  # random initialization
  if (obstimes[1]){
    # compute weights
    index_obs <- index_obs + 1 
    observation <- observations[index_obs, ] # 1 x ydimension 
    logweights <- model$dmeasurement(theta, stepsize[1], xparticles, observation)
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
    if (!is.null(ref_trajectory)){
      xparticles[, nparticles] <- ref_trajectory[, k+1]
      ancestors[nparticles] <- nparticles
    }
    
    # update tree storage
    if (treestorage){
      Tree$update(xparticles, ancestors - 1)    
    } else {
      xtrajectory[k+1, , ] <- xparticles
      ancestries[k, ] <- ancestors
    }
    ancestors <- 1:nparticles
    
    if (obstimes[k+1]){
      # compute weights
      index_obs <- index_obs + 1 
      observation <- observations[index_obs, ] # 1 x ydimension 
      logweights <- logweights + model$dmeasurement(theta, stepsize[k], xparticles, observation)
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
  ancestor <- multinomial_resampling(normweights, 1, runif(1))
  if (treestorage){
    new_trajectory <- Tree$get_path(ancestor - 1)
  } else {
    new_trajectory <- get_path(model, discretization, xtrajectory, ancestries, ancestor)
  }
  return(list(new_trajectory = new_trajectory, ess = ess, log_normconst = log_normconst))
  
}