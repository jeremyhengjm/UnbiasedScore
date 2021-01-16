#' @rdname get_path
#' @title Get a path from the output of a particle filter
#' @description Get a path from the output of a particle filter by tracing an ancestry
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param xtrajectory an array of particle states, of size statelength x xdimension x nparticles
#' @param ancestries a matrix of ancestor indexes, of size statelength x nparticles
#' @param ancestor an ancestor index at the terminal step
#'@return a matrix containing a new trajectory of size xdimension x statelength
#'@export
get_path <- function(model, discretization, xtrajectory, ancestries, ancestor){
  # get model/problem and discretization settings
  xdimension <- model$xdimension
  statelength <- discretization$statelength
  nsteps <- discretization$nsteps
  
  # trace ancestry to get a path
  new_trajectory <- matrix(0, nrow = xdimension, ncol = statelength)
  particle_index <- ancestor
  new_trajectory[, statelength] <- xtrajectory[statelength, , particle_index]
  for (k in nsteps:1){
    particle_index <- ancestries[k, particle_index]
    new_trajectory[, k] <- xtrajectory[k, , particle_index]
  }
  
  return(new_trajectory)
}

#' @rdname get_path_bs
#' @title Get a path from the output of a particle filter using backward sampling (Whiteley, 2010)
#' @description Get a path from the output of a particle filter by sampling an ancestry
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param xtrajectory an array of particle states, of size statelength x xdimension x nparticles
#' @param ancestries a matrix of ancestor indexes, of size statelength x nparticles
#' @param ancestor an ancestor index at the terminal step
#' @param store_logweights a matrix of log-weights, of size nobservations x nparticles
#'@return a matrix containing a new trajectory of size xdimension x statelength
#'@export
get_path_bs <- function(model, discretization, xtrajectory, ancestries, ancestor, store_logweights){
  # get model/problem and discretization settings
  xdimension <- model$xdimension
  statelength <- discretization$statelength
  nsteps <- discretization$nsteps
  obstimes <- discretization$obstimes # vector of length statelength = nsteps + 1
  
  # trace ancestry to get a path
  new_trajectory <- matrix(0, nrow = xdimension, ncol = statelength)
  particle_index <- ancestor
  new_trajectory[, statelength] <- xtrajectory[statelength, , particle_index]
  
  # last observation has already been taken into account if it is at the last time step
  if (obstimes[statelength]){
    index_obs <- sum(obstimes) - 1
  } else {
    index_obs <- sum(obstimes)
  }
  
  for (k in nsteps:1){
    if (obstimes[k]){
      logweights <- store_logweights[index_obs, ]
      bs_logweights <- logweights + model$dtransition(theta, stepsize[k], xtrajectory[k, , ], ref_trajectory[, k+1])
      bs_maxlogweights <- max(bs_logweights)
      bs_weights <- exp(bs_logweights - bs_maxlogweights)
      bs_normweights <- bs_weights / sum(bs_weights)
      index_obs <- index_obs - 1
      
      particle_index <- multinomial_resampling(bs_normweights, 1, runif(1))
    } else {
      particle_index <- ancestries[k, particle_index]
    }
    new_trajectory[, k] <- xtrajectory[k, , particle_index]
  }
  
  return(new_trajectory)
}
