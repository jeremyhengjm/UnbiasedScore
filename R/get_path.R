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

  