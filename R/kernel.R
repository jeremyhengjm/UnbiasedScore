#' @rdname kernel
#' @title Markov kernel
#' @description Runs a Markov kernel that leaves smoothing distribution invariant.
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param ref_trajectory a matrix of reference trajectory, of size xdimension x statelength; 
#' if missing, this function runs a standard particle filter
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013); 
#' if missing, this function store all states and ancestors
#'@return a matrix containing a new trajectory of size xdimension x statelength.
#'@export
kernel <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1,
                   ref_trajectory = NULL, algorithm = "CPF", treestorage = FALSE){
  if (algorithm == "CPF"){
    return(CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory, treestorage))
  }
  
  if (algorithm == "CASPF"){
    return(CASPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory, treestorage))
  }
  
  if (algorithm == "CBSPF"){
    return(CBSPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory))
  }
                     
}
