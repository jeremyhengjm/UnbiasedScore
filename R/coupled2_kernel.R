#' @rdname coupled2_kernel
#' @title Runs a 2-coupled Markov kernel 
#' @description Runs two coupled kernels that leaves smoothing distribution (at each discretization level) invariant. 
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization list containing stepsize, nsteps, statelength and obstimes
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 2-marginal coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param ref_trajectory1 a matrix of first reference trajectory, of size xdimension x statelength
#' @param ref_trajectory2 a matrix of second reference trajectory, of size xdimension x statelength
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
#' if missing, this function store all states and ancestors
#' @return a pair of new trajectories stored as matrices of size xdimension x statelength.
#' @export
coupled2_kernel <- function(model, theta, discretization, observations, nparticles, resampling_threshold = 1, coupled_resampling, 
                         ref_trajectory1, ref_trajectory2, algorithm = "CPF", treestorage = FALSE){
  
  if (algorithm == "CPF"){
    return(coupled2_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, ref_trajectory1, ref_trajectory2, treestorage))
  }
  
  if (algorithm == "CASPF"){
    return(coupled2_CASPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, ref_trajectory1, ref_trajectory2, treestorage))
  }
  
  if (algorithm == "CBSPF"){
    return(coupled2_CBSPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, ref_trajectory1, ref_trajectory2))
  }
  
}