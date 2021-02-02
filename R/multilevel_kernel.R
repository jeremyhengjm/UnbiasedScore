#' @rdname multilevel_kernel
#' @title Runs a multilevel Markov kernel 
#' @description Runs two coupled kernels that leaves the corresponding smoothing distribution (at each discretization level) invariant. 
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param ref_trajectory_coarse a matrix of reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_fine a matrix of reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
#' if missing, this function store all states and ancestors
#' @return two new trajectories stored as matrices of size xdimension x statelength_coarse/fine.
#' @export
multilevel_kernel <- function(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                               ref_trajectory_coarse = NULL, ref_trajectory_fine = NULL, algorithm = "CPF", treestorage = FALSE){
  
  if (algorithm == "CPF"){
    return(multilevel_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                          ref_trajectory_coarse, ref_trajectory_fine, treestorage))
  }
  
  if (algorithm == "CASPF"){
    return(multilevel_CASPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse, ref_trajectory_fine, treestorage))
  }
  
  if (algorithm == "CBSPF"){
    return(multilevel_CBSPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse, ref_trajectory_fine))
  }
  
}