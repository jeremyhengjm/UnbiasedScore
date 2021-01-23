#' @rdname coupled4_kernel
#' @title Runs a 4-coupled Markov kernel 
#' @description Runs four coupled kernels (two at each discretization level) that leaves smoothing distribution at corresponding level invariant
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param discretization lists containing stepsize, nsteps, statelength, obstimes for fine and coarse levels, 
#' and coarsetimes of length statelength_fine indexing time steps of coarse level
#' @param observations a matrix of observations, of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled_resampling a 4-marginal coupled resampling scheme, such as \code{\link{coupled4_maximalchains_maximallevels_independent_residuals}}
#' @param ref_trajectory_coarse1 a matrix of first reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_coarse2 a matrix of second reference trajectory for coarser discretization level, of size xdimension x statelength_coarse
#' @param ref_trajectory_fine1 a matrix of first reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param ref_trajectory_fine2 a matrix of second reference trajectory for finer discretization level, of size xdimension x statelength_fine
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param treestorage logical specifying tree storage of Jacob, Murray and Rubenthaler (2013);
#' if missing, this function store all states and ancestors
#' @return four new trajectories stored as matrices of size xdimension x statelength_coarse/fine
#' @export
coupled4_kernel <- function(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2, algorithm = "CPF", treestorage = FALSE){
  if (algorithm == "CPF"){
    return(coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                        ref_trajectory_coarse1, ref_trajectory_coarse2,
                        ref_trajectory_fine1, ref_trajectory_fine2, treestorage))
    }
  
  if (algorithm == "CASPF"){
    return(coupled4_CASPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                          ref_trajectory_coarse1, ref_trajectory_coarse2,
                          ref_trajectory_fine1, ref_trajectory_fine2, treestorage))
    }
  
  if (algorithm == "CBSPF"){
    return(coupled4_CBSPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                          ref_trajectory_coarse1, ref_trajectory_coarse2, 
                          ref_trajectory_fine1, ref_trajectory_fine2))
  }
  
}
