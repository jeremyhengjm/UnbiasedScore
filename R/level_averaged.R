#' @rdname level_averaged
#' @title Unbiased level-averaged estimator of the gradient of the log-likelihood 
#' @description Estimates the gradient of the log-likelihood using a level-averaged estimator
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param observations a matrix of observations of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled2_resampling a 2-way coupled resampling scheme, such as \code{\link{coupled2_maximal_coupled_residuals}}
#' @param coupled4_resampling a 4-way coupled resampling scheme, such as \code{\link{coupled4_maximal_coupled_residuals}}
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param minimum_level coarsest discretization level 
#' @param maximum_level finest discretization level 
#' @param level_distribution list containing mass_function and tail_function that specify the distribution of levels, 
#' e.g. by calling \code{\link{compute_level_distribution}} 
#' @param lower_level lower level of multilevel estimator
#' @param upper_level upper level of multilevel estimator
#' @return a list with objects such as 
#' random_level is the random level to truncated infinite sum 
#' multilevel_estimator is a multilevel estimator of the gradient of the log-likelihood
#' unbiasedestimator is an unbiased estimator of the gradient of the log-likelihood
#' @export
level_averaged <- function(model, theta, observations, nparticles, resampling_threshold = 1, coupled2_resampling, coupled4_resampling, 
                           k = 0, m = 1, minimum_level, maximum_level, level_distribution, lower_level, upper_level){
  
  # sample random level to truncated infinite sum
  all_levels <- minimum_level:maximum_level
  random_level <- sample(x = all_levels, size = 1, prob = level_distribution$mass_function)
  
  # preallocate
  theta_dimension <- model$theta_dimension
  multilevel_estimator <- rep(0, theta_dimension)
  correction <- rep(0, theta_dimension)
  
  if (random_level > lower_level){
    # bias correction term is non-zero
    
    for (level in (lower_level+1):random_level){
      # unbiased estimation of gradient increment
      discretization <- model$construct_successive_discretization(level)
      increment <- unbiased_gradient_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled4_resampling, 
                                               k = k, m = m, max_iterations = Inf)
      if (level-1 <= upper_level){
        # update multilevel estimator
        gradient_estimator <- increment$unbiasedestimator_coarse
        multilevel_estimator <- multilevel_estimator + gradient_estimator / (upper_level - lower_level + 1)
      }
      
      # update bias correction term
      gradient_increment <- increment$unbiasedestimator
      correction <- correction + gradient_increment * min( 1, (level - lower_level) / (upper_level - lower_level + 1) ) / 
        level_distribution$tail_function[level-minimum_level+1] 
    }
    
    # compute additional gradient terms in multilevel estimator if necessary
    if (random_level <= upper_level){
      for (level in random_level:upper_level){
        discretization <- model$construct_discretization(level)
        estimator <- unbiased_gradient(model, theta, discretization, observations, nparticles, resampling_threshold, 
                                       coupled2_resampling, k = k, m = m, max_iterations = Inf)
        gradient_estimator <- estimator$unbiasedestimator
        multilevel_estimator <- multilevel_estimator + gradient_estimator / (upper_level - lower_level + 1)
      }
    }
    
  } else {
    # bias correction term is zero, return multilevel estimator
    for (level in lower_level:upper_level){
      discretization <- model$construct_discretization(level)
      estimator <- unbiased_gradient(model, theta, discretization, observations, nparticles, resampling_threshold, 
                                     coupled2_resampling, k = k, m = m, max_iterations = Inf)
      gradient_estimator <- estimator$unbiasedestimator
      multilevel_estimator <- multilevel_estimator + gradient_estimator / (upper_level - lower_level + 1)
    }
  }
  
  # compute unbiased gradient estimator
  unbiasedestimator <- multilevel_estimator + correction
  
  return(list(random_level = random_level, multilevel_estimator = multilevel_estimator, unbiasedestimator = unbiasedestimator))
  
}

