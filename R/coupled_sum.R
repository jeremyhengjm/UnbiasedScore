#' @rdname coupled_sum
#' @title Unbiased coupled-sum estimator of the gradient of the log-likelihood 
#' @description Estimates the gradient of the log-likelihood using the coupled-sum estimator of Rhee and Glynn (2015)
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
#' @param is_discrete_observation boolean indicating if observations are discrete
#' @return a list with objects such as 
#' random_level is the random level to truncated infinite sum 
#' unbiasedestimator is an unbiased estimator of the gradient of the log-likelihood
#' cost is the cost to compute the coupled-sum estimator
#' elapsedtime is the time taken to compute the coupled-sum estimator
#' @export
coupled_sum <- function(model, theta, observations, nparticles, resampling_threshold = 1, coupled2_resampling, coupled4_resampling, 
                        k = 0, m = 1, minimum_level, maximum_level, level_distribution, is_discrete_observation = T){
  
  # start timer
  tic()
  
  # for continuous observation process, compute observation sequence corresponding to the discretization level
  if (!is_discrete_observation){
    observations <- model$compute_observations(minimum_level)$observations
  }
  
  # initialize estimate by computing gradient at coarsest discretization level
  discretization <- model$construct_discretization(minimum_level)
  gradient <- unbiased_gradient(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, k = k, m = m, max_iterations = Inf)
  estimator <- gradient$unbiasedestimator
  cost <- nparticles * discretization$nsteps * gradient$cost 
  
  # sample random level to truncated infinite sum
  all_levels <- minimum_level:maximum_level
  random_level <- sample(x = all_levels, size = 1, prob = level_distribution$mass_function)
  
  # increment estimate by computing the difference of the gradient at two successive discretization levels
  if (random_level > minimum_level){
    for (level in (minimum_level+1):random_level){
      
      # for continuous observation process, compute observation sequence corresponding to the discretization level
      if (!is_discrete_observation){
        observations <- model$compute_observations(level)$observations
      }
      
      discretization <- model$construct_successive_discretization(level)
      gradient_increment <- unbiased_gradient_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled4_resampling, 
                                               k = k, m = m, max_iterations = Inf)
      cost <- cost + nparticles * discretization$coarse$nsteps * gradient_increment$cost_coarse
      cost <- cost + nparticles * discretization$fine$nsteps * gradient_increment$cost_fine
      increment <- gradient_increment$unbiasedestimator
      estimator <- estimator + increment / level_distribution$tail_function[level-minimum_level+1]
    }
  }
  
  # end timer and compute elapsed time
  timer <- toc(quiet = TRUE)
  elapsedtime <- timer$toc - timer$tic
  
  return(list(random_level = random_level, unbiasedestimator = estimator, 
              cost = cost, elapsedtime = elapsedtime))
  
}

#' @rdname compute_level_distribution
#' @title Compute distribution of levels 
#' @description Compute distribution of levels using results from Rhee and Glynn (2015) 
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param minimum_level coarsest discretization level 
#' @param maximum_level finest discretization level 
#' @param constant_sigma logical specifying if diffusion coefficient is constant
#' @return a list with objects such as 
#' mass_function is the probability mass function of the level distribution
#' tail_function is the tail probability function of the level distribution
#' @export
compute_level_distribution <- function(model, minimum_level, maximum_level, constant_sigma){
  
  # compute probability mass function of levels
  all_levels <- minimum_level:maximum_level
  num_levels <- length(all_levels)
  max_stepsize <- rep(0, num_levels)
  for (l in 1:num_levels){
    level <- all_levels[l]
    discretization <- model$construct_discretization(level)
    max_stepsize[l] <- max(discretization$stepsize)
  }
  
  pmf_levels <- rep(0, num_levels)
  if (constant_sigma){
    pmf_levels <- max_stepsize * (all_levels + 1) * log2(2 + all_levels)^2
  } else {
    pmf_levels <- sqrt(max_stepsize) * (all_levels + 1) * log2(2 + all_levels)^2
  }
  pmf_levels <- pmf_levels / sum(pmf_levels)
  pmf_tail <- rev(cumsum(rev(pmf_levels)))
  
  return(list(mass_function = pmf_levels, tail_function = pmf_tail))
}