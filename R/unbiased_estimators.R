#' @rdname compute_level_distribution
#' @title Compute distribution of levels 
#' @description Compute distribution of levels using results from Section 4 of Rhee and Glynn (2015) 
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param minimum_level coarsest discretization level 
#' @param maximum_level finest discretization level 
#' @return a list with objects such as 
#' support of the level distribution
#' mass_function is the probability mass function of the level distribution
#' tail_function is the tail probability function of the level distribution
#' @export
compute_level_distribution <- function(model, minimum_level, maximum_level){
  
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
  if (model$constant_sigma){
    pmf_levels <- max_stepsize * all_levels * (log2(1 + all_levels))^2
  } else {
    pmf_levels <- sqrt(max_stepsize) * all_levels * (log2(1 + all_levels))^2
  }
  pmf_levels <- pmf_levels / sum(pmf_levels)
  pmf_tail <- rev(cumsum(rev(pmf_levels)))
  
  return(list(support = all_levels, mass_function = pmf_levels, tail_function = pmf_tail))
}

#' @rdname single_term
#' @title Unbiased single-term estimator of the score function
#' @description Estimates the score function using the single-term estimator of Rhee and Glynn (2015)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param observations a matrix of observations of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled2_resampling a 2-marginal coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param coupled4_resampling a 4-marginal coupled resampling scheme, such as \code{\link{coupled4_maximalchains_maximallevels_independent_residuals}}
#' @param initialization choice of distribution to initialize chains, such as \code{dynamics} or the default \code{particlefilter} 
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param level_distribution list containing mass_function and tail_function that specify the distribution of levels, 
#' e.g. by calling \code{\link{compute_level_distribution}} 
#' @return a list with objects such as 
#' random_level is the random level to truncated infinite sum 
#' unbiasedestimator is an unbiased estimator of the gradient of the log-likelihood
#' cost is the cost to compute the single-term estimator
#' elapsedtime is the time taken to compute the single-term estimator
#' @export
single_term <- function(model, theta, observations, nparticles, resampling_threshold = 1, coupled2_resampling, coupled4_resampling, 
                        initialization = "particlefilter", algorithm = "CPF", k = 0, m = 1, level_distribution){
  
  # start timer
  tic()
  
  # sample random level 
  random_level <- sample(x = level_distribution$support, size = 1, 
                         prob = level_distribution$mass_function)
  # minimum level 
  minimum_level <- min(level_distribution$support)
  
  if (random_level == minimum_level){
    
    # compute score at coarsest discretization level
    discretization <- model$construct_discretization(minimum_level)
    score <- unbiased_discretized_score(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, 
                                        initialization, algorithm, k = k, m = m, max_iterations = Inf)
    estimator <- score$unbiasedestimator / level_distribution$mass_function[1]
    cost <- nparticles * discretization$nsteps * score$cost 
    
  } else {
    
    # compute score increment at random level
    discretization <- model$construct_successive_discretization(random_level)
    score_increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                                initialization, algorithm, k = k, m = m, max_iterations = Inf)
    increment <- score_increment$unbiasedestimator
    estimator <- increment / level_distribution$mass_function[random_level-minimum_level+1]
    cost <- nparticles * discretization$coarse$nsteps * score_increment$cost_coarse
    cost <- cost + nparticles * discretization$fine$nsteps * score_increment$cost_fine
    
  }
  
  # end timer and compute elapsed time
  timer <- toc(quiet = TRUE)
  elapsedtime <- timer$toc - timer$tic
  
  return(list(random_level = random_level, unbiasedestimator = estimator, 
              cost = cost, elapsedtime = elapsedtime))
  
}

#' @rdname independent_sum
#' @title Unbiased independent-sum estimator of the score function
#' @description Estimates the score function using the independent-sum estimator of Rhee and Glynn (2015)
#' @param model a list representing a hidden Markov model, e.g. \code{\link{hmm_ornstein_uhlenbeck}}
#' @param theta a vector of parameters as input to model functions
#' @param observations a matrix of observations of size nobservations x ydimension
#' @param nparticles number of particles
#' @param resampling_threshold ESS proportion below which resampling is triggered (always resample at observation times by default)
#' @param coupled2_resampling a 2-marginal coupled resampling scheme, such as \code{\link{coupled2_maximal_independent_residuals}}
#' @param coupled4_resampling a 4-marginal coupled resampling scheme, such as \code{\link{coupled4_maximalchains_maximallevels_independent_residuals}}
#' @param initialization choice of distribution to initialize chains, such as \code{dynamics} or the default \code{particlefilter} 
#' @param algorithm character specifying type of algorithm desired, i.e. 
#' \code{\link{CPF}} for conditional particle filter, 
#' \code{\link{CASPF}} for conditional ancestor sampling particle filter,
#' \code{\link{CBSPF}} for conditional backward sampling particle filter
#' @param k iteration at which to start averaging (default to 0)
#' @param m iteration at which to stop averaging (default to 1)
#' @param level_distribution list containing mass_function and tail_function that specify the distribution of levels, 
#' e.g. by calling \code{\link{compute_level_distribution}} 
#' @return a list with objects such as 
#' random_level is the random level to truncated infinite sum 
#' unbiasedestimator is an unbiased estimator of the gradient of the log-likelihood
#' cost is the cost to compute the independent-sum estimator
#' elapsedtime is the time taken to compute the independent-sum estimator
#' @export
independent_sum <- function(model, theta, observations, nparticles, resampling_threshold = 1, coupled2_resampling, coupled4_resampling, 
                            initialization = "particlefilter", algorithm = "CPF", k = 0, m = 1, level_distribution){
  
  # start timer
  tic()
  
  # minimum level 
  minimum_level <- min(level_distribution$support)
  
  # initialize estimate by computing gradient at coarsest discretization level
  discretization <- model$construct_discretization(minimum_level)
  cat("running minimum level", "\n")
  score <- unbiased_discretized_score(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, 
                                      initialization, algorithm, k = k, m = m, max_iterations = Inf)
  cat("minimum level completed", "\n")
  estimator <- score$unbiasedestimator
  cost <- nparticles * discretization$nsteps * score$cost 
  
  # sample random level to truncated infinite sum
  all_levels <- minimum_level:maximum_level
  random_level <- sample(x = level_distribution$support, size = 1, prob = level_distribution$mass_function)
  cat("random level:", random_level, "\n")
  
  # increment estimate by computing the difference of the gradient at two successive discretization levels
  if (random_level > minimum_level){
    for (level in (minimum_level+1):random_level){
      cat("running level", level, "\n")
      discretization <- model$construct_successive_discretization(level)
      score_increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                                  initialization, algorithm, k = k, m = m, max_iterations = Inf)
      cat("completed", "\n")
      cost <- cost + nparticles * discretization$coarse$nsteps * score_increment$cost_coarse
      cost <- cost + nparticles * discretization$fine$nsteps * score_increment$cost_fine
      increment <- score_increment$unbiasedestimator
      estimator <- estimator + increment / level_distribution$tail_function[level-minimum_level+1]
    }
  }
  
  # end timer and compute elapsed time
  timer <- toc(quiet = TRUE)
  elapsedtime <- timer$toc - timer$tic
  
  return(list(random_level = random_level, unbiasedestimator = estimator, 
              cost = cost, elapsedtime = elapsedtime))
  
}

