rm(list = ls())
library(UnbiasedScore)
library(tictoc)

file_path <- "inst/logistic_diffusion/results/"
file_sub_name <- 'SGLD'
file_name <- sprintf("%s%s.RData", file_path, file_sub_name)

# load red kangaroo dataset
load("inst/logistic_diffusion/kangaroo.RData")
nobservations <- length(kangaroo$time)
observations <- matrix(0, nrow = nobservations, ncol = 2)
observations[, 1] <- kangaroo$count1
observations[, 2] <- kangaroo$count2

# construct hidden Markov model (inferring diffusivity parameter)
model <- hmm_logistic_diffusion_full(kangaroo$time)

# settings
nparticles <- 2^8
resampling_threshold <- 0.5 # adaptive resampling
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"
niterations <- 5000

theta <- c(2.397, 4.429e-03, 0.840, 17.631)

# estimator settings
k <- 20 # adaptive resampling
a <- rep(1e-2, model$theta_dimension)
a[3] = 1e-4
b <- 100
gam <- 0.6
prior_mean <- c(0, -1, -1, -1)
prior_std <- c(5, 2, 2, 2)

# distribution of levels
minimum_level <- 3
maximum_level <- 13
level_distribution <- compute_level_distribution(model, minimum_level, maximum_level)

sgld.df <- data.frame()
parameters.df <- data.frame(matrix(nrow = niterations, ncol = model$theta_dimension))

# preallocate
colnames(parameters.df) <- model$theta_names

for (i in 1:niterations){
  cat("Iteration:", i, "/", niterations, "\n")
  
  # compute independent sum estimator of score function using simple estimators
  score_estimator <- independent_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                     initialization, algorithm, k = k, m = k, level_distribution)
  
  # compute gradient using chain rule for parameters with positivity
  gradient <- as.numeric(model$theta_positivity) * theta * score_estimator$unbiasedestimator + 
    as.numeric(!model$theta_positivity) * score_estimator$unbiasedestimator
  
  cat(gradient, '\n')
  
  # stochastic gradient algorithm (natural/logarithmic parameterization for unconstrained/constrained parameters)
  learning_rate <- a / (b + i)**gam
  log_theta <- log_transform(model, theta)
  log_theta <- log_theta +
    learning_rate / 2 * ( gradient - sum((log_theta - prior_mean) / prior_std**2 )/2 ) +
    rnorm(model$theta_dimension, mean=0, sd=learning_rate)
  theta <- exp_transform(model, log_theta)
  
  # store output
  sgld.df <- rbind(sgld.df, data.frame(iteration = i,
                                     highest_level = score_estimator$random_level, 
                                     cost = score_estimator$cost,
                                     elapsedtime = score_estimator$elapsedtime))
  parameters.df[i, ] <- theta
  
  # print 
  cat("Highest level:", score_estimator$random_level, "Parameter:", theta, "\n")
  
  # save results
  save("sgld.df", "parameters.df", "i", "theta", file = file_name)
}
