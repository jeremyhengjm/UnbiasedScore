rm(list = ls())
library(UnbiasedScore)
library(tictoc)

# load red kangaroo dataset
load("inst/logistic_diffusion/kangaroo.RData")
nobservations <- length(kangaroo$time)
observations <- matrix(0, nrow = nobservations, ncol = 2)
observations[, 1] <- kangaroo$count1
observations[, 2] <- kangaroo$count2

# construct hidden Markov model (inferring diffusivity parameter)
model <- hmm_logistic_diffusion_full(kangaroo$time)
theta <- c(2.397, 4.429e-03, 0.840, 17.631) 

# settings
nparticles <- 2^8
# resampling_threshold <- 1 # always resampling
resampling_threshold <- 0.5 # adaptive resampling
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"
nrepeats <- 20

# estimator settings
# k <- 23 # always resampling
k <- 20 # adaptive resampling

# distribution of levels
minimum_level <- 3
level_distribution <- compute_level_distribution(model, minimum_level)

# preallocate
estimator.df <- data.frame()

for (i in 1:nrepeats){
  cat("Repetition:", i, "\n")
  
  # compute independent sum estimator of score function using simple estimators
  score_estimator <- independent_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                     initialization, algorithm, k = k, m = k, level_distribution)
  
  estimator.df <- rbind(estimator.df, data.frame(repetition = i,
                                                 highest_level = score_estimator$random_level, 
                                                 cost = score_estimator$cost,
                                                 elapsedtime = score_estimator$elapsedtime,
                                                 theta1 = score_estimator$unbiasedestimator[1],
                                                 theta2 = score_estimator$unbiasedestimator[2],
                                                 theta3 = score_estimator$unbiasedestimator[3],
                                                 theta4 = score_estimator$unbiasedestimator[4]))
    
  cat("Highest level:", score_estimator$random_level, "Estimate:", score_estimator$unbiasedestimator, "\n")
  
  # save results
  save.image(file = "inst/logistic_diffusion/results/estimators_independent.RData")
}



