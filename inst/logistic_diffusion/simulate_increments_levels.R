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

# coupled 2-marginal resampling scheme
coupled2_resampling <- coupled2_maximal_independent_residuals

# coupled 4-marginal resampling scheme
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals

# settings
nparticles <- 2^8
resampling_threshold <- 1
# resampling_threshold <- 0.5
initialization <- "dynamics"
levels <- 3:10
nrepeats <- 100

# estimator settings
k <- 23 # always resampling
# k <- 20 # adaptive resampling

# preallocate
naive.df <- data.frame()
simple.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    cat("Level:", current_level, "Repetition:", i, "\n")
    # compute unbiased increment using naive estimators
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization, algorithm = "CPF", k = 0, m = 0, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    naive.df <- rbind(naive.df, data.frame(level = factor(current_level), 
                                           repetition = i, 
                                           stoppingtime = stopping_time, 
                                           squarednorm = squarednorm, 
                                           cost_coarse = increment$cost_coarse,
                                           cost_fine = increment$cost_fine,
                                           theta1 = increment$unbiasedestimator[1],
                                           theta2 = increment$unbiasedestimator[2],
                                           theta3 = increment$unbiasedestimator[3]))
    cat("Squared norm of naive estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
    
    # compute unbiased increment using simple estimators
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization, algorithm = "CPF", k = k, m = k, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    simple.df <- rbind(simple.df, data.frame(level = factor(current_level), 
                                             repetition = i, 
                                             stoppingtime = stopping_time, 
                                             squarednorm = squarednorm, 
                                             theta1 = increment$unbiasedestimator[1],
                                             theta2 = increment$unbiasedestimator[2],
                                             theta3 = increment$unbiasedestimator[3]))
    cat("Squared norm of simple estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
    
  }
}
save.image(file = "inst/logistic_diffusion/results/increment_levels.RData")
