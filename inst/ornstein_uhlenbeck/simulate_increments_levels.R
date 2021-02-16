rm(list = ls())
library(UnbiasedScore)
library(tictoc)

# load simulated dataset
load("inst/ornstein_uhlenbeck/simulated_data_T25.RData")

# coupled 2-marginal resampling scheme
coupled2_resampling <- coupled2_maximal_independent_residuals
# coupled2_resampling <- coupled2_maximal_coupled_residuals
# coupled2_resampling <- coupled2_resampling

# coupled 4-marginal resampling scheme
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
# coupled4_resampling <- coupled4_maximalchains_maximallevels_coupled_residuals
# coupled4_resampling <- coupled4_resampling
# coupled4_resampling <- coupled4_maximal_independent_residuals
# coupled4_resampling <- coupled4_maximal_coupled_residuals

# settings
nparticles <- 2^7
resampling_threshold <- 1
levels <- 3:8
nrepeats <- 10

# estimator settings
k <- 9 # maximal
k <- 3 # other maximal
# k <- 168 # common uniforms
m <- 10 * k # time average

# preallocate
naive.df <- data.frame()
simple.df <- data.frame()
timeaveraged.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    cat("Level:", current_level, "Repetition:", i, "\n")
    # compute unbiased increment using naive estimators
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = 0, m = 0, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    naive.df <- rbind(naive.df, data.frame(level = factor(current_level), 
                                           repetition = i, 
                                           stoppingtime = stopping_time, 
                                           squarednorm = squarednorm, 
                                           theta1 = increment$unbiasedestimator[1],
                                           theta2 = increment$unbiasedestimator[2],
                                           theta3 = increment$unbiasedestimator[3]))
    cat("Squared norm of naive estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
    
    # compute unbiased increment using simple estimators
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = k, m = k, max_iterations = Inf)
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
    
    # compute unbiased increment using time-averaged estimators
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = k, m = m, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    timeaveraged.df <- rbind(timeaveraged.df, data.frame(level = factor(current_level), 
                                                         repetition = i, 
                                                         stoppingtime = stopping_time, 
                                                         squarednorm = squarednorm, 
                                                         theta1 = increment$unbiasedestimator[1],
                                                         theta2 = increment$unbiasedestimator[2],
                                                         theta3 = increment$unbiasedestimator[3]))
    cat("Squared norm of time-averaged estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
  }
}
save.image(file = "inst/ornstein_uhlenbeck/results/increment_levels_maximal.RData")
