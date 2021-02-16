rm(list = ls())
library(UnbiasedScore)
library(tictoc)

# coupled 2-marginal resampling scheme
coupled2_resampling <- coupled2_maximal_independent_residuals

# coupled 4-marginal resampling scheme
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals

# settings
resampling_threshold <- 1
level <- 6
nrepeats <- 100
grid_nobs <- 25 * 2^(0:4)
length_nobs <- length(grid_nobs)
grid_nparticles <- 2^7 * 2^(0:4) # increase N linearly with T

# preallocate
stoppingtimes.df <- data.frame()

for (i in 1:length_nobs){
  # number of observations
  T <- grid_nobs[i]
  
  # number of particles
  N <- grid_nparticles[i]
  
  # load simulated dataset
  load(paste(paste("inst/ornstein_uhlenbeck/simulated_data_", T, sep = "T"), "RData", sep = "."))
  
  # construct discretization
  discretization <- model$construct_successive_discretization(level)
  for (j in 1:nrepeats){
    # run unbiased estimation with fixed number of particles
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles = grid_nparticles[1], resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = 0, m = 1, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    stoppingtimes.df <- rbind(stoppingtimes.df, data.frame(T = T, repetition = j, stoppingtime = stopping_time, nparticles = factor("fixed")))
    cat("Number of observations:", T, "Repetition:", j, "\n", "Stopping time for fixed N:", stopping_time, "\n")
    
    # run unbiased estimation with increased number of particles
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles = N, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = 0, m = 1, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    stoppingtimes.df <- rbind(stoppingtimes.df, data.frame(T = T, repetition = j, stoppingtime = stopping_time, nparticles = factor("linear")))
    cat("Number of observations:", T, "Repetition:", j, "\n", "Stopping time for linear N:", stopping_time, "\n")
  }
}
save.image(file = "inst/ornstein_uhlenbeck/results/stoppingtimes_nobs_maximal.RData")

