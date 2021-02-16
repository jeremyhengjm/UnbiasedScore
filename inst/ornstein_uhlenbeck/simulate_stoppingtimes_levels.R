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
nrepeats <- 100

# preallocate
stoppingtimes.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    # run unbiased estimation with law of latent dynamics as initialization
    increment <- unbiased_score_increment(model, theta_true, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = "CPF", k = 0, m = 1, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    stoppingtimes.df <- rbind(stoppingtimes.df, data.frame(level = current_level, repetition = i, 
                                                           stoppingtime = stopping_time))
    cat("Level:", current_level, "Repetition:", i, "\n", "Stopping time:", stopping_time, "\n")
    
  }
}
save.image(file = "inst/ornstein_uhlenbeck/results/stoppingtimes_levels_maximal.RData")
