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
levels <- 2:10
nrepeats <- 100

# preallocate
stoppingtimes.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    # run unbiased estimation with law of latent dynamics as initialization
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization, algorithm = "CPF", k = 0, m = 1, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    stoppingtimes.df <- rbind(stoppingtimes.df, data.frame(level = current_level, repetition = i, 
                                                           stoppingtime = stopping_time))
    cat("Level:", current_level, "Repetition:", i, "\n", "Stopping time:", stopping_time, "\n")
    
  }
}
save.image(file = "inst/logistic_diffusion/results/stoppingtimes_levels.RData")
