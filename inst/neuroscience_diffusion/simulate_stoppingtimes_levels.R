rm(list = ls())
library(UnbiasedScore)
library(tictoc)

file_path <- "inst/neuroscience_diffusion/results/"
file_sub_name <- 'stoppingtimes_CPF'
file_name <- sprintf("%s%s.RData", file_path, file_sub_name)

# load dataset
load("inst/neuroscience_diffusion/gridcells.RData")
terminal_time <- 20
spiketimes1 <- neuron1$ts[neuron1$ts < terminal_time]
spiketimes2 <- neuron2$ts[neuron2$ts < terminal_time]

# construct model
level_observation <- 6
model <- hmm_neuroscience_diffusion(spiketimes1, spiketimes2, level_observation, terminal_time)

# plot observation counts for a given discretization 
counts <- model$compute_observations()
observations <- counts$observations
nobservations <- counts$nbins # same as time discretization

theta <- c(1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1)

# settings
nparticles <- 2^8
resampling_threshold <- 0.5 # adaptive resampling
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"

levels <- 11:15
nrepeats <- 10

# preallocate
stoppingtimes.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    # run unbiased estimation with law of latent dynamics as initialization
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = initialization, algorithm = algorithm, k = 0, m = 0, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    i_global <- (igrid - 1) * nrepeats + i
    stoppingtimes.df <- rbind(stoppingtimes.df, data.frame(level = current_level, repetition = i_global, 
                                                           stoppingtime = stopping_time))
    cat("Level:", current_level, "Repetition:", i_global, "\n", "Stopping time:", stopping_time, "\n")
    
  }
  # save results
  save.image(file = file_name)
}



