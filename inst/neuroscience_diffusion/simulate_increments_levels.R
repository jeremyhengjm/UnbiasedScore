rm(list = ls())
library(UnbiasedScore)
library(tictoc)

file_path <- "inst/neuroscience_diffusion/results/"
file_sub_name <- 'increments_CPF'
file_name <- sprintf("%s%s.RData", file_path, file_sub_name)

# load dataset
load("~/Rlocal/neuroscience_diffusion/gridcells.RData")
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
algorithm <- "CPF"

levels <- 11:15
nrepeats <- 10
k = 135 # CPF
# k = 105 # CASPF

# preallocate
naive.df <- data.frame()
simple.df <- data.frame()
for (current_level in levels){
  # construct discretization
  discretization <- model$construct_successive_discretization(current_level)
  for (i in 1:nrepeats){
    i_global <- (igrid - 1) * nrepeats + i
    cat("Level:", current_level, "Repetition:", i, "\n")
    # compute unbiased increment using naive estimators
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = algorithm, k = 0, m = 0, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    naive.df <- rbind(naive.df, data.frame(level = factor(current_level), 
                                           repetition = i_global, 
                                           stoppingtime = stopping_time, 
                                           squarednorm = squarednorm, 
                                           theta1 = increment$unbiasedestimator[1],
                                           theta2 = increment$unbiasedestimator[2],
                                           theta3 = increment$unbiasedestimator[3],
                                           theta4 = increment$unbiasedestimator[4],
                                           theta5 = increment$unbiasedestimator[5],
                                           theta6 = increment$unbiasedestimator[6],
                                           theta7 = increment$unbiasedestimator[7],
                                           theta8 = increment$unbiasedestimator[8],
                                           theta9 = increment$unbiasedestimator[9],
                                           theta10 = increment$unbiasedestimator[10],
                                           theta11 = increment$unbiasedestimator[11],
                                           theta12 = increment$unbiasedestimator[12]))
    cat("Squared norm of naive estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
    
    # compute unbiased increment using simple estimators
    increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                          initialization = "dynamics", algorithm = algorithm, k = k, m = k, max_iterations = Inf)
    stopping_time <- max(increment$meetingtime_coarse, increment$meetingtime_fine)
    squarednorm <- sum(increment$unbiasedestimator^2)
    simple.df <- rbind(simple.df, data.frame(level = factor(current_level), 
                                             repetition = i_global, 
                                             stoppingtime = stopping_time, 
                                             squarednorm = squarednorm, 
                                             theta1 = increment$unbiasedestimator[1],
                                             theta2 = increment$unbiasedestimator[2],
                                             theta3 = increment$unbiasedestimator[3],
                                             theta4 = increment$unbiasedestimator[4],
                                             theta5 = increment$unbiasedestimator[5],
                                             theta6 = increment$unbiasedestimator[6],
                                             theta7 = increment$unbiasedestimator[7],
                                             theta8 = increment$unbiasedestimator[8],
                                             theta9 = increment$unbiasedestimator[9],
                                             theta10 = increment$unbiasedestimator[10],
                                             theta11 = increment$unbiasedestimator[11],
                                             theta12 = increment$unbiasedestimator[12]))
    cat("Squared norm of simple estimator:", squarednorm, "Stopping time:", stopping_time, "\n")
  }
  # save results
  save.image(file = file_name)
}



