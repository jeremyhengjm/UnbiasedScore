rm(list = ls())
library(UnbiasedScore)
library(tictoc)

file_path <- 'inst/neuroscience_diffusion/results/'
file_sub_name <- 'SGA_CASPF'
file_name <- sprintf('%s%s.RData', file_path, file_sub_name)

# load grid cells spike dataset of Hafting et al. (2008)
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

# settings
nparticles <- 2^8
resampling_threshold <- 0.5 # adaptive resampling
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CASPF"
learning_rate <- 1e-3
niterations <- 2000
k <- 100

# initial parameter
theta <- c(1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1)

# distribution of levels
minimum_level <- 11
maximum_level <- 15
level_distribution <- compute_level_distribution(model, minimum_level, maximum_level)

# preallocate
sga.df <- data.frame()
parameters.df <- data.frame(matrix(nrow = niterations, ncol = model$theta_dimension))
colnames(parameters.df) <- model$theta_names

for (i in 1:niterations){
  cat("Iteration:", i, "/", niterations, "\n")
  
  # compute independent sum estimator of score function using simple estimators
  score_estimator <- independent_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                     initialization, algorithm, k = k, m = k, level_distribution)
  
  # compute gradient using chain rule for parameters with positivity
  gradient <- as.numeric(model$theta_positivity) * theta * score_estimator$unbiasedestimator + 
    as.numeric(!model$theta_positivity) * score_estimator$unbiasedestimator
  
  # stochastic gradient algorithm (natural/logarithmic parameterization for unconstrained/constrained parameters)
  log_theta <- log_transform(model, theta) + learning_rate * gradient
  theta <- exp_transform(model, log_theta)
  
  # store output
  sga.df <- rbind(sga.df, data.frame(iteration = i,
                                     highest_level = score_estimator$random_level, 
                                     cost = score_estimator$cost,
                                     elapsedtime = score_estimator$elapsedtime))
  parameters.df[i, ] <- theta
  
  # print 
  cat("Highest level:", score_estimator$random_level, "Parameter:", theta, "\n")
  
  # save results
  save("sga.df", "parameters.df", "i", "theta", file = file_name)
}



