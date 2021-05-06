rm(list = ls())
library(Rcpp)
library(UnbiasedScore)
library(ggplot2)
library(ggthemes)
library(tictoc)

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

# particle filter
level <- 13
discretization <- model$construct_discretization(level)
theta <- rep(1, model$theta_dimension)
nparticles <- 2^8
resampling_threshold <- 0.5
tic()
CPF_output <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory = NULL, treestorage = FALSE)
toc()

# ESS plot
nsteps <- discretization$nsteps
ess.df <- data.frame(step = 0:nsteps, ess = CPF_output$ess * 100 / nparticles)
ggplot(ess.df, aes(x = step, y = ess)) + geom_line() + 
  labs(x = "step", y = "ESS%") + ylim(c(0, 100))

# conditional particle filter (CPF)
xtrajectory <- CPF_output$new_trajectory
CPF_output <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory = xtrajectory)

# 2-CCPF
ref_trajectory1 <- CPF_output$new_trajectory
ref_trajectory2 <- ref_trajectory1
resampling_threshold <- 0.5
coupled_resampling <- coupled2_maximal_independent_residuals
coupled2CPF <- coupled2_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, 
                            ref_trajectory1, ref_trajectory2)
all(coupled2CPF$new_trajectory1 == coupled2CPF$new_trajectory2) # check faithfuless of 2-way coupled CPF

# unbiased score estimation at a level
theta <- rep(1, model$theta_dimension)
level <- 13
discretization <- model$construct_discretization(level)
nparticles <- 2^8
resampling_threshold <- 0.5
coupled_resampling <- coupled2_maximal_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"
discretized_score <- unbiased_discretized_score(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                                                initialization, algorithm, k = 50, m = 50, max_iterations = Inf)
discretized_score

# 4-CCPF
level <- 13
discretization <- model$construct_successive_discretization(level)
theta <- rep(1, model$theta_dimension)
nparticles <- 2^8
resampling_threshold <- 0.5
ref_trajectory_coarse1 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine1 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
coupled_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
coupled4CPF <- coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)

# check 4-CCPF
## chains on coarse level
ref_trajectory_coarse1 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- ref_trajectory_coarse1
ref_trajectory_fine1 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
coupled4CPF <- coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)
all(coupled4CPF$new_trajectory_coarse1 == coupled4CPF$new_trajectory_coarse2)
all(coupled4CPF$new_trajectory_fine1 == coupled4CPF$new_trajectory_fine2)

## chains on fine level
ref_trajectory_coarse1 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine1 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- ref_trajectory_fine1
coupled4CPF <- coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)
all(coupled4CPF$new_trajectory_coarse1 == coupled4CPF$new_trajectory_coarse2)
all(coupled4CPF$new_trajectory_fine1 == coupled4CPF$new_trajectory_fine2)

## chains on both levels
ref_trajectory_coarse1 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- ref_trajectory_coarse1
ref_trajectory_fine1 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- ref_trajectory_fine1
coupled4CPF <- coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, 
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)
all(coupled4CPF$new_trajectory_coarse1 == coupled4CPF$new_trajectory_coarse2)
all(coupled4CPF$new_trajectory_fine1 == coupled4CPF$new_trajectory_fine2)

# unbiased score estimation at a level
theta <- rep(1, model$theta_dimension)
level <- 13
discretization <- model$construct_successive_discretization(level)
nparticles <- 2^8
resampling_threshold <- 0.5
coupled_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"
score_increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                            initialization, algorithm, k = 50, m = 50, max_iterations = Inf)
score_increment

# distribution of levels
minimum_level <- 11
maximum_level <- 15
level_distribution <- compute_level_distribution(model, minimum_level, maximum_level)
pmf.df <- data.frame(support = level_distribution$support, probability = level_distribution$mass_function)
ggplot(pmf.df, aes(x = support, y = probability)) + geom_bar(stat="identity")

# settings
theta <- rep(1, model$theta_dimension)
nparticles <- 2^8
resampling_threshold <- 0.5
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"

# compute single term estimator of score function
score <- single_term(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                     initialization, algorithm, k = 50, m = 50, level_distribution)

# compute independent sum estimator of score function
score <- independent_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                         initialization, algorithm, k = 50, m = 50, level_distribution)

