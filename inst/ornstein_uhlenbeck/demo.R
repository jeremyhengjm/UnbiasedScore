rm(list = ls())
library(Rcpp)
library(UnbiasedScore)
library(ggplot2)
library(ggthemes)
library(tictoc)
setmytheme()
module_tree <<- Module("module_tree", PACKAGE = "UnbiasedGradients")
TreeClass <<- module_tree$Tree

# time interval of interest
terminal_time <- 10
times <- 1:terminal_time
nobservations <- terminal_time # nobservations at unit times

# construct hidden Markov model 
model <- hmm_ornstein_uhlenbeck(times)

# generate states (from true continuous time process) and observations 
set.seed(17)
theta_true <- c(2, 7, 1) # data generating parameter 
X <- rep(0, nobservations)
Y <- rep(0, nobservations)
current_state <- model$rinit(1)
for (k in 1:nobservations){
  # using exact transition on a unit time interval 
  conditional_mean <- theta_true[2] + (current_state - theta_true[2]) * exp(- theta_true[1]) 
  conditional_sd <- model$sigma * sqrt(1 - exp(- 2 * theta_true[1])) / sqrt(2 * theta_true[1])
  new_state <- conditional_mean + conditional_sd * rnorm(1)
  X[k] <- new_state
  Y[k] <- new_state + sqrt(theta_true[3]) * rnorm(1)
  current_state <- new_state
}
observations <- matrix(Y, ncol = 1) # nobservations x 1 

# Kalman smoothing
smoothing <- model$smoothing_moments(theta_true, observations) # needs to be modified

# plot latent and observation process
plot.df <- data.frame(time = times, state = X, process = factor(rep("latent state", nobservations)))
plot.df <- rbind(plot.df, data.frame(time = times, state = Y, process = factor(rep("observation", nobservations))))
plot.df <- rbind(plot.df, data.frame(time = times, state = smoothing$post_mean, process = factor(rep("smoother", nobservations))))
plot.df <- rbind(plot.df, data.frame(time = times, state = smoothing$prior_mean, process = factor(rep("prior", nobservations))))
ggplot(plot.df, aes(x = time, y = state, color = process)) + geom_point() + 
  scale_color_colorblind() + ylab("") + scale_x_continuous(breaks = times)

# particle filter
theta <- theta_true
level <- 4
discretization <- model$construct_discretization(level)
nparticles <- 2^8
resampling_threshold <- 1
CPF_output <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold)

# ESS plot
ess.df <- data.frame(time = seq(0, terminal_time, length.out = discretization$statelength), 
                     ess = CPF_output$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# CPF
xtrajectory <- CPF_output$new_trajectory
CPF_output <- CPF(model, theta, discretization, observations, nparticles, resampling_threshold, ref_trajectory = xtrajectory)

# 2-CCPF
ref_trajectory1 <- CPF_output$new_trajectory
ref_trajectory2 <- ref_trajectory1
coupled_resampling <- coupled2_maximal_coupled_residuals
coupled2CPF <- coupled2_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, 
                            ref_trajectory1, ref_trajectory2)
all(coupled2CPF$new_trajectory1 == coupled2CPF$new_trajectory2) # check faithfuless of 2-CCPF

# unbiased estimation of time-discretized scores
theta <- theta_true
level <- 5
discretization <- model$construct_discretization(level)
nparticles <- 2^8
resampling_threshold <- 1
coupled_resampling <- coupled2_maximal_independent_residuals
# coupled_resampling <- coupled2_maximal_coupled_residual
initialization <- "dynamics"
# initialization <- "particlefilter"
algorithm <- "CPF"
# algorithm <- "CASPF"
# algorithm <- "CBSPF"
avg_meetingtime <- 0
iterations <- 100
for (k in 1:iterations){
  discretized_score <- unbiased_discretized_score(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling, 
                                                  initialization, algorithm, k = 5, m = 10, max_iterations = Inf)
  avg_meetingtime <- avg_meetingtime + discretized_score$meetingtime
}
avg_meetingtime <- avg_meetingtime / iterations
avg_meetingtime

# 4-CCPF
level <- 5
discretization <- model$construct_successive_discretization(level)
theta <- theta_true
nparticles <- 2^8
resampling_threshold <- 1
ref_trajectory_coarse1 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- CPF(model, theta, discretization$coarse, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine1 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- CPF(model, theta, discretization$fine, observations, nparticles, resampling_threshold, ref_trajectory = NULL)$new_trajectory
coupled_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
coupled4CPF <- coupled4_CPF(model, theta, discretization, observations, nparticles, resampling_threshold, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)

# check faithfuless of 4-CCPF
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

# unbiased estimation of score increments
theta <- theta_true
level <- 5
discretization <- model$construct_successive_discretization(level)
nparticles <- 2^8
resampling_threshold <- 1
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
# coupled2_resampling <- coupled2_maximal_coupled_residual
# coupled4_resampling <- coupled4_maximalchains_maximallevels_coupled_residuals
initialization <- "dynamics"
# initialization <- "particlefilter"
algorithm <- "CPF"
#algorithm <- "CASPF"
# algorithm <- "CBSPF"
avg_meetingtime_fine <- 0
avg_meetingtime_coarse <- 0
iterations <- 100
for (k in 1:iterations){
  score_increment <- unbiased_score_increment(model, theta, discretization, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                            initialization, algorithm, k = 5, m = 20, max_iterations = Inf)
  avg_meetingtime_fine <- avg_meetingtime_fine + score_increment$meetingtime_fine
  avg_meetingtime_coarse <- avg_meetingtime_coarse + score_increment$meetingtime_coarse
}
avg_meetingtime_fine <- avg_meetingtime_fine / iterations
avg_meetingtime_coarse <- avg_meetingtime_coarse / iterations
avg_meetingtime_fine
avg_meetingtime_coarse

# compute true score function at DGP
true_score <- model$compute_gradients(theta_true, observations)
true_score

# distribution of levels
minimum_level <- 4
maximum_level <- 20
level_distribution <- compute_level_distribution(model, minimum_level, maximum_level)
pmf.df <- data.frame(support = level_distribution$support, probability = level_distribution$mass_function)
ggplot(pmf.df, aes(x = support, y = probability)) + geom_bar(stat="identity")

# settings
nparticles <- 2^9
resampling_threshold <- 1
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
# coupled2_resampling <- coupled2_maximal_coupled_residual
# coupled4_resampling <- coupled4_maximalchains_maximallevels_coupled_residuals
initialization <- "dynamics"
# initialization <- "particlefilter"
algorithm <- "CPF"
# algorithm <- "CASPF"
# algorithm <- "CBSPF"

# compute single term estimator of score function
score <- single_term(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                     initialization, algorithm, k = 5, m = 10, level_distribution)

# compute independent sum estimator of score function
score <- independent_sum(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                         initialization, algorithm, k = 5, m = 10, level_distribution)

# compute stratified estimators of score function 
score <- stratified_estimator(model, theta, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                              initialization, algorithm, k = 5, m = 10, level_distribution, 
                              nrepeats = 3, stratification = "uniformly-stratified")
score <- stratified_estimator(model, theta_true, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                              initialization, algorithm, k = 5, m = 10, level_distribution, 
                              nrepeats = 3, stratification = "systematic-sampling")