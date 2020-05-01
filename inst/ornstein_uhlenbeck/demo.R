rm(list = ls())
library(Rcpp)
library(UnbiasedGradients)
library(ggplot2)
library(ggthemes)
setmytheme()
module_tree <<- Module("module_tree", PACKAGE = "UnbiasedGradients")
TreeClass <<- module_tree$Tree

# time interval of interest
terminal_time <- 10 # time interval is [0,terminal_time]

# discretization level
current_level <- 5 
stepsize <- 2^(-current_level)
# nofsteps <- 2^max_level 

# construct hidden Markov model for time discretized process
model <- hmm_ornstein_uhlenbeck(terminal_time)
statelength <- model$statelength(current_level)

# generate states from true continuous time process
theta_true <- c(2, 7, 1) # data generating parameter 
X <- rep(0, statelength)
X[1] <- model$rinit(1)
for (k in 2:statelength){
  # using exact transition on a unit time interval 
  conditioanl_mean <- theta_true[2] + (X[k-1] - theta_true[2]) * exp(- theta_true[1]) 
  conditional_sd <- model$sigma * sqrt(1 - exp(- 2 * theta_true[1])) / sqrt(2 * theta_true[1])
  X[k] <- conditioanl_mean + conditional_sd * rnorm(1)
}

# generate observations at unit times
Y <- X[model$obstimes(current_level)] + sqrt(theta_true[3]) * rnorm(terminal_time)
observations <- matrix(Y, ncol = 1) # terminal_time x 1 

# plot latent and observation process
time_interval <- seq(0, terminal_time, by = stepsize)
plot(x = time_interval, y = X, type = "l", xlab = "time", ylab = "")
points(x = 1:terminal_time, y = Y, col = "red")

# CPF
nparticles <- 2^3
CPF_output <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)

# ess plot
ess.df <- data.frame(time = time_interval[-1], ess = CPF_output$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# check faithfulness of 2-way coupled CPF 
coupled_resampling <- coupled2_resampling
# coupled_resampling <- coupled2_maximal_coupled_residuals
ref_trajectory1 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory2 <- ref_trajectory1 
coupled2CPF <- coupled2_CPF(model, theta_true, current_level, observations, nparticles, coupled_resampling, 
                            ref_trajectory1, ref_trajectory2)

# generate 2-way coupled CPF chains
niterations <- 20
meetingtime <- Inf
meet <- FALSE
distance <- rep(0, niterations+1)
functional1 <- matrix(0, niterations+1, model$theta_dimension)
functional2 <- matrix(0, niterations+1, model$theta_dimension)
ref_trajectory1 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory2 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
distance[1] <- sqrt(sum((ref_trajectory1 - ref_trajectory2)^2))
functional1[1, ] <- model$functional(theta_true, current_level, ref_trajectory1, observations)
functional2[1, ] <- model$functional(theta_true, current_level, ref_trajectory2, observations)

for (i in 1:niterations){
  cat("Iteration:", i, "\n")
  coupled2CPF <- coupled2_CPF(model, theta_true, current_level, observations, nparticles, coupled_resampling, 
                              ref_trajectory1, ref_trajectory2)
  ref_trajectory1 <- coupled2CPF$new_trajectory1
  ref_trajectory2 <- coupled2CPF$new_trajectory2
  
  distance[i+1] <- sqrt(sum((ref_trajectory1 - ref_trajectory2)^2))
  functional1[i+1, ] <- model$functional(theta_true, current_level, ref_trajectory1, observations)
  functional2[i+1, ] <- model$functional(theta_true, current_level, ref_trajectory2, observations)
  
  if (all(ref_trajectory1 == ref_trajectory2) & meet == FALSE){
    meetingtime <- i 
    meet <- TRUE
  }
}
cat("Meeting time:", meetingtime, "\n")

# plot distance between chains
distance.df <- data.frame(iteration = 0:niterations, distance = distance)
ggplot(distance.df, aes(x = iteration, y = distance)) + geom_line()

# plot difference in functional values between chains for parameter of interest
parameter <- 3
functional1.df <- data.frame(iteration = 0:niterations, functional = functional1[, parameter])
functional2.df <- data.frame(iteration = 0:niterations, functional = functional2[, parameter])
ggplot() + geom_line(data = functional1.df, aes(x = iteration, y = functional), colour = "blue") + 
  geom_line(data = functional2.df, aes(x = iteration, y = functional), colour = "red")

# check faithfulness of 4-way coupled CPF
# coupled_resampling <- coupled4_resampling
# coupled_resampling <- coupled4_maximal_coupled_residuals
coupled_resampling <- coupled4_maximalchains_maximallevels_coupled_residuals
# coupled_resampling <- coupled4_maximallevels_maximalchains_coupled_residuals
ref_trajectory_coarse1 <- CPF(model, theta_true, current_level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- ref_trajectory_coarse1 # check faithfulness
# ref_trajectory_coarse2 <- CPF(model, theta_true, current_level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory

ref_trajectory_fine1 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- ref_trajectory_fine1 # check faithfulness
# ref_trajectory_fine2 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory

coupled4CPF <- coupled4_CPF(model, theta_true, current_level, observations, nparticles, coupled_resampling,
                            ref_trajectory_coarse1, ref_trajectory_coarse2,
                            ref_trajectory_fine1, ref_trajectory_fine2)
all(coupled4CPF$new_trajectory_coarse1 == coupled4CPF$new_trajectory_coarse2)
all(coupled4CPF$new_trajectory_fine1 == coupled4CPF$new_trajectory_fine2)

# generate 4-way coupled CPF chains
niterations <- 20
meetingtime_coarse <- Inf
meetingtime_fine <- Inf
meet_coarse <- FALSE
meet_fine <- FALSE

distance_coarse <- rep(0, niterations+1)
distance_fine <- rep(0, niterations+1)
functional1 <- matrix(0, niterations+1, model$theta_dimension)
functional2 <- matrix(0, niterations+1, model$theta_dimension)

ref_trajectory_coarse1 <- CPF(model, theta_true, current_level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory_coarse2 <- CPF(model, theta_true, current_level-1, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine1 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory
ref_trajectory_fine2 <- CPF(model, theta_true, current_level, observations, nparticles, ref_trajectory = NULL)$new_trajectory

distance_coarse[1] <- sqrt(sum((ref_trajectory_coarse1 - ref_trajectory_coarse2)^2))
distance_fine[1] <- sqrt(sum((ref_trajectory_fine1 - ref_trajectory_fine2)^2))

functional1[1, ] <- abs(model$functional(theta_true, current_level-1, ref_trajectory_coarse1, observations) - 
  model$functional(theta_true, current_level, ref_trajectory_fine1, observations))
functional2[1, ] <- abs(model$functional(theta_true, current_level-1, ref_trajectory_coarse2, observations) - 
  model$functional(theta_true, current_level, ref_trajectory_fine2, observations))

for (i in 1:niterations){
  cat("Iteration:", i, "\n")
  coupled4CPF <- coupled4_CPF(model, theta_true, current_level, observations, nparticles, coupled_resampling,
                              ref_trajectory_coarse1, ref_trajectory_coarse2,
                              ref_trajectory_fine1, ref_trajectory_fine2)
  ref_trajectory_coarse1 <- coupled4CPF$new_trajectory_coarse1
  ref_trajectory_coarse2 <- coupled4CPF$new_trajectory_coarse2
  ref_trajectory_fine1 <- coupled4CPF$new_trajectory_fine1
  ref_trajectory_fine2 <- coupled4CPF$new_trajectory_fine2
  
  distance_coarse[i+1] <- sqrt(sum((ref_trajectory_coarse1 - ref_trajectory_coarse2)^2))
  distance_fine[i+1] <- sqrt(sum((ref_trajectory_fine1 - ref_trajectory_fine2)^2))
  
  functional1[i+1, ] <- abs(model$functional(theta_true, current_level-1, ref_trajectory_coarse1, observations) - 
    model$functional(theta_true, current_level, ref_trajectory_fine1, observations))
  functional2[i+1, ] <- abs(model$functional(theta_true, current_level-1, ref_trajectory_coarse2, observations) - 
    model$functional(theta_true, current_level, ref_trajectory_fine2, observations))
  
  # check if chains meet
  if (all(ref_trajectory_coarse1 == ref_trajectory_coarse2) & meet_coarse == FALSE){
    meetingtime_coarse <- i 
    meet_coarse <- TRUE
  }
  
  if (all(ref_trajectory_fine1 == ref_trajectory_fine2) & meet_fine == FALSE){
    meetingtime_fine <- i 
    meet_fine <- TRUE
  }
}
cat("Meeting time of coarse chains:", meetingtime_coarse, "\n")
cat("Meeting time of fine chains:", meetingtime_fine, "\n")

# plot distance between chains
distance.df <- data.frame(iteration = 0:niterations, distance = distance_coarse, level = factor(rep("coarse", niterations+1)))
distance.df <- rbind(distance.df, data.frame(iteration = 0:niterations, distance = distance_fine, level = factor(rep("fine", niterations+1))))
ggplot(distance.df, aes(x = iteration, y = distance, color = level)) + geom_line() + scale_color_colorblind()

# plot difference in functional values between chains for parameter of interest
parameter <- 2
functional.df <- data.frame(iteration = 0:niterations, functional = functional1[, parameter], pair = factor(rep("1", niterations+1)))
functional.df <- rbind(functional.df, data.frame(iteration = 0:niterations, functional = functional2[, parameter], pair = factor(rep("2", niterations+1))))
ggplot(functional.df, aes(x = iteration, y = functional, color = pair)) + geom_line() + scale_color_colorblind() + ylab("functional difference")

# compute unbiased gradient estimator
coupled_resampling <- coupled2_resampling
# coupled_resampling <- coupled2_maximal_coupled_residuals
gradient <- unbiased_gradient(model, theta_true, current_level, observations, nparticles, coupled_resampling, 
                              k = 0, m = 1)
gradient

# compute unbiased gradient increment estimator
# coupled_resampling <- coupled4_resampling
# coupled_resampling <- coupled4_maximal_coupled_residuals
coupled_resampling <- coupled4_maximalchains_maximallevels_coupled_residuals
# coupled_resampling <- coupled4_maximallevels_maximalchains_coupled_residuals
increment <- unbiased_gradient_increment(model, theta_true, current_level, observations, nparticles, coupled_resampling, 
                                         k = 0, m = 1)
increment
