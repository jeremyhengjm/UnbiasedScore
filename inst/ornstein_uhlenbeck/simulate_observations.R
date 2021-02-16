rm(list = ls())
library(UnbiasedScore)

# time interval of interest
terminal_time <- 25
times <- 1:terminal_time
nobservations <- length(times) # nobservations at unit times

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

# save simulated dataset
save.image(file = "inst/ornstein_uhlenbeck/simulated_data_T25.RData")
