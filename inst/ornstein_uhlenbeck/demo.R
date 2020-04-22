rm(list = ls())
library(UnbiasedGradients)
module_tree <<- Module("module_tree", PACKAGE = "UnbiasedGradients")
TreeClass <<- module_tree$Tree

# define diffusion model of latent process 
diffusion_model <- list()
diffusion_model$xdimension <- 1
diffusion_model$x_star <- 0 # initial condition
diffusion_model$b_theta <- function(theta, x) - theta * x # drift
diffusion_model$grad_b_theta <- function(theta, x) - x

diffusion_model$sigma <- function(x) 1 # diffusivity coefficient
diffusion_model$a_matrix <- function(x) diffusion_model$sigma(x)^2
diffusion_model$a_inverse <- function(x) 1 / diffusion_model$a_matrix(x)
diffusion_model$nparameters <- 1

# observation process
obs_model <- list()
obs_model$ydimension <- 1
obs_model$h_theta <- function(theta, x) x 
obs_model$grad_h_theta <- function(theta, x) 0 
obs_model$sigma <- 1
obs_model$nparameters <- 0 

# problem setting 
terminal_time <- 2 # time interval is [0,terminal_time]
max_level <- 5 # highest frequency level
nofsteps <- 2^max_level 
stepsize <- 1 / nofsteps 
statelength <- terminal_time * nofsteps + 1 
datalength <- terminal_time + 1 
OU <- get_OU_hmm(terminal_time, diffusion_model, obs_model, datalength)

# generate states from discretized model
true_theta <- 0.8 # data generating parameter 
X <- rep(0, statelength)
X[1] <- diffusion_model$x_star
for (k in 2:statelength){
  X[k] <- OU$rtransition(X[k-1], true_theta, max_level, rnorm(1))
}

# generate observations from discretized model
Y <- rep(0, datalength)
for (k in 1:datalength){
  Y[k] <- obs_model$h_theta(true_theta, X[(k-1) * nofsteps + 1]) + obs_model$sigma * rnorm(1)
}
ytrajectory <- matrix(Y, ncol = 1) # datalength x 1 

# plot latent and observation process
time_interval <- seq(0, terminal_time, by = stepsize)
plot(x = time_interval, y = X, type = "l", xlab = "time", ylab = "")
points(x = 0:terminal_time, y = Y, col = "red")

# CPF
nparticles <- 2^5
output <- CPF(OU, true_theta, max_level, ytrajectory, nparticles, ref_trajectory = NULL)
output$new_trajectory  
output$ess  

# C3PF
nrepeats <- 100
level_grid <- 1:5
level_ngrid <- length(level_grid)
theta <- true_theta
nparticles <- 2^5
coupled_resampling <- CCR_indexmatching #CCR_maxchains_maxlevels_coupled_residuals
k <- 10
m <- 100
cost <- matrix(0, nrepeats, level_ngrid)
diff <- matrix(0, nrepeats, level_ngrid)

for (level in level_grid){
  cat("Level:", level, "\n")
  single_kernel <- function(level, xtrajectory=NULL) CPF_discrete_simple(OU, theta, level, ytrajectory, nparticles,
                                                                         ref_trajectory = xtrajectory)
  
  coupled_coupled_kernel <- function(xtrajectory_coarse1, xtrajectory_coarse2, xtrajectory_fine1, xtrajectory_fine2){
    return(C3PF_discrete_simple(OU, theta, level, ytrajectory, nparticles, coupled_resampling,
                                xtrajectory_coarse1, xtrajectory_coarse2,
                                xtrajectory_fine1, xtrajectory_fine2))
  }
  
  lambda <- function(level, xtrajectory) OU$lambda(theta, level, xtrajectory, ytrajectory)
  
  for (j in 1:nrepeats){
    start_time = Sys.time()
    results <- unbiased_diff_gradient_simple(single_kernel, coupled_coupled_kernel, lambda, level, k = k, m = m)
    cost[j, level] <- Sys.time() - start_time
    diff[j, level] <- results$uestimator
  }
}

mean_cost <- rep(0, level_ngrid)
var_diff <- rep(0, level_ngrid)
for (level in level_grid){
  mean_cost[level] <- mean(cost[,level])
  var_diff[level] <- var(diff[, level])
}


h_l_grid <- 2**(-level_grid)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(x = h_l_grid, y = mean_cost, type = "l", xlab = "", ylab = "", axes=F)
axis(2, ylim=c(0,1),col="black",las=1)
mtext("cost",side=2,line=2.5)
par(new=TRUE)
plot(x = h_l_grid, y = var_diff, type = "l", xlab = "", ylab = "", axes=F, col="red")
mtext("variance",side=4,col="red",line=4) 
axis(4, col="red",col.axis="red",las=1)

axis(1, pretty(range(h_l_grid), 10))
mtext("h_l",side=1,col="black",line=2.5)  

