rm(list = ls())
library(UnbiasedScore)
library(tictoc)

# load simulated dataset
load("inst/ornstein_uhlenbeck/simulated_data_T25.RData")

# settings
resampling_threshold <- 1
coupled2_resampling <- coupled2_maximal_independent_residuals
coupled4_resampling <- coupled4_maximalchains_maximallevels_independent_residuals
initialization <- "dynamics"
algorithm <- "CPF"
nrepeats <- 100

# tuning parameters
grid_nparticles <- c(2^7, 2^9)
k <- 9 # maximal
m <- 10 * k # time average

# distribution of levels
minimum_level <- 3
level_distribution <- compute_level_distribution(model, minimum_level)
pmf.df <- data.frame(support = level_distribution$support, probability = level_distribution$mass_function)
ggplot(pmf.df, aes(x = support, y = probability)) + geom_bar(stat = "identity") + xlab("level") + ylab("probability")

# compute true score function at DGP
score_true <- model$compute_gradients(theta_true, observations)

# preallocate
estimator.df <- data.frame()
for (i in 1:nrepeats){
    cat("Repetition:", i, "\n")
    
    # compute independent sum estimator of score function using naive estimators
    score_estimator <- independent_sum(model, theta_true, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                       initialization, algorithm, k = 0, m = 0, level_distribution)
    squarederror <- sum((score_estimator$unbiasedestimator - score_true)^2)
    estimator.df <- rbind(estimator.df, data.frame(estimator = factor("naive"),
                                                   repetition = i,
                                                   highest_level = score_estimator$random_level, 
                                                   squarederror = squarederror,
                                                   cost = score_estimator$cost,
                                                   elapsedtime = score_estimator$elapsedtime,
                                                   theta1 = score_estimator$unbiasedestimator[1],
                                                   theta2 = score_estimator$unbiasedestimator[2],
                                                   theta3 = score_estimator$unbiasedestimator[3]))
    cat("Squared error of naive estimator:", squarederror, "Highest level:", score_estimator$random_level, "\n")
    
    # compute independent sum estimator of score function using simple estimators
    score_estimator <- independent_sum(model, theta_true, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                       initialization, algorithm, k = k, m = k, level_distribution)
    squarederror <- sum((score_estimator$unbiasedestimator - score_true)^2)
    estimator.df <- rbind(estimator.df, data.frame(estimator = factor("simple"),
                                                   repetition = i,
                                                   highest_level = score_estimator$random_level, 
                                                   squarederror = squarederror,
                                                   cost = score_estimator$cost,
                                                   elapsedtime = score_estimator$elapsedtime,
                                                   theta1 = score_estimator$unbiasedestimator[1],
                                                   theta2 = score_estimator$unbiasedestimator[2],
                                                   theta3 = score_estimator$unbiasedestimator[3]))
    cat("Squared error of simple estimator:", squarederror, "Highest level:", score_estimator$random_level, "\n")
    
    # compute independent sum estimator of score function using time-averaged estimators
    score_estimator <- independent_sum(model, theta_true, observations, nparticles, resampling_threshold, coupled2_resampling, coupled4_resampling, 
                                       initialization, algorithm, k = k, m = m, level_distribution)
    squarederror <- sum((score_estimator$unbiasedestimator - score_true)^2)
    estimator.df <- rbind(estimator.df, data.frame(estimator = factor("time-averaged"),
                                                   repetition = i,
                                                   highest_level = score_estimator$random_level, 
                                                   squarederror = squarederror,
                                                   cost = score_estimator$cost,
                                                   elapsedtime = score_estimator$elapsedtime,
                                                   theta1 = score_estimator$unbiasedestimator[1],
                                                   theta2 = score_estimator$unbiasedestimator[2],
                                                   theta3 = score_estimator$unbiasedestimator[3]))
    cat("Squared error of time-averaged estimator:", squarederror, "Highest level:", score_estimator$random_level, "\n")
    
}
# save results
save.image(file = "inst/ornstein_uhlenbeck/results/estimators_tuning.RData")



