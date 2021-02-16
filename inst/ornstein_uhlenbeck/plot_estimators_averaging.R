library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# load independent repeats
rm(list = ls())
load("inst/ornstein_uhlenbeck/results/estimators_tuning_combined.RData")

# configurations and subset results
maxnrepeats <- 20 # number of independent score estimators to be averaged 
nrealizations <- 100
naive.df <- results.df[results.df$estimator == "naive", ]
simple.df <- results.df[results.df$estimator == "simple", ]
timeaveraged.df <- results.df[results.df$estimator == "time-averaged", ]

# preallocate
estimator.df <- data.frame()

# averaging naive estimators
for (i in 1:nrealizations){
  # subset dataframe
  index <- (i-1) * maxnrepeats + 1:maxnrepeats 
  
  # compute squared error
  squarederror <- (cumsum(naive.df$theta1[index]) / 1:maxnrepeats - score_true[1])^2 
  squarederror <- squarederror + (cumsum(naive.df$theta2[index]) / 1:maxnrepeats - score_true[2])^2 
  squarederror <- squarederror + (cumsum(naive.df$theta3[index]) / 1:maxnrepeats - score_true[3])^2 
  
  # compute cost 
  cost <- cumsum(naive.df$cost[index])
  
  # store results
  estimator.df <- rbind(estimator.df, data.frame(nrepeats = 1:maxnrepeats,
                                                 squarederror = squarederror, 
                                                 cost = cost, 
                                                 estimator = factor(rep("naive", maxnrepeats))))
}

# averaging simple estimators
for (i in 1:nrealizations){
  # subset dataframe
  index <- (i-1) * maxnrepeats + 1:maxnrepeats 
  
  # compute squared error
  squarederror <- (cumsum(simple.df$theta1[index]) / 1:maxnrepeats - score_true[1])^2 
  squarederror <- squarederror + (cumsum(simple.df$theta2[index]) / 1:maxnrepeats - score_true[2])^2 
  squarederror <- squarederror + (cumsum(simple.df$theta3[index]) / 1:maxnrepeats - score_true[3])^2 
  
  # compute cost 
  cost <- cumsum(simple.df$cost[index])
  
  # store results
  estimator.df <- rbind(estimator.df, data.frame(nrepeats = 1:maxnrepeats,
                                                 squarederror = squarederror, 
                                                 cost = cost, 
                                                 estimator = factor(rep("simple", maxnrepeats))))
}

# averaging time-averaged estimators
for (i in 1:nrealizations){
  # subset dataframe
  index <- (i-1) * maxnrepeats + 1:maxnrepeats 
  
  # compute squared error
  squarederror <- (cumsum(timeaveraged.df$theta1[index]) / 1:maxnrepeats - score_true[1])^2 
  squarederror <- squarederror + (cumsum(timeaveraged.df$theta2[index]) / 1:maxnrepeats - score_true[2])^2 
  squarederror <- squarederror + (cumsum(timeaveraged.df$theta3[index]) / 1:maxnrepeats - score_true[3])^2 
  
  # compute cost 
  cost <- cumsum(timeaveraged.df$cost[index])
  
  # store results
  estimator.df <- rbind(estimator.df, data.frame(nrepeats = 1:maxnrepeats,
                                                 squarederror = squarederror, 
                                                 cost = cost, 
                                                 estimator = factor(rep("time-averaged", maxnrepeats))))
}

# compute mean squared error
mse.df <- ddply(estimator.df, c("nrepeats", "estimator"), summarise, mse = mean(squarederror))

# plot mean squared error
g <- ggplot(mse.df, aes(x = nrepeats, y = mse * nrepeats, color = estimator)) + geom_point(size = 3) + 
  scale_y_log10() + xlab("replicates") + ylab("MSE x replicates") + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_meansquarederror_averaging.eps",
       plot = g, device = "eps", width = 8, height = 6)

# compute quantiles of cost
cost.df <- ddply(estimator.df, c("nrepeats", "estimator"), summarise, 
                 lower = quantile(cost, probs = 0.25),
                 median = median(cost),
                 upper = quantile(cost, probs = 0.75))

# plot cost against replicates
g <- ggplot(cost.df, aes(x = nrepeats, y = median, ymin = lower, ymax = upper, color = estimator)) 
g <- g + geom_pointrange() + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log2")
g <- g + xlab("replicates") + ylab("cost") + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_cost_averaging.eps",
       plot = g, device = "eps", width = 8, height = 6)
