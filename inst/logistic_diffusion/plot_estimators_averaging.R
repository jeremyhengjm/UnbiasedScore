library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# load independent repeats
rm(list = ls())
load("inst/logistic_diffusion/results/estimators_combined.RData") 
names(results.df)
dim(results.df)

# approximate true score 
score_true <- colMeans(results.df[, c("theta1", "theta2", "theta3", "theta4")])

# configurations and subset results
maxnrepeats <- 20 # number of independent score estimators to be averaged 
nrealizations <- 1000

# trim tails as empirical averages are not robust 
trim_tails <- function(x){
  for (i in 1:4){
    cat("Parameter", i, "\n")
    mean_x <- mean(x[, i])
    sd_x <- sd(x[, i])
    index <- (x[, i] > (mean_x + 8 * sd_x))
    cat("Number of outliers on right tail", sum(index), "\n")
    x[index, i] <- score_true[i]
    index <- (x[, i] < (mean_x - 8 * sd_x))
    cat("Number of outliers on left tail", sum(index), "\n")
    x[index, i] <- score_true[i]
  }
  return(x)
}
results.df[, c("theta1", "theta2", "theta3", "theta4")] <- trim_tails(results.df[, c("theta1", "theta2", "theta3", "theta4")])

# averaging simple estimators
estimator.df <- data.frame()
for (i in 1:nrealizations){
  # subset dataframe
  index <- (i-1) * maxnrepeats + 1:maxnrepeats 
  
  # compute squared error for each parameter component
  squarederror1 <- (cumsum(results.df$theta1[index]) / 1:maxnrepeats - score_true[1])^2
  squarederror2 <- (cumsum(results.df$theta2[index]) / 1:maxnrepeats - score_true[2])^2
  squarederror3 <- (cumsum(results.df$theta3[index]) / 1:maxnrepeats - score_true[3])^2
  squarederror4 <- (cumsum(results.df$theta4[index]) / 1:maxnrepeats - score_true[4])^2
  
  # compute cost 
  cost <- cumsum(results.df$cost[index])
  
  # store results
  estimator.df <- rbind(estimator.df, data.frame(nrepeats = 1:maxnrepeats,
                                                 squarederror1 = squarederror1, 
                                                 squarederror2 = squarederror2, 
                                                 squarederror3 = squarederror3, 
                                                 squarederror4 = squarederror4, 
                                                 cost = cost))
}

# compute mean squared error
mse.df <- ddply(estimator.df, c("nrepeats"), summarise, 
                mse1 = mean(squarederror1),
                mse2 = mean(squarederror2),
                mse3 = mean(squarederror3),
                mse4 = mean(squarederror4))
                
# create dataframe for plotting
plot.df <- data.frame()
plot.df <- rbind(plot.df, data.frame(nrepeats = mse.df$nrepeats, 
                                     mse = mse.df$mse1,
                                     component = factor(rep(1, maxnrepeats))))
plot.df <- rbind(plot.df, data.frame(nrepeats = mse.df$nrepeats, 
                                     mse = mse.df$mse2,
                                     component = factor(rep(2, maxnrepeats))))
plot.df <- rbind(plot.df, data.frame(nrepeats = mse.df$nrepeats, 
                                     mse = mse.df$mse3,
                                     component = factor(rep(3, maxnrepeats))))
plot.df <- rbind(plot.df, data.frame(nrepeats = mse.df$nrepeats, 
                                     mse = mse.df$mse4,
                                     component = factor(rep(4, maxnrepeats))))

# plot mean squared error
g <- ggplot(plot.df, aes(x = nrepeats, y = mse * nrepeats, colour = component)) + geom_point(size = 3) + 
  scale_y_log10() + xlab("replicates") + ylab("MSE x replicates") + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/logistic_diffusion_meansquarederror_averaging.eps",
       plot = g, device = "eps", width = 8, height = 6)

# compute quantiles of cost
cost.df <- ddply(estimator.df, c("nrepeats"), summarise,
                 lower = quantile(cost, probs = 0.25),
                 median = median(cost),
                 upper = quantile(cost, probs = 0.75))

# plot cost against replicates
g <- ggplot(cost.df, aes(x = nrepeats, y = median, ymin = lower, ymax = upper))
g <- g + geom_pointrange() + scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log2")
g <- g + xlab("replicates") + ylab("cost") + scale_color_colorblind()
g

