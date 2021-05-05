rm(list=ls())
library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# create dataframe for plotting
plot.df <- data.frame()

# compute sum of variance of increments for adaptive resampling
load("inst/logistic_diffusion/results/increments_adaptive_resampling.RData")
var.df <- ddply(results.df, c("level", "estimator"), summarise, 
                variance = log2(var(theta1) + var(theta2) + var(theta3) + var(theta4)))
var.df$level <- levels[var.df$level]
plot.df <- data.frame(level = var.df$level, 
                      estimator = var.df$estimator, 
                      variance = var.df$variance, 
                      resampling = factor(rep("adaptive", nrow(var.df))))

# plot sum of variances against level
g <- ggplot(plot.df, aes(x = level, y = variance, color = estimator)) + geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, size = 1) + 
  labs(x = "level", y = expression(paste("sum of variance ", (log[2]) ))) + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/logistic_diffusion_varincrements_levels.eps",
       plot = g, device = "eps", width = 8, height = 6)

# check slopes
lm(variance ~ level, data = plot.df, subset = (plot.df$estimator == "simple"))
lm(variance ~ level, data = plot.df, subset = (plot.df$estimator == "time-averaged"))
