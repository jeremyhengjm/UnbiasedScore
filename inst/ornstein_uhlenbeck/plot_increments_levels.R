rm(list=ls())
library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# create dataframe for plotting
var.df <- data.frame()

# compute sum of variance of increments for maximal resampling
load("inst/ornstein_uhlenbeck/results/increment_levels_maximal.RData")
var.naive <- ddply(naive.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.simple <- ddply(simple.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.timeaveraged <- ddply(timeaveraged.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
nlevel <- length(levels)
var.df <- rbind(var.df, data.frame(var.naive, estimator = factor(rep("naive", nlevel)), resampling = factor(rep("maximal", nlevel))))
var.df <- rbind(var.df, data.frame(var.simple, estimator = factor(rep("simple", nlevel)), resampling = factor(rep("maximal", nlevel))))
var.df <- rbind(var.df, data.frame(var.timeaveraged, estimator = factor(rep("time-averaged", nlevel)), resampling = factor(rep("maximal", nlevel))))

# compute sum of variance of increments for other maximal resampling
load("inst/ornstein_uhlenbeck/results/increment_levels_othermaximal.RData")
var.naive <- ddply(naive.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.simple <- ddply(simple.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.timeaveraged <- ddply(timeaveraged.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
nlevel <- length(levels)
var.df <- rbind(var.df, data.frame(var.naive, estimator = factor(rep("naive", nlevel)), resampling = factor(rep("other maximal", nlevel))))
var.df <- rbind(var.df, data.frame(var.simple, estimator = factor(rep("simple", nlevel)), resampling = factor(rep("other maximal", nlevel))))
var.df <- rbind(var.df, data.frame(var.timeaveraged, estimator = factor(rep("time-averaged", nlevel)), resampling = factor(rep("other maximal", nlevel))))

# compute sum of variance of increments for common uniforms resampling
load("inst/ornstein_uhlenbeck/results/increment_levels_common.RData")
var.naive <- ddply(naive.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.simple <- ddply(simple.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
var.timeaveraged <- ddply(timeaveraged.df, c("level"), summarise, var = var(theta1) + var(theta2) + var(theta3))
nlevel <- length(levels)
var.df <- rbind(var.df, data.frame(var.naive, estimator = factor(rep("naive", nlevel)), resampling = factor(rep("common uniforms", nlevel))))
var.df <- rbind(var.df, data.frame(var.simple, estimator = factor(rep("simple", nlevel)), resampling = factor(rep("common uniforms", nlevel))))
var.df <- rbind(var.df, data.frame(var.timeaveraged, estimator = factor(rep("time-averaged", nlevel)), resampling = factor(rep("common uniforms", nlevel))))

# plot sum of variances against level
g <- ggplot(var.df, aes(x = level, y = var, color = estimator)) + geom_point(size = 3) + 
  facet_wrap(~resampling, nrow = 3, scales = "free_y") + 
  scale_y_continuous(trans = "log2") + xlab("level") + ylab("sum of variance") + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_varincrements_levels.eps",
       plot = g, device = "eps", width = 8, height = 8)


