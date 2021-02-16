rm(list=ls())
library(UnbiasedScore)
library(ggplot2)
library(plyr)
setmytheme()

# create dataframe for plotting
plot.df <- data.frame()

# maximal resampling
load("inst/ornstein_uhlenbeck/results/stoppingtimes_levels_maximal.RData")
stoppingtimes.df$level <- factor(stoppingtimes.df$level)
plot.df <- rbind(plot.df, data.frame(level = stoppingtimes.df$level, 
                                     stoppingtime = stoppingtimes.df$stoppingtime, 
                                     resampling = factor(rep("maximal", nrepeats))))

# other maximal resampling
load("inst/ornstein_uhlenbeck/results/stoppingtimes_levels_othermaximal.RData")
plot.df <- rbind(plot.df, data.frame(level = stoppingtimes.df$level, 
                                     stoppingtime = stoppingtimes.df$stoppingtime, 
                                     resampling = factor(rep("other maximal", nrepeats))))

# common uniform resampling
load("inst/ornstein_uhlenbeck/results/stoppingtimes_levels_common.RData")
plot.df <- rbind(plot.df, data.frame(level = stoppingtimes.df$level, 
                                     stoppingtime = stoppingtimes.df$stoppingtime, 
                                     resampling = factor(rep("common uniforms", nrepeats))))

# boxplots of stopping times
g <- ggplot(plot.df, aes(x = level, y = stoppingtime)) + geom_boxplot(aes(group = level)) +
  facet_wrap(~resampling, nrow = 3, scales = "free_y") + 
  xlab("level") + ylab("stopping time")
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_stoppingtimes_levels.eps",
plot = g, device = "eps", width = 8, height = 8)

# median and 90% quantile of stoppingtimes at each discretization level
ddply(plot.df, c("level", "resampling"), summarise, 
      median = median(stoppingtime),
      quantile = quantile(stoppingtime, probs = 0.9))

