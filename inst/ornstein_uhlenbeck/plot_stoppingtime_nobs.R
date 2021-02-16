rm(list=ls())
library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# load dataset
load("inst/ornstein_uhlenbeck/results/stoppingtimes_nobs_maximal.RData")

# boxplots of stoppng time
stoppingtimes.df$nparticles <- factor(c(rep("linear N", 500), rep("fixed N", 500)))
g <- ggplot(stoppingtimes.df, aes(x = T, y = stoppingtime)) + geom_boxplot(aes(group = T)) +
  scale_y_continuous(trans = "log2") + scale_x_continuous(trans = "log2", breaks = grid_nobs) +
  facet_wrap(~nparticles, nrow = 2) + xlab("T") + ylab("stopping time") 
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_stoppingtimes_nobs.eps", 
       plot = g, device = "eps", width = 8, height = 6)

# mean and median of stopping time
summarise.df <- ddply(stoppingtimes.df, c("T", "nparticles"), summarise,
                      mean = mean(stoppingtime), median = median(stoppingtime))
ggplot(summarise.df, aes(x = T, y = mean, colour = nparticles)) + geom_point(aes(group = T)) + scale_color_colorblind()
ggplot(summarise.df, aes(x = T, y = median, colour = nparticles)) + geom_point(aes(group = T)) + scale_color_colorblind()



