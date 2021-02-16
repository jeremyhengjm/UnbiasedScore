rm(list=ls())
library(UnbiasedScore)
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# load results
load("inst/ornstein_uhlenbeck/results/estimators_tuning.RData")

# inspect cost vs elapsed time 
ggplot(estimator.df, aes(x = cost, elapsedtime)) + geom_point() + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") + 
  xlab("cost") + ylab("elapsed time")

# squared error against cost
estimator.df$nparticles <- c(rep("N = 128", 3 * nrepeats), rep("N = 512", 3 * nrepeats))
g <- ggplot(estimator.df, aes(x = cost, y = squarederror, color = estimator)) + geom_point(size = 2) + 
  facet_wrap(~nparticles, nrow = 2) + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +  
  xlab("cost") + ylab("squared error") + scale_color_colorblind()
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_squarederror_tuning.eps",
       plot = g, device = "eps", width = 8, height = 6)



