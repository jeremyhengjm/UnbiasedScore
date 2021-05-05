rm(list = ls())
library(Rcpp)
library(UnbiasedScore)
library(ggplot2)
library(ggthemes)
setmytheme()

# load red kangaroo dataset
load("inst/logistic_diffusion/kangaroo.RData")
nobservations <- length(kangaroo$time)
observations <- matrix(0, nrow = nobservations, ncol = 2)
observations[, 1] <- kangaroo$count1
observations[, 2] <- kangaroo$count2
g <- ggplot(kangaroo, aes(x = time)) + geom_segment(aes(x = time, y = count1, xend = time, yend = count2), linetype = "dashed", colour = "black") +
  geom_point(aes(y = count1)) +
  geom_point(aes(y = count2)) +
  xlab("time (years)") + ylab("counts")

ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/logisic_diffusion_observations.eps",
       plot = g, device = "eps", width = 8, height = 6)
