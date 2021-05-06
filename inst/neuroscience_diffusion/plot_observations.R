rm(list = ls())
library(UnbiasedScore)
library(ggplot2)
library(ggthemes)
setmytheme()

# load grid cells spike dataset of Hafting et al. (2008)
load("inst/neuroscience_diffusion/gridcells.RData")
terminal_time <- 20
spiketimes1 <- neuron1$ts[neuron1$ts < terminal_time]
spiketimes2 <- neuron2$ts[neuron2$ts < terminal_time]

# construct model
level_observation <- 6
model <- hmm_neuroscience_diffusion(spiketimes1, spiketimes2, level_observation, terminal_time)

# plot observation counts for a given discretization 
counts <- model$compute_observations()
observations <- counts$observations
nobservations <- counts$nbins # same as time discretization

gridcells <- data.frame(time = counts$time[1:nobservations], 
                        count = observations[, 1], 
                        cell = factor(rep(1, nobservations)))
gridcells <- rbind(gridcells, data.frame(time = counts$time[1:nobservations], 
                                         count = observations[, 2], 
                                         cell = factor(rep(2, nobservations))))
g <- ggplot(gridcells, aes(x = time, y = count, colour = cell)) + geom_point(size = 3) + 
  xlab("time (seconds)") + ylab("counts") + scale_color_colorblind() 
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/neural_network_observations.eps",
       plot = g, device = "eps", width = 9, height = 7)  
