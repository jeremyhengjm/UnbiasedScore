rm(list=ls())
library(ggplot2)
library(ggthemes)
library(plyr)
setmytheme()

# create dataframe for plotting
plot.df <- data.frame()

load("inst/neuroscience_diffusion/results/Stoppingtimes_CPF.RData")
nlevels <- length(as.vector(t(unique(stoppingtimes.df["level"]))))
always.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                   median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$median,
                                     summary = factor(rep("median", nlevels)),
                                     filter = factor(rep('CPF', nlevels))))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevels)),
                                     filter = factor(rep('CPF', nlevels))))

load("inst/neuroscience_diffusion/results/Stoppingtimes_CASPF.RData")
nlevels <- length(as.vector(t(unique(stoppingtimes.df["level"]))))
always.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                   median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$median,
                                     summary = factor(rep("median", nlevels)),
                                     filter = factor(rep('CASPF', nlevels))))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevels)),
                                     filter = factor(rep('CASPF', nlevels))))

# plot stopping times
g <- ggplot(plot.df, aes(x = level, y = stoppingtime, colour = filter)) + 
  facet_wrap(~summary, nrow = 2, scales = "free_y") +
  geom_point(size = 3) + geom_line() + scale_color_colorblind() + 
  xlab("level") + ylab("stopping time")
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/neural_network_stoppingtimes.eps", 
       plot = g, device = "eps", width = 9, height = 7)
