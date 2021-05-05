rm(list=ls())
library(UnbiasedScore)
library(ggplot2)
library(ggthemes)
library(plyr)
setmytheme()

# create dataframe for plotting
plot.df <- data.frame()

# always resampling N = 128
load("inst/logistic_diffusion/results/stoppingtimes_always_resampling_N128.RData")
nlevel <- length(levels)
always.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                   median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# always resampling N = 256
load("inst/logistic_diffusion/results/stoppingtimes_always_resampling_N256.RData")
always.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                   median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# always resampling N = 512
load("inst/logistic_diffusion/results/stoppingtimes_always_resampling_N512.RData")
always.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                   median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = always.df$level, 
                                     stoppingtime = always.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("always", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# adaptive resampling N = 128
load("inst/logistic_diffusion/results/stoppingtimes_adaptive_resampling_N128.RData")
adaptive.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                     median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# adaptive resampling N = 256
load("inst/logistic_diffusion/results/stoppingtimes_adaptive_resampling_N256.RData")
adaptive.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                     median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# adaptive resampling N = 512
load("inst/logistic_diffusion/results/stoppingtimes_adaptive_resampling_N512.RData")
adaptive.df <- ddply(stoppingtimes.df, c("level"), summarise, 
                     median = median(stoppingtime), quantile90 = quantile(stoppingtime, probs = 0.9))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$median,
                                     summary = factor(rep("median", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

plot.df <- rbind(plot.df, data.frame(level = adaptive.df$level, 
                                     stoppingtime = adaptive.df$quantile90,
                                     summary = factor(rep("90%-quantile", nlevel)),
                                     resampling = factor(rep("adaptive", nlevel)), 
                                     N = factor(rep(nparticles, nlevel))))

# plot stopping times
g <- ggplot(plot.df, aes(x = level, y = stoppingtime, colour = N, shape = resampling)) + 
  facet_wrap(~summary, nrow = 2, scales = "free_y") +
  geom_point(size = 3) + geom_line() + scale_color_colorblind() + 
  xlab("level") + ylab("stopping time")
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/logisic_diffusion_stoppingtimes_levels.eps",
       plot = g, device = "eps", width = 9, height = 7)



