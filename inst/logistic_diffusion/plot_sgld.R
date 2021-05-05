rm(list=ls())
library(ggplot2)
#setmytheme()

plot.df <- data.frame()
N <- 5000
comp <- 6
load(sprintf("inst/logistic_diffusion/results/SGLD.RData",comp))
repf <- factor(rep(comp, N))
plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = parameters.df$X1, rep=repf, component = factor(rep(1, N))))
plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = parameters.df$X2, rep=repf, component = factor(rep(2, N))))
plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = parameters.df$X3, rep=repf, component = factor(rep(3, N))))
plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = parameters.df$X4, rep=repf, component = factor(rep(4, N))))

g <- ggplot(plot.df) + geom_line(aes(x = iter, y = theta)) +
  facet_wrap(~component, nrow=2, ncol=2, scales = "free_y") +
  xlab("iteration") + ylab("components") +
  theme_bw() + theme(strip.background=element_rect(fill="white"))
g

#ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_stoppingtimes_levels.eps",
#plot = g, device = "eps", width = 8, height = 8)
