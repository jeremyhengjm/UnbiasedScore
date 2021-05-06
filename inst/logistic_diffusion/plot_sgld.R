rm(list=ls())
library(ggplot2)
#setmytheme()

N <- 5000
comp <- 6
load(sprintf("inst/logistic_diffusion/results/SGLD_%d.RData",comp))

# Parameters of the SGLD
theta_dimension <- 4
a <- rep(1e-2, theta_dimension)
a[3] = 1e-4
b <- 100
gam <- 0.6

delta <- matrix(rep(1:N, theta_dimension), ncol=N, byrow=T)
delta <- a / (b + delta)**gam

plot.df <- data.frame()
repf <- factor(rep(comp, N))
for (dim in 1:theta_dimension){
  trace <- parameters.df[,dim]
  running_avg <- cumsum(trace * delta[dim, ]) / cumsum(delta[dim, ])
  plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = running_avg, rep=repf, component = factor(rep(dim, N))))
}

g <- ggplot(plot.df) + geom_line(aes(x = iter, y = theta)) +
  facet_wrap(~component, nrow=2, ncol=2, scales = "free_y") +
  xlab("iteration") + ylab("components") +
  theme_bw() + theme(strip.background=element_rect(fill="white"))
g

#ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/ou_stoppingtimes_levels.eps",
#plot = g, device = "eps", width = 8, height = 8)
