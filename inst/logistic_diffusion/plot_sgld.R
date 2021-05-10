rm(list=ls())
library(ggplot2)
setmytheme()

N <- 7500
load("inst/logistic_diffusion/results/sgld.RData")

# Parameters of the SGLD
theta_dimension <- 4
a <- rep(1e-2, theta_dimension)
a[3] = 1e-4
b <- 100
gam <- 0.6

delta <- matrix(rep(1:N, theta_dimension), ncol=N, byrow=T)
delta <- a / (b + delta)**gam

plot.df <- data.frame()
for (dim in 1:theta_dimension){
  trace <- parameters.df[,dim]
  running_avg <- cumsum(trace * delta[dim, ]) / cumsum(delta[dim, ])
  # plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = trace, smooth = running_avg, component = factor(rep(dim, N))))
  plot.df <- rbind(plot.df, data.frame(iter = 1:N, theta = trace, smooth = running_avg, component = factor(rep(paste0("theta[", dim, "]"), N))))
}

# Plot
g <- ggplot(plot.df) + geom_line(aes(x = iter, y = smooth)) + 
  facet_wrap(~component, nrow = 2, ncol = 2, scales = "free_y", labeller = label_parsed) + 
  scale_x_continuous(breaks = c(1, 2500, 5000, 7500)) +
  xlab("iteration") + ylab("")
g

ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/logistic_diffusion_sgld.pdf",
      plot = g, device = "pdf", width = 14, height = 6)



