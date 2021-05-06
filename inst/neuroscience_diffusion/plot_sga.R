rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggthemes)
setmytheme()

load("inst/neuroscience_diffusion/results/stochastic_gradient.RData")

# Components to plot
dims <- 1:6
# Maximum number of iterations
N <- 1750

iterations <- as.vector(t(unique(parameters.df[c("j")])))
n <- length(iterations)
dim <- length(dims)

# Compute mean of smoothed SGA paths
theta_mean <- matrix(0, N, dim)
for (comp in 1:dim){
  iloc <- 0
  theta <- matrix(NA, N, n)
  for (i in iterations) {
    iloc <- iloc + 1
    theta_loc <- pull(parameters.df %>% filter(j == i), dims[comp])
    tmp <- cumsum(theta_loc[1:N]) / (1:N)
    theta[, iloc] <- tmp
    theta_mean[, comp] <- theta_mean[, comp] + tmp
  }
  #matplot(theta, type = 'l', col = 'black', lty = 1, xlab = 'iterations', ylab = sprintf('theta[%d]',comp))
}
theta_mean <- theta_mean / n

# Compute variance
theta_var <- matrix(0, N, dim)
for (comp in 1:dim){
  for (i in iterations) {
    theta_loc <- pull(parameters.df %>% filter(j == i), comp)
    tmp <- cumsum(theta_loc[1:N]) / (1:N)
    theta_var[, comp] <- theta_var[, comp] + (tmp - theta_mean[, comp])**2
  }
}
theta_var <- theta_var / n

# Create data frame for plot
plot.df <- data.frame()
for (comp in dims){
  plot.df <- rbind(plot.df, data.frame(iter = 1:N, 
                                       theta = theta_mean[, comp],
                                       var = theta_var[, comp],
                                       parameter = factor(rep(comp, N))))
}

# Plot
g <- ggplot(plot.df, aes(x = iter, y = theta, colour = parameter)) + 
  geom_line(size = 2) +
  geom_ribbon(aes(x = iter, ymax = theta + sqrt(var), ymin = theta - sqrt(var), fill = parameter), alpha = .1, show.legend = FALSE) +
  scale_colour_manual(values = unique(plot.df$parameter),
                      breaks = unique(plot.df$parameter),
                      labels = list(bquote(alpha[1]), bquote(beta[1]), bquote(gamma[1]), bquote(delta[1]), bquote(sigma[1]), bquote(kappa[1]))) +
  xlab("iteration") + ylab("")
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/neural_network_sga.pdf",
       plot = g, device = "pdf", width = 9, height = 7)
