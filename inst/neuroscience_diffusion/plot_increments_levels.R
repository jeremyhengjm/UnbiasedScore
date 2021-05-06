rm(list=ls())
library(plyr)
library(ggplot2)
library(ggthemes)
setmytheme()

# create dataframe for plotting
var.df <- data.frame()

# compute sum of variance of increments for CPF
load("inst/neuroscience_diffusion/results/increments_CPF_naive.RData")
var.naive <- ddply(naive.df, c("level"), summarise, var = log2(var(theta1) + var(theta2) + var(theta3) + var(theta4) + var(theta5) + var(theta6)
                   + var(theta7) + var(theta8) + var(theta9) + var(theta10) + var(theta11) + var(theta12)))
load("inst/neuroscience_diffusion/results/increments_CPF_simple.RData")
var.simple <- ddply(simple.df, c("level"), summarise, var = log2(var(theta1) + var(theta2) + var(theta3) + var(theta4) + var(theta5) + var(theta6)
                    + var(theta7) + var(theta8) + var(theta9) + var(theta10) + var(theta11) + var(theta12)))

nlevel <- nrow(var.naive)
var.df <- rbind(var.df, data.frame(var.naive, estimator = factor(rep("naive", nlevel)), MCMC = factor(rep("CPF", nlevel)), id = factor(rep("naive PF", nlevel))))
nlevel <- nrow(var.simple)
var.df <- rbind(var.df, data.frame(var.simple, estimator = factor(rep("simple", nlevel)), MCMC = factor(rep("CPF", nlevel)), id = factor(rep("simple PF", nlevel))))

# compute sum of variance of increments for CASPF
load("inst/neuroscience_diffusion/results/increments_CASPF_naive.RData")
var.naive <- ddply(naive.df, c("level"), summarise, var = log2(var(theta1) + var(theta2) + var(theta3) + var(theta4) + var(theta5) + var(theta6)
                   + var(theta7) + var(theta8) + var(theta9) + var(theta10) + var(theta11) + var(theta12)))
load("inst/neuroscience_diffusion/results/increments_CASPF_simple.RData")
var.simple <- ddply(simple.df, c("level"), summarise, var = log2(var(theta1) + var(theta2) + var(theta3) + var(theta4) + var(theta5) + var(theta6)
                    + var(theta7) + var(theta8) + var(theta9) + var(theta10) + var(theta11) + var(theta12)))

nlevel <- nrow(var.naive)
var.df <- rbind(var.df, data.frame(var.naive, estimator = factor(rep("naive", nlevel)), MCMC = factor(rep("CASPF", nlevel)), id = factor(rep("naive ASPF", nlevel))))
nlevel <- nrow(var.simple)
var.df <- rbind(var.df, data.frame(var.simple, estimator = factor(rep("simple", nlevel)), MCMC = factor(rep("CASPF", nlevel)), id = factor(rep("simple ASPF", nlevel))))

# plot sum of variances against level
g <- ggplot(var.df, aes(x = level, y = var, color = estimator, shape = MCMC, group = id)) + geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_y_continuous(trans = "log2") + xlab("level") + ylab(expression(paste("sum of variance ", (log[2]) ))) + scale_color_colorblind() 
g
ggsave(filename = "~/Dropbox/UnbiasedGradients/draft/arXiv-v1/neural_network_varincrements_levels.eps",
       plot = g, device = "eps", width = 9, height = 7)




