#===============================================================================
#                         Javier Viviens Martín
#                         Javier.VIVIENS@eui.eu
#                     European University Institute
#===============================================================================

#Clean the workspace:
rm(list = ls())  

#Set your working directory:


#Libraries:
if(!require(ggplot2))install.packages("ggplot2")
library(tidyverse)
if(!require(ggplot2))install.packages("ggplot2")
library(ggplot2)

#-------------------------------------------------------------------------------
#                         Functions
#-------------------------------------------------------------------------------
#Assume prior is a Gamma:
prior <- function(theta){
  prior_density <- dgamma(theta, shape = alpha, scale = sigma)
  return(prior_density)
}

#True DGP is an exponential. We want to estimate its parameter lambda (rate):
likelihood <- function(theta,data){
  likelihood <- prod(dexp(data, rate = theta))
  return(likelihood)
}

#Numerator of the posterior: likelihood * prior:
joint <- function(theta,data){
  joint <- likelihood(theta,data)*prior(theta)
  return(joint)
}

#Integrate over all posible values of theta to get the denominator, this is,
#the marginal distribution of the data.
marginal <- function(data){
  delta <- 0.01
  t <- 0
  density <- 0
  for (i in 1:1000){
    t <- t + delta
    theta <- 1 / t
    density <- density + joint(theta,data)
  }
  density <- density * delta
  return(density)
}

#Posterior is just numerator over denominator.
posterior <- function(theta, data){
  density <- joint(theta,data)/marginal(data)
  return(density)
}

#This function computes the density of the posterior for a given data. We need
#it only for the plot.
posterior_gd <- function(theta){
  density <- map(theta, posterior, data)
  density <- unlist(density)
  return(density)
}

#The following function computes the expected value of the posterior
posterior_mean <- function(data){
  theta <- 0
  mean <- 0
  delta <- 0.01
  for (i in 1:1000){
    theta <- theta + delta
    mean <- mean + theta*posterior(theta,data)
  }
  mean <- mean*delta/integrate(posterior_gd,support_lb,Inf)[[1]]
  return(mean)
}

#The following function normalize the density of the posterior to make it a proper
#distribution and make the graph look nicer. It's not strictly needed.
plot_posterior <- function(theta){
  posterior_gd(theta)/integrate(posterior_gd,support_lb,Inf)[[1]]
}

#-------------------------------------------------------------------------------
#                         Inference small sample
#-------------------------------------------------------------------------------
#Calibrate the parameters of the prior:
alpha <- 3
sigma <- 1
support_lb <- 0
#Load the data:
data <- as.numeric(unlist(read.csv("example_bayes_small.csv")))
#Make it global (for the plot):
#data_gd <<- data
#Compute expectations:
expected_posterior_s <- posterior_mean(data)
prior_m <- function(theta){theta*prior(theta)}
expected_prior <- integrate(prior_m,support_lb,Inf)[[1]]

#Plot
prior_graph <- ggplot() + xlim(0,10) + 
  geom_function(fun = prior, 
                aes(colour = "Prior")) +
  labs(colour = "") +
  scale_color_manual(values = c("mediumseagreen") ) +
  theme(axis.title=element_blank(), axis.text.y = element_blank(),
        legend.position = c(0.8, 0.8), 
        panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent') )  

print(prior_graph)

exponential_small <- ggplot() + xlim(0, 8) + 
  geom_function(fun = plot_posterior, 
                aes(colour = "Posterior")) +
  geom_function(fun = prior, 
                aes(colour = "Prior")) +
  geom_vline(aes(colour = "Posterior", 
                 xintercept = expected_posterior_s),linetype = 2 ) + 
  geom_vline(aes(colour = "Prior", 
                 xintercept = expected_prior), linetype = 2) + 
  labs(colour = "") + 
  scale_color_manual(values = c("dodgerblue4","mediumseagreen") ) +
  theme(axis.title=element_blank(), axis.text.y = element_blank(),
        legend.position = c(0.1, 0.8), 
        panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent') )  


print(exponential_small)
#ggsave(filename = "exponential_small.png", device = "png", dpi = 900)

#-------------------------------------------------------------------------------
#                         Inference larger sample
#-------------------------------------------------------------------------------
#Calibrate the parameters of the prior:
alpha <- 3
sigma <- 1
#Load the data:
data <- as.numeric(unlist(read.csv("example_bayes_large.csv")))
#Make it global (for the plot):
#data_gd <<- data
#Compute expectations:
expected_posterior_l <- posterior_mean(data)
prior_m <- function(theta){theta*prior(theta)}
expected_prior <- integrate(prior_m,support_lb,Inf)[[1]]

exponential_larger <- ggplot() + xlim(0, 8) + 
  geom_function(fun = plot_posterior, 
                aes(colour = "Posterior")) +
  geom_function(fun = prior, 
                aes(colour = "Prior")) +
  geom_vline(aes(colour = "Posterior", 
                 xintercept = expected_posterior_l),linetype = 2 ) + 
  geom_vline(aes(colour = "Prior", 
                 xintercept = expected_prior), linetype = 2) + 
  labs(colour = "") + 
  scale_color_manual(values = c("dodgerblue4","mediumseagreen") ) +
  theme(axis.title=element_blank(), axis.text.y = element_blank(),
        legend.position = c(0.1, 0.8), 
        panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent') )  


print(exponential_larger)
#ggsave(filename = "exponential_larger.png", device = "png", dpi = 900)

#-------------------------------------------------------------------------------
#                         Credible intervals
#-------------------------------------------------------------------------------
# NOTE: The code to compute the credible intervals may take a moment, so don't
#       if you don't get a result immediately.

#The following function computes the cdf of the posterior:
cumulative <- function(p){
  total_density <-  integrate(posterior_gd,support_lb,Inf)[[1]]
  s <- 0
  delta <- 0.01
  for (i in 1:1000){
    s <- s + delta
    temp <- integrate(posterior_gd,support_lb,s)[[1]] / total_density 
    if (temp >= p){
      break
    }
  }
  return(s)
}

credible_interval <- function(alpha){
  lb <- cumulative(alpha/2)
  ub <- cumulative(1-alpha/2)
  return(c(lb,ub))
}


data <- as.numeric(unlist(read.csv("example_bayes_small.csv")))
#data_gd <<- data
credible_interval_10 <- credible_interval(0.1)
credible_plot_s <- exponential_small + 
  geom_vline(aes(colour = "Credible Interval", 
                 xintercept = credible_interval_10[1]), linetype = 1) + 
  geom_vline(aes(colour = "Credible Interval", 
                 xintercept = credible_interval_10[2]), linetype = 1) + 
  scale_color_manual(values = c( "firebrick1","dodgerblue4","mediumseagreen") )

print(credible_plot_s)
#ggsave(filename = "exponential_small_ci.png", device = "png", dpi = 900)

data <- as.numeric(unlist(read.csv("example_bayes_large.csv")))
#data_gd <<- data
credible_interval_10 <- credible_interval(0.1)
credible_plot_l <- exponential_larger + 
  geom_vline(aes(colour = "Credible Interval", 
                 xintercept = credible_interval_10[1]), linetype = 1) + 
  geom_vline(aes(colour = "Credible Interval", 
                 xintercept = credible_interval_10[2]), linetype = 1) + 
  scale_color_manual(values = c( "firebrick1","dodgerblue4","mediumseagreen") )

print(credible_plot_l)
#ggsave(filename = "exponential_larger_ci.png", device = "png", dpi = 900)

