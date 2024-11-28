rm(list = ls())
setwd(here::here())
library(tidyverse)
library(ggplot2)
library(TeachingDemos)
#Import functions
source("Functions.R")


# ------------------------------------------------------------------------------
# 2.1
# ------------------------------------------------------------------------------
# Shape parameters for Gamma Distribution:
alpha <- 2
sigma <- 1
support_lb <- 0
#Load the data:
data <- as.numeric(unlist(read.csv("gasolina_small_sample.csv")))

# Expectations:
expected_posterior_s <- posterior_mean(data)
prior_m <- function(theta){theta*prior(theta)}
expected_prior <- integrate(prior_m,support_lb,Inf)[[1]]
credible_interval_10 <- credible_interval(0.1)
# posterior_a <- plot_posterior(seq(0, 8, length.out = 1000)) # Sanity check

alpha_posterior <- alpha + length(data) # \alpha + n
beta_posterior <- sigma + sum(data) # \beta + \sum_{i=1}^{n} x_i
posterior_b <- dgamma(seq(0, 8, length.out = 1000), shape = alpha_posterior, rate = beta_posterior)
posterior_icdf <- function(p) qgamma(p, shape = alpha_posterior, rate = beta_posterior)

# Compute the HPD interval at 95% confidence level
hpd_interval <- hpd(posterior_icdf, conf = 0.90, tol = 1e-8)
cat("HPD Interval:", hpd_interval, "\n")
cat("Credible Interval:", credible_interval_10, "\n")

gasolina_small <- ggplot() + xlim(0, 8) +
  geom_function(fun = plot_posterior,
                aes(colour = "Posterior")) +
  geom_function(fun = prior,
                aes(colour = "Prior")) +
  geom_vline(aes(colour = "Posterior Expected",
                 xintercept = expected_posterior_s), linetype = 2 ) + # Posterior Expected
  geom_vline(aes(colour = "Prior Expected",
                 xintercept = expected_prior), linetype = 2) +        # Prior Expected
  geom_vline(aes(colour = "Equal Tail Credible Interval",
                 xintercept = credible_interval_10[1]), linetype = 1) +
  geom_vline(aes(colour = "Equal Tail Credible Interval",
                 xintercept = credible_interval_10[2]), linetype = 1) +
  geom_vline(aes(colour = "HPD Interval",
                 xintercept = hpd_interval[1]), linetype = 1) +
  geom_vline(aes(colour = "HPD Interval",
                 xintercept = hpd_interval[2]), linetype = 1) +
  labs(colour = "Legend") + # Add a legend title
  scale_color_manual(values = c(
    "Posterior" = "dodgerblue4",
    "Prior" = "mediumseagreen",
    "Equal Tail Credible Interval" = "purple",
    "HPD Interval" = "red",
    "Posterior Expected" = "blue",
    "Prior Expected" = "green"
  )) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.1, 0.8),
        panel.background = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent'))

ggsave("gasolina_small_plot.png", plot = gasolina_small, width = 8, height = 6, dpi = 300)

print(gasolina_small)


