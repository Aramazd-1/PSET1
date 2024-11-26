rm(list = ls())
setwd(here::here())
library(tidyverse)
library(ggplot2)
library(HDInterval)

#Import functions
source("Functions.R")
# ------------------------------------------------------------------------------
# 2.1
# ------------------------------------------------------------------------------
# Shape parameters of the Gamma prior:
alpha <- 2
sigma <- 1
support_lb <- 0
# Load the data:
data <- as.numeric(unlist(read.csv("gasolina_small_sample.csv")))

# Compute the posterior mean E[λ | Y]:
expected_posterior <- posterior_mean(data)
prior_m <- function(theta){theta * prior(theta)}
expected_prior <- integrate(prior_m, support_lb, Inf)[[1]]

# Print the posterior mean:
cat("Posterior mean E[λ | Y]:", expected_posterior, "\n")

# Compute two different 90% credible intervals:
# First credible interval (using the credible_interval function):
credible_interval_10 <- credible_interval(0.10)
cat("First 90% credible interval (using credible_interval function):\n")
cat("Lower bound:", credible_interval_10[1], "\n")
cat("Upper bound:", credible_interval_10[2], "\n")

# Second credible interval :
# Compute posterior parameters
n <- length(data)
sum_data <- sum(data)
alpha_post <- alpha + n
beta_post <- 1 / sigma + sum_data  # Rate parameter of the posterior Gamma distribution

# Generate samples from the posterior Gamma distribution
posterior_samples <- rgamma(10000, shape = alpha_post, rate = beta_post)
# Compute the 90% highest density interval (HDI)
credible_interval_hdi <- hdi(posterior_samples, credMass = 0.90)

# Print the second credible interval
cat("Second 90% credible interval (using HDInterval package):\n")
cat("Lower bound:", credible_interval_hdi[1], "\n")
cat("Upper bound:", credible_interval_hdi[2], "\n")

# ------------------------------------------------------------------------------
# Plotting the prior and posterior distributions with credible intervals
# ------------------------------------------------------------------------------

# Define a range of lambda values for plotting
lambda_values <- seq(0, expected_posterior * 3, length.out = 1000)
prior_density_values <- dgamma(lambda_values, shape = alpha, scale = sigma)

# Compute the posterior density values
# Posterior density using the analytical form of the Gamma distribution
posterior_density_values <- dgamma(lambda_values, shape = alpha_post, rate = beta_post)
df_plot <- data.frame(lambda = lambda_values,
                      Prior = prior_density_values,
                      Posterior = posterior_density_values)
df_melted <- df_plot %>% gather(key = "Distribution", value = "Density", -lambda)

plotto <- ggplot(df_melted, aes(x = lambda, y = Density, color = Distribution)) +
  geom_line(size = 1) +
  geom_vline(xintercept = credible_interval_10, color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(credible_interval_hdi[1], credible_interval_hdi[2]), color = "green", linetype = "dotted") +
  labs(title = "Prior and Posterior Distributions of Lambda with Credible Intervals",
       x = expression(lambda),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "darkorange")) +
  annotate("text", x = credible_interval_10[1], y = max(df_melted$Density) * 0.8, label = "First CI", color = "red", hjust = -0.1, vjust = 1) +
  annotate("text", x = credible_interval_hdi[1], y = max(df_melted$Density) * 0.7, label = "Second CI", color = "green", hjust = -0.1, vjust = 1)
print(plotto)

# 2.2
# Load the new data:
# Shape parameters of the Gamma prior:
alpha <- 2
sigma <- 1
support_lb <- 0
# Load the data:
data <- as.numeric(unlist(read.csv("gasolina_large_sample.csv")))

# Compute the posterior mean E[λ | Y]:
expected_posterior <- posterior_mean(data)
prior_m <- function(theta){theta * prior(theta)}
expected_prior <- integrate(prior_m, support_lb, Inf)[[1]]

# Print the posterior mean:
cat("Posterior mean E[λ | Y]:", expected_posterior, "\n")

# Credible interval (using the credible_interval function):
credible_interval_10 <- credible_interval(0.10)
cat("First 90% credible interval (using credible_interval function):\n")
cat("Lower bound:", credible_interval_10[1], "\n")
cat("Upper bound:", credible_interval_10[2], "\n")

# Define a range of lambda values for plotting
lambda_values <- seq(0, expected_posterior * 3, length.out = 1000)
prior_density_values <- dgamma(lambda_values, shape = alpha, scale = sigma)

# Compute the posterior density values
# Posterior density using the analytical form of the Gamma distribution
posterior_density_values <- dgamma(lambda_values, shape = alpha_post, rate = beta_post)

df_plot <- data.frame(
  lambda = lambda_values,
  Prior = prior_density_values,
  Posterior = posterior_density_values
)
df_melted <- tidyr::pivot_longer(df_plot, cols = c("Prior", "Posterior"), names_to = "Distribution", values_to = "Density")

plotto2 <- ggplot(df_melted, aes(x = lambda, y = Density, color = Distribution)) +
  geom_line(size = 1) +
  geom_vline(xintercept = credible_interval_10, color = "red", linetype = "dashed") +
  labs(title = "Prior and Posterior Distributions of Lambda with Credible Interval",
       x = expression(lambda),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "darkorange")) +
  annotate("text", x = credible_interval_10[1], y = max(df_plot$Posterior) * 0.8,
           label = "Lower Bound", color = "red", hjust = -0.1, vjust = 1) +
  annotate("text", x = credible_interval_10[2], y = max(df_plot$Posterior) * 0.8,
           label = "Upper Bound", color = "red", hjust = 1.1, vjust = 1)

print(plotto2)

# 2.3
# Shape parameters of the Gamma prior:
alpha <- 1
sigma <- 1
support_lb <- 0
# Load the data:
data <- as.numeric(unlist(read.csv("gasolina_large_sample.csv")))

# Compute the posterior mean E[λ | Y]:
expected_posterior <- posterior_mean(data)
prior_m <- function(theta){theta * prior(theta)}
expected_prior <- integrate(prior_m, support_lb, Inf)[[1]]

# Print the posterior mean:
cat("Posterior mean E[λ | Y]:", expected_posterior, "\n")

# Credible interval (using the credible_interval function):
credible_interval_10 <- credible_interval(0.10)
cat("First 90% credible interval (using credible_interval function):\n")
cat("Lower bound:", credible_interval_10[1], "\n")
cat("Upper bound:", credible_interval_10[2], "\n")

# Define a range of lambda values for plotting
lambda_values <- seq(0, expected_posterior * 3, length.out = 1000)
prior_density_values <- dgamma(lambda_values, shape = alpha, scale = sigma)

# Compute the posterior density values
# Posterior density using the analytical form of the Gamma distribution
posterior_density_values <- dgamma(lambda_values, shape = alpha_post, rate = beta_post)

df_plot <- data.frame(
  lambda = lambda_values,
  Prior = prior_density_values,
  Posterior = posterior_density_values
)
df_melted <- tidyr::pivot_longer(df_plot, cols = c("Prior", "Posterior"), names_to = "Distribution", values_to = "Density")

plotto3 <- ggplot(df_melted, aes(x = lambda, y = Density, color = Distribution)) +
  geom_line(size = 1) +
  geom_vline(xintercept = credible_interval_10, color = "red", linetype = "dashed") +
  labs(title = "Prior and Posterior Distributions of Lambda with Credible Interval",
       x = expression(lambda),
       y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "darkorange")) +
  annotate("text", x = credible_interval_10[1], y = max(df_plot$Posterior) * 0.8,
           label = "Lower Bound", color = "red", hjust = -0.1, vjust = 1) +
  annotate("text", x = credible_interval_10[2], y = max(df_plot$Posterior) * 0.8,
           label = "Upper Bound", color = "red", hjust = 1.1, vjust = 1)

print(plotto3)