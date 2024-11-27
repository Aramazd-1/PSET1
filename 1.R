setwd(here::here())
library(tidyverse)
getwd()

# ------------------------------------------------------------------------------
# 1.1 Exact Fisher P-values for Absolute Difference in Sample Averages
# ------------------------------------------------------------------------------
# Read the data and  Check the first few rows of the data
data <- read.csv("data_ex1.csv")
head(data)

obs_diff <- abs(mean(data$Y[data$W == 1]) - mean(data$Y[data$W == 0]))
N <- nrow(data)
N_treated <- sum(data$W)  # Number of treated units
N_control <- N - N_treated  # Number of control units


# Function to generate all possible treatment assignment vectors (W)
# Arguments:
#   N - Total number of units
#   N_treated - Number of treated units
# Uses:
#   combn(N, N_treated): Generates all possible combinations of indices for treated units.
#   function(x): For each combination of treated indices (x), creates a binary vector W
#                where W[i] = 1 if unit i is treated, and W[i] = 0 otherwise.
# Returns:
#   A list of all possible treatment assignment vectors.
generate_assignments <- function(N, N_treated) {
  combn(N, N_treated, function(x) {
    W <- rep(0, N)
    W[x] <- 1
    return(W)
  }, simplify = FALSE)
}
all_assignments <- generate_assignments(N, N_treated)

# Calculate the test statistic for each possible assignment
test_stats <- sapply(all_assignments, function(W_perm) {
  diff_mean <- abs(mean(data$Y[W_perm == 1]) - mean(data$Y[W_perm == 0]))
  return(diff_mean)
})

p_value <- mean(test_stats >= obs_diff)
cat("Exact Fisher p-value:", p_value, "\n")
# ------------------------------------------------------------------------------
# 1.2 Estimation of Finite Sample Average Treatment Effect (ATE)
# ------------------------------------------------------------------------------

# Separate the outcome variable by treatment status
Y_treated <- data$Y[data$W == 1]
Y_control <- data$Y[data$W == 0]

# Sample sizes and means
n_treated <- length(Y_treated)
n_control <- length(Y_control)

mean_treated <- mean(Y_treated)
mean_control <- mean(Y_control)

# Estimate of the ATE (difference in means)
ATE_hat <- mean_treated - mean_control

# Sample variances
s2_treated <- var(Y_treated)
s2_control <- var(Y_control)

# Estimate of the variance of the ATE estimator (design-based variance)
Var_ATE_hat <- s2_treated / n_treated + s2_control / n_control

# Standard error of the ATE estimator
SE_ATE_hat <- sqrt(Var_ATE_hat)

# 95% confidence interval
alpha_95 <- 0.05
t_alpha_95 <- qt(1 - alpha_95 / 2, df = 9)
CI_lower_95 <- ATE_hat - t_alpha_95 * SE_ATE_hat
CI_upper_95 <- ATE_hat + t_alpha_95 * SE_ATE_hat

# 99% confidence interval
alpha_99 <- 0.01
t_alpha_99 <- qt(1 - alpha_99/2, df = 9)
CI_lower_99 <- ATE_hat - t_alpha_99 * SE_ATE_hat
CI_upper_99 <- ATE_hat + t_alpha_99 * SE_ATE_hat

# Output the results
cat("Estimated ATE:", ATE_hat, "\n")
cat("Estimated variance of ATE estimator:", Var_ATE_hat, "\n")
cat("Standard error:", SE_ATE_hat, "\n")
cat("95% confidence interval: [", CI_lower_95, ",", CI_upper_95, "]\n")
cat("99% confidence interval: [", CI_lower_99, ",", CI_upper_99, "]\n")