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


# CI

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

