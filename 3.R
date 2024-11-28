#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#       EXERCISE 3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Housekeeping
rm(list = ls())

if(!require(tidyverse))install.packages("tidyverse"); library(tidyverse)
if(!require(combinat))install.packages("combinat"); library(combinat)
if(!require(xtable))install.packages("xtable"); library(xtable)
if(!require(dplyr))install.packages("dplyr"); library(dplyr)
if(!require(ggplot2))install.packages("ggplot2"); library(ggplot2)
if(!require(MCMCpack))install.packages("MCMCpack"); library(MCMCpack)
if(!require(gridExtra))install.packages("gridExtra"); library(gridExtra)
if(!require(tidyr))install.packages("tidyr"); library(tidyr)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Load an explore the data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

Path = setwd("working//directory")
setwd(Path)
rm(Path)

data_orig <- read.table("Hbse.txt", header = T, sep = "")

# Rename columns to common notation of assignment, treatment & outcome
data <- data_orig
colnames(data) <- c("Z", "W", "Y")

# Here we are using functions from another file. You can ignore this part.
source("do_not_look_at_this.R")

## Explore the data a bit
data %>%
  dplyr::select(Z,W) %>% #Note: It seems that there is another function called "select()" in another package we used. To avoid ambiguity, we specify that we want the select() function from package dplyr
  ftable()
print(head(data_orig))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    3.1. Intention to treat effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## ITT of Z on Y
# Computed as a ratio 
ITT_Y <- with(data,cov(Y,Z)/var(Z))

# Computed from regression
ITT_Y_reg <- lm(Y~Z,data=data) # Note that this yields exactly the same estimate as the ITT computed as a ratio + it gives us the std. errors etc.

ITT_Y_table <- ITT_Y_reg %>%
  broom::tidy() %>%
  filter(term=="Z") %>% 
  dplyr::select(estimate, std.error)

stargazer::stargazer(ITT_Y_reg,
                     type = "latex", 
                     header=F,
                     keep.stat = "n",
                     digits=4,
                     title = "ITT Causal Effect of Assignment on Outcome")


## ITT of Z on W
# Computed as a ratio
ITT_W <- with(data,cov(W,Z)/var(Z))

# Computed from regression
ITT_W_reg <- lm(W~Z,data=data)

ITT_W_table <- ITT_W_reg %>%
  broom::tidy() %>%
  filter(term=="Z") %>% 
  dplyr::select(estimate, std.error)

stargazer::stargazer(ITT_W_reg,
                     type = "latex", 
                     header=F,
                     keep.stat = "n",
                     digits=4,
                     title = "ITT Causal Effect of Assignment on Treatment")


## Alternatively, you can compute ITTs "by hand"
summary_SZ <- data %>%
  summarise(count000 = sum(Z==0 & W==0 & Y==0),
            count001 = sum(Z==0 & W==0 & Y==1),
            count100 = sum(Z==1 & W==0 & Y==0),
            count101 = sum(Z==1 & W==0 & Y==1),
            count110 = sum(Z==1 & W==1 & Y==0),
            count111 = sum(Z==1 & W==1 & Y==1))
  
ITT_W_hand <- (summary_SZ$count111 + summary_SZ$count110) / (summary_SZ$count111 + summary_SZ$count110 + summary_SZ$count100 + summary_SZ$count101)
ITT_Y_hand <- (summary_SZ$count111 + summary_SZ$count101)/(summary_SZ$count100 + summary_SZ$count101 + summary_SZ$count110 + summary_SZ$count111) - (summary_SZ$count001)/(summary_SZ$count000 + summary_SZ$count001)
  
ITT_W_hand
ITT_Y_hand


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   3.2. Compute CACE / LATE
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## CACE 
# Computed as a ratio
CACE_IV_rat <- ITT_Y/ITT_W 
CACE_IV_rat


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   3.4.1 Run Gibbs Sampler & Plot CACE Histogram 
#                       with Exclusion Restriction
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set a random seed for reproducibility
set.seed(1900)

## Define the prior: assume flat priors
prior <- list(a0c=1,b0c=1, # beta distribution for untreated compliers 
              a1c=1,b1c=1, # beta distribution for treated compliers
              a0n=1,b0n=1, # beta distribution for untreated noncompliers
              a1n=1,b1n=1, # beta distribution for treated noncompliers
              ap=1,bp=1    # beta distribution for prob. of compliance
)

# Run Gibbs sampler under exclusion restriction
res_ER <- gibbs(data,prior,2000,excl_restr=T)

## Plot histogram of the CACE
ggER <- ggplot(res_ER, aes(x = cace_sp, fill = "CACE")) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.8, bins = 150) +
  scale_fill_manual(values = c("CACE" = "steelblue2")) +
  labs(x = "Value", y = "Density", fill="")

CACE_mean_ER <- mean(res_ER$cace_sp)
CACE_sd_ER <- sd(res_ER$cace_sp)

ggER + 
  stat_function(
    fun = dnorm,
    args = list(mean = CACE_mean_ER, sd = CACE_sd_ER),
    aes(color = "Normal Distribution (Super Population)"),
    size = 0.5, 
    color = "royalblue") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(y="density", x="",fill="") 
#ggsave("ex3.4.1.png")



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   3.4.2. Run Gibbs Sampler & Plot CACE Histogram 
#                           without Exclusion Restriction
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## Define the prior (same as before, not necessary but for completeness)
prior <- list(a0c=1,b0c=1, # beta distribution for untreated compliers 
              a1c=1,b1c=1, # beta distribution for treated compliers
              a0n=1,b0n=1, # beta distribution for untreated noncompliers
              a1n=1,b1n=1, # beta distribution for treated noncompliers
              ap=1,bp=1    # beta distribution for prob. of compliance
)

## Run Gibbs Sampler for relaxed exclusion restriction
res_NER <- gibbs(data,prior,2000,excl_restr=F)

## Plot histogram of the CACE
ggNER <- ggplot(res_NER, aes(x = cace_sp, fill = "CACE")) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.8, bins = 150) +
  scale_fill_manual(values = c("CACE" = "steelblue2")) +
  labs(x = "Value", y = "Density", fill="")

CACE_mean_NER <- mean(res_NER$cace_sp)
CACE_sd_NER <- sd(res_NER$cace_sp)

ggNER + 
  stat_function(
    fun = dnorm,
    args = list(mean = CACE_mean_NER, sd = CACE_sd_NER),
    aes(color = "Normal Distribution (Super Population)"),
    size = 0.5, 
    color = "royalblue") +
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(y="density", x="",fill="") 
#ggsave("ex3.4.2.png")


## Compare Histogram of CACE under ER and no ER 
# Adjust axis for better comparability
x_limits <- layer_scales(ggNER)$x$dimension()
y_limits <- layer_scales(ggER)$y$dimension()

# Apply the axis limits to ggER & ggNER 
ggNER <- ggNER +
  scale_y_continuous(limits = y_limits)

ggER <- ggER +
  scale_x_continuous(limits = x_limits)

# Display both histograms in one graph
grid.arrange(ggER, ggNER, nrow = 2)
#ggsave("ex3.4.png")

