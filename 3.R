rm(list = ls())
setwd(here::here())

library(tidyverse)
library(combinat)
library(xtable)
library(dplyr)
library(ggplot2)
library(MCMCpack)
library(gridExtra)
library(tidyr)
library(ivreg)

# Load the data:
data_texto <- read.table("Hbse.txt", header = TRUE, sep = "")
head(data_texto)

data <- data_texto
colnames(data) <- c("Z", "W", "Y")

data %>%
  dplyr::select(Z,W) %>%
  ftable()

## ITT of Z on Y
ITT_Y <- with(data,cov(Y,Z)/var(Z))

# Computed from regression
ITT_Y_reg <- lm(Y~Z,data=data) # Note that this yields exactly the same estimate as the ITT computed as a ratio + it gives us the std. errors etc.

ITT_Y_table <- ITT_Y_reg %>%
  broom::tidy() %>%
  filter(term=="Z") %>%
  dplyr::select(estimate, std.error)

stargazer::stargazer(ITT_Y_reg,
                     type = "text",
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
                     type = "text",
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

## CACE
# Computed as a ratio
CACE_IV_rat <- ITT_Y/ITT_W

# Computed from IV regression (taking the assignment as an instrument)
CACE_IV_reg <- ivreg(Y~W|Z,data=data)

CACE_IV_table <- CACE_IV_reg %>%
  broom::tidy() %>%
  filter(term=="W") %>%
  dplyr::select(estimate, std.error)

stargazer::stargazer(CACE_IV_reg,
                     type = "text",
                     header=F,
                     keep.stat = "n",
                     digits=4,
                     title = "Compliers Average Causal Effect")
