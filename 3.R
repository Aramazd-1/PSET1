rm(list = ls())
setwd(here::here())

# Load the data:
data <- read.table("Hbse.txt", header = TRUE, sep = "\t")  # Assuming tab-separated