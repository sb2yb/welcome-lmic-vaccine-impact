# this script is combining the mal-ed and gems subclincal pathogen data for 0-2 years 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

load("Scripts - Aim 3/Output/subclinical_gems_0_2.Rdata")
load("Scripts - Aim 3/Output/subclinical_maled - 14March.Rdata")

# bind together
subclinical_gems_maled_0_2 <- cbind(subclinical_gems_0_2, subclinical_maled)

#save
save(subclinical_gems_maled_0_2,file="Scripts - Aim 3/Output/subclinical_gems_maled_0_2.Rdata")
