# Modifying DHS data to be what we want it to be

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

# load country list from Joe
load("Scripts - Aim 3/Joe data/diarAbx.Rdata")

# check array formatting
str(diarAbx)

# change name of diarAbx
DHS_diar <- diarAbx

# change bootstraps from 5,000 to 1,000
DHS_diar <- DHS_diar[1:1000,]
str(DHS_diar)

# save the IHME dataset 
save(DHS_diar,file="Scripts - Aim 3/Output/DHS_diar_matrix.Rdata")
