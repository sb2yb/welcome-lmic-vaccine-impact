# Modifying DHS data to be what we want it to be

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

# load country list from Joe
load("Scripts - Aim 3/Joe data/ariAbx.Rdata")

# check array formatting
str(ariAbx)

# change name of diarAbx
DHS_ari <- ariAbx

# change bootstraps from 5,000 to 1,000
DHS_ari <- DHS_ari[1:1000,]
str(DHS_ari)

# save the IHME dataset 
save(DHS_ari,file="Scripts - Aim 3/Output/DHS_ari_matrix.Rdata")
