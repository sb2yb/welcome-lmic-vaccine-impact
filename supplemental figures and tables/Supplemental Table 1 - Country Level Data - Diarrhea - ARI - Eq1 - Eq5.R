# combining the ARI and diarrhea supplementary tables for Eq1 and Eq5 into one table 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

# load scripts 
diar_eq1_eq5 <- read.csv("Scripts - Aim 3/Tables/SuppTable_DIAR_Eq1_Eq5_NOPATHOGENS_18Dec2025.csv")
ari_eq1_eq5 <- read.csv("Scripts - Aim 3/Tables/SuppTable_ARI_Eq1_Eq5_18Dec2025.csv")

# rename columns
diar_eq1_eq5 <- diar_eq1_eq5 %>%
  rename(diarrhea_Eq1 = Eq1, 
         diarrhea_Eq5 = Eq5)

ari_eq1_eq5 <- ari_eq1_eq5 %>%
  rename(ARI_Eq1 = Eq1, 
         ARI_Eq5 = Eq5)

# bind diarrhea and ARI together 
diar_ari_Eq1_Eq5 <- cbind(diar_eq1_eq5, ari_eq1_eq5)

# keep only the columns we want
diar_ari_Eq1_Eq5 <- diar_ari_Eq1_Eq5[ , -c(1, 5, 6)]

# save the table
write.csv(diar_ari_Eq1_Eq5, file="Scripts - Aim 3/Tables/Supp Table 1 - DIAR_ARI_Eq1_Eq5.csv")

