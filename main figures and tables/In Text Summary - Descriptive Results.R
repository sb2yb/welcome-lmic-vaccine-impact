# this script has the results that will be presented in-text

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

# Eq`1 and Eq5
Diarrhea_Eq1_Eq5 <- read.csv("Scripts - Aim 3/Tables/SuppTable_DIAR_Eq1_Eq5_NOPATHOGENS_forINTEXTresults_18Dec2025.csv")
ARI_Eq1_Eq5 <- read.csv("Scripts - Aim 3/Tables/SuppTable_ARI_Eq1_Eq5_forINTEXTresults_18Dec2025.csv")


# Eq 1 and 5 - Diarrhea
Summary_Diarrhea_Eq1_Eq5 <- Diarrhea_Eq1_Eq5 %>%
  group_by(WHO_REGION, equation) %>%
  summarise(
    median_value = round(median(median, na.rm = TRUE), 2),
    iqr_lower = round(quantile(median, 0.25, na.rm = TRUE), 2),
    iqr_upper = round(quantile(median, 0.75, na.rm = TRUE), 2),
    .groups = "drop")


# Eq 1 and 5 - Respiratory
Summary_ARI_Eq1_Eq5 <- ARI_Eq1_Eq5 %>%
  group_by(WHO_REGION, equation) %>%
  summarise(
    median_value = round(median(median, na.rm = TRUE), 2),
    iqr_lower = round(quantile(median, 0.25, na.rm = TRUE), 2),
    iqr_upper = round(quantile(median, 0.75, na.rm = TRUE), 2),
    .groups = "drop")

write.csv(Summary_Diarrhea_Eq1_Eq5, file="Scripts - Aim 3/Tables/IN TEXT MEDIAN SUMMARY - Diarrhea_Eq1_Eq5_18Dec2025.csv")
write.csv(Summary_ARI_Eq1_Eq5, file="Scripts - Aim 3/Tables/IN TEXT MEDIAN SUMMARY - ARI_Eq1_Eq5_18Dec2025.csv")
