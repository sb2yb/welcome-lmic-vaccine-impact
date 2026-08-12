# This script combines RSV and shigella into one plot 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(png)
library(cowplot)

# Load the saved ggplot objects
Eq3_rsv <- readRDS("Scripts - Aim 3/Output/Eq3_rsv_DHS.rds")
Eq7_rsv <- readRDS("Scripts - Aim 3/Output/Eq7_rsv_DHS.rds")

Eq3_shig <- readRDS("Scripts - Aim 3/Output/Eq3_shigella_DHS.rds")
Eq7_shig <- readRDS("Scripts - Aim 3/Output/Eq7_shigella_DHS.rds")

# Combine them vertically
Eq3_combined <- plot_grid(Eq3_rsv, Eq3_shig, ncol = 1, labels = c("A", "B"))
Eq7_combined <- plot_grid(Eq7_rsv, Eq7_shig, ncol = 1, labels = c("A", "B"))

# View in RStudio
print(Eq3_combined)
print(Eq7_combined)

# Save to PDF
ggsave("Scripts - Aim 3/Maps/Eq3_RSV_Shigella.pdf", Eq3_combined, width = 12, height = 8) 
ggsave("Scripts - Aim 3/Maps/Eq7_RSV_Shigella.pdf", Eq7_combined, width = 12, height = 8) 

