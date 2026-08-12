# this script takes this dataframe: bystander_region_estimates_0_5_years_final and this country list: country_region_list
# and applies these region specific estimates to each country within the region

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(janitor)

load("Scripts - Aim 3/Output/bystander_region_estimates_0_5_years_final.Rdata")
load("Scripts - Aim 3/country_region_list.Rdata")

# rename columns
bystander_region_estimates_0_5_years_final <- rename(bystander_region_estimates_0_5_years_final, AFR = AFRO_0_5, AMR = AMRO_0_5, EMR = EMRO_0_5, EUR = EURO_0_5, SEAR = SEARO_0_5, WPR = WPRO_0_5)

bystander <- t(bystander_region_estimates_0_5_years_final)

bystander1 <- bystander
WHO_REGION_CODE <- rownames(bystander1)
rownames(bystander1) <- NULL
bystander2 <- cbind(WHO_REGION_CODE,bystander1)
bystander2 <- as.data.frame(bystander2)

# slim down country list 
country <- country_region_list %>%
  select(countriesSub, WHO_REGION_CODE)

# merge
bystander_country_estimates <- left_join(country, bystander2, by="WHO_REGION_CODE") # left_join prevents it from sorting
bystander_country_estimates <- bystander_country_estimates[,c(2,1,3:1002)]

bystander_country_estimates_final <- t(bystander_country_estimates)

# remove first 2 rows
bystander_country_estimates_final <- bystander_country_estimates_final[-c(1, 2), ]

# save matrix
save(bystander_country_estimates_final,file="Scripts - Aim 3/Output/bystander_country_estimates_final_matrix.Rdata")
