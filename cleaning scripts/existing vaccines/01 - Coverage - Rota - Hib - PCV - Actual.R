# this script takes the actual data for rota, hib, and PCV

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(openxlsx)

# vaccine coverage provided by joe - NO LONGER USING VACCINE COVERAGE DATA FROM JOE 
# coverage <- read.csv("Scripts - Aim 3/Joe data/averaged.csv")

# WUENIC 2024 coverage estimates 
coverage <- read.xlsx("Scripts - Aim 3/James data/james - vaccination coverage data.xlsx", sheet=1)

# load country region list
load("Scripts - Aim 3/country_region_list.Rdata")


# clean coverage data 
# Nicaragua doesn't have WUENIC estimates - taking the ADMIN estimates for them 
coverage1 <- coverage %>%
  filter(YEAR == 2024 & ((COVERAGE_CATEGORY == "WUENIC" & (ANTIGEN == "HIB3" | ANTIGEN == "PCV3" | ANTIGEN == "ROTAC") & CODE != "NIC") |
          (COVERAGE_CATEGORY == "ADMIN" & (ANTIGEN == "HIB3" | ANTIGEN == "PCV3" | ANTIGEN == "ROTAC") & CODE == "NIC")))

coverage2 <- coverage1 %>%
  select(CODE, ANTIGEN, COVERAGE) %>%
  pivot_wider(
    names_from = ANTIGEN,
    values_from = COVERAGE)

# remove the 6 countries that we do not need: 
# Democratic Peoples Republic of Korea, Palau, Argentina, Occupied Palestinian territory, Syria, duplicate of Georgia
country_region_list <- country_region_list %>%
  filter(DisplayString!="Democratic People's Republic of Korea" & DisplayString!="Palau" & DisplayString!="Argentina" & 
           DisplayString!="occupied Palestinian territory, including east Jerusalem" & DisplayString!="Syrian Arab Republic") %>%
  select(countriesSub, DisplayString) %>%
  distinct() # this removes one instance of Georgia (doesn't matter which one since this is only a list and not tied to actual data)

# coverage data needs to be merged into country data so that we can extract the data we want on our 129 countries 
# left_join to preserve the country order for country_region_list
coverage_final <- country_region_list %>% 
  left_join(coverage2, by = c("countriesSub" = "CODE"))

summary(coverage_final$ROTAC)
summary(coverage_final$HIB3)
summary(coverage_final$PCV3)

# make NAs = 0 
coverage_final <- coverage_final %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# coverage needs to be expressed from 0-1 - dividing all by 100
coverage_final1 <- coverage_final %>%
  mutate(
    HIB3 = HIB3 / 100,
    PCV3 = PCV3 / 100,
    ROTAC = ROTAC / 100)


# save dataframes 
write.csv(coverage_final1, file="Scripts - Aim 3/vaccine_coverage_Rota_Hib_PCV_ACTUAL.csv")
save(coverage_final1,file="Scripts - Aim 3/vaccine_coverage_Rota_Hib_PCV_ACTUAL.Rdata")

# For the manuscript, i need to include the proportion of countries which have introduced each vaccines 
# Countries with <5% coverage were not considered as having "introduced" the vaccine 

sum(coverage_final1$HIB3 > 0.05, na.rm = TRUE) 
128/129*100
sum(coverage_final1$PCV3 > 0.05, na.rm = TRUE) 
102/129*100
sum(coverage_final1$ROTAC > 0.05, na.rm = TRUE) 
87/129*100











