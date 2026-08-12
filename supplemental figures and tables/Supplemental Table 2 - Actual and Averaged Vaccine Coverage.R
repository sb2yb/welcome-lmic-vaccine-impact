# this script is going to combine the actual coverage data with the averaged coverage data and organize it by WHO region 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(stringr)

# load country region list
load("Scripts - Aim 3/country_region_list.Rdata")

# actual coverage data
coverage_actual <- read.csv("Scripts - Aim 3/vaccine_coverage_Rota_Hib_PCV_ACTUAL.csv")

# averaged coverage data 
coverage_average <- read.csv("Scripts - Aim 3/vaccine_coverage.csv")

# remove the 6 countries that we do not need: 
# Democratic Peoples Republic of Korea, Palau, Argentina, Occupied Palestinian territory, Syria, duplicate of Georgia
country_region_list <- country_region_list %>%
  filter(DisplayString!="Democratic People's Republic of Korea" & DisplayString!="Palau" & DisplayString!="Argentina" & 
           DisplayString!="occupied Palestinian territory, including east Jerusalem" & DisplayString!="Syrian Arab Republic") %>%
  select(DisplayString, WHO_REGION) %>%
  distinct() # this removes one instance of Georgia (doesn't matter which one since this is only a list and not tied to actual data)

###################################
# bind the coverage data together #
###################################

coverage <- cbind(coverage_actual, coverage_average)

# keep what we need and rename variables
coverage <- coverage[, c(2:6, 10:13)]

coverage <- coverage %>%
  rename(hib_actual = HIB3,
         pcv_actual = PCV3,
         rota_actual= ROTAC,
         hib_est = HIB3.1,
         pcv_est = PCV3.1, 
         rota_est = ROTAC.1)

# merge WHO region data into the coverage data 
country_coverage <- left_join(country_region_list, coverage, by = "DisplayString")

# clean the name of 14 countries to match the world country list
country_coverage <- country_coverage %>%
  mutate(DisplayString = str_replace(DisplayString, "Bolivia \\(Plurinational State of\\)", "Bolivia")) %>%
  mutate(DisplayString = ifelse(DisplayString== "Congo", "Republic of the Congo",DisplayString)) %>%
  mutate(DisplayString = str_replace(DisplayString, "Cote d'Ivoire", "Ivory Coast")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Eswatini", "eSwatini")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Iran \\(Islamic Republic of\\)", "Iran")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Lao People's Democratic Republic", "Laos")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Micronesia \\(Federated States of\\)", "Federated States of Micronesia")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Republic of Moldova", "Moldova")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Russian Federation", "Russia")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Sao Tome and Principe", "São Tomé and Principe")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Serbia", "Republic of Serbia")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Timor-Leste", "East Timor")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Turkiye", "Turkey")) %>%
  mutate(DisplayString = str_replace(DisplayString, "Viet Nam", "Vietnam"))


########################################
# get the data into the desired format #
########################################

# set order of regions
region_levels <- c("Africa", "Americas", "Eastern Mediterranean", "Europe", "South-East Asia", "Western Pacific")  # whatever order you prefer

country_coverage <- country_coverage %>%
  mutate(WHO_REGION = factor(WHO_REGION, levels = region_levels))

# Sort alphabetically by WHO_REGION and then by country name
country_coverage <- country_coverage %>%
  arrange(WHO_REGION, DisplayString)

# round to 2 decimal places
country_coverage[, 4:10] <- round(country_coverage[, 4:10], 2)

# save dataframes 
write.csv(country_coverage, file="Scripts - Aim 3/Tables/Supp Table 2 - Actual and Averaged Vaccine Coverage.csv")

# vaccine coverage 
table(coverage_actual$HIB3>0.05)
table(coverage_actual$PCV3>0.05)
table(coverage_actual$ROTAC>0.05)

