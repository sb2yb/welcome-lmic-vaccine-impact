# this script will look at the different vaccine coverages 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(openxlsx)

# WUENIC 2024 coverage estimates 
coverage <- read.xlsx("Scripts - Aim 3/James data/james - vaccination coverage data.xlsx", sheet=1)

# needed for china below 
# unicef/who vaccine coverage for the last 3 vaccines china has introduced 
hepb <- read.xlsx("Scripts - Aim 3/Hepatitis B vaccination coverage 2025-06-10 10-19 UTC.xlsx", sheet = 1)
polio <- read.xlsx("Scripts - Aim 3/Poliomyelitis vaccination coverage 2025-06-10 10-39 UTC.xlsx", sheet = 1)
rubella <- read.xlsx("Scripts - Aim 3/Rubella vaccination coverage 2025-06-10 10-30 UTC.xlsx", sheet = 1)

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

# DEAL WITH IRAN 
# iran has 2% coverage for rotavirus making this a 0 
coverage_final <- coverage_final %>%
  mutate(ROTAC = ifelse(countriesSub == "IRN", 0, ROTAC))

# create a column that averages the values for rota, hib, and pcv for each country,excluding the 0s
coverage_final <- coverage_final %>%
  rowwise() %>%
  mutate(avg_vaxcov = mean(c(ROTAC, HIB3, PCV3)[c(ROTAC, HIB3, PCV3) != 0], na.rm = TRUE)) %>%
  ungroup()

# DEAL WITH NO COVERAGE IN CHINA 
hepb <- hepb %>%
  filter(COVERAGE_CATEGORY=="WUENIC" & ANTIGEN=="HEPB_BD")

polio <- polio %>%
  filter(COVERAGE_CATEGORY=="WUENIC" & (ANTIGEN=="IPV1" | ANTIGEN== "IPV2"))

rubella <- rubella %>%
  filter(COVERAGE_CATEGORY=="WUENIC")

hepb_polio_rubella <- rbind(hepb, polio, rubella)
hepb_polio_rubella <- hepb_polio_rubella %>%
  select(NAME, YEAR, ANTIGEN, COVERAGE)

# find the mean of the COVERAGE column
summary(hepb_polio_rubella$COVERAGE)

# add china's coverage to the coverage_final dataset 
coverage_final_c <- coverage_final %>%
  mutate(avg_vaxcov = case_when(
    countriesSub=="CHN" ~ 94.56,
    TRUE ~ avg_vaxcov))

# any country with missing data for a vaccine will get the averaged coverage data in that place 
coverage_final1 <- coverage_final_c %>%
  mutate(
    HIB3  = ifelse(is.na(HIB3)  | HIB3  == 0, avg_vaxcov, HIB3),
    PCV3  = ifelse(is.na(PCV3)  | PCV3  == 0, avg_vaxcov, PCV3),
    ROTAC = ifelse(is.na(ROTAC) | ROTAC == 0, avg_vaxcov, ROTAC))

# coverage needs to be expressed from 0-1 - dividing all by 100
coverage_final1 <- coverage_final1 %>%
  mutate(
    HIB3 = HIB3 / 100,
    PCV3 = PCV3 / 100,
    ROTAC = ROTAC / 100,
    avg_vaxcov = avg_vaxcov / 100)


# save dataframes 
write.csv(coverage_final1, file="Scripts - Aim 3/vaccine_coverage.csv")
save(coverage_final1,file="Scripts - Aim 3/vaccine_coverage.Rdata")

