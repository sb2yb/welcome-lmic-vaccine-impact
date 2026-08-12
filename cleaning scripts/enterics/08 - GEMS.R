# The "gemstac" R data file has 4 datasets within it. The only we need is "data". This script is saving "data" as an Rda file. 
      # stddata is the same dataset but with standard curves applied, we ultimately went away from this approach.
      # gemsconv is the full gems original dataset. If you had need additional variables, that would be the source but sounds like you don’t.
      # shigdata is updated Shigella typing data that have already been incorporated in “data”

# This will also serve as the codebook for the GEMS data file. 
      # If you filter to case.control==0 and valid_sample==1 you should get 5447 controls that you can use for this.
      #“country” will get you site
      #“age” will get you the age stratum for GEMS where 1 = <12 months, 2 = 12-23 months, and 3 = 24-59 months
      # If you need age more granularly, “agemonths” is age in months but I think is not derived for controls. 
              # So you can either derived from the matching case (using “Case.ID”) but this wasn’t exact, within 2-4 months depending on the stratum, 
              # so we may need to dig for DOB/enrollment date and re-derive. Depends what you need.
      # TAC pathogens should look reasonably familiar, start on column 6 and end on column 80.

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

load("Scripts - Aim 3/James data/gemstac.Rdata")
# save the "data" df as is
# save(data,file="Scripts - Aim 3/James data/data.Rda")

# filter to controls and valid samples
GEMS_data <- data %>%
  filter(case.control==0 & valid_sample==1) 

# keep only the needed variables 
# this will have all ages
GEMS_data <- GEMS_data %>%
  select(Case.ID, age, case.control, country, adenovirus:STEC)

# this will have ages 0-2
GEMS_data_0_2 <- GEMS_data %>%
  filter(age==1 | age==2)

# this will have ages 3-5
GEMS_data_3_5 <- GEMS_data %>%
  filter(age==3)

# save the GEMS dataset that we will use for the bystander piece
save(GEMS_data,file="Scripts - Aim 3/James data/GEMS_data.Rda")
save(GEMS_data_0_2,file="Scripts - Aim 3/James data/GEMS_data_0_2.Rda")
save(GEMS_data_3_5,file="Scripts - Aim 3/James data/GEMS_data_3_5.Rda")
