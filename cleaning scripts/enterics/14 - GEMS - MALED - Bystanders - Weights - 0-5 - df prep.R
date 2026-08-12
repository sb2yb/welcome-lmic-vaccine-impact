# this script is to set up the weights calculations
# where the bystander pathogen exposures will be calculated by region 
# first, calculate weights for kids 0-2 in GEMS and MALED
# second, calculate weights for kids 3-5 in GEMS
# third, extract the number of rows that should be pulled from each site
# finally, combine the weights and end up with 1 weight for all kids 0-5 = the average number of bystanders per abx course by region

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(janitor)

load("Scripts - Aim 3/Output/subclinical_gems_maled_0_2.Rdata")
load("Scripts - Aim 3/Output/subclinical_gems_3_5.Rdata")
load("Scripts - Aim 3/James data/WHO region draw weights for Wellcome AMR.Rdata")

################
# for ages 0-2 #
################

# multiply GPDS weights by 1000 and round to whole numbers
weights_maled_gems_1000 <- L2Weights %>%
  mutate(across(where(is.numeric), ~ round(.x * 1000,0)))

# sum the rows - we will use this a check. 
#   if there is a row with less than 1000 (due to rounding), then we need to randomly sample additional bootstraps
#   if there is a row with more than 1000 (due to rounding), then we need to randomly remove bootstraps
#   AFRO is 1001, WPRO is 1001, and EMRO is 998
weights_maled_gems_1000_sum <- weights_maled_gems_1000 %>%
  mutate(sum = rowSums(across(where(is.numeric)))) 

# random.org will decide which site will have a bootstrap added or removed
#   AFRO: random.org = 4 = bg, which means 1 bootstrap from bg will be removed - instead of 86 it will now be 85
#   WPRO: random.org = 5 = SA, which means 1 bootstrap from SA will be removed - instead of 62 it will now be 61
#   EMRO: random.org = 2 = mz, which means 1 bootstrap from mz will be added - instead of 42 it will now be 43
#   EMRO: random.org = 1 = gm, which means 1 bootstrap from gm will be added - instead of 171 it will now be 172

# build in the above to the final dataset for ages 0-2
# remove sums column 
weights_maled_gems_1000_final <- weights_maled_gems_1000_sum[,-17]
weights_maled_gems_1000_final <- weights_maled_gems_1000_final %>%
  mutate(bg = ifelse(bg==86,85,bg)) %>%
  mutate(SA = ifelse(SA==62,61,SA)) %>%
  mutate(mz = ifelse(mz==42,43,mz)) %>%
  mutate(gm = ifelse(gm==171,172,gm)) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) # recalculate sum (as a check)

# remove sums column  
# this is the final dataset for ages 0-2 years
weights_maled_gems_1000_final <- weights_maled_gems_1000_final[,-17]


################
# for ages 3-5 #
################

# we need to:
# 1) remove maled data 
weights_gems_1000 <- L2Weights[,c(1:8)]

# 2) make the sum of the row without maled equal to 1 (because that's what it was with MAL-ED + GEMS)
# we will sum up the row and divide each number in that row by the sum of the row
weights_gems_1000 <-weights_gems_1000 %>%
  mutate(sum = rowSums(across(where(is.numeric)))) # sum across rows

weights_gems_1000 <- weights_gems_1000 %>% # override the varibles in each column where the new value is now the old value/row sum
  mutate(across(where(is.numeric), ~ .x / sum)) # the "sum" column should now total 1
  
# multiply weights by 1000 and round to whole numbers
weights_gems_1000 <- weights_gems_1000 %>%
  mutate(across(where(is.numeric), ~ round(.x * 1000,0)))  

# remove sum column  
weights_gems_1000 <- weights_gems_1000[,-9]

# sum the rows - we will use this a check. 
#   if there is a row with less than 1000 (due to rounding), then we need to randomly sample additional bootstraps
#   if there is a row with more than 1000 (due to rounding), then we need to randomly remove bootstraps
#   AFRO is 999, EURO is 999, and WPRO is 1001
weights_gems_1000_sum <- weights_gems_1000 %>%
  mutate(sum = rowSums(across(where(is.numeric)))) 

# random.org will decide which site will have a bootstrap added or removed
#   AFRO: random.org = 1 = gm, which means 1 bootstrap from gm will be added - instead of 21 it will now be 22
#   EURO: random.org = 3 = mz, which means 1 bootstrap from mz will be added - instead of 54 it will now be 55
#   WPRO: random.org = 1 = gm, which means 1 bootstrap from gm will be removed - instead of 4 it will now be 3

# build in the above to the final dataset for ages 3-5
# remove sums column 
weights_gems_1000_sum_final <- weights_gems_1000_sum[,-9]
weights_gems_1000_sum_final <- weights_gems_1000_sum_final %>%
  mutate(gm = ifelse(gm==21,22,gm)) %>%
  mutate(mz = ifelse(mz==54,55,mz)) %>%
  mutate(gm = ifelse(gm==4,3,gm)) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) # recalculate sum (as a check)

# remove sums column  
# this is the final dataset for ages 3-5 years
weights_gems_1000_final <- weights_gems_1000_sum_final[,-9]

######################
# save both datasets #
######################

save(weights_maled_gems_1000_final,file="Scripts - Aim 3/weights_maled_gems_1000_final.Rdata")
save(weights_gems_1000_final,file="Scripts - Aim 3/weights_gems_1000_final.Rdata")

############################################################################################################
# THE NEXT DATAFRAME SHOULD FIGURE OUT HOW TO USE THESE NUMBERS TO PULL THE ROWS FROM THE SUBCLINICAL DATA #
############################################################################################################

# The _final datasets indicate how many rows should be pulled from each site in the subclinical_gems_maled_0_2.Rdata and subclinical_gems_3_5.Rdata datasets
#     weights_maled_gems_1000_sum_final indicates how many rows for each country from each region should be pulled from subclinical_gems_maled_0_2
# so that we will have 1000 rows for each region (as the draws will total 1000 from all the countries for that region)

# Function where I specify the unique num_rows_to_select for each region - 
# the code below is pulling the number of rows specified in weights_maled_gems_1000_final from subclinical_gems_maled_0_2 and then pasting them together in 1 column


select_rows_from_columns <- function(data, num_rows_to_select, region) {
  # Set the seed for reproducibility
  set.seed(42)
  
  # Initialize a list to store selected data
  selected_data <- list()
  
  # Randomly select unique rows from each column and bind them into a single column
  for (i in seq_along(num_rows_to_select)) {
    selected_data[[i]] <- sample(data[[i]], num_rows_to_select[i], replace = FALSE)
  }
  
  # Combine the selected data into a single column dataframe
  selected_data <- as.data.frame(do.call(c, selected_data))
  
  # Rename the column if needed
  colnames(selected_data) <- region
  
  # Return the selected data
  return(selected_data)
}

# num_rows_to_select is overwritten each time. not sure that's the best way to go about this, but it works if done in the proper sequence.

########
# AFRO #
########

# 0-2 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(4,1,84,0,0,85,0,0,238,50,144,0,29,280,85)  # Different number of rows for each column
# Call the function
AFRO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "AFRO_0_2")

# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(22,5,480,0,0,493,0)  # Different number of rows for each column
# Call the function
AFRO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "AFRO_3_5")

########
# AMRO #
########

# 0-2 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(0,0,0,0,0,87,24,0,671,0,28,0,0,177,13)  # Different number of rows for each column
# Call the function
AMRO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "AMRO_0_2")

# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(0,0,0,0,0,783,217)  # Different number of rows for each column
# Call the function
AMRO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "AMRO_3_5")

########
# EURO #
########

# 0-2 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(0,0,0,0,0,0,0,0,886,0,110,0,0,4,0)  # Different number of rows for each column
# Call the function
EURO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "EURO_0_2")

# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(144,83,55,218,117,182,201)  # Different number of rows for each column
# Call the function
EURO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "EURO_3_5")

#########
# SEARO #
#########

# 0-2 years #
num_rows_to_select <- c(339,0,0,0,0,262,12,0,284,6,0,0,0,97,0)  # Different number of rows for each column
# Call the function
SEARO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "SEARO_0_2")

# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(552,0,0,0,0,428,20)  # Different number of rows for each column
# Call the function
SEARO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "SEARO_3_5")

########
# WPRO #
########

# 0-2 years #
num_rows_to_select <- c(0,0,0,0,0,0,0,0,883,9,13,0,34,61,0)  # Different number of rows for each column
# Call the function
WPRO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "WPRO_0_2")
  
# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(3,1,2,1,1,991,1)  # Different number of rows for each column
# Call the function
WPRO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "WPRO_3_5")

########
# EMRO #
########

# 0-2 years #
num_rows_to_select <- c(172,0,43,0,0,174,6,0,261,28,72,0,14,188,42)  # Different number of rows for each column
# Call the function
EMRO_0_2 <- select_rows_from_columns(subclinical_gems_maled_0_2, num_rows_to_select, "EMRO_0_2")

# 3-5 years #
# Define the number of rows to randomly select for each column
num_rows_to_select <- c(435,1,106,0,0,442,16)  # Different number of rows for each column
# Call the function
EMRO_3_5 <- select_rows_from_columns(subclinical_gems_3_5, num_rows_to_select, "EMRO_3_5")

################################
# Bind the dataframes together #
################################

region_estimates_0_2_years <- cbind(AFRO_0_2, AMRO_0_2, EMRO_0_2, EURO_0_2, SEARO_0_2, WPRO_0_2)
region_estimates_3_5_years <- cbind(AFRO_3_5, AMRO_3_5, EMRO_3_5, EURO_3_5, SEARO_3_5, WPRO_3_5)

##########################################################################################
# Draws need to be randomly re-ordered as they are currently grouped together by country #
##########################################################################################

bystander_region_estimates_0_2_years_final <- region_estimates_0_2_years[sample(nrow(region_estimates_0_2_years)), ]
bystander_region_estimates_3_5_years_final <- region_estimates_3_5_years[sample(nrow(region_estimates_3_5_years)), ]

######################################
# Save final age specific dataframes # 
######################################

save(bystander_region_estimates_0_2_years_final,file="Scripts - Aim 3/Output/bystander_region_estimates_0_2_years_final.Rdata")
save(bystander_region_estimates_3_5_years_final,file="Scripts - Aim 3/Output/bystander_region_estimates_3_5_years_final.Rdata")

########################################################################################
# Combine age specific dataframes into one FINAL dataframe for all kids ages 0-5 years #
########################################################################################

# first, bind the two dataframes together
bystander_region_estimates_combined <- cbind(bystander_region_estimates_0_2_years_final,bystander_region_estimates_3_5_years_final)

# make entire dataset numeric
bystander_region_estimates_combined <- mutate_all(bystander_region_estimates_combined, as.numeric)

# 40% (0-2 years) and 60% weighting (>2-5 years)
bystander_region_estimates_0_5_years <- bystander_region_estimates_combined %>%
  mutate(AFRO_0_5 = (AFRO_0_2*0.4) + (AFRO_3_5*0.6)) %>%
  mutate(AMRO_0_5 = (AMRO_0_2*0.4) + (AMRO_3_5*0.6)) %>%
  mutate(EURO_0_5 = (EURO_0_2*0.4) + (EURO_3_5*0.6)) %>%
  mutate(SEARO_0_5 = (SEARO_0_2*0.4) + (SEARO_3_5*0.6)) %>%
  mutate(WPRO_0_5 = (WPRO_0_2*0.4) + (WPRO_3_5*0.6)) %>%
  mutate(EMRO_0_5 = (EMRO_0_2*0.4) + (EMRO_3_5*0.6))

# the FINAL bystander pathogen data for kids 0-5 years
# this is the average number of bystanders per abx course by region
bystander_region_estimates_0_5_years_final <- bystander_region_estimates_0_5_years %>%
  select(AFRO_0_5:EMRO_0_5)

########################
# Save final dataframe # 
########################

save(bystander_region_estimates_0_5_years_final,file="Scripts - Aim 3/Output/bystander_region_estimates_0_5_years_final.Rdata")
