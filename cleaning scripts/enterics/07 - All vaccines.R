# IGNORE ROTAVIRUS MATRIX

# Combining the rotavirus vaccine + all other vaccine output into 1 dataframe with 1000 bootstrapped rows, but all vaccines included
# All countries from all regions are getting the same estimates - There is not enough data to make this country/region specific

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(abind)

# both datasets belowo were originally named the same so we have to re-name
# bring in datasets
load("abx output/aim3_Combo_Vaccines_NoRota_boot_df_7Mar.Rda")
NoRota <- aim2_boot_df

load("abx output/aim3_Combo_Vaccines_Rota_boot_df_7Mar.Rda")
Rota <- aim2_boot_df

######################################
# Shigella, Campy, ETEC, Noro, Adeno #
######################################

# keep only relative differences
NoRota1 <-NoRota[,c(4,10,16,22,28)]

# Setting up  dataframe for arrays
# The  dataframe needs to be broken up by site and then transformed where the 1000 bootstraps are columns instead of rows
# renaming columns

NoRota1 <- NoRota1 %>% rename(rel_diff_Eq2_shig = rel_diff_Eq2)
NoRota1 <- NoRota1 %>% rename(rel_diff_Eq2_campy = rel_diff_Eq2.1)
NoRota1 <- NoRota1 %>% rename(rel_diff_Eq2_ETEC = rel_diff_Eq2.2)
NoRota1 <- NoRota1 %>% rename(rel_diff_Eq2_noro = rel_diff_Eq2.3)
NoRota1 <- NoRota1 %>% rename(rel_diff_Eq2_adeno = rel_diff_Eq2.4)

#############
# Rotavirus #
#############

# keep only relative differences
Rota1 <- Rota %>% 
  select(rel_diff_Eq2)

# Setting up  dataframe for arrays
# The  dataframe needs to be broken up by site and then transformed where the 1000 bootstraps are columns instead of rows
# renaming columns

Rota1 <- Rota1 %>% rename(rel_diff_Eq2_rota = rel_diff_Eq2)

##########################################################
# binding all pathogens together to make the final array #
##########################################################

all_vaccines <- cbind(NoRota1, Rota1)

# re-order columns
all_vaccines <- all_vaccines[, c(1,2,3,4,6,5)]

# convert to matrix
all_vaccines <- as.matrix(all_vaccines)

# divide up into individual vaccine matrices 
# then duplicate to 135 columns for the 135 countries
shig_vaccine_matrix <-all_vaccines[,c(1)]
shig_vaccine_matrix <- replicate(135, shig_vaccine_matrix)

campy_vaccine_matrix <-all_vaccines[,c(2)]
campy_vaccine_matrix <- replicate(135, campy_vaccine_matrix)

ETEC_vaccine_matrix <-all_vaccines[,c(3)]
ETEC_vaccine_matrix <- replicate(135, ETEC_vaccine_matrix)

noro_vaccine_matrix <-all_vaccines[,c(4)]
noro_vaccine_matrix <- replicate(135, noro_vaccine_matrix)

rota_vaccine_matrix <-all_vaccines[,c(5)]
rota_vaccine_matrix <- replicate(135, rota_vaccine_matrix)

adeno_vaccine_matrix <-all_vaccines[,c(6)]
adeno_vaccine_matrix <- replicate(135, adeno_vaccine_matrix)

# unclear if we are going to do the vaccine piece as individual matrices or if they should be in an array
# this will have 1000 rows, 135 columns, and they are grouped by vaccine
all_vaccines_array <- array(data = c(shig_vaccine_matrix, campy_vaccine_matrix, ETEC_vaccine_matrix, noro_vaccine_matrix, rota_vaccine_matrix, adeno_vaccine_matrix),
                    dim = c(1000, 135, 6))

# double checking that the above array is correct 
subset <- all_vaccines_array[1:1000, 1:135, 1]

########
# save #
########

save(shig_vaccine_matrix,file="Scripts - Aim 3/Output/shig_vaccine_matrix.Rdata")
save(campy_vaccine_matrix,file="Scripts - Aim 3/Output/campy_vaccine_matrix.Rdata")
save(ETEC_vaccine_matrix,file="Scripts - Aim 3/Output/ETEC_vaccine_matrix.Rdata")
save(noro_vaccine_matrix,file="Scripts - Aim 3/Output/noro_vaccine_matrix.Rdata")
save(rota_vaccine_matrix,file="Scripts - Aim 3/Output/rota_vaccine_matrix.Rdata")
save(adeno_vaccine_matrix,file="Scripts - Aim 3/Output/adeno_vaccine_matrix.Rdata")


save(all_vaccines_array,file="Scripts - Aim 3/Output/all_vaccines_array.Rdata")
