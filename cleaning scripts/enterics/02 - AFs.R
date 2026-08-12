# this script cleans up the GPDS AF data sent by James 
# AF.region.error.weighted.draws.20172018 – this contains a list of dataframes that have 1000 draws for each pathogen from GPDS 2017-2018. 
#        It is a list, with one item for each region. 
#        In each item, there is a data frame that has 16 rows (one for each pathogen) and 1000 columns (one for each draw)

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(janitor)

load("Scripts - Aim 3/James data/GPDS AF region draws for Wellcome AMR.Rdata")
load("Scripts - Aim 3/country_region_list.Rdata")

# making the dataframe into a "flat" matrix
AF_matrix <- rbind(AF.region.error.weighted.draws.20172018$AFRO,AF.region.error.weighted.draws.20172018$AMRO,AF.region.error.weighted.draws.20172018$EURO,AF.region.error.weighted.draws.20172018$SEARO,AF.region.error.weighted.draws.20172018$WPRO)
# 1000 columns and 80 rows (16 pathogens x 5 regions)
str(AF_matrix)

# Define row names to be removed
rows_to_remove <- c("aeromonas", "astrovirus", "cryptosporidium", "cyclospora", "e_histolytica", "isospora", "salmonella", "sapovirus", "v_cholerae", "TEPEC")

# Remove rows based on row names
# now there are 30 rows (6 pathogens x 5 regions)
AF_matrix <- AF_matrix[!(rownames(AF_matrix) %in% rows_to_remove), ]
str(AF_matrix)

# Get the number of columns (to rename column names to match the other matrices)
num_cols <- ncol(AF_matrix)

# Generate new column names starting with "1" and continuing with consecutive numbers
new_col_names <- paste0(seq_len(num_cols))

# Assign new column names to the data frame
colnames(AF_matrix) <- new_col_names

# separate matrix out into 6 matrices (one for each pathogen)
# this will create 6 matrices  with 1000 columns and 5 rows - where each row is a WHO region in this order (AFRO, AMRO, EURO, SEARO, WPRO)
AF_adeno <- AF_matrix[c(1,7,13,19,25),]
AF_campy <- AF_matrix[c(2,8,14,20,26),]
AF_noro <- AF_matrix[c(3,9,15,21,27),]
AF_rota <- AF_matrix[c(4,10,16,22,28),]
AF_shigella <- AF_matrix[c(5,11,17,23,29),]
AF_ETEC <- AF_matrix[c(6,12,18,24,30),]

# rename rows
rownames(AF_adeno) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")
rownames(AF_campy) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")
rownames(AF_noro) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")
rownames(AF_rota) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")
rownames(AF_shigella) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")
rownames(AF_ETEC) <- c("AFR", "AMR", "EUR", "SEAR", "WPR")

# now that we have the region estimates in the format we need - we now need each country to get the data for the region in which they are part of
# we will make this a function, but start with one pathogen first 

country_path <- function(AF_path) {
  
AF_path1 <- AF_path
WHO_REGION_CODE <- rownames(AF_path1)
rownames(AF_path1) <- NULL
AF_path2 <- cbind(WHO_REGION_CODE,AF_path1)
AF_path2 <- as.data.frame(AF_path2)

# slim down country list 
country <- country_region_list %>%
  select(countriesSub, WHO_REGION_CODE)

# merge
AF_country_estimates <- left_join(country, AF_path2, by="WHO_REGION_CODE") # left_join prevents it from sorting
AF_country_estimates <- AF_country_estimates[,c(2,1,3:1002)]

AF_country_estimates_final <- t(AF_country_estimates)

# remove first 2 rows
AF_country_estimates_final <- AF_country_estimates_final[-c(1, 2), ]

# Return the selected data
return(AF_country_estimates_final)
}

AF_country_estimates_adeno_matrix <- country_path(AF_adeno)
AF_country_estimates_campy_matrix <- country_path(AF_campy)
AF_country_estimates_noro_matrix <- country_path(AF_noro)
AF_country_estimates_rota_matrix <- country_path(AF_rota)
AF_country_estimates_shigella_matrix <- country_path(AF_shigella)
AF_country_estimates_ETEC_matrix <- country_path(AF_ETEC)

#####################
# ADDED 3 JUNE 2024 #
#####################

# Fill in data for the 14 EMRO countries: #66 Afghanistan (use V93), #52 Djibouti (use V64), #68 Egypt (use V64), #69 Iran (use V93), #70 Iraq (use V93), 
                                          #71 Jordan (use V93), #72 Lebanon (use V93), #73 Libya (use V64), #74 Morocco (use V64), #95 Pakistan (use V93), 
                                          #60 Somalia (use V64), #62 Sudan (use V64), #77 Tunisia (use V64), #79 Yemen (use V93)
# For those now assigned to SEARO use column: V93 India
# For those now assigned to AFRO use column: V64 Tanzania

# make matrices a df so that mutate works
AF_country_estimates_adeno_matrix <- as.data.frame(AF_country_estimates_adeno_matrix)
AF_country_estimates_campy_matrix <- as.data.frame(AF_country_estimates_campy_matrix)
AF_country_estimates_noro_matrix <- as.data.frame(AF_country_estimates_noro_matrix)
AF_country_estimates_rota_matrix <- as.data.frame(AF_country_estimates_rota_matrix)
AF_country_estimates_shigella_matrix <- as.data.frame(AF_country_estimates_shigella_matrix)
AF_country_estimates_ETEC_matrix <- as.data.frame(AF_country_estimates_ETEC_matrix)

fill <- function(df) {
# create function to loop through the below for all dfs
df <- df %>%
  mutate(V66 = V93) %>% #66 Afghanistan (use V93)
  mutate(V52 = V64) %>% #52 Djibouti (use V64)
  mutate(V68 = V64) %>% #68 Egypt (use V64)
  mutate(V69 = V93) %>% #69 Iran (use V93)
  mutate(V70 = V93) %>% #70 Iraq (use V93)
  mutate(V71 = V93) %>% #71 Jordan (use V93)
  mutate(V72 = V93) %>% #72 Lebanon (use V93)
  mutate(V73 = V64) %>% #73 Libya (use V64)
  mutate(V74 = V64) %>% #74 Morocco (use V64)
  mutate(V95 = V93) %>% #95 Pakistan (use V93)
  mutate(V60 = V64) %>% #60 Somalia (use V64)
  mutate(V62 = V64) %>% #62 Sudan (use V64)
  mutate(V77 = V64) %>% #77 Tunisia (use V64)
  mutate(V79 = V93) #79 Yemen (use V93)

return(df)
}

# call the function and back to a matrix
AF_country_estimates_adeno_matrix <- as.matrix(fill(AF_country_estimates_adeno_matrix))
AF_country_estimates_campy_matrix <- as.matrix(fill(AF_country_estimates_campy_matrix))
AF_country_estimates_noro_matrix <- as.matrix(fill(AF_country_estimates_noro_matrix))
AF_country_estimates_rota_matrix <- as.matrix(fill(AF_country_estimates_rota_matrix))
AF_country_estimates_shigella_matrix <- as.matrix(fill(AF_country_estimates_shigella_matrix))
AF_country_estimates_ETEC_matrix <- as.matrix(fill(AF_country_estimates_ETEC_matrix))

# end of addition
#####################
#####################

# this will have 1000 rows, 135 columns, and they are grouped by vaccine
# pathogen order is: shigella, campy, etec, noro, rota, adeno
AF_country_estimates_array <- array(data = c(AF_country_estimates_shigella_matrix, AF_country_estimates_campy_matrix, AF_country_estimates_ETEC_matrix, 
                                             AF_country_estimates_noro_matrix, AF_country_estimates_rota_matrix, AF_country_estimates_adeno_matrix),
                            dim = c(1000, 135, 6))

# double checking that the above array is correct 
subset <- AF_country_estimates_array[1:1000, 1:135, 1]

########
# save #
########

save(AF_country_estimates_shigella_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_shigella_matrix.Rdata")
save(AF_country_estimates_campy_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_campy_matrix.Rdata")
save(AF_country_estimates_ETEC_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_ETEC_matrix.Rdata")
save(AF_country_estimates_noro_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_noro_matrix.Rdata")
save(AF_country_estimates_rota_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_rota_matrix.Rdata")
save(AF_country_estimates_adeno_matrix,file="Scripts - Aim 3/Output/AF_country_estimates_adeno_matrix.Rdata")


save(AF_country_estimates_array,file="Scripts - Aim 3/Output/AF_country_estimates_array.Rdata")













