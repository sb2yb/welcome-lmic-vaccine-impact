# this script will go through diarrhea pathogens and all equations for Aim 3
# Producing region level, not country level estimates 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)

# bring in matrices 
# incidence
load("Scripts - Aim 3/Output/DHS_diar_matrix.Rdata")

# AF - these are character at the moment because of we have NAs
load("Scripts - Aim 3/Output/AF_country_estimates_shigella_matrix.Rdata")
load("Scripts - Aim 3/Output/AF_country_estimates_campy_matrix.Rdata")
load("Scripts - Aim 3/Output/AF_country_estimates_ETEC_matrix.Rdata")
load("Scripts - Aim 3/Output/AF_country_estimates_noro_matrix.Rdata")
load("Scripts - Aim 3/Output/AF_country_estimates_adeno_matrix.Rdata")

# vaccine efficacy
load("Scripts - Aim 3/Output/shig_vaccine_matrix.Rdata")
load("Scripts - Aim 3/Output/campy_vaccine_matrix.Rdata")
load("Scripts - Aim 3/Output/ETEC_vaccine_matrix.Rdata")
load("Scripts - Aim 3/Output/noro_vaccine_matrix.Rdata")
load("Scripts - Aim 3/Output/adeno_vaccine_matrix.Rdata")

# rotavirus vaccine efficacy - DONE SEPARATELY FROM THE REST AS THERE IS REAL WORLD DATA 
load("Scripts - Aim 3/Joe data/veRota.Rdata")

# vaccine coverage - updated 
coverage <- read.csv("Scripts - Aim 3/vaccine_coverage.csv")

# bystander pathogens
load("Scripts - Aim 3/Output/bystander_country_estimates_final_matrix.Rdata")

# load country region list
load("Scripts - Aim 3/country_region_list.Rdata")

# bring in curU5 for population level estimates
pop_estimates<-read.csv(file="Scripts - Aim 3/Joe data/curU5.csv", h=T)

##################
# Clean the data #
##################

# FOR ALL LOADED SCRIPTS: REMOVE THE FOLLOWING #

# 1) Remove: #45 Democratic Peoples Republic of Korea, #84 Palau, and #108 Argentina as they go back and forth on being an LMIC (data not in DHS at all)
# 2) Remove: #75 Occupied Palestinian territory, including east Jerusalem as we don't have any data to extrapolate to this
# 3) Remove: #76 Syria (due to conflict)
# 4) Remove: #18 Georgia (the one located in Position #17 is the correct one as confirmed by Joe)

# make dataframes to make it easy to remove columns (this can be done as a matrix but then the column names get re-numbered and for tracking reasons i'd rather that not happen)
DHS_diar <- as.data.frame(DHS_diar)

AF_country_estimates_shigella_matrix <- as.data.frame(AF_country_estimates_shigella_matrix)
AF_country_estimates_campy_matrix <- as.data.frame(AF_country_estimates_campy_matrix)
AF_country_estimates_ETEC_matrix <- as.data.frame(AF_country_estimates_ETEC_matrix)
AF_country_estimates_noro_matrix <- as.data.frame(AF_country_estimates_noro_matrix)
AF_country_estimates_adeno_matrix <- as.data.frame(AF_country_estimates_adeno_matrix)

shig_vaccine_matrix <- as.data.frame(shig_vaccine_matrix)
campy_vaccine_matrix <- as.data.frame(campy_vaccine_matrix)
ETEC_vaccine_matrix <- as.data.frame(ETEC_vaccine_matrix)
noro_vaccine_matrix <- as.data.frame(noro_vaccine_matrix)
adeno_vaccine_matrix <- as.data.frame(adeno_vaccine_matrix)

bystander_country_estimates_final <- as.data.frame(bystander_country_estimates_final)

# remove countries not needed
DHS_diar <- DHS_diar %>% select(-c(V18, V45, V84, V108, V75, V76))

AF_country_estimates_shigella_matrix <- AF_country_estimates_shigella_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
AF_country_estimates_campy_matrix <- AF_country_estimates_campy_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
AF_country_estimates_ETEC_matrix <- AF_country_estimates_ETEC_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
AF_country_estimates_noro_matrix <- AF_country_estimates_noro_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
AF_country_estimates_adeno_matrix <- AF_country_estimates_adeno_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))

shig_vaccine_matrix <- shig_vaccine_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
campy_vaccine_matrix <- campy_vaccine_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
ETEC_vaccine_matrix <- ETEC_vaccine_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
noro_vaccine_matrix <- noro_vaccine_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))
adeno_vaccine_matrix <- adeno_vaccine_matrix %>% select(-c(V18, V45, V84, V108, V75, V76))

bystander_country_estimates_final <- bystander_country_estimates_final %>% select(-c(V18, V45, V84, V108, V75, V76))

#########################################
# Creating vaccine matrix for rotavirus #
#########################################

# the vaccine matrices only have 1 column. we need to make 129 identical columns for all countries. 
veRota <- replicate(129, veRota)

# the vaccine matrices have 10,000 rows. we only need 1000
veRota <- veRota[1:1000,]

rota_vaccine_matrix <- as.data.frame(veRota)

###################
# ADDING COVERAGE #
###################

# creating individual vaccine coverage dataframes 
rota_coverage <- coverage %>% select(ROTAC)
shigella_coverage <- coverage %>% select(avg_vaxcov) # using average as shigella doesn't have coverage data
campy_coverage <- coverage %>% select(avg_vaxcov) # using average as campy doesn't have coverage data
ETEC_coverage <- coverage %>% select(avg_vaxcov) # using average as ETEC doesn't have coverage data
noro_coverage <- coverage %>% select(avg_vaxcov) # using average as noro doesn't have coverage data
adeno_coverage <- coverage %>% select(avg_vaxcov) # using average as adeno doesn't have coverage data

# coverage has to be transformed so that it has 129 columns and 1000 rows. 
rotacov_wide <- as.data.frame(t(rota_coverage$ROTAC))
shigellacov_wide <- as.data.frame(t(shigella_coverage$avg_vaxcov)) #using average as shigella doesn't have coverage data
campycov_wide <- as.data.frame(t(campy_coverage$avg_vaxcov)) #using average as campy doesn't have coverage data
ETECcov_wide <- as.data.frame(t(ETEC_coverage$avg_vaxcov)) #using average as ETEC doesn't have coverage data
norocov_wide <- as.data.frame(t(noro_coverage$avg_vaxcov)) #using average as noro doesn't have coverage data
adenocov_wide <- as.data.frame(t(adeno_coverage$avg_vaxcov)) #using average as adeno doesn't have coverage data

# the _wide datasets need to have 1000 rows
rotacov_wide <- rotacov_wide[rep(1, 1000), ]
shigellacov_wide <- shigellacov_wide[rep(1, 1000), ]
campycov_wide <- campycov_wide[rep(1, 1000), ]
ETECcov_wide <- ETECcov_wide[rep(1, 1000), ]
norocov_wide <- norocov_wide[rep(1, 1000), ]
adenocov_wide <- adenocov_wide[rep(1, 1000), ]


#############
# Equations #
#############

#####################################################################
# Equation 1: Incidence of pathogen specific treatment for diarrhea #
#####################################################################
# Incidence of antibiotic use for diarrhea*AF of antibiotic treated diarrhea

# as of now, the AFs are not numeric (due to missing EMRO so we have to make numeric here - won't be needed once we have all the data)
AF_country_estimates_shigella_matrix <- AF_country_estimates_shigella_matrix %>% mutate_all(as.numeric)
AF_country_estimates_campy_matrix <- AF_country_estimates_campy_matrix %>% mutate_all(as.numeric)
AF_country_estimates_ETEC_matrix <- AF_country_estimates_ETEC_matrix %>% mutate_all(as.numeric)
AF_country_estimates_noro_matrix <- AF_country_estimates_noro_matrix %>% mutate_all(as.numeric)
AF_country_estimates_adeno_matrix <- AF_country_estimates_adeno_matrix %>% mutate_all(as.numeric)

Eq1_shigella <- DHS_diar * AF_country_estimates_shigella_matrix
Eq1_campy <- DHS_diar * AF_country_estimates_campy_matrix
Eq1_ETEC <- DHS_diar * AF_country_estimates_ETEC_matrix
Eq1_noro <- DHS_diar * AF_country_estimates_noro_matrix
Eq1_rota <- DHS_diar 
Eq1_adeno <- DHS_diar * AF_country_estimates_adeno_matrix

####################################################################################
# Equation 2: Incidence of pathogen specific treatment for diarrhea with a vaccine #
####################################################################################
# Equation 1 * Vaccine efficacy

# vaccine efficacies need to be numeric
shig_vaccine_matrix <- shig_vaccine_matrix %>% mutate_all(as.numeric)
campy_vaccine_matrix <- campy_vaccine_matrix %>% mutate_all(as.numeric)
ETEC_vaccine_matrix <- ETEC_vaccine_matrix %>% mutate_all(as.numeric)
noro_vaccine_matrix <- noro_vaccine_matrix %>% mutate_all(as.numeric)
rota_vaccine_matrix <- rota_vaccine_matrix %>% mutate_all(as.numeric)
adeno_vaccine_matrix <- adeno_vaccine_matrix %>% mutate_all(as.numeric)

Eq2_shigella <- Eq1_shigella * (1- shig_vaccine_matrix * shigellacov_wide)
Eq2_campy <- Eq1_campy * (1- campy_vaccine_matrix * campycov_wide)
Eq2_ETEC <- Eq1_ETEC * (1- ETEC_vaccine_matrix * ETECcov_wide)
Eq2_noro <- Eq1_noro * (1- noro_vaccine_matrix * norocov_wide)
Eq2_rota <- Eq1_rota * (1- rota_vaccine_matrix * rotacov_wide)
Eq2_adeno <- Eq1_adeno * (1- adeno_vaccine_matrix * adenocov_wide)

#######################################################################
# Equation 3: The absolute reduction in antibiotic use with a vaccine #
#######################################################################
# Equation 2 - Equation 1
# multiplying by 100 to make this 100 child years (DHS is per 1 child year)
Eq3_shigella <- (Eq2_shigella - Eq1_shigella) *100
Eq3_campy <- (Eq2_campy - Eq1_campy) *100
Eq3_ETEC <- (Eq2_ETEC - Eq1_ETEC) *100
Eq3_noro <- (Eq2_noro - Eq1_noro) *100
Eq3_rota <- (Eq2_rota - Eq1_rota) *100
Eq3_adeno <- (Eq2_adeno - Eq1_adeno) *100

######################################################################
# Equation 4: The percent reduction in antibiotic use with a vaccine #
######################################################################
# (1 - (Equation 2 / Equation 1)) * 100

Eq4_shigella <- (1-(Eq2_shigella / Eq1_shigella))*100
Eq4_campy <- (1-(Eq2_campy / Eq1_campy))*100
Eq4_ETEC <- (1-(Eq2_ETEC / Eq1_ETEC))*100
Eq4_noro <- (1-(Eq2_noro / Eq1_noro))*100
Eq4_rota <- (1-(Eq2_rota / Eq1_rota))*100
Eq4_adeno <- (1-(Eq2_adeno / Eq1_adeno))*100

###################################################################
# Equation 5: The incidence of bystander exposures to antibiotics #
###################################################################
# Equation 1 * average number of asymptomatic pathogens

# bystander pathogens need to be numeric
bystander_country_estimates_final <- bystander_country_estimates_final %>% mutate_all(as.numeric)

Eq5_shigella <- Eq1_shigella * bystander_country_estimates_final
Eq5_campy <- Eq1_campy * bystander_country_estimates_final
Eq5_ETEC <- Eq1_ETEC * bystander_country_estimates_final
Eq5_noro <- Eq1_noro * bystander_country_estimates_final
Eq5_rota <- Eq1_rota * bystander_country_estimates_final
Eq5_adeno <- Eq1_adeno * bystander_country_estimates_final

#############################################################################
#Equation 6: Incidence of bystander exposures to antibiotics with a vaccine #
#############################################################################
#Equation 2 * average number of asymptomatic pathogens

Eq6_shigella <- Eq2_shigella * bystander_country_estimates_final
Eq6_campy <- Eq2_campy * bystander_country_estimates_final
Eq6_ETEC <- Eq2_ETEC * bystander_country_estimates_final
Eq6_noro <- Eq2_noro * bystander_country_estimates_final
Eq6_rota <- Eq2_rota * bystander_country_estimates_final
Eq6_adeno <- Eq2_adeno * bystander_country_estimates_final

###########################################################################################
# Equation 7: The absolute reduction of bystander exposures to antibiotics with a vaccine #
###########################################################################################
# Equation 6 - Equation 5
# multiplying by 100 to make this 100 child years (DHS is per 1 child year)
Eq7_shigella <- (Eq6_shigella - Eq5_shigella) *100
Eq7_campy <- (Eq6_campy - Eq5_campy) *100
Eq7_ETEC <- (Eq6_ETEC - Eq5_ETEC) *100
Eq7_noro <- (Eq6_noro - Eq5_noro) *100
Eq7_rota <- (Eq6_rota - Eq5_rota) *100
Eq7_adeno <- (Eq6_adeno - Eq5_adeno) *100

###########################################################################################
# Equation 8: The absolute reduction of bystander exposures to antibiotics with a vaccine #
###########################################################################################
# (1-(Equation 6 / Equation 5))*100

Eq8_shigella <- (1-(Eq6_shigella / Eq5_shigella))*100
Eq8_campy <- (1-(Eq6_campy / Eq5_campy))*100
Eq8_ETEC <- (1-(Eq6_ETEC / Eq5_ETEC))*100
Eq8_noro <- (1-(Eq6_noro / Eq5_noro))*100
Eq8_rota <- (1-(Eq6_rota / Eq5_rota))*100
Eq8_adeno <- (1-(Eq6_adeno / Eq5_adeno))*100

## bind dataframes together

bind <- cbind(Eq1_shigella, Eq1_campy, Eq1_ETEC, Eq1_noro, Eq1_rota, Eq1_adeno, 
          Eq2_shigella, Eq2_campy, Eq2_ETEC, Eq2_noro, Eq2_rota, Eq2_adeno, 
          Eq3_shigella, Eq3_campy, Eq3_ETEC, Eq3_noro, Eq3_rota, Eq3_adeno, 
          Eq4_shigella, Eq4_campy, Eq4_ETEC, Eq4_noro, Eq4_rota, Eq4_adeno, 
          Eq5_shigella, Eq5_campy, Eq5_ETEC, Eq5_noro, Eq5_rota, Eq5_adeno, 
          Eq6_shigella, Eq6_campy, Eq6_ETEC, Eq6_noro, Eq6_rota, Eq6_adeno,
          Eq7_shigella, Eq7_campy, Eq7_ETEC, Eq7_noro, Eq7_rota, Eq7_adeno,
          Eq8_shigella, Eq8_campy, Eq8_ETEC, Eq8_noro, Eq8_rota, Eq8_adeno)

################################################
# Make one dataframe per pathogen per equation #
################################################

# Having R make a new dataframe every 129 columns for the 6 pathogens and 8 equations = 48 dataframes (keeps me from having to manually do this)
total_columns <- 6192 # total number of columns in bind 
chunk_size <- 129 # total number of countries 
num_chunks <- total_columns / chunk_size

# Create a list to store the smaller dataframes
dfs <- list()

# Split the dataframe into chunks
for (i in 1:num_chunks) {
  start_col <- (i - 1) * chunk_size + 1
  end_col <- i * chunk_size
  chunk_df <- bind[, start_col:end_col]
  dfs[[paste0("df_chunk_", i)]] <- chunk_df
}

# Optionally assign the smaller dataframes to variables in the global environment
# This will create 36 smaller dataframes named df_chunk_1 to df_chunk_36
for (name in names(dfs)) {
  assign(name, dfs[[name]], envir = .GlobalEnv)
}


####################
# ADD WEIGHTS HERE #
####################

# use WHO region codes in country_region_list

# remove last 3 columns
country_region_list <- country_region_list %>% select(countriesSub, DisplayString, WHO_REGION_CODE, WHO_REGION)
country_region_list <- country_region_list[-c(18, 45, 84, 108, 75, 76), ] # remove the countries we don't want to keep 


## below code is trying to find a way to specifically tell R which columns belong to which regions for proper addition across rows ##
# it needs to be done from country_region_list because the columns are in the format that matches the chunks where they kept their original row number going to 135
# create column that adds "V" to the row number 
country_region_list$V <- paste0("V", rownames(country_region_list)) # IMPORTANT that it's the row number and simply counting from 1 (needs to match original row #)

# collapse all Vs for each region into a single string
country_region_table <- country_region_list %>%
  group_by(WHO_REGION_CODE) %>%
  summarise(Vs = paste(V, collapse = ", "), .groups = "drop")

# print without truncation
options(width = 10000)  # set a very wide console width
print(country_region_table, row.names = FALSE)

## end code ## 


# keep columns from pop_estimates that we want 
pop_estimates <- pop_estimates %>% select(Country.Code, curU5)

# add population estimates to country region list 
country_pop_est <- left_join(country_region_list, pop_estimates, by=c("countriesSub" = "Country.Code"))

# keep only curU5 and transform so that it's 129 columns and 1000 rows (of identical data)
# we will then use this and multiply it against all the chunks 
country_pop_est_wide <- country_pop_est %>%
  select(curU5) %>%
  t() %>%        
  as.data.frame() 

# duplicate row 1000 times
country_pop_est_wide <- country_pop_est_wide[rep(1, 1000), ]
rownames(country_pop_est_wide) <- NULL # unnaming the rows

# define which columns belong to which regions 
AFR_cols <- c("V38", "V39", "V40", "V41", "V42", "V43", "V50", "V51", "V53", "V54",
              "V55", "V56", "V57", "V58", "V59", "V61", "V63", "V64", "V65", "V67",
              "V101","V109","V110","V111","V112","V113","V114","V117","V118","V119",
              "V120","V121","V122","V123","V124","V125","V126","V127","V128","V129",
              "V130","V131","V132","V133","V134","V135")

AMR_cols <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13",
              "V14","V31","V32","V33","V34","V35","V36","V37","V115","V116")

EMR_cols <- c("V52","V60","V62","V66","V68","V69","V70","V71","V72","V73","V74","V77",
              "V79","V95")

EUR_cols <- c("V15","V16","V17","V19","V20","V22","V23","V24","V25","V26","V27","V28",
              "V29","V30","V46","V47","V48","V49","V78")

SEAR_cols <- c("V91","V92","V93","V94","V97","V100","V102","V104","V105","V106")

WPR_cols <- c("V21","V44","V80","V81","V82","V83","V85","V86","V87","V88","V89","V90",
              "V96","V98","V99","V103","V107")

# calculate the total population in each region 
region_pop_df <- country_pop_est %>%
  group_by(WHO_REGION_CODE) %>%
  summarise(region_pop = sum(curU5, na.rm = TRUE), .groups = "drop")

# define population denominator by region
region_pop <- c(
  AFR = 176555118,
  AMR = 42319195,
  EMR = 82833169,
  EUR = 28516443,
  SEAR = 166759589,
  WPR = 106391504)


# eventually make this a function, but we will try it out with 1 datachunk to see if we can get this up and running

# step 1:
# multiply the dataframe with vaccine estimates against the population estimates for each country 
# df_chunk_1_mult <- df_chunk_1*country_pop_est_wide

# AFR: V38, V39, V40, V41, V42, V43, V50, V51, V53, V54, V55, V56, V57, V58, V59, V61, V63, V64, V65, V67, V101, V109, V110, V111, V112, V113, V114, V117, V118, V119, V120, V121, V122, V123, V124, V125, V126, V127, V128, V129, V130, V131, V132, V133, V134, V135
# AMR: V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V31, V32, V33, V34, V35, V36, V37, V115, V116  
# EMR: V52, V60, V62, V66, V68, V69, V70, V71, V72, V73, V74, V77, V79, V95
# EUR: V15, V16, V17, V19, V20, V22, V23, V24, V25, V26, V27, V28, V29, V30, V46, V47, V48, V49, V78
# SEAR: V91, V92, V93, V94, V97, V100, V102, V104, V105, V106
# WPR: V21, V44, V80, V81, V82, V83, V85, V86, V87, V88, V89, V90, V96, V98, V99, V103, V107 

# step 2:
# sum the region specific rows and create a column to hold those sums. then divide by the region population denominator to weight the estimates
# df_chunk_1_mult_weighted <- df_chunk_1_mult %>%
#  mutate(
#    AFR = rowSums(across(all_of(AFR_cols)), na.rm = TRUE) / region_pop["AFR"],
#    AMR = rowSums(across(all_of(AMR_cols)), na.rm = TRUE) / region_pop["AMR"],
#    EMR = rowSums(across(all_of(EMR_cols)), na.rm = TRUE) / region_pop["EMR"],
#    EUR = rowSums(across(all_of(EUR_cols)), na.rm = TRUE) / region_pop["EUR"],
#    SEAR = rowSums(across(all_of(SEAR_cols)), na.rm = TRUE) / region_pop["SEAR"],
#    WPR = rowSums(across(all_of(WPR_cols)), na.rm = TRUE) / region_pop["WPR"])

# making sure that the above worked
#EMR_check <- df_chunk_1_mult_add %>%
#  select(V52, V60, V62, V66, V68, V69, V70, V71, V72, V73, V74, V77, V79, V95, EMR)

## FUNCTION TO WEIGHT EACH OF THE 48 CHUNKS ## 
# (keep above as a way of "going inside" the function)

weight_by_region <- function(df_chunks) {
  
  # Step 1: multiply by population estimates
  df_mult <- df_chunks * country_pop_est_wide
  
  # Step 2: sum by region and divide by region population
  df_chunks_weighted <- df_mult %>%
    mutate(
      AFR = rowSums(across(all_of(AFR_cols)), na.rm = TRUE) / region_pop["AFR"],
      AMR = rowSums(across(all_of(AMR_cols)), na.rm = TRUE) / region_pop["AMR"],
      EMR = rowSums(across(all_of(EMR_cols)), na.rm = TRUE) / region_pop["EMR"],
      EUR = rowSums(across(all_of(EUR_cols)), na.rm = TRUE) / region_pop["EUR"],
      SEAR = rowSums(across(all_of(SEAR_cols)), na.rm = TRUE) / region_pop["SEAR"],
      WPR = rowSums(across(all_of(WPR_cols)), na.rm = TRUE) / region_pop["WPR"]
    )
  
  return(df_chunks_weighted)
}

# calling the function for all 48 chunks
df_chunk_1_mult_weight <- weight_by_region(df_chunk_1)
df_chunk_2_mult_weight <- weight_by_region(df_chunk_2)
df_chunk_3_mult_weight <- weight_by_region(df_chunk_3)
df_chunk_4_mult_weight <- weight_by_region(df_chunk_4)
df_chunk_5_mult_weight <- weight_by_region(df_chunk_5)
df_chunk_6_mult_weight <- weight_by_region(df_chunk_6)
df_chunk_7_mult_weight <- weight_by_region(df_chunk_7)
df_chunk_8_mult_weight <- weight_by_region(df_chunk_8)
df_chunk_9_mult_weight <- weight_by_region(df_chunk_9)
df_chunk_10_mult_weight <- weight_by_region(df_chunk_10)
df_chunk_11_mult_weight <- weight_by_region(df_chunk_11)
df_chunk_12_mult_weight <- weight_by_region(df_chunk_12)
df_chunk_13_mult_weight <- weight_by_region(df_chunk_13)
df_chunk_14_mult_weight <- weight_by_region(df_chunk_14)
df_chunk_15_mult_weight <- weight_by_region(df_chunk_15)
df_chunk_16_mult_weight <- weight_by_region(df_chunk_16)
df_chunk_17_mult_weight <- weight_by_region(df_chunk_17)
df_chunk_18_mult_weight <- weight_by_region(df_chunk_18)
df_chunk_19_mult_weight <- weight_by_region(df_chunk_19)
df_chunk_20_mult_weight <- weight_by_region(df_chunk_20)
df_chunk_21_mult_weight <- weight_by_region(df_chunk_21)
df_chunk_22_mult_weight <- weight_by_region(df_chunk_22)
df_chunk_23_mult_weight <- weight_by_region(df_chunk_23)
df_chunk_24_mult_weight <- weight_by_region(df_chunk_24)
df_chunk_25_mult_weight <- weight_by_region(df_chunk_25)
df_chunk_26_mult_weight <- weight_by_region(df_chunk_26)
df_chunk_27_mult_weight <- weight_by_region(df_chunk_27)
df_chunk_28_mult_weight <- weight_by_region(df_chunk_28)
df_chunk_29_mult_weight <- weight_by_region(df_chunk_29)
df_chunk_30_mult_weight <- weight_by_region(df_chunk_30)
df_chunk_31_mult_weight <- weight_by_region(df_chunk_31)
df_chunk_32_mult_weight <- weight_by_region(df_chunk_32)
df_chunk_33_mult_weight <- weight_by_region(df_chunk_33)
df_chunk_34_mult_weight <- weight_by_region(df_chunk_34)
df_chunk_35_mult_weight <- weight_by_region(df_chunk_35)
df_chunk_36_mult_weight <- weight_by_region(df_chunk_36)
df_chunk_37_mult_weight <- weight_by_region(df_chunk_37)
df_chunk_38_mult_weight <- weight_by_region(df_chunk_38)
df_chunk_39_mult_weight <- weight_by_region(df_chunk_39)
df_chunk_40_mult_weight <- weight_by_region(df_chunk_40)
df_chunk_41_mult_weight <- weight_by_region(df_chunk_41)
df_chunk_42_mult_weight <- weight_by_region(df_chunk_42)
df_chunk_43_mult_weight <- weight_by_region(df_chunk_43)
df_chunk_44_mult_weight <- weight_by_region(df_chunk_44)
df_chunk_45_mult_weight <- weight_by_region(df_chunk_45)
df_chunk_46_mult_weight <- weight_by_region(df_chunk_46)
df_chunk_47_mult_weight <- weight_by_region(df_chunk_47)
df_chunk_48_mult_weight <- weight_by_region(df_chunk_48)

# we only want to keep the last 6 columns of each chunk (because those are where our weighted estimates are)
# Loop over all 48 weighted dataframes
for(i in 1:48) {
  
  # Construct the name of the dataframe as a string.
  # For example, when i = 1, df_name becomes "df_chunk_1_mult_weight"
  df_name <- paste0("df_chunk_", i, "_mult_weight")
  
  # Use get() to retrieve the dataframe object using the string name
  # Then select only the last 6 columns with select(tail(names(.), 6))
  # Finally, assign() writes the result back to the same dataframe name
  assign(df_name, get(df_name) %>% select(tail(names(.), 6)))
  
  # Summary:
  # - get(df_name) fetches the existing dataframe
  # - tail(names(.), 6) grabs the last 6 column names
  # - select(...) keeps only those columns
  # - assign(...) replaces the old dataframe with the trimmed one
}


################################
# Getting the median + 95% CIs #
################################

# Creating a function to loop through the below with all 36 dataframes
CI_med_function <- function(chunk_df, lower, upper, median) {

# grabbing the 2.5% and 97.5%
chunk_CI <- as.data.frame(t(apply(chunk_df, 2, quantile, c(0.025, 0.975), na.rm=TRUE)))
names(chunk_CI)[1] <- lower
names(chunk_CI)[2] <- upper

# grabbing the median
chunk_med <- as.data.frame(apply(chunk_df, 2, quantile, 0.5, na.rm=TRUE))
names(chunk_med)[1] <- median

# bind together
chunk_med_CIs <- cbind(chunk_med, chunk_CI)

# Return the selected data
return(chunk_med_CIs)

}


Eq1_shigella_med_CI <- CI_med_function(df_chunk_1_mult_weight, "Eq1_shigella_lower", "Eq1_shigella_upper", "Eq1_shigella_median")
Eq1_campy_med_CI <- CI_med_function(df_chunk_2_mult_weight, "Eq1_campy_lower", "Eq1_campy_upper", "Eq1_campy_median")
Eq1_ETEC_med_CI <- CI_med_function(df_chunk_3_mult_weight, "Eq1_ETEC_lower", "Eq1_ETEC_upper", "Eq1_ETEC_median")
Eq1_noro_med_CI <- CI_med_function(df_chunk_4_mult_weight, "Eq1_noro_lower", "Eq1_noro_upper", "Eq1_noro_median")
Eq1_rota_med_CI<- CI_med_function(df_chunk_5_mult_weight, "Eq1_rota_lower", "Eq1_rota_upper", "Eq1_rota_median")
Eq1_adeno_med_CI <- CI_med_function(df_chunk_6_mult_weight, "Eq1_adeno_lower", "Eq1_adeno_upper", "Eq1_adeno_median")

Eq2_shigella_med_CI <- CI_med_function(df_chunk_7_mult_weight, "Eq2_shigella_lower", "Eq2_shigella_upper", "Eq2_shigella_median")
Eq2_campy_med_CI <- CI_med_function(df_chunk_8_mult_weight, "Eq2_campy_lower", "Eq2_campy_upper", "Eq2_campy_median")
Eq2_ETEC_med_CI <- CI_med_function(df_chunk_9_mult_weight, "Eq2_ETEC_lower", "Eq2_ETEC_upper", "Eq2_ETEC_median")
Eq2_noro_med_CI <- CI_med_function(df_chunk_10_mult_weight, "Eq2_noro_lower", "Eq2_noro_upper", "Eq2_noro_median")
Eq2_rota_med_CI <- CI_med_function(df_chunk_11_mult_weight, "Eq2_rota_lower", "Eq2_rota_upper", "Eq2_rota_median")
Eq2_adeno_med_CI<- CI_med_function(df_chunk_12_mult_weight, "Eq2_adeno_lower", "Eq2_adeno_upper", "Eq2_adeno_median")

Eq3_shigella_med_CI <- CI_med_function(df_chunk_13_mult_weight, "Eq3_shigella_lower", "Eq3_shigella_upper", "Eq3_shigella_median")
Eq3_campy_med_CI <- CI_med_function(df_chunk_14_mult_weight, "Eq3_campy_lower", "Eq3_campy_upper", "Eq3_campy_median")
Eq3_ETEC_med_CI <- CI_med_function(df_chunk_15_mult_weight, "Eq3_ETEC_lower", "Eq3_ETEC_upper", "Eq3_ETEC_median")
Eq3_noro_med_CI <- CI_med_function(df_chunk_16_mult_weight, "Eq3_noro_lower", "Eq3_noro_upper", "Eq3_noro_median")
Eq3_rota_med_CI <- CI_med_function(df_chunk_17_mult_weight, "Eq3_rota_lower", "Eq3_rota_upper", "Eq3_rota_median")
Eq3_adeno_med_CI <- CI_med_function(df_chunk_18_mult_weight, "Eq3_adeno_lower", "Eq3_adeno_upper", "Eq3_adeno_median")

Eq4_shigella_med_CI <- CI_med_function(df_chunk_19_mult_weight, "Eq4_shigella_lower", "Eq4_shigella_upper", "Eq4_shigella_median")
Eq4_campy_med_CI <- CI_med_function(df_chunk_20_mult_weight, "Eq4_campy_lower", "Eq4_campy_upper", "Eq4_campy_median")
Eq4_ETEC_med_CI <- CI_med_function(df_chunk_21_mult_weight, "Eq4_ETEC_lower", "Eq4_ETEC_upper", "Eq4_ETEC_median")
Eq4_noro_med_CI <- CI_med_function(df_chunk_22_mult_weight, "Eq4_noro_lower", "Eq4_noro_upper", "Eq4_noro_median")
Eq4_rota_med_CI <- CI_med_function(df_chunk_23_mult_weight, "Eq4_rota_lower", "Eq4_rota_upper", "Eq4_rota_median")
Eq4_adeno_med_CI <- CI_med_function(df_chunk_24_mult_weight, "Eq4_adeno_lower", "Eq4_adeno_upper", "Eq4_adeno_median")

Eq5_shigella_med_CI <- CI_med_function(df_chunk_25_mult_weight, "Eq5_shigella_lower", "Eq5_shigella_upper", "Eq5_shigella_median")
Eq5_campy_med_CI <- CI_med_function(df_chunk_26_mult_weight, "Eq5_campy_lower", "Eq5_campy_upper", "Eq5_campy_median")
Eq5_ETEC_med_CI <- CI_med_function(df_chunk_27_mult_weight, "Eq5_ETEC_lower", "Eq5_ETEC_upper", "Eq5_ETEC_median")
Eq5_noro_med_CI <- CI_med_function(df_chunk_28_mult_weight, "Eq5_noro_lower", "Eq5_noro_upper", "Eq5_noro_median")
Eq5_rota_med_CI <- CI_med_function(df_chunk_29_mult_weight, "Eq5_rota_lower", "Eq5_rota_upper", "Eq5_rota_median")
Eq5_adeno_med_CI <- CI_med_function(df_chunk_30_mult_weight, "Eq5_adeno_lower", "Eq5_adeno_upper", "Eq5_adeno_median")

Eq6_shigella_med_CI <- CI_med_function(df_chunk_31_mult_weight, "Eq6_shigella_lower", "Eq6_shigella_upper", "Eq6_shigella_median")
Eq6_campy_med_CI <- CI_med_function(df_chunk_32_mult_weight, "Eq6_campy_lower", "Eq6_campy_upper", "Eq6_campy_median")
Eq6_ETEC_med_CI <- CI_med_function(df_chunk_33_mult_weight, "Eq6_ETEC_lower", "Eq6_ETEC_upper", "Eq6_ETEC_median")
Eq6_noro_med_CI <- CI_med_function(df_chunk_34_mult_weight, "Eq6_noro_lower", "Eq6_noro_upper", "Eq6_noro_median")
Eq6_rota_med_CI <- CI_med_function(df_chunk_35_mult_weight, "Eq6_rota_lower", "Eq6_rota_upper", "Eq6_rota_median")
Eq6_adeno_med_CI <- CI_med_function(df_chunk_36_mult_weight, "Eq6_adeno_lower", "Eq6_adeno_upper", "Eq6_adeno_median")

Eq7_shigella_med_CI <- CI_med_function(df_chunk_37_mult_weight, "Eq7_shigella_lower", "Eq7_shigella_upper", "Eq7_shigella_median")
Eq7_campy_med_CI <- CI_med_function(df_chunk_38_mult_weight, "Eq7_campy_lower", "Eq7_campy_upper", "Eq7_campy_median")
Eq7_ETEC_med_CI <- CI_med_function(df_chunk_39_mult_weight, "Eq7_ETEC_lower", "Eq7_ETEC_upper", "Eq7_ETEC_median")
Eq7_noro_med_CI <- CI_med_function(df_chunk_40_mult_weight, "Eq7_noro_lower", "Eq7_noro_upper", "Eq7_noro_median")
Eq7_rota_med_CI <- CI_med_function(df_chunk_41_mult_weight, "Eq7_rota_lower", "Eq7_rota_upper", "Eq7_rota_median")
Eq7_adeno_med_CI <- CI_med_function(df_chunk_42_mult_weight, "Eq7_adeno_lower", "Eq7_adeno_upper", "Eq7_adeno_median")

Eq8_shigella_med_CI <- CI_med_function(df_chunk_43_mult_weight, "Eq8_shigella_lower", "Eq8_shigella_upper", "Eq8_shigella_median")
Eq8_campy_med_CI <- CI_med_function(df_chunk_44_mult_weight, "Eq8_campy_lower", "Eq8_campy_upper", "Eq8_campy_median")
Eq8_ETEC_med_CI <- CI_med_function(df_chunk_45_mult_weight, "Eq8_ETEC_lower", "Eq8_ETEC_upper", "Eq8_ETEC_median")
Eq8_noro_med_CI <- CI_med_function(df_chunk_46_mult_weight, "Eq8_noro_lower", "Eq8_noro_upper", "Eq8_noro_median")
Eq8_rota_med_CI <- CI_med_function(df_chunk_47_mult_weight, "Eq8_rota_lower", "Eq8_rota_upper", "Eq8_rota_median")
Eq8_adeno_med_CI <- CI_med_function(df_chunk_48_mult_weight, "Eq8_adeno_lower", "Eq8_adeno_upper", "Eq8_adeno_median")

# Bind all of the dfs together (keep all in case they are wanted later, but really only need 1,4,3,6)

diarrhea_all_eq_med_CI_weighted <- cbind(Eq1_shigella_med_CI, Eq1_campy_med_CI, Eq1_ETEC_med_CI, Eq1_noro_med_CI, Eq1_rota_med_CI, Eq1_adeno_med_CI, 
                                     Eq2_shigella_med_CI, Eq2_campy_med_CI, Eq2_ETEC_med_CI, Eq2_noro_med_CI, Eq2_rota_med_CI, Eq2_adeno_med_CI, 
                                     Eq3_shigella_med_CI, Eq3_campy_med_CI, Eq3_ETEC_med_CI, Eq3_noro_med_CI, Eq3_rota_med_CI, Eq3_adeno_med_CI, 
                                     Eq4_shigella_med_CI, Eq4_campy_med_CI, Eq4_ETEC_med_CI, Eq4_noro_med_CI, Eq4_rota_med_CI, Eq4_adeno_med_CI, 
                                     Eq5_shigella_med_CI, Eq5_campy_med_CI, Eq5_ETEC_med_CI, Eq5_noro_med_CI, Eq5_rota_med_CI, Eq5_adeno_med_CI, 
                                     Eq6_shigella_med_CI, Eq6_campy_med_CI, Eq6_ETEC_med_CI, Eq6_noro_med_CI, Eq6_rota_med_CI, Eq6_adeno_med_CI,
                                     Eq7_shigella_med_CI, Eq7_campy_med_CI, Eq7_ETEC_med_CI, Eq7_noro_med_CI, Eq7_rota_med_CI, Eq7_adeno_med_CI,
                                     Eq8_shigella_med_CI, Eq8_campy_med_CI, Eq8_ETEC_med_CI, Eq8_noro_med_CI, Eq8_rota_med_CI, Eq8_adeno_med_CI)


# save 
write.csv(diarrhea_all_eq_med_CI_weighted, file="Scripts - Aim 3/Output/Diarrhea_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.csv")
save(diarrhea_all_eq_med_CI_weighted,file="Scripts - Aim 3/Output/Diarrhea_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")
