# this script will go through respiratory pathogens and all equations for Aim 3
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
load("Scripts - Aim 3/Output/DHS_ari_matrix.Rdata")

# vaccine efficacy
load("Scripts - Aim 3/Joe data/veHib.Rdata")
load("Scripts - Aim 3/Joe data/vePcv.Rdata")
load("Scripts - Aim 3/Joe data/veRsv.Rdata")

# vaccine coverage 
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

# the vaccine matrices only have 1 column. we need to make 129 identical columns for all countries. 
veHib <- replicate(129, veHib)
vePcv <- replicate(129, vePcv)
veRsv <- replicate(129, veRsv)

# the vaccine matrices have 10,000 rows. we only need 1000
veHib <- veHib[1:1000,]
vePcv <- vePcv[1:1000,]
veRsv <- veRsv[1:1000,]

# make dataframes to make it easy to remove columns (this can be done as a matrix but then the column names get re-numbered and for tracking reasons i'd rather that not happen)
DHS_ari <- as.data.frame(DHS_ari)

hib_vaccine_matrix <- as.data.frame(veHib)
pcv_vaccine_matrix <- as.data.frame(vePcv)
rsv_vaccine_matrix <- as.data.frame(veRsv)

bystander_country_estimates_final <- as.data.frame(bystander_country_estimates_final)

# FOR ALL LOADED SCRIPTS: REMOVE THE FOLLOWING #

# 1) Remove: #45 Democratic Peoples Republic of Korea, #84 Palau, and #108 Argentina as they go back and forth on being an LMIC (data not in DHS at all)
# 2) Remove: #75 Occupied Palestinian territory, including east Jerusalem as we don't have any data to extrapolate to this
# 3) Remove: #76 Syria (due to conflict)
# 4) Remove: #18 Georgia (the one located in Position #17 is the correct one as confirmed by Joe)
DHS_ari <- DHS_ari %>% select(-c(V18, V45, V84, V108, V75, V76))
bystander_country_estimates_final <- bystander_country_estimates_final %>% select(-c(V18, V45, V84, V108, V75, V76))

# creating individual vaccine coverage dataframes 
hib_coverage <- coverage %>% select(HIB3)
pcv_coverage <- coverage %>% select(PCV3)
rsv_coverage <- coverage %>% select(avg_vaxcov) #using average as rsv doesn't have coverage data

# coverage has to be transformed so that it has 129 columns and 1000 rows. 
hibcov_wide <- as.data.frame(t(hib_coverage$HIB3))
pcvcov_wide <- as.data.frame(t(pcv_coverage$PCV3))
rsvcov_wide <- as.data.frame(t(rsv_coverage$avg_vaxcov)) #using average as rsv doesn't have coverage data

# the _wide datasets need to have 1000 rows
hibcov_wide <- hibcov_wide[rep(1, 1000), ]
pcvcov_wide <- pcvcov_wide[rep(1, 1000), ]
rsvcov_wide <- rsvcov_wide[rep(1, 1000), ]

#############
# Equations #
#############

##################################################################
# Equation 1: Incidence of treatment for respiratory infections  #
##################################################################
# Incidence of antibiotic use for respiratory infections 

Eq1_hib <- DHS_ari 
Eq1_pcv <- DHS_ari 
Eq1_rsv <- DHS_ari 

##################################################################################################
# Equation 2: Incidence of pathogen specific treatment for respiratory infections with a vaccine #
##################################################################################################
# Equation 1 * (1-vaccine efficacy * vaccine coverage)

# vaccine efficacies need to be numeric
hib_vaccine_matrix <- hib_vaccine_matrix %>% mutate_all(as.numeric)
pcv_vaccine_matrix <- pcv_vaccine_matrix %>% mutate_all(as.numeric)
rsv_vaccine_matrix <- rsv_vaccine_matrix %>% mutate_all(as.numeric)


Eq2_hib <- Eq1_hib * (1- hib_vaccine_matrix * hibcov_wide)
Eq2_pcv <- Eq1_pcv * (1- pcv_vaccine_matrix * pcvcov_wide)
Eq2_rsv <- Eq1_rsv * (1- rsv_vaccine_matrix * rsvcov_wide)


#######################################################################
# Equation 3: The absolute reduction in antibiotic use with a vaccine #
#######################################################################
# Equation 2 - Equation 1
# multiplying by 100 to make this 100 child years (DHS is per 1 child year)
Eq3_hib <- (Eq2_hib - Eq1_hib) * 100
Eq3_pcv <- (Eq2_pcv - Eq1_pcv) * 100
Eq3_rsv <- (Eq2_rsv - Eq1_rsv) * 100

######################################################################
# Equation 4: The percent reduction in antibiotic use with a vaccine #
######################################################################
# (1 - (Equation 2 / Equation 1)) * 100

Eq4_hib <- (1-(Eq2_hib / Eq1_hib))*100
Eq4_pcv <- (1-(Eq2_pcv / Eq1_pcv))*100
Eq4_rsv <- (1-(Eq2_rsv / Eq1_rsv))*100

###################################################################
# Equation 5: The incidence of bystander exposures to antibiotics #
###################################################################
# Equation 1 * average number of asymptomatic pathogens

# bystander pathogens need to be numeric
bystander_country_estimates_final <- bystander_country_estimates_final %>% mutate_all(as.numeric)

Eq5_hib <- Eq1_hib * bystander_country_estimates_final
Eq5_pcv <- Eq1_pcv * bystander_country_estimates_final
Eq5_rsv <- Eq1_rsv * bystander_country_estimates_final

#############################################################################
#Equation 6: Incidence of bystander exposures to antibiotics with a vaccine #
#############################################################################
#Equation 2 * average number of asymptomatic pathogens

Eq6_hib <- Eq2_hib * bystander_country_estimates_final
Eq6_pcv <- Eq2_pcv * bystander_country_estimates_final
Eq6_rsv <- Eq2_rsv * bystander_country_estimates_final

###########################################################################################
# Equation 7: The absolute reduction of bystander exposures to antibiotics with a vaccine #
###########################################################################################
# Equation 6 - Equation 5
# multiplying by 100 to make this 100 child years (DHS is per 1 child year)
Eq7_hib <- (Eq6_hib - Eq5_hib) * 100
Eq7_pcv <- (Eq6_pcv - Eq5_pcv) * 100
Eq7_rsv <- (Eq6_rsv - Eq5_rsv) * 100

###########################################################################################
# Equation 8: The absolute reduction of bystander exposures to antibiotics with a vaccine #
###########################################################################################
# (1-(Equation 6 / Equation 5))*100

Eq8_hib <- (1-(Eq6_hib / Eq5_hib))*100
Eq8_pcv <- (1-(Eq6_pcv / Eq5_pcv))*100
Eq8_rsv <- (1-(Eq6_rsv / Eq5_rsv))*100

## bind dataframes together

bind <- cbind(Eq1_hib, Eq1_pcv, Eq1_rsv,
              Eq2_hib, Eq2_pcv, Eq2_rsv,
              Eq3_hib, Eq3_pcv, Eq3_rsv,
              Eq4_hib, Eq4_pcv, Eq4_rsv,
              Eq5_hib, Eq5_pcv, Eq5_rsv,
              Eq6_hib, Eq6_pcv, Eq6_rsv,
              Eq7_hib, Eq7_pcv, Eq7_rsv,
              Eq8_hib, Eq8_pcv, Eq8_rsv)

################################################
# Make one dataframe per pathogen per equation #
################################################

# Having R make a new dataframe every 129 columns for the 3 pathogens and 8 equations = 24 dataframes (keeps me from having to manually do this)
total_columns <- 3096 # total number of columns in bind 
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

# define which columns belong to which regions (we need this for the population)
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

# CHECK ON THIS
 # in the function df_chunks (V1-V135) is multiplied by country_pop_est_wide (V1-V129)
CHECK <- df_chunk_1 * country_pop_est_wide

CHECK1 <- CHECK %>%
  mutate(
    AFR = rowSums(across(all_of(AFR_cols)), na.rm = TRUE) / region_pop["AFR"],
    AMR = rowSums(across(all_of(AMR_cols)), na.rm = TRUE) / region_pop["AMR"],
    EMR = rowSums(across(all_of(EMR_cols)), na.rm = TRUE) / region_pop["EMR"],
    EUR = rowSums(across(all_of(EUR_cols)), na.rm = TRUE) / region_pop["EUR"],
    SEAR = rowSums(across(all_of(SEAR_cols)), na.rm = TRUE) / region_pop["SEAR"],
    WPR = rowSums(across(all_of(WPR_cols)), na.rm = TRUE) / region_pop["WPR"])

  
# calling the function for all 24 chunks
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

# we only want to keep the last 6 columns of each chunk (because those are where our weighted estimates are)
# Loop over all 48 weighted dataframes
for(i in 1:24) {
  
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

Eq1_hib_med_CI <- CI_med_function(df_chunk_1_mult_weight, "Eq1_hib_lower", "Eq1_hib_upper", "Eq1_hib_median")
Eq1_pcv_med_CI <- CI_med_function(df_chunk_2_mult_weight, "Eq1_pcv_lower", "Eq1_pcv_upper", "Eq1_pcv_median")
Eq1_rsv_med_CI <- CI_med_function(df_chunk_3_mult_weight, "Eq1_rsv_lower", "Eq1_rsv_upper", "Eq1_rsv_median")

Eq2_hib_med_CI <- CI_med_function(df_chunk_4_mult_weight, "Eq2_hib_lower", "Eq2_hib_upper", "Eq2_hib_median")
Eq2_pcv_med_CI <- CI_med_function(df_chunk_5_mult_weight, "Eq2_pcv_lower", "Eq2_pcv_upper", "Eq2_pcv_median")
Eq2_rsv_med_CI <- CI_med_function(df_chunk_6_mult_weight, "Eq2_rsv_lower", "Eq2_rsv_upper", "Eq2_rsv_median")

Eq3_hib_med_CI <- CI_med_function(df_chunk_7_mult_weight, "Eq3_hib_lower", "Eq3_hib_upper", "Eq3_hib_median") 
Eq3_pcv_med_CI <- CI_med_function(df_chunk_8_mult_weight, "Eq3_pcv_lower", "Eq3_pcv_upper", "Eq3_pcv_median")
Eq3_rsv_med_CI <- CI_med_function(df_chunk_9_mult_weight, "Eq3_rsv_lower", "Eq3_rsv_upper", "Eq3_rsv_median")

Eq4_hib_med_CI <- CI_med_function(df_chunk_10_mult_weight, "Eq4_hib_lower", "Eq4_hib_upper", "Eq4_hib_median")
Eq4_pcv_med_CI <- CI_med_function(df_chunk_11_mult_weight, "Eq4_pcv_lower", "Eq4_pcv_upper", "Eq4_pcv_median")
Eq4_rsv_med_CI <- CI_med_function(df_chunk_12_mult_weight, "Eq4_rsv_lower", "Eq4_rsv_upper", "Eq4_rsv_median")

Eq5_hib_med_CI <- CI_med_function(df_chunk_13_mult_weight, "Eq5_hib_lower", "Eq5_hib_upper", "Eq5_hib_median")
Eq5_pcv_med_CI <- CI_med_function(df_chunk_14_mult_weight, "Eq5_pcv_lower", "Eq5_pcv_upper", "Eq5_pcv_median")
Eq5_rsv_med_CI <- CI_med_function(df_chunk_15_mult_weight, "Eq5_rsv_lower", "Eq5_rsv_upper", "Eq5_rsv_median")

Eq6_hib_med_CI <- CI_med_function(df_chunk_16_mult_weight, "Eq6_hib_lower", "Eq6_hib_upper", "Eq6_hib_median")
Eq6_pcv_med_CI <- CI_med_function(df_chunk_17_mult_weight, "Eq6_pcv_lower", "Eq6_pcv_upper", "Eq6_pcv_median")
Eq6_rsv_med_CI <- CI_med_function(df_chunk_18_mult_weight, "Eq6_rsv_lower", "Eq6_rsv_upper", "Eq6_rsv_median")

Eq7_hib_med_CI <- CI_med_function(df_chunk_19_mult_weight, "Eq7_hib_lower", "Eq7_hib_upper", "Eq7_hib_median") 
Eq7_pcv_med_CI <- CI_med_function(df_chunk_20_mult_weight, "Eq7_pcv_lower", "Eq7_pcv_upper", "Eq7_pcv_median")
Eq7_rsv_med_CI <- CI_med_function(df_chunk_21_mult_weight, "Eq7_rsv_lower", "Eq7_rsv_upper", "Eq7_rsv_median")

Eq8_hib_med_CI <- CI_med_function(df_chunk_22_mult_weight, "Eq8_hib_lower", "Eq8_hib_upper", "Eq8_hib_median")
Eq8_pcv_med_CI <- CI_med_function(df_chunk_23_mult_weight, "Eq8_pcv_lower", "Eq8_pcv_upper", "Eq8_pcv_median")
Eq8_rsv_med_CI <- CI_med_function(df_chunk_24_mult_weight, "Eq8_rsv_lower", "Eq8_rsv_upper", "Eq8_rsv_median")

# Bind all of the dfs together (keep all in case they are wanted later, but really only need 1,4,3,6)

ari_all_eq_med_CI_weighted <-  cbind(Eq1_hib_med_CI, Eq1_pcv_med_CI, Eq1_rsv_med_CI, 
                                     Eq2_hib_med_CI, Eq2_pcv_med_CI, Eq2_rsv_med_CI, 
                                     Eq3_hib_med_CI, Eq3_pcv_med_CI, Eq3_rsv_med_CI,
                                     Eq4_hib_med_CI, Eq4_pcv_med_CI, Eq4_rsv_med_CI,
                                     Eq5_hib_med_CI, Eq5_pcv_med_CI, Eq5_rsv_med_CI, 
                                     Eq6_hib_med_CI, Eq6_pcv_med_CI, Eq6_rsv_med_CI,
                                     Eq7_hib_med_CI, Eq7_pcv_med_CI, Eq7_rsv_med_CI, 
                                     Eq8_hib_med_CI, Eq8_pcv_med_CI, Eq8_rsv_med_CI)


# save dataframes 
write.csv(ari_all_eq_med_CI_weighted, file="Scripts - Aim 3/Output/ARI_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.csv")
save(ari_all_eq_med_CI_weighted,file="Scripts - Aim 3/Output/ARI_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")


