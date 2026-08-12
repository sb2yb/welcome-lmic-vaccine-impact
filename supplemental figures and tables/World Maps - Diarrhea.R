# this script will create world maps for the 6 equations in aim 3

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)

# bring in matrices 
# incidence
load("Scripts - Aim 3/Output/All_Diarrhea_Pathogens_AllEqs_Med_CI_DHS_24Nov2025.Rdata")

# bring in world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Remove Antarctica
world <- world %>% 
  filter(admin != "Antarctica")

# load country region list
load("Scripts - Aim 3/country_region_list.Rdata")

###############################
# clean country_region script #
###############################

# remove the 6 countries that we do not need: 
# Democratic Peoples Republic of Korea, Palau, Argentina, Occupied Palestinian territory, Syria, duplicate of Georgia
country_region_list <- country_region_list %>%
  filter(DisplayString!="Democratic People's Republic of Korea" & DisplayString!="Palau" & DisplayString!="Argentina" & 
           DisplayString!="occupied Palestinian territory, including east Jerusalem" & DisplayString!="Syrian Arab Republic") %>%
  select(countriesSub, DisplayString) %>%
  distinct() # this removes one instance of Georgia (doesn't matter which one since this is only a list and not tied to actual data)

# clean the name of 14 countries to match the world country list
country_region_list <- country_region_list %>%
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

##############################################
# THIS FUNCTION IS FOR EQUATIONS 1,2,4,5,6,8 #
##############################################

# make a function where it will loop through the DHS and IHME datasets
df_function <- function(df, eq_path, eq_path_name, title) {

# select all median values
eq <- df %>%
  select(Eq1_shigella_median, Eq1_campy_median, Eq1_ETEC_median, Eq1_noro_median, Eq1_rota_median, Eq1_adeno_median,
         Eq2_shigella_median, Eq2_campy_median, Eq2_ETEC_median, Eq2_noro_median, Eq2_rota_median, Eq2_adeno_median,
         Eq3_shigella_median, Eq3_campy_median, Eq3_ETEC_median, Eq3_noro_median, Eq3_rota_median, Eq3_adeno_median,
         Eq4_shigella_median, Eq4_campy_median, Eq4_ETEC_median, Eq4_noro_median, Eq4_rota_median, Eq4_adeno_median,
         Eq5_shigella_median, Eq5_campy_median, Eq5_ETEC_median, Eq5_noro_median, Eq5_rota_median, Eq5_adeno_median,
         Eq6_shigella_median, Eq6_campy_median, Eq6_ETEC_median, Eq6_noro_median, Eq6_rota_median, Eq6_adeno_median,
         Eq7_shigella_median, Eq7_campy_median, Eq7_ETEC_median, Eq7_noro_median, Eq7_rota_median, Eq7_adeno_median,
         Eq8_shigella_median, Eq8_campy_median, Eq8_ETEC_median, Eq8_noro_median, Eq8_rota_median, Eq8_adeno_median)

# add the country_region_list to the all_path_all_eq_med_CI_IHME
eq_countries <- cbind(country_region_list, eq)

# number 1-129 so we can ensure that the order does not change with the merge
eq_countries$order <- seq(1, 129)

# add eq1 countries to world
merge <- merge(eq_countries, world, by.x = "DisplayString", by.y="admin", all.y=T)

# sort by order
merge <- merge[order(merge$order), ]


#########
# Plots #
#########

# Make a function over this so we can swap out pathogens
Eq <- ggplot(data = merge) 
Eq = Eq + geom_sf(aes(geometry=geometry, fill= !!rlang::sym(eq_path)))
Eq = Eq + scale_fill_gradient(name = eq_path_name,
                                          low = "lightpink", high = "deeppink", na.value = "gray", guide = "colorbar")
Eq = Eq + theme_void() 
Eq = Eq + theme(legend.position = "bottom")
Eq = Eq +  labs(title = title)
Eq 

return(Eq)
}

#######################################################################################################
#                            THIS FUNCTION IS FOR EQUATIONS 3 AND 7                                   #
# SINCE THE VALUES ARE NEGATIVE THE COLORING NEEDS TO BE FLIPPED TO MAKE THE INTERPRETATION INTUITIVE #
#######################################################################################################

# make a function where it will loop through the DHS and IHME datasets
df_function_Eq3_7 <- function(df, eq_path, eq_path_name, title) {
  
  # select all median values
  eq <- df %>%
    select(Eq1_shigella_median, Eq1_campy_median, Eq1_ETEC_median, Eq1_noro_median, Eq1_rota_median, Eq1_adeno_median,
           Eq2_shigella_median, Eq2_campy_median, Eq2_ETEC_median, Eq2_noro_median, Eq2_rota_median, Eq2_adeno_median,
           Eq3_shigella_median, Eq3_campy_median, Eq3_ETEC_median, Eq3_noro_median, Eq3_rota_median, Eq3_adeno_median,
           Eq4_shigella_median, Eq4_campy_median, Eq4_ETEC_median, Eq4_noro_median, Eq4_rota_median, Eq4_adeno_median,
           Eq5_shigella_median, Eq5_campy_median, Eq5_ETEC_median, Eq5_noro_median, Eq5_rota_median, Eq5_adeno_median,
           Eq6_shigella_median, Eq6_campy_median, Eq6_ETEC_median, Eq6_noro_median, Eq6_rota_median, Eq6_adeno_median,
           Eq7_shigella_median, Eq7_campy_median, Eq7_ETEC_median, Eq7_noro_median, Eq7_rota_median, Eq7_adeno_median,
           Eq8_shigella_median, Eq8_campy_median, Eq8_ETEC_median, Eq8_noro_median, Eq8_rota_median, Eq8_adeno_median)
  
  # add the country_region_list to the all_path_all_eq_med_CI_IHME
  eq_countries <- cbind(country_region_list, eq)
  
  # number 1-129 so we can ensure that the order does not change with the merge
  eq_countries$order <- seq(1, 129)
  
  # add eq1 countries to world
  merge <- merge(eq_countries, world, by.x = "DisplayString", by.y="admin", all.y=T)
  
  # sort by order
  merge <- merge[order(merge$order), ]
  
  
  #########
  # Plots #
  #########
  
  # Make a function over this so we can swap out pathogens
  Eq <- ggplot(data = merge) 
  Eq = Eq + geom_sf(aes(geometry=geometry, fill= !!rlang::sym(eq_path)))
  Eq = Eq + scale_fill_gradient(name = eq_path_name,
                                high = "lightpink", low = "deeppink", na.value = "gray", guide = "colorbar")
  Eq = Eq + theme_void() 
  Eq = Eq + theme(legend.position = "bottom")
  Eq = Eq +  labs(title = title)
  Eq 
  
  return(Eq)
}

#######################################################################################################
#                            THIS FUNCTION IS FOR EQUATIONS 3 - SHIGELLA & RSV                        #
# SINCE THE VALUES ARE NEGATIVE THE COLORING NEEDS TO BE FLIPPED TO MAKE THE INTERPRETATION INTUITIVE #
#######################################################################################################

# make a function where it will loop through the DHS and IHME datasets
df_function_Eq3_SHIG_RSV <- function(df, eq_path, eq_path_name, title) {
  
  # select all median values
  eq <- df %>%
    select(Eq1_shigella_median, Eq1_campy_median, Eq1_ETEC_median, Eq1_noro_median, Eq1_rota_median, Eq1_adeno_median,
           Eq2_shigella_median, Eq2_campy_median, Eq2_ETEC_median, Eq2_noro_median, Eq2_rota_median, Eq2_adeno_median,
           Eq3_shigella_median, Eq3_campy_median, Eq3_ETEC_median, Eq3_noro_median, Eq3_rota_median, Eq3_adeno_median,
           Eq4_shigella_median, Eq4_campy_median, Eq4_ETEC_median, Eq4_noro_median, Eq4_rota_median, Eq4_adeno_median,
           Eq5_shigella_median, Eq5_campy_median, Eq5_ETEC_median, Eq5_noro_median, Eq5_rota_median, Eq5_adeno_median,
           Eq6_shigella_median, Eq6_campy_median, Eq6_ETEC_median, Eq6_noro_median, Eq6_rota_median, Eq6_adeno_median,
           Eq7_shigella_median, Eq7_campy_median, Eq7_ETEC_median, Eq7_noro_median, Eq7_rota_median, Eq7_adeno_median,
           Eq8_shigella_median, Eq8_campy_median, Eq8_ETEC_median, Eq8_noro_median, Eq8_rota_median, Eq8_adeno_median)
  
  # add the country_region_list to the all_path_all_eq_med_CI_IHME
  eq_countries <- cbind(country_region_list, eq)
  
  # number 1-129 so we can ensure that the order does not change with the merge
  eq_countries$order <- seq(1, 129)
  
  # add eq1 countries to world
  merge <- merge(eq_countries, world, by.x = "DisplayString", by.y="admin", all.y=T)
  
  # sort by order
  merge <- merge[order(merge$order), ]
  
  
  #########
  # Plots #
  #########
  
  # Make a function over this so we can swap out pathogens
  Eq <- ggplot(data = merge) 
  Eq = Eq + geom_sf(aes(geometry=geometry, fill= !!rlang::sym(eq_path)))
  Eq = Eq + scale_fill_gradient(name = eq_path_name,
                                high = "lightpink", low = "deeppink", na.value = "gray", guide = "colorbar", limits = c(-10, 0))
  Eq = Eq + theme_void() 
  Eq = Eq + theme(legend.position = "bottom")
  Eq = Eq +  labs(title = title)
  Eq 
  
  return(Eq)
}

#######################################################################################################
#                            THIS FUNCTION IS FOR EQUATIONS 7 - SHIGELLA & RSV                        #
# SINCE THE VALUES ARE NEGATIVE THE COLORING NEEDS TO BE FLIPPED TO MAKE THE INTERPRETATION INTUITIVE #
#######################################################################################################

# make a function where it will loop through the DHS and IHME datasets
df_function_Eq7_SHIG_RSV <- function(df, eq_path, eq_path_name, title) {
  
  # select all median values
  eq <- df %>%
    select(Eq1_shigella_median, Eq1_campy_median, Eq1_ETEC_median, Eq1_noro_median, Eq1_rota_median, Eq1_adeno_median,
           Eq2_shigella_median, Eq2_campy_median, Eq2_ETEC_median, Eq2_noro_median, Eq2_rota_median, Eq2_adeno_median,
           Eq3_shigella_median, Eq3_campy_median, Eq3_ETEC_median, Eq3_noro_median, Eq3_rota_median, Eq3_adeno_median,
           Eq4_shigella_median, Eq4_campy_median, Eq4_ETEC_median, Eq4_noro_median, Eq4_rota_median, Eq4_adeno_median,
           Eq5_shigella_median, Eq5_campy_median, Eq5_ETEC_median, Eq5_noro_median, Eq5_rota_median, Eq5_adeno_median,
           Eq6_shigella_median, Eq6_campy_median, Eq6_ETEC_median, Eq6_noro_median, Eq6_rota_median, Eq6_adeno_median,
           Eq7_shigella_median, Eq7_campy_median, Eq7_ETEC_median, Eq7_noro_median, Eq7_rota_median, Eq7_adeno_median,
           Eq8_shigella_median, Eq8_campy_median, Eq8_ETEC_median, Eq8_noro_median, Eq8_rota_median, Eq8_adeno_median)
  
  # add the country_region_list to the all_path_all_eq_med_CI_IHME
  eq_countries <- cbind(country_region_list, eq)
  
  # number 1-129 so we can ensure that the order does not change with the merge
  eq_countries$order <- seq(1, 129)
  
  # add eq1 countries to world
  merge <- merge(eq_countries, world, by.x = "DisplayString", by.y="admin", all.y=T)
  
  # sort by order
  merge <- merge[order(merge$order), ]
  
  
  #########
  # Plots #
  #########
  
  # Make a function over this so we can swap out pathogens
  Eq <- ggplot(data = merge) 
  Eq = Eq + geom_sf(aes(geometry=geometry, fill= !!rlang::sym(eq_path)))
  Eq = Eq + scale_fill_gradient(name = eq_path_name,
                                high = "lightpink", low = "deeppink", na.value = "gray", guide = "colorbar", limits = c(-15, 0))
  Eq = Eq + theme_void() 
  Eq = Eq + theme(legend.position = "bottom")
  Eq = Eq +  labs(title = title)
  Eq 
  
  return(Eq)
}


##############
# EQUATION 1 #
##############

Eq1_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_shigella_median", "Eq1_shigella", "DHS incidence of Shigella specific treatment for diarrhea by country")
Eq1_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_campy_median", "Eq1_campy", "DHS incidence of Campylobacter specific treatment for diarrhea by country")
Eq1_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_ETEC_median", "Eq1_etec", "DHS incidence of ETEC specific treatment for diarrhea by country")
Eq1_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_noro_median", "Eq1_noro", "DHS incidence of Norovirus specific treatment for diarrhea by country")
Eq1_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_rota_median", "Eq1_rota", "DHS incidence of Rotavirus specific treatment for diarrhea by country")
Eq1_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq1_adeno_median", "Eq1_adeno", "DHS incidence of Adenovirus specific treatment for diarrhea by country")

##############
# EQUATION 2 #
##############
 
Eq2_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_shigella_median", "Eq2_shigella", "DHS incidence of Shigella specific treatment for diarrhea with a Shigella vaccine by country")
Eq2_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_campy_median", "Eq2_campy", "DHS incidence of Campylobacter specific treatment for diarrhea with a Campylobacter vaccine by country")
Eq2_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_ETEC_median", "Eq2_etec", "DHS incidence of ETEC specific treatment for diarrhea with an ETEC vaccine by country")
Eq2_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_noro_median", "Eq2_noro", "DHS incidence of Norovirus specific treatment for diarrhea with a Norovirus vaccine by country")
Eq2_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_rota_median", "Eq2_rota", "DHS incidence of Rotavirus specific treatment for diarrhea with a Rotavirus vaccine by country")
Eq2_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq2_adeno_median", "Eq2_adeno", "DHS incidence of Adenovirus specific treatment for diarrhea with an Adenovirus vaccine by country")

##############
# EQUATION 3 #
##############

Eq3_shigella_DHS<- df_function_Eq3_SHIG_RSV(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_shigella_median", "Eq3_shigella", "DHS absolute reductions in antibiotic use with a Shigella vaccine by country")
Eq3_campy_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_campy_median", "Eq3_campy", "DHS absolute reductions in antibiotic use with a Campylobacter vaccine by country")
Eq3_ETEC_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_ETEC_median", "Eq3_etec", "DHS absolute reductions in antibiotic use with an ETEC vaccine by country")
Eq3_noro_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_noro_median", "Eq3_noro", "DHS absolute reductions in antibiotic use with a Norovirus vaccine by country")
Eq3_rota_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_rota_median", "Eq3_rota", "DHS absolute reductions in antibiotic use with a Rotavirus vaccine by country")
Eq3_adeno_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq3_adeno_median", "Eq3_adeno", "DHS absolute reductions in antibiotic use with an Adenovirus vaccine by country")

##############
# EQUATION 4 #
##############

Eq4_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_shigella_median", "Eq4_shigella", "DHS percent reductions in antibiotic use with a Shigella vaccine by country")
Eq4_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_campy_median", "Eq4_campy", "DHS percent reductions in antibiotic use with a Campylobacter vaccine by country")
Eq4_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_ETEC_median", "Eq4_etec", "DHS percent reductions in antibiotic use with an ETEC vaccine by country")
Eq4_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_noro_median", "Eq4_noro", "DHS percent reductions in antibiotic use with a Norovirus vaccine by country")
Eq4_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_rota_median", "Eq4_rota", "DHS percent reductions in antibiotic use with a Rotavirus vaccine by country")
Eq4_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq4_adeno_median", "Eq4_adeno", "DHS percent reductions in antibiotic use with an Adenovirus vaccine by country")

##############
# EQUATION 5 #  
##############

Eq5_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_shigella_median", "Eq5_shigella", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of Shigella by country")
Eq5_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_campy_median", "Eq5_campy", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of Campylobacter by country")
Eq5_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_ETEC_median", "Eq5_etec", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of ETEC by country")
Eq5_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_noro_median", "Eq5_noro", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of Norovirus by country")
Eq5_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_rota_median", "Eq5_rota", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of Rotavirus by country")
Eq5_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq5_adeno_median", "Eq5_adeno", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of Adenovirus by country")

##############
# EQUATION 6 # 
##############

Eq6_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_shigella_median", "Eq6_shigella", "DHS incidence of bystander pathogen exposures to antibiotics with a Shigella vaccine by country")
Eq6_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_campy_median", "Eq6_campy", "DHS incidence of bystander pathogen exposures to antibiotics with a Campylobacter vaccine by country")
Eq6_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_ETEC_median", "Eq6_etec", "DHS incidence of bystander pathogen exposures to antibiotics with an ETEC vaccine by country")
Eq6_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_noro_median", "Eq6_noro", "DHS incidence of bystander pathogen exposures to antibiotics with a Norovirus vaccine by country")
Eq6_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_rota_median", "Eq6_rota", "DHS incidence of bystander pathogen exposures to antibiotics with a Rotavirus vaccine by country")
Eq6_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq6_adeno_median", "Eq6_adeno", "DHS incidence of bystander pathogen exposures to antibiotics with an Adenovirus vaccine by country")

##############
# EQUATION 7 #
##############

Eq7_shigella_DHS<- df_function_Eq7_SHIG_RSV(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_shigella_median", "Eq7_shigella", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a Shigella vaccine by country")
Eq7_campy_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_campy_median", "Eq7_campy", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a Campylobacter vaccine by country")
Eq7_ETEC_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_ETEC_median", "Eq7_etec", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a ETEC vaccine by country")
Eq7_noro_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_noro_median", "Eq7_noro", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a Norovirus vaccine by country")
Eq7_rota_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_rota_median", "Eq7_rota", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a Rotavirus vaccine by country")
Eq7_adeno_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_DIAR, "Eq7_adeno_median", "Eq7_adeno", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a Adenovirus vaccine by country")

##############
# EQUATION 8 #
##############

Eq8_shigella_DHS<- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_shigella_median", "Eq8_shigella", "DHS percent reductions of bystander pathogen exposures to antibiotics with a Shigella vaccine by country")
Eq8_campy_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_campy_median", "Eq8_campy", "DHS percent reductions of bystander pathogen exposures to antibiotics with a Campylobacter vaccine by country")
Eq8_ETEC_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_ETEC_median", "Eq8_etec", "DHS percent reductions of bystander pathogen exposures to antibiotics with a ETEC vaccine by country")
Eq8_noro_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_noro_median", "Eq8_noro", "DHS percent reductions of bystander pathogen exposures to antibiotics with a Norovirus vaccine by country")
Eq8_rota_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_rota_median", "Eq8_rota", "DHS percent reductions of bystander pathogen exposures to antibiotics with a Rotavirus vaccine by country")
Eq8_adeno_DHS <- df_function(all_path_all_eq_med_CI_DHS_DIAR, "Eq8_adeno_median", "Eq8_adeno", "DHS percent reductions of bystander pathogen exposures to antibiotics with a Adenovirus vaccine by country")


###############
# Saving Maps #
###############

# EQUATION 1
ggsave(file="Scripts - Aim 3/Maps/Eq1_Shigella_Map_DHS.png",Eq1_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_Adenovirus_Map_DHS.png",Eq1_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_Campy_Map_DHS.png",Eq1_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_ETEC_Map_DHS.png",Eq1_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_Norovirus_Map_DHS.png",Eq1_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_Rotavirus_Map_DHS.png",Eq1_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 2
ggsave(file="Scripts - Aim 3/Maps/Eq2_Shigella_Map_DHS.png",Eq2_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_Adenovirus_Map_DHS.png",Eq2_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_Campy_Map_DHS.png",Eq2_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_ETEC_Map_DHS.png",Eq2_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_Norovirus_Map_DHS.png",Eq2_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_Rotavirus_Map_DHS.png",Eq2_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 3
ggsave(file="Scripts - Aim 3/Maps/Eq3_Shigella_Map_DHS_24Nov2025.png",Eq3_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Adenovirus_Map_DHS_24Nov2025.png",Eq3_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Campy_Map_DHS_24Nov2025.png",Eq3_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_ETEC_Map_DHS_24Nov2025.png",Eq3_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Norovirus_Map_DHS_24Nov2025.png",Eq3_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Rotavirus_Map_DHS_24Nov2025.png",Eq3_rota_DHS,width=22,height=10,dpi=300)

ggsave(file="Scripts - Aim 3/Maps/Eq3_Shigella_Map_DHS_24Nov2025.pdf",Eq3_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Adenovirus_Map_DHS_24Nov2025.pdf",Eq3_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Campy_Map_DHS_24Nov2025.pdf",Eq3_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_ETEC_Map_DHS_24Nov2025.pdf",Eq3_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Norovirus_Map_DHS_24Nov2025.pdf",Eq3_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_Rotavirus_Map_DHS_24Nov2025.pdf",Eq3_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 4
ggsave(file="Scripts - Aim 3/Maps/Eq4_Shigella_Map_DHS.png",Eq4_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_Adenovirus_Map_DHS.png",Eq4_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_Campy_Map_DHS.png",Eq4_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_ETEC_Map_DHS.png",Eq4_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_Norovirus_Map_DHS.png",Eq4_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_Rotavirus_Map_DHS.png",Eq4_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 5
ggsave(file="Scripts - Aim 3/Maps/Eq5_Shigella_Map_DHS.png",Eq5_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_Adenovirus_Map_DHS.png",Eq5_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_Campy_Map_DHS.png",Eq5_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_ETEC_Map_DHS.png",Eq5_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_Norovirus_Map_DHS.png",Eq5_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_Rotavirus_Map_DHS.png",Eq5_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 6
ggsave(file="Scripts - Aim 3/Maps/Eq6_Shigella_Map_DHS.png",Eq6_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_Adenovirus_Map_DHS.png",Eq6_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_Campy_Map_DHS.png",Eq6_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_ETEC_Map_DHS.png",Eq6_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_Norovirus_Map_DHS.png",Eq6_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_Rotavirus_Map_DHS.png",Eq6_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 7
ggsave(file="Scripts - Aim 3/Maps/Eq7_Shigella_Map_DHS_24Nov2025.png",Eq7_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Adenovirus_Map_DHS_24Nov2025.png",Eq7_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Campy_Map_DHS_24Nov2025.png",Eq7_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_ETEC_Map_DHS_24Nov2025.png",Eq7_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Norovirus_Map_DHS_24Nov2025.png",Eq7_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Rotavirus_Map_DHS_24Nov2025.png",Eq7_rota_DHS,width=22,height=10,dpi=300)

ggsave(file="Scripts - Aim 3/Maps/Eq7_Shigella_Map_DHS_24Nov2025.pdf",Eq7_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Adenovirus_Map_DHS_24Nov2025.pdf",Eq7_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Campy_Map_DHS_24Nov2025.pdf",Eq7_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_ETEC_Map_DHS_24Nov2025.pdf",Eq7_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Norovirus_Map_DHS_24Nov2025.pdf",Eq7_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_Rotavirus_Map_DHS_24Nov2025.pdf",Eq7_rota_DHS,width=22,height=10,dpi=300)

# EQUATION 8
ggsave(file="Scripts - Aim 3/Maps/Eq8_Shigella_Map_DHS.png",Eq8_shigella_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_Adenovirus_Map_DHS.png",Eq8_adeno_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_Campy_Map_DHS.png",Eq8_campy_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_ETEC_Map_DHS.png",Eq8_ETEC_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_Norovirus_Map_DHS.png",Eq8_noro_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_Rotavirus_Map_DHS.png",Eq8_rota_DHS,width=22,height=10,dpi=300)

saveRDS(Eq3_shigella_DHS, "Scripts - Aim 3/Output/Eq3_shigella_DHS.rds")
saveRDS(Eq7_shigella_DHS, "Scripts - Aim 3/Output/Eq7_shigella_DHS.rds")
