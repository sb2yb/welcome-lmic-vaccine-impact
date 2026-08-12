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
load("Scripts - Aim 3/Output/All_ARI_Pathogens_AllEqs_Med_CI_DHS_24Nov2025.Rdata")

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
  select(Eq1_hib_median, Eq1_pcv_median, Eq1_rsv_median, # !! ALL RSV NEEDS TO BE UPDATED WITH ACTUAL RSV DATA - THIS IS MCV AS A PLACEHOLDER !!
         Eq2_hib_median, Eq2_pcv_median, Eq2_rsv_median,
         Eq3_hib_median, Eq3_pcv_median, Eq3_rsv_median,
         Eq4_hib_median, Eq4_pcv_median, Eq4_rsv_median,
         Eq5_hib_median, Eq5_pcv_median, Eq5_rsv_median, 
         Eq6_hib_median, Eq6_pcv_median, Eq6_rsv_median,
         Eq7_hib_median, Eq7_pcv_median, Eq7_rsv_median, 
         Eq8_hib_median, Eq8_pcv_median, Eq8_rsv_median)

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
    select(Eq1_hib_median, Eq1_pcv_median, Eq1_rsv_median, 
           Eq2_hib_median, Eq2_pcv_median, Eq2_rsv_median,
           Eq3_hib_median, Eq3_pcv_median, Eq3_rsv_median,
           Eq4_hib_median, Eq4_pcv_median, Eq4_rsv_median,
           Eq5_hib_median, Eq5_pcv_median, Eq5_rsv_median, 
           Eq6_hib_median, Eq6_pcv_median, Eq6_rsv_median,
           Eq7_hib_median, Eq7_pcv_median, Eq7_rsv_median, 
           Eq8_hib_median, Eq8_pcv_median, Eq8_rsv_median)
  
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
    select(Eq1_hib_median, Eq1_pcv_median, Eq1_rsv_median, 
           Eq2_hib_median, Eq2_pcv_median, Eq2_rsv_median,
           Eq3_hib_median, Eq3_pcv_median, Eq3_rsv_median,
           Eq4_hib_median, Eq4_pcv_median, Eq4_rsv_median,
           Eq5_hib_median, Eq5_pcv_median, Eq5_rsv_median, 
           Eq6_hib_median, Eq6_pcv_median, Eq6_rsv_median,
           Eq7_hib_median, Eq7_pcv_median, Eq7_rsv_median, 
           Eq8_hib_median, Eq8_pcv_median, Eq8_rsv_median)
  
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
    select(Eq1_hib_median, Eq1_pcv_median, Eq1_rsv_median, 
           Eq2_hib_median, Eq2_pcv_median, Eq2_rsv_median,
           Eq3_hib_median, Eq3_pcv_median, Eq3_rsv_median,
           Eq4_hib_median, Eq4_pcv_median, Eq4_rsv_median,
           Eq5_hib_median, Eq5_pcv_median, Eq5_rsv_median, 
           Eq6_hib_median, Eq6_pcv_median, Eq6_rsv_median,
           Eq7_hib_median, Eq7_pcv_median, Eq7_rsv_median, 
           Eq8_hib_median, Eq8_pcv_median, Eq8_rsv_median)
  
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

Eq1_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq1_hib_median", "Eq1_hib", "DHS incidence of hib specific treatment for diarrhea by country")
Eq1_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq1_pcv_median", "Eq1_pcv", "DHS incidence of pcv specific treatment for diarrhea by country")
Eq1_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq1_rsv_median", "Eq1_rsv", "DHS incidence of rsv specific treatment for diarrhea by country")

##############
# EQUATION 2 #
##############

Eq2_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq2_hib_median", "Eq2_hib", "DHS incidence of hib specific treatment for diarrhea with a hib vaccine by country")
Eq2_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq2_pcv_median", "Eq2_pcv", "DHS incidence of pcv specific treatment for diarrhea with a pcv vaccine by country")
Eq2_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq2_rsv_median", "Eq2_rsv", "DHS incidence of rsv specific treatment for diarrhea with an rsv vaccine by country")

##############
# EQUATION 3 #
##############

Eq3_hib_DHS<- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_ARI, "Eq3_hib_median", "Eq3_hib", "DHS absolute reductions in antibiotic use with a hib vaccine by country")
Eq3_pcv_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_ARI, "Eq3_pcv_median", "Eq3_pcv", "DHS absolute reductions in antibiotic use with a pcv vaccine by country")
Eq3_rsv_DHS <- df_function_Eq3_SHIG_RSV(all_path_all_eq_med_CI_DHS_ARI, "Eq3_rsv_median", "Eq3_rsv", "DHS absolute reductions in antibiotic use with an rsv vaccine by country")

##############
# EQUATION 4 #
##############

Eq4_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq4_hib_median", "Eq4_hib", "DHS percent reductions in antibiotic use with a hib vaccine by country")
Eq4_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq4_pcv_median", "Eq4_pcv", "DHS percent reductions in antibiotic use with a pcv vaccine by country")
Eq4_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq4_rsv_median", "Eq4_rsv", "DHS percent reductions in antibiotic use with an rsv vaccine by country")

##############
# EQUATION 5 #  
##############

Eq5_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq5_hib_median", "Eq5_hib", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of hib by country")
Eq5_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq5_pcv_median", "Eq5_pcv", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of pcv by country")
Eq5_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq5_rsv_median", "Eq5_rsv", "DHS incidence of bystander pathogen exposures to antibiotics due to treatment of rsv by country")

##############
# EQUATION 6 # 
##############

Eq6_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq6_hib_median", "Eq6_hib", "DHS incidence of bystander pathogen exposures to antibiotics with a hib vaccine by country")
Eq6_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq6_pcv_median", "Eq6_pcv", "DHS incidence of bystander pathogen exposures to antibiotics with a pcv vaccine by country")
Eq6_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq6_rsv_median", "Eq6_rsv", "DHS incidence of bystander pathogen exposures to antibiotics with an rsv vaccine by country")

##############
# EQUATION 7 #
##############

Eq7_hib_DHS<- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_ARI, "Eq7_hib_median", "Eq7_hib", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a hib vaccine by country")
Eq7_pcv_DHS <- df_function_Eq3_7(all_path_all_eq_med_CI_DHS_ARI, "Eq7_pcv_median", "Eq7_pcv", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a pcv vaccine by country")
Eq7_rsv_DHS <- df_function_Eq7_SHIG_RSV(all_path_all_eq_med_CI_DHS_ARI, "Eq7_rsv_median", "Eq7_rsv", "DHS absolute reductions of bystander pathogen exposures to antibiotics with a rsv vaccine by country")

##############
# EQUATION 8 #
##############

Eq8_hib_DHS<- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq8_hib_median", "Eq8_hib", "DHS percent reductions of bystander pathogen exposures to antibiotics with a hib vaccine by country")
Eq8_pcv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq8_pcv_median", "Eq8_pcv", "DHS percent reductions of bystander pathogen exposures to antibiotics with a pcv vaccine by country")
Eq8_rsv_DHS <- df_function(all_path_all_eq_med_CI_DHS_ARI, "Eq8_rsv_median", "Eq8_rsv", "DHS percent reductions of bystander pathogen exposures to antibiotics with a rsv vaccine by country")

###############
# Saving Maps #
###############

# EQUATION 1
ggsave(file="Scripts - Aim 3/Maps/Eq1_hib_Map_DHS.png",Eq1_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_pcv_Map_DHS.png",Eq1_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq1_rsv_Map_DHS.png",Eq1_rsv_DHS,width=22,height=10,dpi=300) 

# EQUATION 2
ggsave(file="Scripts - Aim 3/Maps/Eq2_hib_Map_DHS.png",Eq2_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_pcv_Map_DHS.png",Eq2_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq2_rsv_Map_DHS.png",Eq2_rsv_DHS,width=22,height=10,dpi=300) 

# EQUATION 3
ggsave(file="Scripts - Aim 3/Maps/Eq3_hib_Map_DHS_24Nov2025.png",Eq3_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_pcv_Map_DHS_24Nov2025.png",Eq3_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_rsv_Map_DHS_24Nov2025.png",Eq3_rsv_DHS,width=22,height=10,dpi=300)

ggsave(file="Scripts - Aim 3/Maps/Eq3_hib_Map_DHS_24Nov2025.pdf",Eq3_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_pcv_Map_DHS_24Nov2025.pdf",Eq3_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq3_rsv_Map_DHS_24Nov2025.pdf",Eq3_rsv_DHS,width=22,height=10,dpi=300)

# EQUATION 4
ggsave(file="Scripts - Aim 3/Maps/Eq4_hib_Map_DHS.png",Eq4_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_pcv_Map_DHS.png",Eq4_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq4_rsv_Map_DHS.png",Eq4_rsv_DHS,width=22,height=10,dpi=300) 

# EQUATION 5
ggsave(file="Scripts - Aim 3/Maps/Eq5_hib_Map_DHS.png",Eq5_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_pcv_Map_DHS.png",Eq5_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq5_rsv_Map_DHS.png",Eq5_rsv_DHS,width=22,height=10,dpi=300)

# EQUATION 6
ggsave(file="Scripts - Aim 3/Maps/Eq6_hib_Map_DHS.png",Eq6_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_pcv_Map_DHS.png",Eq6_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq6_rsv_Map_DHS.png",Eq6_rsv_DHS,width=22,height=10,dpi=300) 

# EQUATION 7
ggsave(file="Scripts - Aim 3/Maps/Eq7_hib_Map_DHS_24Nov2025.png",Eq7_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_pcv_Map_DHS_24Nov2025.png",Eq7_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_rsv_Map_DHS_24Nov2025.png",Eq7_rsv_DHS,width=22,height=10,dpi=300) 

ggsave(file="Scripts - Aim 3/Maps/Eq7_hib_Map_DHS_24Nov2025.pdf",Eq7_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_pcv_Map_DHS_24Nov2025.pdf",Eq7_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq7_rsv_Map_DHS_24Nov2025.pdf",Eq7_rsv_DHS,width=22,height=10,dpi=300) 

# EQUATION 8
ggsave(file="Scripts - Aim 3/Maps/Eq8_hib_Map_DHS.png",Eq8_hib_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_pcv_Map_DHS.png",Eq8_pcv_DHS,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Maps/Eq8_rsv_Map_DHS.png",Eq8_rsv_DHS,width=22,height=10,dpi=300) 


saveRDS(Eq3_rsv_DHS, "Scripts - Aim 3/Output/Eq3_rsv_DHS.rds")
saveRDS(Eq7_rsv_DHS, "Scripts - Aim 3/Output/Eq7_rsv_DHS.rds")

