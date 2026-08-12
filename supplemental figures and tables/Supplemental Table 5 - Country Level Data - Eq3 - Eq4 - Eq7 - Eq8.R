# this script will create the supplemental table which lists the data from the world maps
# for all the diarrhea pathogens and all the respiratory pathogens for Eq 3,4,7,8

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(stringr)

# diarrhea 
load("Scripts - Aim 3/Output/All_Diarrhea_Pathogens_AllEqs_Med_CI_DHS_24Nov2025.Rdata")

# ari 
load("Scripts - Aim 3/Output/All_ARI_Pathogens_AllEqs_Med_CI_DHS_24Nov2025.Rdata")

#combined - REMOVED THIS - YOU CANNOT ADD INCIDENCE TOGETHER
# load("Scripts - Aim 3/Output/CombinedImpact_Eq3_Eq4_Eq7_Eq8_hib_pcv_rota_Med_CI_DHS_16Oct2025.Rdata")

#country list
load("Scripts - Aim 3/country_region_list.Rdata")

########################
# clean country_region #
########################

# remove the 6 countries that we do not need: 
# Democratic Peoples Republic of Korea, Palau, Argentina, Occupied Palestinian territory, Syria, duplicate of Georgia
country_region_list <- country_region_list %>%
  filter(DisplayString!="Democratic People's Republic of Korea" & DisplayString!="Palau" & DisplayString!="Argentina" & 
           DisplayString!="occupied Palestinian territory, including east Jerusalem" & DisplayString!="Syrian Arab Republic") %>%
  select(DisplayString, WHO_REGION) %>%
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


###################################
# combine the 4 datasets together #
###################################

# we need to add the country data as a column
country_diar <- cbind(country_region_list, all_path_all_eq_med_CI_DHS_DIAR) 
country_diar_ari <- cbind(country_diar, all_path_all_eq_med_CI_DHS_ARI)

###############################
# keep only equations 3,4,7,8 #
###############################

country_diar_ari <- country_diar_ari %>%
  select(-matches("Eq1|Eq2|Eq5|Eq6")) 

########################################
# get the data into the desired format #
########################################

# set order of regions
region_levels <- c("Africa", "Americas", "Eastern Mediterranean", "Europe", "South-East Asia", "Western Pacific")  # whatever order you prefer

country_diar_ari <- country_diar_ari %>%
  mutate(WHO_REGION = factor(WHO_REGION, levels = region_levels))

# pivot data long 
country_diar_ari_long <- country_diar_ari %>%
  pivot_longer(
    cols = -c(DisplayString, WHO_REGION),
    names_to = "colname",
    values_to = "value")

# Pathogen vector (your exact order)
pathogen_vector <- c(
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv")

# Stat vector repeats: median, lb, ub
stat_vector <- rep(c("median", "lb", "ub"), length.out = length(pathogen_vector))

# Repeat for every country
country_diar_ari_long <- country_diar_ari_long %>%
  group_by(DisplayString) %>%
  mutate(
    pathogen = pathogen_vector[1:n()],
    stat     = stat_vector[1:n()],
    equation = str_extract(colname, "^Eq\\d+")) %>%  # keep equation extraction from colname
  ungroup()


# pivot wider - one row per country/pathogen/equation
country_diar_ari_wide <- country_diar_ari_long %>%
  select(DisplayString, WHO_REGION, equation, pathogen, stat, value) %>%
  pivot_wider(names_from = stat, values_from = value)


# round to 2 decimal places - %.2f ensures that if numbers rounded to 0 are closer to "negative 0" that they will appear as such
country_diar_ari_wide <- country_diar_ari_wide %>%
  mutate(
    # copy original lb and ub
    lb_orig = lb,
    ub_orig = ub,
    # negate median, and swap+negate lb/ub for Eq4 and Eq8
    median = ifelse(equation %in% c("Eq4","Eq8"), -median, median),
    lb     = ifelse(equation %in% c("Eq4","Eq8"), -ub_orig, lb),
    ub     = ifelse(equation %in% c("Eq4","Eq8"), -lb_orig, ub)) %>%
  select(-lb_orig, -ub_orig) %>%
  mutate(
    est_CI = paste0(
      sprintf("%.2f", median),
      " (", sprintf("%.2f", lb), ", ", sprintf("%.2f", ub), ")" ))

# Pivot wider to have one column per equation
country_diar_ari_wide_final <- country_diar_ari_wide %>%
  select(DisplayString, WHO_REGION, pathogen, equation, est_CI) %>%
  pivot_wider(names_from = equation, values_from = est_CI)


#####################################
# making a table for the manuscript #
#####################################

country_diar_ari_wide_final1 <- lapply(region_levels, function(region) {
  
  # WHO_REGION header row
  region_header <- data.frame(
    DisplayString = region,
    WHO_REGION = "",
    pathogen = "",
    Eq3 = "", Eq4 = "", Eq7 = "", Eq8 = "",
    stringsAsFactors = FALSE
  )
  
  
  # Countries in this region, alphabetically
  countries <- sort(unique(country_diar_ari_wide_final$DisplayString[
    country_diar_ari_wide_final$WHO_REGION == region
  ]))
  
  # For each country, create header + pathogen rows
  country_rows <- do.call(rbind, lapply(countries, function(ctry) {
    # country header
    header <- data.frame(
      DisplayString = ctry,
      WHO_REGION = "",
      pathogen = "",
      Eq3 = "", Eq4 = "", Eq7 = "", Eq8 = "",
      stringsAsFactors = FALSE
    )
    
    # pathogen rows
    data <- country_diar_ari_wide_final %>%
      filter(DisplayString == ctry) %>%
      mutate(DisplayString = "", WHO_REGION = "")  # blank for manuscript
    
    bind_rows(header, data)
  }))
  
  # Combine region header + country rows
  bind_rows(region_header, country_rows)
  
}) %>% bind_rows()


# Manuscript-friendly pathogen names
pathogen_map <- c(
  shigella = "Shigella",
  campy = "Campylobacter",
  ETEC = "ETEC",
  noro = "Norovirus",
  rota = "Rotavirus",
  adeno = "Adenovirus 40/41",
  hib = "Hib",
  pcv = "PCV",
  rsv = "RSV")

# Apply renaming only to pathogen rows (DisplayString != country)
country_diar_ari_wide_final1 <- country_diar_ari_wide_final1 %>%
  mutate(
    pathogen = ifelse(
      pathogen != "",
      pathogen_map[pathogen],
      pathogen)) %>%
  select(-WHO_REGION)

# remove Eq4 because it's identical to Eq8 and rename column Eq4 and Eq8
country_diar_ari_wide_final2 <- country_diar_ari_wide_final1 %>%
  select(-Eq4) %>%
  rename(Eq4_Eq8 = Eq8)

# save the table
write.csv(country_diar_ari_wide_final2, file="Scripts - Aim 3/Tables/Supp Table 4 - Diar_ARI_Combined_Eq3_Eq4_Eq7_Eq8.csv")


