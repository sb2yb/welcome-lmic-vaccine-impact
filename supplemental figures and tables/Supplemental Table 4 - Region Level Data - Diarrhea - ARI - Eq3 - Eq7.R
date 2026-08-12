# this script will create the REGIONAL supplemental table for Eq 3,7 (DIARRHEA & ARI)
# this table will be pathogen & REGION specific 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# diarrhea 
load("Scripts - Aim 3/Output/Diarrhea_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")

# ARI
load("Scripts - Aim 3/Output/ARI_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")

# THIS IS NOT REGIONAL DATA THAT IS PRODUCED

# diarrhea 
# make rownames a column - regions are now a column 
diarrhea_all_eq_med_CI_weighted <- rownames_to_column(diarrhea_all_eq_med_CI_weighted, var = "region")

# keep only equation3
eq3_7_diarrhea <- diarrhea_all_eq_med_CI_weighted %>%
  select(matches("region|Eq3|Eq7"))

# ari 
# make rownames a column - regions are now a column 
ari_all_eq_med_CI_weighted <- rownames_to_column(ari_all_eq_med_CI_weighted, var = "region")

# keep only equation3
eq3_7_ari <- ari_all_eq_med_CI_weighted %>%
  select(matches("Eq3|Eq7"))

# bind eq3_diarrhea and eq3_ari 
eq3_7_combined <- cbind(eq3_7_diarrhea, eq3_7_ari)


########################################
# get the data into the desired format #
########################################

# pivot data long 
eq3_7_combined_long <- eq3_7_combined %>%
  pivot_longer(
    cols = -c(region),
    names_to = "colname",
    values_to = "value")

# Pathogen vector (your exact order)
pathogen_vector <- c(
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "shigella","shigella","shigella","campy","campy","campy",
  "ETEC","ETEC","ETEC","noro","noro","noro","rota","rota","rota","adeno","adeno","adeno",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv",
  "hib","hib","hib","pcv","pcv","pcv","rsv","rsv","rsv")

# Stat vector repeats: median, lb, ub
stat_vector <- rep(c("median", "lb", "ub"), length.out = length(pathogen_vector))

# Repeat for every country
eq3_7_combined_long <- eq3_7_combined_long %>%
  group_by(region) %>%
  mutate(
    pathogen = pathogen_vector[1:n()],
    stat     = stat_vector[1:n()],
    equation = str_extract(colname, "^Eq\\d+")) %>%  # keep equation extraction from colname
  ungroup()


# pivot wider - one row per region/pathogen/equation
eq3_7_combined_wide <- eq3_7_combined_long %>%
  select(region, equation, pathogen, stat, value) %>%
  pivot_wider(names_from = stat, values_from = value)


# round to 2 decimal places - %.2f ensures that if numbers rounded to 0 are closer to "negative 0" that they will appear as such
eq3_7_combined_wide <- eq3_7_combined_wide %>%
  mutate(
    est_CI = paste0(
      sprintf("%.2f", median),
      " (", sprintf("%.2f", lb), ", ", sprintf("%.2f", ub), ")"))

# Pivot wider to have one column per equation
eq3_7_combined_wide_final <- eq3_7_combined_wide %>%
  select(region, pathogen, equation, est_CI) %>%
  pivot_wider(names_from = equation, values_from = est_CI)


#####################################
# making a table for the manuscript #
#####################################

eq3_7_combined_wide_final1 <- lapply(unique(eq3_7_combined_wide_final$region), function(region) {
  
  # WHO_REGION header row
  region_header <- data.frame(
    region = region,
    pathogen = "",
    Eq3 = "", Eq7 = "",
    stringsAsFactors = FALSE)
  
  # Pathogen rows for that region from your data
  pathogen_rows <- eq3_7_combined_wide_final %>%
    filter(region == !!region) %>%
    select(region, pathogen, Eq3, Eq7) %>%
    mutate(region = "")  # blank out for formatting under header
  
  # Combine region header and pathogen rows
  bind_rows(region_header, pathogen_rows)}) %>% bind_rows()


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
eq3_7_combined_wide_final1 <- eq3_7_combined_wide_final1 %>%
  mutate(
    pathogen = ifelse(
      pathogen != "",
      pathogen_map[pathogen],
      pathogen))


# save the table 
write.csv(eq3_7_combined_wide_final1, file="Scripts - Aim 3/Tables/Supp Table 3 - Weighted_Regional_Diar_ARI_Eq3_Eq7.csv")


