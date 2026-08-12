# This script is for Aim 3 analyses
# The below script is modified from the comparison paper, which only included the 5 sites that hadn't introduced rotavirus vaccine estimates
# This will be for all vaccines EXCEPT: rotavirus, norovirus + rotavirus, adenovirus + norovirus + rotavirus
# INFLATION FACTOR NEEDS TO BE WHAT WE USED FOR SHIGELLA PAPER SINCE THAT IS FOR ALL CHILDREN FROM 8 SITES
# We are only interestined in 1 outcome: The outcome is pathogen specific diarrhea episodes for the pathogen included in the vaccine 
#        (i.e., we do not care about abx treated shigella diarrhea episodes for a norovirus vaccine). 
#        Essentially, just Eq #2 that we used in the Shigella paper. 
#       We are going to scrap combination vaccines because combination vaccines with rotavirus are complicated and limited.

## Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(tidyr)
library(splitstackshape)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(lme4)
library(nloptr)
library(CMAverse)
library(doRNG)
library(Rlab)

# load abx_all_courses1 - needed for all analyses - converted NAs for pathogens to 0
load("abx output/abx_all_courses1_Aim3_VE.Rda")

# set seed (random.org from 1-1000)
set.seed(229) # new set seed as of 7Mar24

# NAs for score need to be converted to 0
abx_all_courses1$score[is.na(abx_all_courses1$score)] <- 0

## load list function for returning list from for each loop
source('list.R')

R = 1000

## run parallel with foreach/dopar
# This allows the loop to run faster
#nocores<-detectCores()
nocores<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))- 1
cl<-makeCluster(nocores-1)
registerDoParallel(cl)

## create blank log file
writeLines(c(""), paste0("amr_log.txt"))

lst <- foreach(s = 1:R,.packages=c("lme4","nloptr","dplyr", "CMAverse")) %dorng% {
  
  source('list.R')
  
  ## monitor progress
  # SB: Kicks out a log file in the AMR box folder
  cat(paste("\n","Bootstrap Iteration",s,"\n"),file=paste0("amr_log.txt"),append=TRUE)
  print(s)
  
  
  # All variables and equations below need to be based off 50000 simulated dataframe
  
  # 1. Sample 50,000 with replacement
  # grabbing the unique pids from siteagedata and creating a pids list (1595 kids)
  pids <- unique(abx_all_courses1$pid)
  # grabbing a pid from pids then putting it back in the drawing to be picked again and creating a pidsstar list (sample pids, same length with replacement)
  pidsstar2 <- sample(pids, length(pids), replace = T) # bootstrap step 
  pidsstar <- sample(pidsstar2, 50000, replace = T) # simulation step 
  # lapply performing the function(pid) on every item in pidsstar - identify rows in old dataset that correspond with the new order of data in pidsstar
  # do.call"c" concatenates them into a vector
  rows <- do.call("c", lapply(pidsstar, function(pid) which(abx_all_courses1$pid == pid)))
  # rowsPerPid says how many rows per pid in pidsstar but pids are in the order in which they were selected 
  rowsPerPid <- do.call("c", lapply(pidsstar, function(pid) sum(abx_all_courses1$pid == pid)))

  boot <- abx_all_courses1[rows, ]
  # adding new var (newpid) to boot df - assigning a new pid per child to be able to group them 
  boot$newpid <- rep(1:50000, times = rowsPerPid)
  
  # 61  = 2 months
  # 75  = 2 months + 14 days
  # 122 = 4 months 
  # 136 = 4 months + 14 days
  # 183 = 6 months
  # 197 = 6 months + 14 days
  # 274 = 9 months
  # 288 = 9 months + 14 days

  
  # function starts here
  scenario <- function(path_afe, afe, age1, age2, dose, efficacy, mult, est2) { 
    # dose = 1 for 1 dose and 2 for 2 dose 
    # rbeta(1, 6, 4)  #60% efficacy 
    # rbeta(1, 6, 1.5) #80% efficacy 
    # mult will be 2/3 for the 60% VE non-severe episodes, and 0.75 for the 80% VE non-severe episodes
    
    #------------------------#
    # Monte Carlo Simulation # 
    #------------------------#
    
    # 2. Calculate probability 
    boot <- boot %>% mutate(prevented=case_when(
      !!rlang::sym(path_afe) > afe & agedays >=age1 & agedays <=age2 & score>6 ~ rbinom(n(),1,efficacy/dose),  # severe
      !!rlang::sym(path_afe) > afe & agedays >=age1 & agedays <=age2 & score<=6 ~ rbinom(n(),1,efficacy*mult/dose),  # non-severe
      
      !!rlang::sym(path_afe) > afe & agedays > age2 & score>6~ rbinom(n(),1,efficacy), # severe
      !!rlang::sym(path_afe) > afe & agedays > age2 & score<=6~ rbinom(n(),1,efficacy*mult),# non-severe
      
      TRUE~ 0L)) 
    
    #-----------------------------------------------------------#
    # 2. Antibiotic treated pathogen specific diarrhea episodes
    #-----------------------------------------------------------#
    # (((N antibiotic treated diarrhea episodes attributed to Shigella)*(total number of diarrhea associated courses)/(Total number of linked diarrhea associated courses))/
    # (total child years of follow up))*100
    
    # The est2 variables are now created in prep 1 df
    # create variable where shigella_eiec_afe >0.5 (should go in separate script since there will be more cases of this)
    # boot <- boot %>% mutate(est2_shig= ifelse(shigella_eiec_afe>0.5 & diarassociated==1, 1, 0)) # SUBSETTING TO DIARASSOCIATED DID NOT ADD ANYTHING
    
    # create equation for: total number of diarrhea associated courses/Total number of linked diarrhea associated courses = 3029/2199
    diarnew_diar <- sum(boot$diarassociatednew==1)/sum(boot$diarassociated==1)
    
    # equation 
    # eq2 <- ((colSums(dplyr::select(filter(boot,est2_shig==1 & prevented==0), starts_with("est2_shig")))*diarnew_diar)/107523.5)*100
    
    # equation 
    eq2 <- ((colSums(dplyr::select(filter(boot, !!rlang::sym(est2)==1 & prevented==0), !!rlang::sym(est2)))*diarnew_diar)/107523.5)*100
    
    
    
    
  #  eq2 <- ((colSums(dplyr::select(filter(boot, est2_shig==1), est2_shig))*diarnew_diar)/107523.5)*100
    
    return(c(eq2))
    
  } # func scenario bracket
  
  ############################################################
  #                     Scenarios                            #
  #   path_afe, afe, age1, age2, dose, efficacy, mult, est2  #
  ############################################################

  # need the following:
  # Direct protection only
  
  # Shigella, Campy, ETEC - 2 doses at 6 and 9 months, 60/40% efficacy
  # Norovirus - 2 doses at 2 and 4 months - 60/40% efficacy
  
  
  # 75  = 2 months + 14 days
  # 136 = 4 months + 14 days
  # 197 = 6 months + 14 days
  # 288 = 9 months + 14 days
  
  ##############
  # NO VACCINE #
  ##############
  scenario_0_shig <-  as.data.frame(scenario("shigella_eiec_afe",0,0,0,0,0,0, "est2_shig"))
  scenario_0_campy <-  as.data.frame(scenario("campylobacter_jejuni_coli_afe",0,0,0,0,0,0, "est2_campy")) 
  scenario_0_ETEC <-  as.data.frame(scenario("ETEC_afe",0,0,0,0,0,0, "est2_etec"))
  scenario_0_noro <-  as.data.frame(scenario("norovirus_gii_afe",0,0,0,0,0,0, "est2_noro")) 
  
  scenario_0_adeno <-  as.data.frame(scenario("adenovirus_40_41_afe",0.5,75,136,2,0,2/3, "est2_adeno"))

  ##########################
  # DIRECT PROTECTION ONLY #
  ##########################
  
  # pathogen, pathogen1, pathogen2, afe, age1, age2, dose, efficacy, mult
  # SIMULATED scenario - AFe > 0.5, 30% VE 14 days after 1st dose, 60% VE 14 days after 2nd dose
  
  # rbeta(1, 6, 4)  #60% efficacy 
  # mult 2/3
  
  scenario_shig <- as.data.frame(cbind(data.frame(pathogen="shigella",scenario_est=scenario("shigella_eiec_afe", 0.5,197,288,2,rbeta(n(), 6, 4),2/3, "est2_shig")),scenario_0_shig)) # Schedule: 6 & 9 months 
  rownames(scenario_shig) <- c("Eq2") 
  colnames(scenario_shig)[3] <- "scenario_0"
  
  scenario_campy <- as.data.frame(cbind(data.frame(pathogen="campy",scenario_est=scenario("campylobacter_jejuni_coli_afe", 0.5,197,288,2,rbeta(n(), 6, 4),2/3, "est2_campy")),scenario_0_campy)) # Schedule: 6 & 9 months 
  rownames(scenario_campy) <- c("Eq2") 
  colnames(scenario_campy)[3] <- "scenario_0"
  
  scenario_ETEC <- as.data.frame(cbind(data.frame(pathogen="ETEC",scenario_est=scenario("ETEC_afe",0.5,197,288,2,rbeta(n(), 6, 4),2/3, "est2_etec")),scenario_0_ETEC)) # Schedule: 6 & 9 months 
  rownames(scenario_ETEC) <- c("Eq2") 
  colnames(scenario_ETEC)[3] <- "scenario_0"
  
  scenario_noro <- as.data.frame(cbind(data.frame(pathogen="noro",scenario_est=scenario("norovirus_gii_afe", 0.5,75,136,2,rbeta(n(), 6, 4),2/3, "est2_noro")),scenario_0_noro)) # Schedule: 2 & 4 months 
  rownames(scenario_noro) <- c("Eq2") 
  colnames(scenario_noro)[3] <- "scenario_0"
  

  scenario_adeno <- as.data.frame(cbind(data.frame(pathogen="adeno",scenario_est=scenario("adenovirus_40_41_afe",0.5,75,136,2,rbeta(n(), 6, 4),2/3, "est2_adeno")),scenario_0_adeno)) # Schedule: 2 & 4 months 
  rownames(scenario_adeno) <- c("Eq2") 
  colnames(scenario_adeno)[3] <- "scenario_0"

  
  ####################
  # COMBINING OUTPUT #
  ####################
  
  # combine into one df
  all_scenarios <- list(scenario_shig, scenario_campy, scenario_ETEC, scenario_noro, scenario_adeno)  # direct 
  
  ########################################################
  # Relative and absolute differences and percent change #
  ########################################################
  
  #######################
  # Relative difference #
  #######################
  
  # function for relative difference 
  rel_diff <- function(path_scenario) {
    model <- path_scenario%>%
      mutate(rel_diff=scenario_est/scenario_0)
    return(model)
  }
  
  rel_diff_output <- lapply(all_scenarios, rel_diff) # compares path specific scenario 0 to scenario 1
  
  #######################
  # Absolute difference #
  #######################
  
  # function for absolute difference 
  ab_diff <- function(path_scenario) {
    model <- path_scenario %>%
      mutate(abs_diff= abs(scenario_est - scenario_0))
    return(model)
  }
  
  ab_diff_output <- lapply(rel_diff_output, ab_diff) # add abs diff to relative diff
  
  ##################
  # Percent Change #
  ##################
  
  # function for percent change
  perc_change <- function(path_scenario) {
    model <- path_scenario %>%
      mutate(perc_change= ((scenario_est - scenario_0)/scenario_0)*100) %>%
      tibble::rownames_to_column(var="Equation")
    
    return(model)
  }
  
  perc_change_output <- lapply(ab_diff_output, perc_change)
  
  
  ################################
  # Bind all dataframes together #
  ################################
  
  all_dfs <- do.call("rbind", perc_change_output)
  
  all_dfs <- all_dfs %>% tidyr::pivot_wider(names_from = Equation, values_from=c(scenario_est:perc_change))
  
  ###############################################
  # Each equation needs to be its own dataframe #
  ###############################################
  
  
  shig <- all_dfs[1,]
  campy <- all_dfs[2,]
  ETEC <- all_dfs[3,]
  noro <- all_dfs[4,]
  
  adeno <- all_dfs[5,]
  
  return(c(shig, campy, ETEC, noro, adeno))
  
}


#########################
#       Formatting      #
#########################

# saving the bootstrapped CIs as a df
aim3_boot_df <- as.data.frame(do.call("rbind", lst))

# WORK ON FORMATTING BELOW. 

# data frame columns are lists, not vectors - need to transform back to normal format
aim3_boot_df <- sapply(aim3_boot_df, unlist)
aim3_boot_df <- as.data.frame(aim3_boot_df)

write.csv(aim3_boot_df, "abx output/aim3_Combo_Vaccines_NoRota_boot_df_7Mar.csv")
save(aim3_boot_df, file=paste0("abx output/aim3_Combo_Vaccines_NoRota_boot_df_7Mar.Rda"))

# Uncomment if needing to load partb_boot_df to work on formatting 
# aim3_boot_df <- read.csv(file = 'abx output/aim3_Combo_Vaccines_NoRota_boot_df_22Feb.csv')

# remove pathogen columns
aim3_boot_df <- aim3_boot_df[,-c(1,7,13,19,25)]
aim3_boot_df <- sapply(aim3_boot_df, as.numeric)
aim3_boot_df <- as.data.frame(aim3_boot_df)

# grabbing the 2.5% and 97.5%
aim3_boot_CI <- as.data.frame(t(apply(aim3_boot_df, 2, quantile, c(0.025, 0.975), na.rm=TRUE)))
# grabbing the median
aim3_boot_med <- as.data.frame(apply(aim3_boot_df, 2, quantile, 0.5, na.rm=TRUE))
names(aim3_boot_med)[1] <- 'median'

# bind together
aim3_boot_med_CIs <- cbind(aim3_boot_med, aim3_boot_CI)
# remove extra row
#aim3_boot_med_CIs <- aim3_boot_med_CIs[-1,]

# save estimates with CIs
write.csv(aim3_boot_med_CIs, "abx output/aim3_Combo_Vaccines_NoRota_boot_med_CIs_7Mar.csv")
save(aim3_boot_med_CIs, file=paste0("abx output/aim3_Combo_Vaccines_NoRota_boot_med_CIs_7Mar.Rda"))

# round to 2 decimal places
aim3_boot_med_CIs <- aim3_boot_med_CIs %>% mutate_if(is.numeric, ~round(., 3))

# combine estimates and CIs into one column 
aim3_boot_med_CIs$CI_bind <- paste0("(", paste(aim3_boot_med_CIs$"2.5%", aim3_boot_med_CIs$"97.5%", sep= ", "),paste0(")"))
aim3_boot_med_CIs$est_CI_bind <- paste(aim3_boot_med_CIs$"median", aim3_boot_med_CIs$CI_bind, sep= " ")

# keep only combined column
aim3_boot_med_CIs <- aim3_boot_med_CIs %>% select(est_CI_bind)
# make row names a column
aim3_boot_med_CIs <- aim3_boot_med_CIs %>% tibble::rownames_to_column()

# create equation column
aim3_boot_med_CIs$Eq <- c("Eq2")

# create pathogen column
aim3_boot_med_CIs$pathogen <- c("Shigella","Shigella","Shigella","Shigella","Shigella", "Campy", "Campy","Campy","Campy","Campy","ETEC","ETEC","ETEC","ETEC","ETEC", "Norovirus","Norovirus","Norovirus","Norovirus","Norovirus","Adenovirus","Adenovirus","Adenovirus","Adenovirus","Adenovirus")                                                                                                         

all_eq <- aim3_boot_med_CIs
# remove characters starting at the second underscore
cut<- as.data.frame(stringr::str_extract(all_eq$rowname, "^[^_]*_[^_]*"))
colnames(cut) <- 'scenario'
all_eq <- cbind(cut, all_eq)
all_eq <- all_eq[,-2]
# reorder columns
all_eq <- all_eq[, c(3,4,1,2)]

# format to more closely resemble how we formatted for shigella
all_eq <- all_eq %>% tidyr::pivot_wider(names_from = Eq, values_from=est_CI_bind)

# separate out each pathogen into its own dataframe
#shig <- as.data.frame(aim3_boot_med_CIs[1:40,])
#campy <- as.data.frame(aim3_boot_med_CIs[41:80,])
#etec <- as.data.frame(aim3_boot_med_CIs[81:120,])
#noro <- as.data.frame(aim3_boot_med_CIs[121:160,])
#rota <- as.data.frame(aim3_boot_med_CIs[161:200,])

# bind together
#aim3_med_CIs <- cbind(shig, campy, etec, noro, rota)
# keep only wanted columns
#aim3_med_CIs <- aim3_med_CIs[,c(1,2,4,6,8,10)]

##################
# Save dataframe #
##################

write.csv(all_eq, "abx output/aim3_Combo_Vaccines_NoRota_med_CIs_FINAL_7Mar.csv")

