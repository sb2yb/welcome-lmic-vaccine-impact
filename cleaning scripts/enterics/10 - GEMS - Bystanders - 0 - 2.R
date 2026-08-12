# This script is to calculate the average number of asymptomatic pathogens in MAL-ED by site 

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

# load the MALED tac data
load("Scripts - Aim 3/James data/GEMS_data_0_2.Rda")

# filter and order stools (stooltype monthlies, remove 60 month stools, excluded samples with bad ms2
# select only the needed variables
GEMS_data1 <- GEMS_data_0_2 %>%
  select(Case.ID, country, aeromonas, aEPEC, campylobacter_pan, EAEC, LT_ETEC, ST_ETEC, h_pylori, salmonella, STEC, shigella_eiec, TEPEC, v_cholerae)

# create dichotomous variables 
GEMS_data2 <- GEMS_data1 %>%
  mutate(aeromonas_sub = ifelse(aeromonas<35, 1, 0)) %>%
  mutate(aEPEC_sub = ifelse(aEPEC<35, 1, 0)) %>%
  mutate(campylobacter_pan_sub = ifelse(campylobacter_pan<35, 1, 0)) %>%
  mutate(EAEC_sub = ifelse(EAEC<35, 1, 0)) %>%
  mutate(h_pylori_sub = ifelse(h_pylori<35, 1, 0)) %>%
  mutate(salmonella_sub = ifelse(salmonella<35, 1, 0)) %>%
  mutate(STEC_sub = ifelse(STEC<35, 1, 0)) %>%
  mutate(shigella_eiec_sub = ifelse(shigella_eiec<35, 1, 0)) %>%
  mutate(tEPEC_sub = ifelse(TEPEC<35, 1, 0)) %>%
  mutate(v_cholerae_sub = ifelse(v_cholerae<35, 1, 0))

# create ETEC variable 
GEMS_data2 <- GEMS_data2 %>%
  mutate(ETEC_sub = ifelse(LT_ETEC<35 | ST_ETEC<35,1,0))

# set seed (random.org from 1-1000)
 set.seed(151) # seed set 14 March 24


#############################
#         BOOTSTRAPS        #
#############################

## load list function for returning list from for each loop
source('list.R')

R = 1000

## run parallel with foreach/dopar
# This allows the loop to run faster
#nocores<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))- 1
nocores<-detectCores()
cl<-makeCluster(nocores-1)
registerDoParallel(cl)

## create blank log file
writeLines(c(""), paste0("med_log.txt"))

lst <- foreach(s = 1:R,.packages=c("lme4","nloptr","dplyr", "CMAverse")) %dorng% {
  
  source('list.R')
  
  ## monitor progress
  # SB: Kicks out a log file in the AMR box folder
  cat(paste("\n","Bootstrap Iteration",s,"\n"),file=paste0("med_log.txt"),append=TRUE)
  print(s)
  
  
  # grabbing the unique pids from siteagedata and creating a pids list (1260 kids)
  pids <- unique(GEMS_data2$Case.ID)
  # grabbing a Case.ID from pids then putting it back in the drawing to be picked again and creating a pidsstar list (sample pids, same length with replacement)
  pidsstar <- sample(pids, length(pids), replace = T)
  # lapply performing the function(Case.ID) on every item in pidsstar - identify rows in old dataset that correspond with the new order of data in pidsstar
  # do.call"c" concatenates them into a vector
  rows <- do.call("c", lapply(pidsstar, function(Case.ID) which(GEMS_data2$Case.ID == Case.ID)))
  # rowsPerPid says how many rows per Case.ID in pidsstar but pids are in the order in which they were selected 
  rowsPerPid <- do.call("c", lapply(pidsstar, function(Case.ID) sum(GEMS_data2$Case.ID == Case.ID)))
  # 
  boot <- GEMS_data2[rows, ]
  # adding new var (newpid) to boot df - assigning a new Case.ID per child to be able to group them 
  boot$newpid <- rep(1:length(pids), times = rowsPerPid)
  
  GEMS_data3 <- boot %>%
    rowwise() %>%
    mutate(sub_sum = sum(c_across(c(aeromonas_sub:ETEC_sub)), na.rm = T))
  
  GEMS_data4 <- GEMS_data3 %>%
    group_by(country) %>%
    summarise(mean= mean(sub_sum))
  
  GEMS_data5 <- t(GEMS_data4)
  colnames(GEMS_data5) = GEMS_data5[1, ]  # make first row column names

  GEMS_data5 <- as.data.frame(GEMS_data5)
  GEMS_data5 <- GEMS_data5[-c(1), ]

return(GEMS_data5)

}

subclinical_gems_0_2 <- as.data.frame(do.call("rbind", lst))

# re-order so that the gems sites match the order for the weights
subclinical_gems_0_2 <- subclinical_gems_0_2[,c(2,5,6,4,3,1,7)]

# save the dataset 
save(subclinical_gems_0_2,file="Scripts - Aim 3/Output/subclinical_gems_0_2.Rdata")
