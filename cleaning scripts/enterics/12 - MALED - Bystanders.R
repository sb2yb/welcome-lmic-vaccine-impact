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
load("microtac.Rda")

# filter and order stools (stooltype monthlies, remove 60 month stools, excluded samples with bad ms2
# select only the needed variables
microtac1 <- microtac %>%
  filter(stooltype == "M1" & month_ss <=24 & tac_ms2_exclude==0) %>%
  select(pid, sid, country_id, aeromonas, aEPEC, campylobacter_pan, EAEC, ETEC, h_pylori, salmonella, STEC, shigella_eiec, tEPEC, v_cholerae)
microtac1 = microtac1[order(microtac1$pid,microtac1$sid),]

# create dichotomous variables 
microtac1 <- microtac1 %>%
  mutate(aeromonas_sub = ifelse(aeromonas<35, 1, 0)) %>%
  mutate(aEPEC_sub = ifelse(aEPEC<35, 1, 0)) %>%
  mutate(campylobacter_pan_sub = ifelse(campylobacter_pan<35, 1, 0)) %>%
  mutate(EAEC_sub = ifelse(EAEC<35, 1, 0)) %>%
  mutate(ETEC_sub = ifelse(ETEC<35, 1, 0)) %>%
  mutate(h_pylori_sub = ifelse(h_pylori<35, 1, 0)) %>%
  mutate(salmonella_sub = ifelse(salmonella<35, 1, 0)) %>%
  mutate(STEC_sub = ifelse(STEC<35, 1, 0)) %>%
  mutate(shigella_eiec_sub = ifelse(shigella_eiec<35, 1, 0)) %>%
  mutate(tEPEC_sub = ifelse(tEPEC<35, 1, 0)) %>%
  mutate(v_cholerae_sub = ifelse(v_cholerae<35, 1, 0))
 
# set seed (random.org from 1-1000)
 set.seed(571) # seed set 14 March 24


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
  pids <- unique(microtac1$pid)
  # grabbing a pid from pids then putting it back in the drawing to be picked again and creating a pidsstar list (sample pids, same length with replacement)
  pidsstar <- sample(pids, length(pids), replace = T)
  # lapply performing the function(pid) on every item in pidsstar - identify rows in old dataset that correspond with the new order of data in pidsstar
  # do.call"c" concatenates them into a vector
  rows <- do.call("c", lapply(pidsstar, function(pid) which(microtac1$pid == pid)))
  # rowsPerPid says how many rows per pid in pidsstar but pids are in the order in which they were selected 
  rowsPerPid <- do.call("c", lapply(pidsstar, function(pid) sum(microtac1$pid == pid)))
  # 
  boot <- microtac1[rows, ]
  # adding new var (newpid) to boot df - assigning a new pid per child to be able to group them 
  boot$newpid <- rep(1:length(pids), times = rowsPerPid)
  
  

  microtac2 <- boot %>%
    rowwise() %>%
    mutate(sub_sum = sum(c_across(c(aeromonas_sub:v_cholerae_sub)), na.rm = T))
  
  microtac3 <- microtac2 %>%
    group_by(country_id) %>%
    summarise(mean= mean(sub_sum))
  
  microtac4 <- t(microtac3)
  colnames(microtac4) = microtac4[1, ]  # make first row column names

  microtac4 <- as.data.frame(microtac4)
  microtac4 <- microtac4[-c(1), ]

return(microtac4)

}

subclinical_maled <- as.data.frame(do.call("rbind", lst))

# save the dataset 
save(subclinical_maled,file="Scripts - Aim 3/Output/subclinical_maled - 14March.Rdata")
