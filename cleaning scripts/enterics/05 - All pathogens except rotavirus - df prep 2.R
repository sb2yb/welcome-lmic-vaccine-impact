# This was the comparison paper prep 2 script that additionally filtered to the 5 sites without rotavirus vaccination (amongst other things)
# This script is removing that filter statement and the rotavirus derived variables 

## Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(splitstackshape)
library(tidyr)
library(naniar)

# load scripts created in "AMR - Aim 1 Part B - Code for Part B - df prep step 1"
# load("abx output/abx_all_courses.Rda")
# while looks identical to abx_all_courses, this one does not have NAs for the pathogens 
load("abx output/abx_all_courses1_Aim3_VE.Rda")

abxpathogens = c("adenovirus_40_41","astrovirus","campylobacter_jejuni_coli","cryptosporidium","norovirus_gii","rotavirus","sapovirus","shigella_eiec","tEPEC","ETEC")
subclinicalpathogens = c("EAEC", "giardia", "campylobacter_pan", "aEPEC", "e_bieneusi", "ETEC", "tEPEC", "shigella_eiec", "cryptosporidium")

# not sure they should be named this way? not sure if i should pick the subclinical or diarrheal var name or the cause var i created 
causeandsubpathogens = c("campy_cause", "crypto_cause", "ETEC_cause", "shigella_cause", 'tEPEC_cause')
onlysubclinicalpathogens = c("EAEC", "giardia", "aEPEC", "e_bieneusi", "tEPEC")

#################################################################################
#       PREPARING THE DATAFRAME FOR THE 4 CALCULATION STEPS TO FOLLOW           #
#################################################################################



################ 
#    FOR #1    # 
################ 

# We need to identify each instance where each of the 9 subclinical pathogens are truly subclinical 
abx_all_courses2 <- abx_all_courses1 %>%
  mutate(aEPEC_sub=ifelse(aEPEC >0 & aEPEC_afrep <0.5, 1, 0)) %>%
  mutate(campy_sub=ifelse(campylobacter_pan >0 & campylobacter_jejuni_coli_afe <0.5, 1, 0)) %>%
  mutate(crypto_sub=ifelse(cryptosporidium >0 & cryptosporidium_afe <0.5, 1, 0)) %>%
  mutate(EAEC_sub=ifelse(EAEC >0 & EAEC_afrep <0.5, 1, 0)) %>%
  mutate(e_bieneusi_sub=ifelse(e_bieneusi >0, 1, 0)) %>% 
  mutate(ETEC_sub=ifelse(ETEC >0 & ETEC_afe <0.5, 1, 0)) %>%
  mutate(giardia_sub=ifelse(giardia >0, 1, 0)) %>%
  mutate(shigella_sub=ifelse(shigella_eiec >0 & shigella_eiec_afe <0.5, 1, 0)) %>%
  mutate(tEPEC_sub=ifelse(tEPEC >0 & tEPEC_afe < 0.5, 1, 0))

# summing all the instances of subclinical bacterial pathogens
abx_all_courses2$bacterialpath_sub <- rowSums(abx_all_courses2[c("aEPEC_sub",'campy_sub','EAEC_sub','ETEC_sub','shigella_sub','tEPEC_sub')])

# summing all the instances of subclinical parasitic pathogens
abx_all_courses2$parasiticpath_sub <- rowSums(abx_all_courses2[c("crypto_sub","e_bieneusi_sub","giardia_sub")])  
  
# The 4 additional versions of combined bacterial pathogens due to those that can be cause + subclinical: 
abx_all_courses2$bacterialpath_sub_no_campy <- rowSums(abx_all_courses2[c("aEPEC_sub",'EAEC_sub','ETEC_sub','shigella_sub','tEPEC_sub')])
abx_all_courses2$bacterialpath_sub_no_ETEC <- rowSums(abx_all_courses2[c("aEPEC_sub",'campy_sub','EAEC_sub','shigella_sub','tEPEC_sub')])
abx_all_courses2$bacterialpath_sub_no_shigella <- rowSums(abx_all_courses2[c("aEPEC_sub",'campy_sub','EAEC_sub','ETEC_sub','tEPEC_sub')])
abx_all_courses2$bacterialpath_sub_no_tEPEC <- rowSums(abx_all_courses2[c("aEPEC_sub",'campy_sub','EAEC_sub','ETEC_sub','shigella_sub')])

  #mutate(bacterialpath_sub_no_campy=ifelse(EAEC_sub >0 | aEPEC_sub >0 | ETEC_sub >0 | tEPEC_sub >0 | shigella_sub >0, 1, 0))
  #mutate(bacterialpath_sub_no_ETEC=ifelse(EAEC_sub >0 | campy_sub >0 | aEPEC_sub >0 | tEPEC_sub >0 | shigella_sub >0, 1, 0)) 
  #mutate(bacterialpath_sub_no_shigella=ifelse(EAEC_sub >0 | campy_sub >0 | aEPEC_sub >0 | ETEC_sub >0 | tEPEC_sub >0, 1, 0)) 
  #mutate(bacterialpath_sub_no_tEPEC=ifelse(EAEC_sub >0 | campy_sub >0 | aEPEC_sub >0 | ETEC_sub >0 | shigella_sub >0, 1, 0)) 

# There is 1 additional version of combined parasitic pathogens due to crypot being a cause + subclinical:
abx_all_courses2$parasiticpath_sub_no_crypto <- rowSums(abx_all_courses2[c("e_bieneusi_sub","giardia_sub")])   

#mutate(parasiticpath_sub_no_crypto=ifelse(giardia_sub >0 | e_bieneusi_sub >0, 1, 0))
  
#check <- abx_all_courses2 %>%
  #dplyr::select(ends_with("sub"), bacterialpath_sub, parasiticpath_sub)

## add a variable that is "other abx"
abx_all_courses2 <- abx_all_courses2 %>%
  mutate(otherany=ifelse(cephaloany!=1 & fluoroany!=1 & macroany!=1 & sulfonany!=1,1,0))

#############
#  FOR #5   #
#############

# ALL ABX
abx_all_courses2$aEPEC_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$aEPEC_sub==1),6])
abx_all_courses2$campy_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$campy_sub==1),6])
abx_all_courses2$EAEC_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$EAEC_sub==1),6])
abx_all_courses2$ETEC_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$ETEC_sub==1),6])
abx_all_courses2$shigella_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$shigella_sub==1),6])
abx_all_courses2$tEPEC_abxdays <- sum(abx_all_courses2[which(abx_all_courses2$tEPEC_sub==1),6])

#CEPHALO
abx_all_courses2$aEPEC_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$aEPEC_sub==1 & abx_all_courses2$cephaloany==1),6])
abx_all_courses2$campy_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$campy_sub==1 & abx_all_courses2$cephaloany==1),6])
abx_all_courses2$EAEC_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$EAEC_sub==1 & abx_all_courses2$cephaloany==1),6])
abx_all_courses2$ETEC_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$ETEC_sub==1 & abx_all_courses2$cephaloany==1),6])
abx_all_courses2$shigella_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$shigella_sub==1 & abx_all_courses2$cephaloany==1),6])
abx_all_courses2$tEPEC_abxdaysc <- sum(abx_all_courses2[which(abx_all_courses2$tEPEC_sub==1 & abx_all_courses2$cephaloany==1),6])

#FLUORO
abx_all_courses2$aEPEC_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$aEPEC_sub==1 & abx_all_courses2$fluoroany==1),6])
abx_all_courses2$campy_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$campy_sub==1 & abx_all_courses2$fluoroany==1),6])
abx_all_courses2$EAEC_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$EAEC_sub==1 & abx_all_courses2$fluoroany==1),6])
abx_all_courses2$ETEC_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$ETEC_sub==1 & abx_all_courses2$fluoroany==1),6])
abx_all_courses2$shigella_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$shigella_sub==1 & abx_all_courses2$fluoroany==1),6])
abx_all_courses2$tEPEC_abxdaysf <- sum(abx_all_courses2[which(abx_all_courses2$tEPEC_sub==1 & abx_all_courses2$fluoroany==1),6])

#MACRO
abx_all_courses2$aEPEC_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$aEPEC_sub==1 & abx_all_courses2$macroany==1),6])
abx_all_courses2$campy_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$campy_sub==1 & abx_all_courses2$macroany==1),6])
abx_all_courses2$EAEC_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$EAEC_sub==1 & abx_all_courses2$macroany==1),6])
abx_all_courses2$ETEC_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$ETEC_sub==1 & abx_all_courses2$macroany==1),6])
abx_all_courses2$shigella_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$shigella_sub==1 & abx_all_courses2$macroany==1),6])
abx_all_courses2$tEPEC_abxdaysm <- sum(abx_all_courses2[which(abx_all_courses2$tEPEC_sub==1 & abx_all_courses2$macroany==1),6])

#SULFON
abx_all_courses2$aEPEC_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$aEPEC_sub==1 & abx_all_courses2$sulfonany==1),6])
abx_all_courses2$campy_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$campy_sub==1 & abx_all_courses2$sulfonany==1),6])
abx_all_courses2$EAEC_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$EAEC_sub==1 & abx_all_courses2$sulfonany==1),6])
abx_all_courses2$ETEC_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$ETEC_sub==1 & abx_all_courses2$sulfonany==1),6])
abx_all_courses2$shigella_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$shigella_sub==1 & abx_all_courses2$sulfonany==1),6])
abx_all_courses2$tEPEC_abxdayss <- sum(abx_all_courses2[which(abx_all_courses2$tEPEC_sub==1 & abx_all_courses2$sulfonany==1),6])

################ 
#    FOR #2    # 
################ 

# We need an indicator for each of the 4 pathogens that can be causative and subclinical for when they are CAUSATIVE
# campy, crypto, ETEC, shigella 
abx_all_courses2 <- abx_all_courses2 %>%
  mutate(campy_cause= ifelse(campylobacter_jejuni_coli_afe >0.5, 1, 0)) %>%
  mutate(crypto_cause= ifelse(cryptosporidium_afe >0.5, 1, 0)) %>%
  mutate(ETEC_cause= ifelse(ETEC_afe >0.5, 1, 0)) %>%
  mutate(shigella_cause= ifelse(shigella_eiec_afe >0.5, 1, 0)) %>%
  mutate(tEPEC_cause=ifelse(tEPEC_afe >0.5, 1, 0)) %>%
# ADDED 9/10/21: We are now adding aEPEC and EAEC to fig 2 and I need a "cause" pathogen to do so
  mutate(aEPEC_cause=ifelse(aEPEC_afrep >0.5, 1, 0)) %>%
  mutate(EAEC_cause=ifelse(EAEC_afrep >0.5, 1, 0))
  
################ 
#    FOR #3    # 
################ 

# Within the same row, for each subclinical path (9), indicate when each diarrheal path (10) is the cause of diarrhea
# e.g. In the same row, shigella is present as subclinical infection & rotavirus afe >0.5 (but want to make sure shigella is not the cause)

# These subclinical pathogens (aEPEC, EAEC, e_bieneusi, giardia, and tEPEC) can have the same code format where their Ct must be >0 and the AFes
# for all 10 diarrhea causing pathogens are >0.5 
# The remain subclinical pathogens need a contingency that they are not also the cause of the diarrhea (that they are TRULy subclinical)


###############################################
# 6/30/21
# TRYING TO ADD MORE COLUMNS WHERE EAEC, CAMPY, ETEC, AEPEC, TEPEC, AND SHIGELLA WERE EXPOSED TO ABX DURING DIARRHEA DUE TO *ANY* CAUSE - "ANY PATHOGEN OR NO PATHOGEN DETECTED"
# NUMERATOR: N times pathogen X exposed to abx during diarrhea due to any cause 
# Just count instances in which a diarrheal stool was linked (should be N=13629) and substract out the 10 diarrheal causing pathogens

# anycause: if all 10 pathogens = 0 then they're assigned a 1 for "anycause" other than the 10 pathogens 
abx_all_courses2 <- abx_all_courses2 %>%
  mutate(anycause_afe = ifelse(adenovirus_40_41_afe<0.5 & astrovirus_afe<0.5 & campylobacter_jejuni_coli_afe<0.5 & cryptosporidium_afe<0.5 & 
                             norovirus_gii_afe<0.5 & rotavirus_afe<0.5 & sapovirus_afe<0.5 & shigella_eiec_afe<0.5 & tEPEC_afe<0.5 
                           & ETEC_afe<0.5 & diarassociatednew==1, 1, 0))
###############################################

### TIM'S HANDY WORK ###

filterscp <- function(scp) {
  # Check for cause and subclinical
  # grep is a character search function - try and match the first 4 letters of scp anywhere in causeandsubpathogens and then tells you were it found a match
  causepathidx <- grep(substr(scp, 1, 5), causeandsubpathogens)
  # if it found a match, the length of the causematchidx >0 
  if(length(causepathidx) > 0) {
    if(length(causepathidx) > 1) {
      # this is a check on grep in case there are 2 or more pathogens that have the same 4 letters in a row- this will throw off the matching
      stop("MULTIPLE PATHOGENS MATCHED")
    }
    # all the causeandsubpathogens pathogens are getting stored in causefilter. e.g., for shigella, it says if causepathidx=4 it then finds shigella_cause in abx_all_courses2
    CauseFilter <- abx_all_courses2[,causeandsubpathogens[causepathidx]]
    
  } else {  
    # so make 0 for the other subclinical pathogens that are never the cause of diarrhea, their causefilter = 0 (e.g. giardia)
    CauseFilter <- 0
  }
  
  # selecting all diarrheal pathogens by using "ends with afe". get rid of ST_ETEC. 2 means columns. 
  # function (col) is going through the diarrheal pathogens. col is essentially going through all 10 diarrheal pathogens column by column. swapping them out one by one. 
  # scp is subclinical pathogen holder that we call below. col refers to diarrheal pathogens to get afe >0.5. scp is essentially all 5 only sub clin pathogens (1 at a time)
  # ifelse is looking row by row between the diarrheal pathogens and subclinical pathogens and assigning 0/1 
  ret <- apply(dplyr::select(abx_all_courses2, ends_with("afe"), -ST_ETEC_afe), 2, 
    function(col) ifelse(abx_all_courses2[,scp] > 0 & CauseFilter != 1 & col > 0.5, 1, 0))
  # dpath is essentially grabbing the name of the diarrhea pathogens - splitting before the _ (for naming them below)
  dpath <- sapply(strsplit(colnames(ret), "_"), "[", 1)
  # how we want the new variable names returned to us - subclinicalpath_exp_diarrhealpath. ret is a matrix. one row for every stool and one column for each pathogen.
  # naming the columns for all the variables we've made 
  colnames(ret) <- paste(scp, "exp", dpath, sep = "_")
  
  # so this is where we drop when a subclinical pathogen is the same as a diarrheal pathogen (drop the 4 column for the ones that overlap (e.g. shigella as cause and sub))
  if(length(causepathidx) > 0) { #Drop a column
    dropidx <- grep(substr(scp, 1, 5), dpath)
    ret <- ret[,-dropidx]
  }
  
  return(ret)
}

# looping through the onlysubclinicalpathogens through the function above 
# all_scp is what is saving it as a bunch of different matrices saved in a list. so you could call up each subclinical pathogen's data (10 columns each) on it's own
# onlysubclinicalpathogens is only able to link to scp because I put the exact columns names in the vector (onlysubclinicalpathogens)
# lapply works because it is saying "onlysubclinicalpathogens" = scp in the filterscp function
# lappy returns a list of 5 onlysubclinicalpathogens
all_scp <- lapply(subclinicalpathogens, filterscp)
# this names the 5 lists that lapply makes - names all of the 5 elements that live in all_scp - how to view it: str(all_scp)
names(all_scp) <- subclinicalpathogens
# swap about pathogen in quotes to see what is in all_scp for that pathogen
# sanity checks
head(all_scp[["tEPEC"]])

# combining the individuals matrices into one dataframe instead of each subclinical pathogen having their own matrix/df
all_scp_df <- data.frame(do.call("cbind", all_scp))

# adding newly created variables to abx_all_courses2 and saving as a new dataframe 
abx_all_courses3 <- as.data.frame(cbind(abx_all_courses2, all_scp_df))


#### ADDING IN ADDITIONAL 20 VARIABLES FOR "OVERALL BACTERIAL" AND "OVERALL PARASITIC" PATHOGENS MATCHED TO DIARRHEAL PATHOGENS ####

# BACTERIAL 
# REVISED Sept 23, 2021
# This makes sure that each of the subclinical pathogens is counted and not simply 0/1
retbac <- apply(dplyr::select(abx_all_courses2, ends_with("afe"), -ST_ETEC_afe), 2, 
              function(col) ifelse(col > 0.5, paste(abx_all_courses2$bacterialpath_sub), 0))

dpath <- sapply(strsplit(colnames(retbac), "_"), "[", 1)
colnames(retbac) <- paste("bacterialpath_sub", "exp", dpath, sep = "_")

retbac <- as.data.frame(retbac)
retbac <- retbac %>%
  select(-c(bacterialpath_sub_exp_campylobacter, bacterialpath_sub_exp_shigella, bacterialpath_sub_exp_ETEC, bacterialpath_sub_exp_tEPEC))

# Converts the character "NAs" to be actual NAs
retbac <- sapply(retbac,as.numeric)
retbac <- as.data.frame(retbac)

# PARASITIC
# REVISED Sept 23, 2021
# This makes sure that each of the subclinical pathogens is counted and not simply 0/1
retpara <- apply(dplyr::select(abx_all_courses2, ends_with("afe"), -ST_ETEC_afe), 2, 
                 function(col) ifelse(col > 0.5, paste(abx_all_courses2$parasiticpath_sub), 0))
dpath <- sapply(strsplit(colnames(retpara), "_"), "[", 1)
colnames(retpara) <- paste("parasiticpath_sub", "exp", dpath, sep = "_")

retpara <- as.data.frame(retpara)
retpara <- retpara %>%
  select(-c(parasiticpath_sub_exp_cryptosporidium)) 

# Converts the character "NAs" to be actual NAs
retpara <- sapply(retpara,as.numeric)
retpara <- as.data.frame(retpara)

# combine bacterial and parasitic dataframes into 1
bacpara <- cbind(retbac, retpara)

# Adding the bacterial and parasitic grouped pathogens to the abx_all_courses_3 df
abx_all_courses3 <- cbind(abx_all_courses3, bacpara)

# Adding the "special cases" for the 5 bacterial/parasitic pathogens that are subclinical and causative
abx_all_courses3 <- abx_all_courses3 %>%
  mutate(bacterialpath_sub_no_campy_exp_campylobacter= ifelse(bacterialpath_sub_no_campy> 0 & campylobacter_jejuni_coli_afe >0.5, paste(abx_all_courses3$bacterialpath_sub_no_campy), 0)) %>%
  mutate(bacterialpath_sub_no_ETEC_exp_ETEC= ifelse(bacterialpath_sub_no_ETEC> 0 & ETEC_afe  >0.5, paste(abx_all_courses3$bacterialpath_sub_no_ETEC), 0)) %>%
  mutate(bacterialpath_sub_no_shigella_exp_shigella= ifelse(bacterialpath_sub_no_shigella> 0 & shigella_eiec_afe  >0.5, paste(abx_all_courses3$bacterialpath_sub_no_shigella), 0))%>%
  mutate(bacterialpath_sub_no_tEPEC_exp_tEPEC= ifelse(bacterialpath_sub_no_tEPEC>0 & tEPEC_afe > 0.5,  paste(abx_all_courses3$bacterialpath_sub_no_tEPEC), 0)) %>%
  mutate(parasiticpath_sub_no_cryptosporidium_exp_cryptosporidium=ifelse(parasiticpath_sub_no_crypto> 0 & cryptosporidium_afe  >0.5, paste(abx_all_courses3$parasiticpath_sub_no_crypto), 0))

# Converting to numeric
abx_all_courses3$bacterialpath_sub_no_campy_exp_campylobacter <- as.numeric(abx_all_courses3$bacterialpath_sub_no_campy_exp_campylobacter)
abx_all_courses3$bacterialpath_sub_no_ETEC_exp_ETEC <- as.numeric(abx_all_courses3$bacterialpath_sub_no_ETEC_exp_ETEC)
abx_all_courses3$bacterialpath_sub_no_shigella_exp_shigella <- as.numeric(abx_all_courses3$bacterialpath_sub_no_shigella_exp_shigella)
abx_all_courses3$bacterialpath_sub_no_tEPEC_exp_tEPEC <- as.numeric(abx_all_courses3$bacterialpath_sub_no_tEPEC_exp_tEPEC)
abx_all_courses3$parasiticpath_sub_no_cryptosporidium_exp_cryptosporidium <- as.numeric(abx_all_courses3$parasiticpath_sub_no_cryptosporidium_exp_cryptosporidium)

##################################################################
## MAKING VARIABLES THAT WILL GO INTO FIG 3 (RESP INFECTIONS)
# need 6 subclinical pathogens (aEPEC, campy, EAEC, ETEC, shigella, tEPEC) and hierarchy (alri, uri, diarrhea, dysentery, other)

abx_all_courses3 <- abx_all_courses3 %>%
  mutate(aEPEC_sub_and_alri=ifelse(aEPEC_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(aEPEC_sub_and_diarrhea=ifelse(aEPEC_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(aEPEC_sub_and_dysentery=ifelse(aEPEC_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(aEPEC_sub_and_uri=ifelse(aEPEC_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(aEPEC_sub_and_other=ifelse(aEPEC_sub==1 & hierarchy=="Other", 1,0)) %>%
  
  mutate(campy_sub_and_alri=ifelse(campy_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(campy_sub_and_diarrhea=ifelse(campy_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(campy_sub_and_dysentery=ifelse(campy_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(campy_sub_and_uri=ifelse(campy_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(campy_sub_and_other=ifelse(campy_sub==1 & hierarchy=="Other", 1,0)) %>%
  
  mutate(EAEC_sub_and_alri=ifelse(EAEC_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(EAEC_sub_and_diarrhea=ifelse(EAEC_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(EAEC_sub_and_dysentery=ifelse(EAEC_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(EAEC_sub_and_uri=ifelse(EAEC_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(EAEC_sub_and_other=ifelse(EAEC_sub==1 & hierarchy=="Other", 1,0)) %>%
  
  mutate(ETEC_sub_and_alri=ifelse(ETEC_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(ETEC_sub_and_diarrhea=ifelse(ETEC_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(ETEC_sub_and_dysentery=ifelse(ETEC_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(ETEC_sub_and_uri=ifelse(ETEC_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(ETEC_sub_and_other=ifelse(ETEC_sub==1 & hierarchy=="Other", 1,0)) %>%
  
  mutate(shigella_sub_and_alri=ifelse(shigella_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(shigella_sub_and_diarrhea=ifelse(shigella_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(shigella_sub_and_dysentery=ifelse(shigella_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(shigella_sub_and_uri=ifelse(shigella_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(shigella_sub_and_other=ifelse(shigella_sub==1 & hierarchy=="Other", 1,0)) %>%
  
  mutate(tEPEC_sub_and_alri=ifelse(tEPEC_sub==1 & hierarchy=="ALRI", 1,0)) %>%
  mutate(tEPEC_sub_and_diarrhea=ifelse(tEPEC_sub==1 & hierarchy=="Diarrhea", 1,0)) %>%
  mutate(tEPEC_sub_and_dysentery=ifelse(tEPEC_sub==1 & hierarchy=="Dysentery", 1,0)) %>% 
  mutate(tEPEC_sub_and_uri=ifelse(tEPEC_sub==1 & hierarchy=="URI", 1,0)) %>%
  mutate(tEPEC_sub_and_other=ifelse(tEPEC_sub==1 & hierarchy=="Other", 1,0)) %>%

  # adding in "any sub pathogen"
  mutate(any_sub_and_alri=ifelse(bacterialpath_sub>=1 & hierarchy=="ALRI", bacterialpath_sub,0)) %>%
  mutate(any_sub_and_diarrhea=ifelse(bacterialpath_sub>=1 & hierarchy=="Diarrhea", bacterialpath_sub,0)) %>%
  mutate(any_sub_and_dysentery=ifelse(bacterialpath_sub>=1 & hierarchy=="Dysentery", bacterialpath_sub,0)) %>% 
  mutate(any_sub_and_uri=ifelse(bacterialpath_sub>=1 & hierarchy=="URI", bacterialpath_sub,0)) %>%
  mutate(any_sub_and_other=ifelse(bacterialpath_sub>=1 & hierarchy=="Other", bacterialpath_sub,0))


test <- abx_all_courses3 %>%
  select(pid, hierdi_alri, hierdi_diarrhea, hierdi_dysentery, hierdi_uri, hierdi_other, any_sub_and_alri,  any_sub_and_diarrhea, any_sub_and_dysentery,  any_sub_and_uri, any_sub_and_other, bacterialpath_sub)

#############################################  
# FINAL DATAFRAME FOR AIM 1 PART B ANALYSES #  
#############################################  
  
# After all indicator variables have been made, need to subset dataframe to be linked=1 ONLY
# Note that none of the variables made above will have correctly populated for the unlinked cases, but I wanted the variables in 1 dataframe

abx_all_courses4 <- abx_all_courses3 %>%
  filter(linked==1)

###############################################

# to check on subpathogen_exp_anycause
table(abx_all_courses4$EAEC_exp_anycause)
table(abx_all_courses4$campylobacter_pan_exp_anycause)
table(abx_all_courses4$ETEC_exp_anycause)
table(abx_all_courses4$aEPEC_exp_anycause)
table(abx_all_courses4$tEPEC_exp_anycause)
table(abx_all_courses4$shigella_eiec_exp_anycause)

###############################################
# Saving the 3 dataframes made in this script # 
###############################################

save(abx_all_courses2,file="abx output/abx_all_courses2_Aim3_VE.Rda")
save(abx_all_courses3,file="abx output/abx_all_courses3_Aim3_VE.Rda")
save(abx_all_courses4,file="abx output/abx_all_courses4_Aim3_VE.Rda")
