# This was the comparison paper prep 1 script that additionally filtered to the 5 sites without rotavirus vaccination (amongst other things)
# This script is removing that filter statement and the rotavirus derived variables 

## Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# load libraries
library(dplyr)
library(splitstackshape)
library(tidyr)

# load the MALED tac data
load("microtac.Rda")

# load in abxcourses df that Liz made (TO BE MODIFIED BELOW)
abxcourses<-read.csv(file="Abx courses envelopes datasets/abxcourses.csv", h=T)
#save(abxcourses, file="Abx courses envelopes datasets/abxcourses.Rda")
#load("Abx courses envelopes datasets/abxcourses.Rda")

# creating variable that indicates if a child had the path prior to 6 (183 days), 9 (274 days), 12 (365 days) and 15 (456 days) months  (for aim 2)
# THE EDIT NEEDS TO BE IN MICROTAC - which is capturing diarrheal episodes

# need to group by PID to do this - capture first instance of path afe > 0.5 and pull that, if never above 0.5, then pull first instance afe >0?

# make all _afe NAs = 0 (this will help with aim 2 equations later)
microtac$shigella_eiec_afe[is.na(microtac$shigella_eiec_afe)] <- 0
microtac$campylobacter_jejuni_coli_afe[is.na(microtac$campylobacter_jejuni_coli_afe)] <- 0
microtac$ETEC_afe[is.na(microtac$ETEC_afe)] <- 0
microtac$norovirus_afe[is.na(microtac$norovirus_afe)] <- 0
microtac$adenovirus_40_41_afe[is.na(microtac$adenovirus_40_41_afe)] <- 0

# we need to remove these instances so that an instance where Ct=35 is chosen for that participant instead below
micro_test <- microtac %>%
  filter(!is.na(agedays))

# Function to loop over all 5 pathogens
path_loop <- function(path_ct) {

micro_test <- micro_test %>%
  group_by(pid) %>%
  filter(!!rlang::sym(path_ct)<35) %>%
  slice(1) # grabs first instance where path <35, problem is we still need a variable to be grabbed for those who never reach that threshold

micro_test1 <- microtac %>%
group_by(pid) %>%
  filter(!!rlang::sym(path_ct)==35) %>% # pull instances where path = 35 for each child - not perfect as it will still pull from pids who have afe > 0.5
  slice(1)

micro_test1 <- micro_test1 %>% select(pid, sid, agedays, stooltype, !!rlang::sym(path_ct)) # create small df of wanted variables (1 obs per kid)

# merge them together and then keep the pid with the lowest ct
test_merge <- rbind(micro_test, micro_test1)

# one obs per child where the lowest Ct for each child was kept (for some kids this means Ct 35)
test_merge1 <- test_merge %>%
  group_by(pid) %>%
  slice(1)

small <- test_merge1 %>% select(pid, sid, agedays, stooltype, !!rlang::sym(path_ct))

return(small)
}

shigella <- path_loop("shigella_eiec")
campy <- path_loop("campylobacter_jejuni_coli")
ETEC <- path_loop("ETEC")
norovirus <- path_loop("norovirus")
adenovirus_40_41 <- path_loop("adenovirus_40_41")

# Create the path_6mo, etc variables within each of the 5 dataframes and merge back in - 
# Keep in mind, this won't make sense when merged back into the larger dataset as the "agedays" will be off when merged back in

age_loop <- function(path_ct, dat){

  dat <- dat %>%
  mutate(mo2 = ifelse(!!rlang::sym(path_ct)<35 & agedays <61, 1, 0)) %>%
  mutate(mo4 = ifelse(!!rlang::sym(path_ct)<35 & agedays <122, 1, 0)) %>%  
  mutate(mo6 = ifelse(!!rlang::sym(path_ct)<35 & agedays <183, 1, 0)) %>%
  mutate(mo9 = ifelse(!!rlang::sym(path_ct)<35 & agedays  <274, 1, 0))

return(dat)
}

shigella1 <- age_loop("shigella_eiec", shigella)
shigella1 <- shigella1 %>% dplyr::rename_with(~ paste0("shig_", .x), matches("mo")) # paste the "shig" in front of the months

campy1 <- age_loop("campylobacter_jejuni_coli", campy)
campy1 <- campy1 %>% dplyr::rename_with(~ paste0("campy_", .x), matches("mo"))

ETEC1 <- age_loop("ETEC", ETEC)
ETEC1 <- ETEC1 %>% dplyr::rename_with(~ paste0("ETEC_", .x), matches("mo"))

norovirus1 <- age_loop("norovirus", norovirus)
norovirus1 <- norovirus1 %>% dplyr::rename_with(~ paste0("noro_", .x), matches("mo"))

adenovirus1 <- age_loop("adenovirus_40_41", adenovirus_40_41)
adenovirus1 <- adenovirus1 %>% dplyr::rename_with(~ paste0("adeno_", .x), matches("mo"))

# only keep "months" columns - merge into abxcourses
shigella1 = subset(shigella1, select = c(pid, shig_mo2:shig_mo9))
campy1 = subset(campy1, select = c(pid, campy_mo2:campy_mo9))
ETEC1 = subset(ETEC1, select = c(pid, ETEC_mo2:ETEC_mo9))
norovirus1 = subset(norovirus1, select = c(pid, noro_mo2:noro_mo9))
adenovirus1 = subset(adenovirus1, select = c(pid, adeno_mo2:adeno_mo9))

test_merge1 <- merge(shigella1, campy1, by="pid")
test_merge1 <- merge(test_merge1, ETEC1, by="pid")
test_merge1 <- merge(test_merge1, norovirus1, by="pid")
test_merge1 <- merge(test_merge1, adenovirus1, by="pid")

# Removing agedays var
abxcourses_mod <- abxcourses %>%
  dplyr::select(-agedays)

# Renaming Pid and age to be consistent with other datasets
# Rename EAEC_afe and aEPEC_afe to distinguish from diarrheal causing pathogens afe 
abxcourses_mod1 <- abxcourses_mod %>%
  dplyr::rename(pid = Pid) %>%
  dplyr::rename(agedays = age) %>%
  dplyr::rename(aEPEC_afrep = aEPEC_afe) %>%
  dplyr::rename(EAEC_afrep = EAEC_afe) %>%
# need to recode resp = 1 if ALRI = 1  
  mutate(resp=ifelse(resp==1 | resp==0 & alri==1, 1, 0)) %>%
# make an uri variable 
  mutate(uri=ifelse(resp==1 & alri==0, 1, 0))

# Hierarchy needed for respiratory fig 3:
abxcourses_mod1 <- abxcourses_mod1 %>%
mutate(hierarchy = case_when(
  alri == 1 ~ "ALRI",
  maxbnew==1 & alri==0 ~ "Dysentery",
  diarassociatednew==1 & maxbnew==0 ~ "Diarrhea",
  uri==1 & maxbnew==0 & diarassociatednew==0 ~ "URI",
  alri==0 & maxbnew==0 & uri==0 & diarassociatednew==0 ~ "Other"))

# Might be easier to adjust the current code if each hierarchy option was its own variable
abxcourses_mod1 <- abxcourses_mod1 %>%
  mutate(hierdi_alri = ifelse(hierarchy=="ALRI", 1,0)) %>%
  mutate(hierdi_dysentery = ifelse(hierarchy=="Dysentery", 1,0)) %>%
  mutate(hierdi_uri = ifelse(hierarchy=="URI", 1,0)) %>%
  mutate(hierdi_diarrhea = ifelse(hierarchy=="Diarrhea", 1,0)) %>%
  mutate(hierdi_other = ifelse(hierarchy=="Other", 1,0)) 

################
# END ADDITION #
################

# create list for abxpathogens - pathogens known to cause diarrhea (10)
# NOTICED - ETEC is in for ST-ETEC, this is okay. In manuscript refer to ETEC as ST-ETEC when referring to causal (since LT will never pop up). Refer to ETEC when subclinical.
abxpathogens = c("adenovirus_40_41","astrovirus","campylobacter_jejuni_coli","cryptosporidium","norovirus_gii","rotavirus","sapovirus","shigella_eiec","tEPEC","ETEC")
for(i in abxpathogens) { microtac[,i] = ifelse(microtac[,i]>=35,35,microtac[,i]); microtac[,i] = (35-microtac[,i])/3.322 }

# create list for subclinical pathogens (Fig 2 Liz's growth paper) - pathogens not known to cause diarrhea 
# do not add viruses as we do not care about sub-clinical pressure for viruses
# NEEDED TO REMOVE 4 PATHOGENS THAT BELONG HERE BECAUSE THEY'VE ALREADY BEEN TRANSFORMED ("ETEC", "tEPEC", "shigella_eiec", "cryptosporidium")
subclinicalpathogens = c("EAEC", "giardia", "campylobacter_pan", "aEPEC", "e_bieneusi")
for(i in subclinicalpathogens) { microtac[,i] = ifelse(microtac[,i]>=35,35,microtac[,i]); microtac[,i] = (35-microtac[,i])/3.322 }

# Updating subclinical pathogens list (9)
subclinicalpathogens = c("EAEC", "giardia", "campylobacter_pan", "aEPEC", "e_bieneusi", "ETEC", "tEPEC", "shigella_eiec", "cryptosporidium")

# from siteagedata keep: pid, sid, agedays, ct values for those in abxpathogens + subclinical pathogens
# transform ct values 
## filter and order stools (any diarrhea stool + stooltype monthlies, remove 60 month stools, excluded samples with bad ms2
siteagedata = filter(microtac,stooltype == "D1" | (stooltype == "M1" & month_ss <=24))
siteagedata = filter(siteagedata,tac_ms2_exclude==0)
siteagedata = siteagedata[order(siteagedata$pid,siteagedata$sid),]

##################################################

siteagedata1 <- siteagedata %>%
  dplyr::select(pid, sid, agedays, stooltype, adenovirus_40_41, astrovirus, campylobacter_jejuni_coli, cryptosporidium, norovirus_gii,
                rotavirus, sapovirus, shigella_eiec, tEPEC, EAEC, giardia, campylobacter_pan, aEPEC, e_bieneusi, ETEC) %>%
  dplyr::rename(agedaysnew = agedays)

# modifying abxcourses to be only 1 row per course
abxcourses_mod2 <- abxcourses_mod1 %>%
  filter(fstab==1)

# duplicate each row of the abxcourse dataset and link the stools (siteagedata) to that
# So the abxcourse dataset is the “primary” and we’re linking stools in.
# I think we should stick with 21 days, and if we are unable to link stools for a substantial proportion of abx courses, 
# then we can think about expanding to 30 days

# modifying david's code

# Creating a variable to be able to group by later on
abxcourses_mod2$rowid <- as.numeric(row.names(abxcourses_mod2))

# Using the expandRows function from splitstackshape I can create 22 duplicate rows for 0-21 days. 
# I'm not sure how to make this so that it doesn't include day 0 (the day the course was started) - Interested in 1-21 (days prior to abx course)
# THIS HAS BEEN UPDATED TO NOW EXTEND OUT 30 DAYS 
abxcourses_mod3 <- expandRows(abxcourses_mod2, count = 31, count.is.col = FALSE)

# I created a way to identify each individual row within an ageday/pid group so that the first row in the group would get a tag of 1, 
# the next would get a 2 and so on up to 21. I then used that tag to subtract from the age in days so that we could get a daily age. 
# Then I merged this over to the siteagedata1 dataframe
# siteagedata1 is being joined INTO abxcourses3 so abxcourses on the "left" and siteagedata1 on the "right"

joined_courses_siteage <- abxcourses_mod3 %>% 
  group_by(pid, agedays) %>% 
  mutate(id = (row_number(pid) - 1)) %>% 
  mutate(agedaysnew = agedays - id) %>%
  left_join(siteagedata1, by = c("pid", "agedaysnew")) %>%
  dplyr::rename(
    sid_diarrhea = sid.x,
    sid_subclinical = sid.y)

testing <- joined_courses_siteage %>% select(pid, agedays, agedaysnew, id, rowid, sid_diarrhea, sid_subclinical)

# Lowest ID (days) per row_id (Course)
# 13629 courses could be linked to a stool out of 15697 courses (87% match)
# ONLY NEED TO DROP MISSING CTS FOR SUBCLINICAL PATHOGENS
unique_courses <- joined_courses_siteage %>% 
  drop_na(sid_subclinical, campylobacter_pan, cryptosporidium, shigella_eiec, tEPEC, EAEC, giardia, aEPEC, e_bieneusi, ETEC)  %>%  # Dropping missing sid_subclinical & missing Cts
  filter(id != 0)    %>%  # dropping missing id=0
  group_by(rowid)    %>% 
  arrange(sid_subclinical, id) %>% 
  slice_tail()       %>%  # Grabs the most recent stool 
  mutate(linked =1)  %>%  # Creating indicator variable that identifies all courses that could be linked
  ungroup()

# Merge the linked data into the original dataframe where no columns are duplicated 
abx_all_courses <- merge(unique_courses, abxcourses_mod2, all.y=TRUE) # so we want all columns from unique courses but all rows from abx_courses_mod2 (linked and unlinked dataset)

# If linked is "NA" making it 0
abx_all_courses$linked[is.na(abx_all_courses$linked)] <- 0

# the creation of abx_all_courses results in sid being added in from abxcourses_mod2 so now we have sid, sid_diarrhea, and sid_subclinical
# sid and sid_diarrhea are the same except sid_diarrhea only covers stools linked to subclinical pathogens (because that's what we need)
# however sid has linked and unlinked so it's a more complete version of the variable. Either one will work, just might be nice to have the more complete one
abx_all_courses <- abx_all_courses %>%
dplyr::rename(sid_diarrhea_complete = sid)

# Save dataframe for analyses 
save(abx_all_courses,file="abx output/abx_all_courses_Aim3_VE.Rda")

#################################################################################
# Convert all NAs for pathogen AFes to be 0
# Again, ETEC is in here in place of ST_ETEC. Any AFe that pops up for ETEC is actually ST_ETEC. 

# diarassociated:	1/0 indicator for whether the course was associated with diarrhea and could be linked to a TAC'd stool 
#                 (note - if the course was associated with diarrhea but a stool wasn’t collected/tested, this will be 0)
# diarassociatednew: 	1/0 indicator for whether the course was associated with diarrhea regardless of whether it could be linked to a TAC'd stool 
#                 (note - if the course was associated with diarrhea but a stool wasn’t collected/tested, this will be 1)

# Sept 14, 2022 - need to see how much missing TAC data there is for the 10 diarrheal pathogens
check1 <- abx_all_courses %>% select(pid, sid_diarrhea_complete, diarassociated, diarassociatednew, hierarchy, linked, adenovirus_40_41_afe, astrovirus_afe, campylobacter_jejuni_coli_afe, cryptosporidium_afe, norovirus_gii_afe,
                                     rotavirus_afe, sapovirus_afe, shigella_eiec_afe, tEPEC_afe, ETEC_afe) %>% filter(diarassociated==1)

table(check1$adenovirus_40_41_afe>0, useNA = "ifany")
table(check1$astrovirus_afe>0, useNA = "ifany")
table(check1$campylobacter_jejuni_coli_afe>0, useNA = "ifany")
table(check1$cryptosporidium_afe>0, useNA = "ifany")
table(check1$norovirus_gii_afe>0, useNA = "ifany")
table(check1$rotavirus_afe>0, useNA = "ifany")
table(check1$sapovirus_afe>0, useNA = "ifany")
table(check1$shigella_eiec_afe>0, useNA = "ifany")
table(check1$tEPEC_afe>0, useNA = "ifany")
table(check1$ETEC_afe>0, useNA = "ifany")

test <- abx_all_courses %>%
  filter(diarassociated==1 & is.na(adenovirus_40_41_afe) & is.na(astrovirus_afe) & is.na(campylobacter_jejuni_coli_afe) &
                                   is.na(cryptosporidium_afe) & is.na(norovirus_gii_afe) & is.na(rotavirus_afe) &
                                   is.na(sapovirus_afe) & is.na(shigella_eiec_afe) & is.na(tEPEC_afe)) %>%
  select(pid, sid_diarrhea_complete, diarassociated, diarassociatednew, hierarchy, linked, adenovirus_40_41_afe, astrovirus_afe, campylobacter_jejuni_coli_afe, cryptosporidium_afe, norovirus_gii_afe,
         rotavirus_afe, sapovirus_afe, shigella_eiec_afe, tEPEC_afe) 

table(test$linked)


# Sept 14, 2022 - if all pathogen AFes (except ETEC) are NA then make diarassociated= 0.
abx_all_courses1 <- abx_all_courses %>%
  mutate(diarassociated = ifelse(diarassociated==1 & is.na(adenovirus_40_41_afe) & is.na(astrovirus_afe) & is.na(campylobacter_jejuni_coli_afe) &
                                   is.na(cryptosporidium_afe) & is.na(norovirus_gii_afe) & is.na(rotavirus_afe) &
                                   is.na(sapovirus_afe) & is.na(shigella_eiec_afe) & is.na(tEPEC_afe),0,diarassociated))


# convert afes that are NA to 0
abx_all_courses1 <- abx_all_courses1 %>%
  dplyr::rename(country_id = Country_ID)
abx_all_courses1$adenovirus_40_41_afe[is.na(abx_all_courses1$adenovirus_40_41_afe)] <- 0
abx_all_courses1$astrovirus_afe[is.na(abx_all_courses1$astrovirus_afe)] <- 0
abx_all_courses1$campylobacter_jejuni_coli_afe[is.na(abx_all_courses1$campylobacter_jejuni_coli_afe)] <- 0
abx_all_courses1$cryptosporidium_afe[is.na(abx_all_courses1$cryptosporidium_afe)] <- 0
abx_all_courses1$norovirus_gii_afe[is.na(abx_all_courses1$norovirus_gii_afe)] <- 0
abx_all_courses1$rotavirus_afe[is.na(abx_all_courses1$rotavirus_afe)] <- 0
abx_all_courses1$sapovirus_afe[is.na(abx_all_courses1$sapovirus_afe)] <- 0
abx_all_courses1$shigella_eiec_afe[is.na(abx_all_courses1$shigella_eiec_afe)] <- 0
abx_all_courses1$ETEC_afe[is.na(abx_all_courses1$ETEC_afe)] <- 0 
abx_all_courses1$tEPEC_afe[is.na(abx_all_courses1$tEPEC_afe)] <- 0

abx_all_courses1$aEPEC_afrep[is.na(abx_all_courses1$aEPEC_afrep)] <- 0
abx_all_courses1$EAEC_afrep[is.na(abx_all_courses1$EAEC_afrep)] <- 0


table(abx_all_courses1$adenovirus_40_41_afe>0, useNA="ifany")
table(abx_all_courses1$astrovirus_afe>0, useNA="ifany")
table(abx_all_courses1$campylobacter_jejuni_coli_afe>0, useNA="ifany")
table(abx_all_courses1$cryptosporidium_afe>0, useNA="ifany")
table(abx_all_courses1$norovirus_gii_afe>0, useNA="ifany")
table(abx_all_courses1$rotavirus_afe>0, useNA="ifany")
table(abx_all_courses1$sapovirus_afe>0, useNA="ifany")
table(abx_all_courses1$shigella_eiec_afe>0, useNA="ifany")
table(abx_all_courses1$ETEC_afe>0, useNA="ifany")
table(abx_all_courses1$tEPEC_afe>0, useNA="ifany")

### MERGING PATH MONTH VARIABLES IN 
abx_all_courses1 <- merge(abx_all_courses1, test_merge1, by="pid")

# SEPTEMBER 30, 2022 #
# We are merging in score for DIARRHEA episodes for Aim 2 analyses 
diarrheascores<-read.csv(file="diarrheascores.csv", h=T)
abx_all_courses1 <- merge(abx_all_courses1, diarrheascores, by.x = c("pid", "agestart"), by.y = c("Pid", "agestart"), all.x = T)


#################################################################################
# NEW AS OF 22 FEB 24
# need to create this statement here:  boot <- boot %>% mutate(est2_shig= ifelse(shigella_eiec_afe>0.5 & diarassociated==1, 1, 0))
# since we have 5 pathogens, we should create the est2 variables here instead of in the boot statement in the analysis script

abx_all_courses1 <- abx_all_courses1 %>%
  mutate(est2_shig= ifelse(shigella_eiec_afe>0.5 & diarassociated==1, 1, 0)) %>%
  mutate(est2_adeno= ifelse(adenovirus_40_41_afe>0.5 & diarassociated==1, 1, 0)) %>%
  mutate(est2_campy= ifelse(campylobacter_jejuni_coli_afe>0.5 & diarassociated==1, 1, 0)) %>%
  mutate(est2_etec= ifelse(ETEC_afe>0.5 & diarassociated==1, 1, 0)) %>%
  mutate(est2_noro= ifelse(norovirus_gii_afe>0.5 & diarassociated==1, 1, 0)) 


################################################
# Use this Rda for the next step in formatting #
################################################

save(abx_all_courses1,file="abx output/abx_all_courses1_Aim3_VE.Rda")
write.csv(abx_all_courses1,"abx output/abx_all_courses1_Aim3_VE.csv")

#################################################################################
#################################################################################

