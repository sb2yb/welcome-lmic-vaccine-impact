# This script is for a figure of Eq3 (diarrhea and ARI)

# Tidy up the workspace
rm(list=ls())

library(data.table) 
library(dplyr)
library(gridExtra)
library(forcats) 
library(ggplot2)
library(stringr)
library(tidyr)
library(grid)
library(tidyverse)

# diarrhea 
load("Scripts - Aim 3/Output/Diarrhea_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")

# ARI
load("Scripts - Aim 3/Output/ARI_AllEqs_Med_CI_DHS - Weighted by Region_24Nov2025.Rdata")

# diarrhea 
# make rownames a column - regions are now a column 
diarrhea_all_eq_med_CI_weighted <- rownames_to_column(diarrhea_all_eq_med_CI_weighted, var = "region")

# keep only equation3
eq3_diarrhea <- diarrhea_all_eq_med_CI_weighted %>%
  select(matches("region|Eq3"))

# ari 
# make rownames a column - regions are now a column 
ari_all_eq_med_CI_weighted <- rownames_to_column(ari_all_eq_med_CI_weighted, var = "region")

# keep only equation3
eq3_ari <- ari_all_eq_med_CI_weighted %>%
  select(matches("Eq3"))

# bind eq3_diarrhea and eq3_ari 
eq3_combined <- cbind(eq3_diarrhea, eq3_ari)

# make long 
eq3_combined_long <- eq3_combined %>%
  pivot_longer(
    cols = -region,
    names_to = c("pathogen", ".value"),
    names_pattern = "Eq3_(.*?)_(median|lower|upper)")

eq3_combined_long$pathogen <- factor(eq3_combined_long$pathogen, levels = c("rsv", "pcv", "hib","adeno","rota","noro","ETEC","campy","shigella"),
                                  labels = c("RSV", "PCV", "Hib","Adenovirus 40/41","Rotavirus","Norovirus","ETEC","Campylobacter","Shigella"))

eq3_combined_long$region <- factor(eq3_combined_long$region, levels = c("WPR","SEAR","EUR","EMR","AMR","AFR"),
                                     labels = c("WPR","SEAR","EUR","EMR","AMR","AFR"))

# make median, lower, and upper positive
eq3_combined_long <- eq3_combined_long %>%
  mutate(median = abs(median)) %>%
  mutate(upper_pos = abs(lower)) %>% # since making negative, need to swap lower and upper
  mutate(lower_pos = abs(upper)) # since making negative, need to swap lower and upper

# sort by pathogen (forces the above level to kick in)
eq3_combined_long <- eq3_combined_long[order(eq3_combined_long$pathogen ),]


q = ggplot(eq3_combined_long, aes(fill=pathogen, x=region, y=median,)) 
q = q + geom_bar(position="dodge" , stat="identity", colour="black") 
q = q + geom_errorbar(aes(ymin=lower_pos, ymax=upper_pos), width=.3, position=position_dodge(.9))
# q = q + geom_text(aes(label=round(median, digits=2)), position=position_dodge(width=0.9), vjust=0.4, hjust=-0.6)
q = q + theme_bw() 

q = q + scale_fill_manual(
  name="",
  breaks=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  labels=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  values=c("#e3342f", "#f6993f", "#ffed4a", "#38c172", "#4dc0b5", "#3490dc", "#6574cd", "#9561e2", "#f66d9b"),
  guide = guide_legend(reverse = FALSE))
q = q + theme(legend.position="bottom", legend.title = element_text(size=14), legend.text = element_text(size= 14)) + guides(fill=guide_legend(nrow=1, byrow=TRUE))
q = q + 
  labs(
    x = "",
    y = "Number prevented per 100 child years (95% CI)")
# add space before the axis title 
q = q + theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
q = q + theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
# remove panel boarder
q = q + theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.border = element_blank(),plot.margin=unit(c(1, 1, 1, 1),"cm"))
# keeps axis lines 
q = q + theme(axis.line = element_line(color = 'black'))
q = q + scale_y_continuous(expand = c(0, 0), limits = c(0, 25))
q = q + scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
q = q + theme(axis.text=element_text(size=16),
              axis.title=element_text(size=16))
q = q + theme(panel.spacing = unit(2, "lines"))
q = q + coord_flip()
q


ggsave(file="Scripts - Aim 3/Figures/Figure - Weighted Regional Data - Diarrhea - ARI - Eq3_25Nov2025.pdf",q,width=22,height=10,dpi=300)
ggsave(file="Scripts - Aim 3/Figures/Figure - Weighted Regional Data - Diarrhea - ARI - Eq3_25Nov2025.png",q,width=22,height=10,dpi=300)

# save eq3_combined_long so that it can be combined with the eq7 data
save(eq3_combined_long, file = "Scripts - Aim 3/Output/Eq3_Figure_output_25Nov2025.Rdata")
