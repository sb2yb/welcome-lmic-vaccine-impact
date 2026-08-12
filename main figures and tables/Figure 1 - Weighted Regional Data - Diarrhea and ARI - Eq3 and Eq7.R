# This script is for a figure of Eq3 (diarrhea and ARI) AND figure of Eq7 (diarrhea and ARI)

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

# load scripts
load("Scripts - Aim 3/Output/Eq3_Figure_output_25Nov2025.Rdata")
load("Scripts - Aim 3/Output/Eq7_Figure_output_25Nov2025.Rdata")


# equation 3
q3 = ggplot(eq3_combined_long, aes(fill=pathogen, x=region, y=median,)) 
q3 = q3 + geom_bar(position="dodge" , stat="identity", colour="black") 
q3 = q3 + geom_errorbar(aes(ymin=lower_pos, ymax=upper_pos), width=.3, position=position_dodge(.9))
# q3 = q3 + geom_text(aes(label=round(median, digits=2)), position=position_dodge(width=0.9), vjust=0.4, hjust=-0.6)
q3 = q3 + theme_bw() 

q3 = q3 + scale_fill_manual(
  name="",
  breaks=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  labels=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  values=c("#e3342f", "#f6993f", "#ffed4a", "#38c172", "#4dc0b5", "#3490dc", "#6574cd", "#9561e2", "#f66d9b"),
  guide = guide_legend(reverse = FALSE))
q3 = q3 + theme(legend.position="bottom", legend.title = element_text(size=14), legend.text = element_text(size= 14)) + guides(fill=guide_legend(nrow=1, byrow=TRUE))
q3 = q3 + 
  labs(
    x = "",
    y = "Number prevented per 100 child years (95% CI)")
# add space before the axis title 
q3 = q3 + theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
q3 = q3 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
# remove panel boarder
q3 = q3 + theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.border = element_blank(),plot.margin=unit(c(1, 1, 1, 1),"cm"))
# keeps axis lines 
q3 = q3 + theme(axis.line = element_line(color = 'black'))
q3 = q3 + scale_y_continuous(expand = c(0, 0), limits = c(0, 35))
q3 = q3 + scale_x_discrete(labels = function(x) str_wrap(x, width = 35))
q3 = q3 + theme(axis.text=element_text(size=16),
              axis.title=element_text(size=16))
q3 = q3 + theme(panel.spacing = unit(2, "lines"))
q3 = q3 + coord_flip()
q3


# equation 7
q7 = ggplot(eq7_combined_long, aes(fill=pathogen, x=region, y=median,)) 
q7 = q7 + geom_bar(position="dodge" , stat="identity", colour="black") 
q7 = q7 + geom_errorbar(aes(ymin=lower_pos, ymax=upper_pos), width=.3, position=position_dodge(.9))
# q7 = q7 + geom_text(aes(label=round(median, digits=2)), position=position_dodge(width=0.9), vjust=0.4, hjust=-0.6)
q7 = q7 + theme_bw() 

q7 = q7 + scale_fill_manual(
  name="",
  breaks=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  labels=c("Shigella","Campylobacter","ETEC","Norovirus","Rotavirus","Adenovirus 40/41", "Hib", "PCV", "RSV"),
  values=c("#e3342f", "#f6993f", "#ffed4a", "#38c172", "#4dc0b5", "#3490dc", "#6574cd", "#9561e2", "#f66d9b"),
  guide = guide_legend(reverse = FALSE))
q7 = q7 + theme(legend.position="bottom", legend.title = element_text(size=14), legend.text = element_text(size= 14)) + guides(fill=guide_legend(nrow=1, byrow=TRUE))
q7 = q7 + 
  labs(
    x = "",
    y = "Number prevented per 100 child years (95% CI)")
# add space before the axis title 
q7 = q7 + theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
q7 = q7 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
# remove panel boarder
q7 = q7 + theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.border = element_blank(),plot.margin=unit(c(1, 1, 1, 1),"cm"))
# keeps axis lines 
q7 = q7 + theme(axis.line = element_line(color = 'black'))
q7 = q7 + scale_y_continuous(expand = c(0, 0), limits = c(0, 35))
q7 = q7 + scale_x_discrete(labels = function(x) str_wrap(x, width = 35))
q7 = q7 + theme(axis.text=element_text(size=16),
              axis.title=element_text(size=16))
q7 = q7 + theme(panel.spacing = unit(2, "lines"))
q7 = q7 + coord_flip()
q7

# arrange side by side 
q3_q7 <- grid.arrange(q3, q7, ncol = 1)
q3_q7

ggsave(file="Scripts - Aim 3/Figures/Figure - Weighted Regional Data - Diarrhea - ARI - Eq3 - Eq7_25Nov2025.pdf",q3_q7,width=15,height=22,dpi=300)
ggsave(file="Scripts - Aim 3/Figures/Figure - Weighted Regional Data - Diarrhea - ARI - Eq3 - Eq7_25Nov2025.png",q3_q7,width=15,height=22,dpi=300)
