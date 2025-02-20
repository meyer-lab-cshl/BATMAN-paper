rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)
library(viridis)
library(RColorBrewer)



############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is


distances <- read.csv('peptide_distances_multihamming.csv')

#rename
distances$activation[distances$activation=='0']='No'
distances$activation[distances$activation=='1']='Weak'
distances$activation[distances$activation=='2']='Strong'

# Remove index peptides and 1 HD peptides
distances <- distances[!distances$d_hamming %in% c(0,1),]

# Subset for TCRs with 1HD and >1HD data found from same paper or >1 HD activation
# class clearly mentioned
distances <- distances[distances$tcr %in% c('a3a','A3-05','A3-10',
                                            'TIL1383I','c259','868Z11'),]

############
## Plot ####
############

p <- ggplot(distances, aes(y=d,x=factor(d_hamming))) +
     geom_boxplot(aes(color=factor(activation,
                                       levels=c("Strong","Weak","No")),
                      fill=factor(activation,
                               levels=c("Strong","Weak","No"))),
                      outlier.colour = NA,varwidth = FALSE,
                      show.legend = FALSE, alpha=1,lwd=0.2,
                  position=position_dodge2(0.1,preserve = "single")) +
  scale_color_manual(values=c("#bf212f","#f9a73e","#264b96")) +
  scale_fill_manual(values=c("#FFD6D7","#FFD6D7","#d2ebff")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(activation,
                              levels=c("Strong","Weak","No")),
                 group = factor(activation,levels=c("Strong","Weak","No"))), 
              alpha=0.7, 
              position = position_jitterdodge(jitter.width=0.25,
                                              dodge.width = 0.65),size=0.1) +
               scale_color_manual(values=c("#bf212f","#f9a73e","#264b96")) +
  labs(y="BATMAN distance",
       x="Hamming distance",
       color="") +
  theme_cowplot() +
  facet_wrap(~factor(tcr,levels=c('a3a','A3-05','A3-10',
                                  'TIL1383I','c259','868Z11')),
             nrow=2,scales = "free") +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  theme(axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.y = element_text(size=10))

# Save plot
ggsave(plot=p, "multihamming_distances_for_5a.pdf",
       width=8, height=2.5)
# 
# # Save raw data
write.csv(distances, "raw_data_fig_5a.csv")



