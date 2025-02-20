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


distances <- read.csv('../../figure_4/4d/mutant_distances_a3a.csv')

#rename
distances$activation[distances$activation=='0']='No'
distances$activation[distances$activation=='1']='Weak'
distances$activation[distances$activation=='2']='Strong'

# Remove index peptides
distances <- distances[!distances$d_hamming==0,]

# Separate 1HD and >1 HD
distances$HD <- "1 Hamming mutants"
distances$HD[distances$d_hamming>1] <- ">1 Hamming mutants"

# Change to long form for 2 distance types
distances <- pivot_longer(distances,
                          cols=c("d_pan","d_al"),
                          names_to = "d_type")
distances$d_type[distances$d_type=="d_pan"] <- "Pan-TCR"
distances$d_type[distances$d_type=="d_al"] <- "AL"

# Add labels to well-known mutants
distances$label <- NA

distances$label[distances$peptide=='ESDPIVAQY'] <- "TITIN"
distances$label[distances$peptide=='EVDPIGHVY'] <- "MAGE-A6"
distances$label[distances$peptide=='ETDPVNHMV'] <- "FAT2"



############
## Plot ####
############

p <- ggplot(distances, aes(y=value,x=factor(d_type,levels=c("Pan-TCR","AL")),
                           )) +
     geom_boxplot(aes(color=factor(activation,
                                       levels=c("Strong","Weak","No")),
                      fill=factor(activation,
                               levels=c("Strong","Weak","No"))),
                      outlier.colour = NA,varwidth = FALSE,
                      show.legend = FALSE, alpha=1,
                      position=position_dodge2(0.85,preserve = "single")) +
  scale_color_manual(values=c("#bf212f","#f9a73e","#264b96")) +
  scale_fill_manual(values=c("#FFD6D7","#FFD6D7","#d2ebff")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(d_hamming),
                 group = factor(activation,levels=c("Strong","Weak","No"))), 
              alpha=0.8, 
              position = position_jitterdodge(0.6),size=0.6) +
  scale_color_manual(values = colorRampPalette(brewer.pal(6, "Set1"))(11)[2:10]) +
  labs(y="Peptide-to-index distance",
       x="",
       color="") +
  theme_cowplot() +
  facet_wrap(~factor(HD, levels=c("1 Hamming mutants",
                                  ">1 Hamming mutants")),nrow=1,scales = "free_y") +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  theme(axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12))

# Save plot
ggsave(plot=p, "a3a_distances_for_5b.pdf",
       width=5, height=3)
# 
# # Save raw data
write.csv(distances, "raw_data_fig_5b.csv")



