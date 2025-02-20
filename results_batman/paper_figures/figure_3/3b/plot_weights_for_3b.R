rm(list=ls())#Clear environment

# Packages
library(readxl)
library(pROC)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

###########  Read NLV TCR weight Data ##########
weights_nlv <- read.csv('inferred_weights_NLV_TCRs.csv')

colnames(weights_nlv) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')
weights_nlv <- pivot_longer(weights_nlv, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))

colnames(weights_nlv) <- c('tcr','position','weight')

###################################
######## Plotting #################
###################################

#MHCI
p<- ggplot(weights_nlv, aes(x=factor(position), y=weight)) +
    geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
           alpha=1,color="#3A85A8") +
     labs(y="Positional weight",
       x="Position",
       color="") +
      theme_cowplot() +
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        legend.text =element_text(size=20)) +
    facet_wrap(vars(tcr),nrow = 6)+
  theme(strip.text = element_text(
    size = 10, color = "black")) +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  theme(legend.position = "bottom")

# Save plot
 ggsave(plot=p, "weights_NLV_TCRs_for_3b.pdf",
        width=6, height=6)
# 
# # Save raw data
 write.csv(weights_nlv, "raw_data_fig_3b.csv")
 
 
