rm(list=ls())#Clear environment

# Packages
library(readxl)
library(pROC)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyverse)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# selected TCRs for plotting

###########  Read within TCR weights ##########
weights_9mer <- read.csv('../../figure_3/3a-b/inferred_parameters/inferred_weights_9_mers.csv')[,1:10]

colnames(weights_9mer) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')
weights_9mer <- pivot_longer(weights_9mer, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))

colnames(weights_9mer) <- c('tcr','position','weight')

weights_a3a <- weights_9mer[weights_9mer$tcr=='a3a',]

# Normalize
weights_a3a$weight <- (weights_a3a$weight - min(weights_a3a$weight))/(
                            max(weights_a3a$weight) - min(weights_a3a$weight))

weights_a3a$tcr <- 'Within TCR'

# Read Pan-TCR and AL weights
weights_al <- read.csv('a3a_weights.csv')

colnames(weights_al) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')
weights_al <- pivot_longer(weights_al, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))

colnames(weights_al) <- c('tcr','position','weight')

# Normalize Pan-TCR weights
weights_al$weight[weights_al$tcr=='Pan-TCR'] <- 
  (weights_al$weight[weights_al$tcr=='Pan-TCR'] - min(weights_al$weight[weights_al$tcr=='Pan-TCR']))/(
  max(weights_al$weight[weights_al$tcr=='Pan-TCR']) - min(weights_al$weight[weights_al$tcr=='Pan-TCR']))


# Add all data
weights_a3a <- rbind(weights_a3a,weights_al)

# Add sequence
weights_a3a$position <- replace(weights_a3a$position,
                                weights_a3a$position %in% c('P1','P2','P3','P4',
                                                            'P5','P6','P7','P8',
                                                            'P9'),
                                c('E','V','D','P','I','G','H','L','Y'))


###################################
######## Plotting #################
###################################

p<- ggplot(weights_a3a, aes(x=factor(position,levels=c('E','V','D','P','I','G',
                                                       'H','L','Y')), 
                            y=weight)) +
    geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
           alpha=1,fill="#E41A1C") +
     labs(y="Positional weight",
       x=" ",
       color="") +
    ylim(0, 1) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.title.x = element_blank(),
        legend.text =element_text(size=10)) +
  facet_wrap(vars(factor(tcr,levels=c("Pan-TCR","AL","Within TCR"))),ncol = 1)+
  theme(strip.text = element_text(
    size = 12, color = "black")) +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  theme(legend.position = "bottom")

# Save plot
 ggsave(plot=p, "a3a_weights_for_4d.pdf",
        width=1.6, height=3)
# 
# # Save raw data
 write.csv(weights_a3a, "raw_data_fig_4d.csv")
 
 
