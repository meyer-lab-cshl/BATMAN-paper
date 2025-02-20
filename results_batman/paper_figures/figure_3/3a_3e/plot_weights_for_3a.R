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
all_TCRs <- c('pan-TCR','a3a','NYE-S1','TCR7',
              'T1','FLT3DY','TCR2-T','TCR-F5','TCR-3598-2')

###########  Read MHCI Data ##########
weights_9mer <- read.csv('inferred_parameters/inferred_weights_9_mers.csv')[,1:10]

colnames(weights_9mer) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')
weights_9mer <- pivot_longer(weights_9mer, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))

colnames(weights_9mer) <- c('tcr','position','weight')

weights_9mer <- weights_9mer[weights_9mer$tcr %in% all_TCRs,]

###########  Read MHCII Data ##########
weights_mhcii <- rbind(read.csv('inferred_parameters/inferred_weights_TCR-F5.csv')[,1:21],
                       read.csv('inferred_parameters/inferred_weights_TCR-3598-2.csv')[,1:21])

colnames(weights_mhcii) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10',
                             'P11','P12','P13','P14','P15','P16','P17','P18',
                             'P19','P20')
weights_mhcii <- pivot_longer(weights_mhcii, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10',
                                      'P11','P12','P13','P14','P15','P16','P17','P18',
                                      'P19','P20'))

colnames(weights_mhcii) <- c('tcr','position','weight')

weights_mhcii <- weights_mhcii[weights_mhcii$tcr %in% all_TCRs,]

###################################
######## Plotting #################
###################################

#MHCI
p<- ggplot(weights_9mer, aes(x=factor(position), y=weight, 
                            fill=factor(tcr,levels=all_TCRs))) +
    geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
           alpha=1) +
     labs(y="Positional weight",
       x="Position",
       color="") +
  scale_fill_manual(values = c('black',"#E41A1C", "#3A85A8", "#46A169", "#629363",
                                        "#77777C", "#A6548B")) +
    ylim(0, 1) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        legend.text =element_text(size=20)) +
  facet_wrap(vars(factor(tcr,levels=all_TCRs)),nrow = 1)+
  theme(strip.text = element_text(
    size = 20, color = "black")) +
  theme(legend.position = "bottom")

# Save plot
 ggsave(plot=p, "weights_9mer_for_3a.pdf",
        width=8, height=2)
# 
# # Save raw data
 write.csv(weights_9mer, "raw_data_mhci_fig_3a.csv")
 
 
#MHCI
 p<- ggplot(weights_mhcii, aes(x=factor(position,levels=c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10',
                                                          'P11','P12','P13','P14','P15','P16','P17','P18',
                                                          'P19','P20')), y=weight, 
                              fill=factor(tcr,levels=all_TCRs))) +
   geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
            alpha=1) +
   labs(y="Positional weight",
        x="Position",
        color="") +
   scale_fill_manual(values = c("#F2E631","#D8B62E")) +
                                         ylim(0, 1) +
   theme_cowplot() +
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(size=20),
         axis.title.y = element_text(size=20),
         axis.title.x = element_text(size=20),
         legend.text =element_text(size=20)) +
   facet_wrap(vars(factor(tcr,levels=all_TCRs)),nrow = 1)+
   theme(strip.text = element_text(
     size = 20, color = "black")) +
   theme(legend.position = "bottom")
 
 # Save plot
 ggsave(plot=p, "weights_mhcii_for_3a.pdf",
        width=8, height=2)
 # 
 # # Save raw data
 write.csv(weights_9mer, "raw_data_mhcii_fig_3a.csv") 





  







