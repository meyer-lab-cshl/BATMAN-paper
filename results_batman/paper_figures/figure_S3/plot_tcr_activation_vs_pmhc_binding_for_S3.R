rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)
library(viridis)
library(loo)
library(pROC)
library(RColorBrewer)


############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to current

# Get TCR activation data
tcr_data <- read_excel(paste0('../../tcr_epitope_datasets/mutational_scan_datasets',
                              '/train_test_data_folds.xlsx'))

# Discard index peptides with incomplete Hamming scan and MHCII
tcr_data <- tcr_data[!tcr_data$index_peptide %in% c('YMDGTMSQV','VPSVWRSSL',
                                                    "FRDYVDRFYKTLRAEQASQE",
                                                    "LPVPGVLLKEFTVSGNILTI",
                                                    "ASQKRPSQRSK",
                                                   'FEAQKAKANKAVD'),]

# Rename TCR activation classes
tcr_data$activation[tcr_data$activation==0] <- 'No'
tcr_data$activation[tcr_data$activation==1] <- 'Weak'
tcr_data$activation[tcr_data$activation==2] <- 'Strong'

# Remove WT peptides
tcr_data <- tcr_data[tcr_data$mutation_position != 0,]

# Count number of data entries for different pMHC binding and TCR activation
counts <- tcr_data %>% dplyr::count(mhc, BindLevel, activation)

############
## Plot ####
############

p <- ggplot(tcr_data,aes(y=factor(activation,levels=c('No','Weak','Strong')),
                   x=factor(BindLevel,levels=c('NB','WB','SB')))) + 
  geom_jitter(aes(color=factor(mutation_position)),size=0.7,alpha=0.75)+
  
#  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(15)) +
                               
                                         labs(x="Predicted MHC binding",
                                              y="TCR activation",
                                              color="")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        axis.text.y = element_text(size=25)) +
  facet_wrap(vars(factor(mhc)))+
  theme(strip.text = element_text(
    size = 15, color = "black")) +
  theme(strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "cm")))+
  theme(legend.text = element_text(size=20)) +
  labs(color="Peptide mutation position") +
  theme(legend.title = element_text(angle = -90,size=20))

ggsave(plot=p, "figS3.pdf",
       width=15, height=10)




