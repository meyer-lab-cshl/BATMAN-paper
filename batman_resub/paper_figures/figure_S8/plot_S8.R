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


########## Plot weight profiles #################################
# selected TCRs for plotting
all_TCRs <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
              'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11')

# Read and add to weights data

### Inferred Symmetric
weights_all <- read.csv('inferred_weights/inferred_weights_symm.csv')
colnames(weights_all) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')

weights_all <- pivot_longer(weights_all, 
                             cols = c('P1','P2','P3','P4','P5','P6','P7','P8','P9'))

colnames(weights_all) <- c('tcr','position','weight')
weights_all$aa_matrix <- 'Inferred Symmetric'

for (aa in c('full','Hamming','BLOSUM100','PAM10','DAYHOFF','GONNET')){
  
  weights <- read.csv(paste0('inferred_weights/inferred_weights_',aa,'.csv'))
  colnames(weights) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8','P9')
  
  weights <- pivot_longer(weights, 
                              cols = c('P1','P2','P3','P4','P5',
                                       'P6','P7','P8','P9'))
  
  colnames(weights) <- c('tcr','position','weight')
  weights$aa_matrix <- aa
  
  weights_all <- rbind(weights_all,weights)
}

# Normalize to 1 for each TCR and method
weights_all <- weights_all %>% 
  group_by(tcr,aa_matrix) %>% 
  mutate(w_max=max(weight))

weights_all$weight <- weights_all$weight/weights_all$w_max


weights_all$aa_matrix[weights_all$aa_matrix=='full'] <- 'Inferred Full'

# Take out P from position
weights_all$position <- gsub("P","",as.character(weights_all$position)) 

###################################
######## Plotting #################
###################################


p <- ggplot(weights_all, aes(x=factor(position), y=weight, 
                            fill=factor(tcr,levels=all_TCRs))) +
    geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
           alpha=1) +
     labs(y="Positional weight",
       x="Position",
       color="") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)[1:14]) +
    ylim(0, 1) +
  
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        legend.text =element_text(size=20)) +
  facet_grid(vars(factor(tcr,levels=all_TCRs)),vars(factor(aa_matrix,
                                                           levels = c(c('Inferred Full',
                                                                        'Inferred Symmetric',
                                'Hamming','BLOSUM100','PAM10','DAYHOFF','GONNET')))))+
  theme(strip.text = element_text(
    size = 12, color = "black"))

# Save plot
ggsave(plot=p, "figS8.pdf",
         width=15, height=15)
# 
# # Save raw data
write.csv(weights_all, "raw_data_fig_S8.csv")



 
