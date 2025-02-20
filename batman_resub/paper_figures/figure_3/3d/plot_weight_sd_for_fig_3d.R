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

###########  Read Ova TCR weight Data ##########

# Read mean inferred weights
weights_ova <- read.csv('inferred_weights_SIINFEKL.csv')

colnames(weights_ova) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8')
weights_ova <- pivot_longer(weights_ova, 
                            cols = c('P1','P2','P3','P4','P5','P6','P7','P8'))

colnames(weights_ova) <- c('tcr','position','weight')

# Read HDR of inferred weights
weights_sd_ova <- read.csv('inferred_weights_sd_SIINFEKL.csv')

colnames(weights_sd_ova) <- c('tcr','P1','P2','P3','P4','P5','P6','P7','P8')
weights_sd_ova <- pivot_longer(weights_sd_ova, 
                            cols = c('P1','P2','P3','P4','P5','P6','P7','P8'))

colnames(weights_sd_ova) <- c('tcr','position','weight_sd')

# Merge data 
weights_data <- merge(weights_ova,weights_sd_ova,by=c('tcr','position'))

# Normalize weights and HDR to max for each TCR
# for (tcr in  unique(weights_data$tcr)){
#   max_weight <- 1#max(weights_data$weight[weights_data$tcr==tcr])
#   weights_data$weight[weights_data$tcr==tcr] <- weights_data$weight[weights_data$tcr==tcr]/max_weight
#   weights_data$weight_sd[weights_data$tcr==tcr] <- weights_data$weight_sd[weights_data$tcr==tcr]/max_weight
#}

# Add ratio column
weights_data$ratio <- weights_data$weight_sd/weights_data$weight

# Add TCR education label
weights_data$tcr_type <- "Naive TCRs" #Initialize

educated_tcr_list <- weights_data %>% 
  filter(grepl('Ed', tcr) | 
           tcr == 'OT1')

weights_data$tcr_type[weights_data$tcr %in% unique(educated_tcr_list$tcr)] <- "Educated TCRs"

# Use only non-anchor positions
weights_data <- weights_data[!weights_data$position %in% c('P2','P3','P5','P8'),]

# Rename positions by residue
weights_data$position <- replace(weights_data$position,
                                 weights_data$position %in% c('P1','P4','P6','P7'),
                                 c('S1','N4','E6','K7'))


###########################
#### Plot data ############
###########################
p <- ggplot(weights_data, aes(y=weight_sd,x=tcr_type)) + 
  geom_boxplot(aes(color = tcr_type,
                   fill = tcr_type),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1) +
  scale_color_manual(values=c('#1f78b4'	,'#e41a1c')) +
  scale_fill_manual(values=c("#a6cee3","#fbb4ae")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(position,levels= c('S1','N4','E6','K7')),
                 group=factor(position,levels= c('S1','N4','E6','K7'))),
             alpha=0.85,
             size=3,
             position = position_jitterdodge(jitter.width=0.1)) +
  # geom_point(aes(color=tcr_type), 
  #                          alpha=0.85,
  #                          size=3) +
  labs(y="Posterior HDI length",
       x=" ",
       color="") +
  theme_cowplot() +
  #  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=20,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20)) +
  ylim(0,1) +
  theme(legend.position="right") 

# Save plot
ggsave(plot=p, "weight_posterior_for_3d.pdf",
       width=6, height=7)

# Save raw data
write.csv(weights_data, "raw_data_fig_3d.csv")













