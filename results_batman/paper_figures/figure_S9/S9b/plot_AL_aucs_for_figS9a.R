rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# Subset selected TCRs and columns
all_TCRs <- c('a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
              'A6','T1','FLT3DY','A23','TCR2-T','868Z11')

# Collect AUCs for different AL modes

# AL: No random
auc_all <- read.csv('auc_AL_no_random_all_TCRs.csv')
colnames(auc_all) <- c("al_round",all_TCRs)
auc_all <- pivot_longer(auc_all,
                     cols = all_of(all_TCRs),
                     names_to = 'tcr',
                     values_to = 'auc')
auc_all$al_type <- 'No random'
auc_all$al_round[auc_all$al_round==0] <- 'LOO TCR'

# AL: Even steps random

#Average over random realizations
for (tcr in all_TCRs){
  auc_files <- list.files(path = "../../figure_4/4a-b/BATMAN_AL_aucs/",
                          pattern = paste0("batman_AL_auc_",
                                           tcr,"_.*.csv"),
                          full.names = TRUE)
  
  auc_df <- data.frame()
  for (file in auc_files){
    # Read data
    file_data <- read.csv(file) 
    auc_df <- rbind(auc_df, file_data)
  }
  
  colnames(auc_df) <- c("al_round","auc")
  
  # Average over random seed
  # Specify data frame
  auc_mean <- auc_df %>%
    # Specify group indicator, column, function
    group_by(al_round) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(auc),
                 list(auc = mean))
  
  auc_mean$tcr <- tcr
  auc_mean$al_type <- 'Even steps random'
  
  auc_all <- rbind(auc_all,auc_mean)
}

# Load within-TCR no-mhc 5 fold data
within_aucs <- read.csv('../../figure_2/2d/BATMAN_no_mhc_AUCs.csv')[,1:2]
colnames(within_aucs) <- c('tcr','auc')

within_aucs <- within_aucs[within_aucs$tcr %in% all_TCRs,]

within_aucs$al_round <- 'Within TCR'
within_aucs$al_type <- "No AL"
auc_all <- rbind(auc_all,within_aucs)

# AL: All Random

# load BATMAN random learning AUCs
for (tcr in all_TCRs){
  auc_files <- list.files(path = "../../figure_4/4a-b/BATMAN_RL_aucs/",
                          pattern = paste0("batman_RL_auc_",
                                           tcr,"_.*.csv"),
                          full.names = TRUE)
  
  auc_df <- data.frame()
  for (file in auc_files){
    # Read data
    file_data <- read.csv(file) 
    auc_df <- rbind(auc_df, file_data)
  }
  
  colnames(auc_df) <- c("al_round","auc")
  
  # Average over random seed
  # Specify data frame
  auc_mean <- auc_df %>%
    # Specify group indicator, column, function
    group_by(al_round) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(auc),
                 list(auc = mean))
  
  auc_mean$tcr <- tcr
  # discard LOO TCR data (already collected above)
  auc_mean <- auc_mean[!auc_mean$al_round==0,]
  
  auc_mean$al_type <- "All random"
  
  # Add to main datafile
  auc_all <- rbind(auc_all,auc_mean)
}

# Remove multiple occurrences of LOO TCR
auc_all <- auc_all[auc_all$al_round!=0,]




###################################
######## Plotting       ###########
###################################

p <- ggplot(auc_all, aes(y=auc,x=factor(al_round,levels=c("LOO TCR","1","2","3",
                                                          "4","5","6","7","8","9","10",
                                                          "Within TCR")))) + 
  geom_boxplot(aes(color = al_type),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1,
               position=position_dodge2(1,preserve = "single")) +
#  scale_color_manual(values=c("#FAA0A0",'#89CFF0','#AA4A44')) +
#  scale_fill_manual(show.legend = TRUE) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels = all_TCRs), 
                 group=al_type), 
             alpha=0.85,
             size=1,
             position = position_jitterdodge(jitter.width=0.7,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)[
    c(1,3,4,5,6,7,8,9,10,11,12,14)]) +
  labs(y="TCR activation prediction AUC",
       x="Experimental round",
       color="") +
  theme_cowplot() +
  #  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15))
#  ylim(0.46,0.93)
#  theme(legend.position="right") +
#  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "figS9.pdf",
       width=7, height=6)


# Save raw data
write.csv(auc_all, "raw_data_fig_S9b.csv")

