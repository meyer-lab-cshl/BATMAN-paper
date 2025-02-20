rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)


# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# Subset selected TCRs and columns
all_TCRs <- c('a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                   'A6','T1','FLT3DY','A23','TCR2-T','868Z11')

###############################
## AL with BATMAN #############
###############################

# open empty df for storing AUC
auc_all <- data.frame(matrix(ncol=3,nrow=0, 
                                     dimnames=list(NULL, c("tcr", "al_round", 
                                                           "auc"))))
#Average over random realizations
for (tcr in all_TCRs){
auc_files <- list.files(path = "../4a-b/BATMAN_AL_aucs/",
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

auc_all <- rbind(auc_all,auc_mean)
}


# Subset for first AUC rounds only
auc_active <- auc_all[auc_all$al_round==1,2:3] 
colnames(auc_active) <- c('auc_al','tcr')

auc_active_loo <- auc_all[auc_all$al_round==0,2:3] 
colnames(auc_active_loo) <- c('auc_loo','tcr')

auc_active <- merge(auc_active,auc_active_loo,by="tcr")

# Get median mutant2index distance
distances <- read.csv('median_mutant_distances.csv')
colnames(distances) <- c('tcr','d_mean','d_median','d_std')

# Merge dfs
all_data <- merge(auc_active,distances,by="tcr")

###################################
######## Plotting: MHCI ###########
###################################

p <- ggplot(all_data, aes(y=auc_al-auc_loo,x=d_median,label=tcr)) + 
  geom_smooth(data = all_data, aes(y=auc_al-auc_loo,x=d_median),
              colour="gray",
              method = "lm") +
   geom_point(aes(color=factor(tcr,levels = all_TCRs)), 
             alpha=0.85,show.legend = FALSE,
             size=5) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)[
    c(1,3,4,5,6,7,8,9,10,11,12,14)]) +
  geom_text_repel(size = 4,nudge_y=0.02,nudge_x=0.02)+
  labs(y="AUC increase by AL with 9 peptides",
       x="Median peptide-to-index distance",
       color="") +
  theme_cowplot() +
  geom_hline(yintercept=0, linetype='dashed', col = '#966919') +
  geom_vline(xintercept=0.8, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
#  ylim(0.6,0.95) +
  theme(legend.position="right")
#  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "al_aucs_vs_distance_for_4c.pdf",
       width=5, height=5)


# Save raw data
write.csv(all_data, "raw_data_fig_4c.csv")

