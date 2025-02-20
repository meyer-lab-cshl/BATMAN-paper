rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(viridis)
library(forcats)
library(tidyr)
library(RColorBrewer)
library(patchwork)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

###############################
## BATMAN regression with AA matrices#####
###############################
# List all files
corr_files <- list.files(path = "S5a/Corrs_TCR_AA_matrix/",
                         pattern = "corr_.*.csv",
                         full.names = TRUE)
# open empty df for storing
corr_df <- data.frame()
for (file in corr_files){
  # Read data
  file_data <- read.csv(file) 
  corr_df <- rbind(corr_df, file_data)
}

colnames(corr_df) <- c("X","corr",'aa','tcr')

# Make AA matrix names in all caps
corr_df$aa <- toupper(corr_df$aa)


###################################
######## Plotting P1###############
###################################

p1 <- ggplot(corr_df,aes(y=corr,x=fct_inorder((factor(aa))))) + 
  geom_boxplot(outlier.colour = NA,
               show.legend = FALSE) +
  geom_jitter(aes(color=tcr,group=aa),
              alpha=0.5, width=0.1, height=0,size=0.1, show.legend = FALSE) +
  #  scale_colour_manual(values = palette(viridis(151))) +
  labs(y="Average Spearman r",
       x="",
       color="") +
  #  xlim(NA, 0.9) +
  theme_cowplot() +
  geom_hline(aes(yintercept=0.25),colour="gray")+
  geom_hline(aes(yintercept=0.50),colour="gray")+
  geom_hline(aes(yintercept=0.75),colour="gray")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=5),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank())

##############################################
## BATMAN classification with AA matrices#####
##############################################
# List all files
auc_files <- list.files(path = "S5b-c/aucs_TCR_AA_matrix/",
                        pattern = "auc_.*.csv",
                        full.names = TRUE)
# open empty df for storing
auc_df <- data.frame()
for (file in auc_files){
  # Read data
  file_data <- read.csv(file) 
  auc_df <- rbind(auc_df, file_data)
}

auc_df$mean_auc <- (auc_df$WS + auc_df$NS +auc_df$NW)/3

# Make AA matrix names in all caps
auc_df$aa_matrix <- toupper(auc_df$aa_matrix)

###################################
######## Plotting P2###############
###################################

p2 <- ggplot(auc_df,aes(y=mean_auc,x=fct_inorder(factor(aa_matrix)))) + 
  geom_boxplot(outlier.colour = NA,
               show.legend = FALSE) +
  geom_jitter(aes(color=tcr,group=aa_matrix),
              alpha=0.5, width=0.1, height=0,size=0.1, show.legend = FALSE) +
  #  scale_colour_manual(values = palette(viridis(151))) +
  labs(y="Average AUC",
       x="",
       color="") +
  #  xlim(NA, 0.9) +
  theme_cowplot() +
  geom_hline(aes(yintercept=0.4),colour="gray")+
  geom_hline(aes(yintercept=0.6),colour="gray")+
  geom_hline(aes(yintercept=0.8),colour="gray")+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=5),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank())

#### Combine all plots #####
multi <- (p1 + p2) + 
  plot_layout(ncol=1) +
  plot_annotation(tag_levels = 'A') #add figure labels

ggsave(plot=multi, "figS5.pdf",width=7, height=6)

# Save raw data
write.csv(corr_df, "raw_data_fig_S5a.csv")
write.csv(auc_df, "raw_data_fig_S5b.csv")
