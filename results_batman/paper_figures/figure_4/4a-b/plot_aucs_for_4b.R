rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# Subset selected TCRs and columns
all_TCRs <- c('a3a','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                   'A6','T1','FLT3DY','A23','TCR2-T','868Z11')

###############################
## AL with BATMAN #############
###############################

# open empty df for storing AUC
auc_all <- data.frame(matrix(ncol=4,nrow=0, 
                                     dimnames=list(NULL, c("tcr", "npeptide","nstep", 
                                                           "auc"))))
npeptide_list <- c(4,8,12,2,4,6,1,2,3)
nstep_list <- c(1,1,1,2,2,2,4,4,4)

#average AL AUCs over random seeds
for (al_params in 1:9){
  npeptide <- npeptide_list[al_params]
  nstep <- nstep_list[al_params]
  
  # List files
  auc_files <- list.files(path = "BATMAN_AL_nstep_npeptide/",
                          pattern = paste0("auc_",npeptide,"_",nstep,
                                           "_.*.csv"),
                          full.names = TRUE)
  
  auc_df <- data.frame()
  for (file in auc_files){
    # Read data
    file_data <- read.csv(file)[,c(1,1+nstep)] #Rad only AUC from last step
    auc_df <- rbind(auc_df, file_data)
  }
  
  colnames(auc_df) <- c("tcr","auc")
  
  # Average over random seed
  # Specify data frame
  auc_mean <- auc_df %>%
    # Specify group indicator, column, function
    group_by(tcr) %>%
    # Calculate the mean of each group
    summarise_at(vars(auc),
                 list(auc = mean))
  
  auc_mean$npeptide <- npeptide*nstep*9 #Total number of peptides sampled
  auc_mean$nstep <- nstep
  auc_all <- rbind(auc_all,auc_mean)
}



# Read LOO TCR data ############################
# open empty df for storing AUC
auc_loo <- data.frame(matrix(ncol=3,nrow=0, 
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
  auc_loo <- rbind(auc_loo,auc_mean)
}


# Subset for LOO AUC only
auc_loo <- auc_loo[auc_loo$al_round==0,2:3] 
colnames(auc_loo) <- c('auc','tcr')
auc_loo$npeptide <- 'LOO TCR'
auc_loo$nstep <- 0

auc_all <- rbind(auc_all,auc_loo)

# Load within-TCR no-mhc 5 fold data ##########
within_aucs <- read.csv('../../figure_2/2d/BATMAN_no_mhc_AUCs.csv')[,1:2]
colnames(within_aucs) <- c('tcr','auc')

within_aucs <- within_aucs[within_aucs$tcr %in% all_TCRs,]

within_aucs$npeptide <- 'Within TCR'
within_aucs$nstep <- 0

auc_all <- rbind(auc_all,within_aucs)


###################################
######## Plotting       ###########
###################################

p <- ggplot(auc_all, aes(y=auc,x=factor(npeptide,levels=c("LOO TCR","36",
                                                          "72","108",
                                                          "Within TCR")))) + 
  geom_boxplot(aes(color = factor(nstep),
    fill = factor(nstep)),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1,
               position=position_dodge2(1,preserve = "single")) +
  scale_color_manual(values=c('#AA4A44',"#FFD6D7",'#FF999C','#FF7074')) +
  scale_fill_manual(values=c("#FFD6D7","#FFEAEB","#FFEAEB","#FFEAEB")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels = all_TCRs), 
                 group=factor(nstep)), 
             alpha=0.85,
             size=1,
             position = position_jitterdodge(jitter.width=0.7,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)[
    c(1,3,4,5,6,7,8,9,10,11,12,14)]) +
  labs(y="TCR activation prediction AUC",
       x="Peptides sampled",
       color="") +
  theme_cowplot() +
#  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=20,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20))
#ylim(0.46,0.93)
#  theme(legend.position="right") +
#  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "al_aucs_for_4b.pdf",
       width=6, height=6)


# Save raw data
write.csv(auc_all, "raw_data_fig_4b.csv")

