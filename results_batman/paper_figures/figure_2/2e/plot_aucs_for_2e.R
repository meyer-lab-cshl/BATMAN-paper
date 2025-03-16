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
TCR_list_9mer <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                   'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11')

TCR_list_10mer <- c('TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va')

TCR_list_mhcii <- c('TCR-F5','TCR-3598-2','MBP-TCR','B3K508')


all_TCRs <- c(TCR_list_9mer,TCR_list_10mer)

###############################
## Pooled BATMAN ##############
###############################

# open empty df for storing AUC
auc_all <- data.frame(matrix(ncol=3,nrow=0, 
                                     dimnames=list(NULL, c("tcr", "fpeptide", 
                                                           "auc"))))

for (fpeptide in c(1,2,3,4,5,7,9)){
auc_files <- list.files(path = "pooled_batman_outputs/",
                           pattern = paste0("aucs_batman_pooled_.*_",
                                            fpeptide,"_.*.csv"),
                           full.names = TRUE)

auc_df <- data.frame()
for (file in auc_files){
  # Read data
  file_data <- read.csv(file) 
  auc_df <- rbind(auc_df, file_data)
}

colnames(auc_df) <- c("tcr","auc")

# Average over folds and random seed
# Specify data frame
auc_mean <- auc_df %>%
  # Specify group indicator, column, function
  group_by(tcr) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(auc),
               list(auc = mean))

auc_mean$fpeptide <- fpeptide/10

auc_all <- rbind(auc_all,auc_mean)
}

# Load within-TCR 5 fold data
within_aucs <- read.csv('../2d/BATMAN_AUCs.csv')[,1:2]
colnames(within_aucs) <- c('tcr','auc')

within_aucs <- within_aucs[within_aucs$tcr %in% all_TCRs,]

within_aucs$fpeptide <- 1
auc_all <- rbind(auc_all,within_aucs)

auc_all$pooling <- 'BATMAN pooled'

###############################
## Unpooled BATMAN ############
###############################
# open empty df for storing AUC
auc_unpooled <- data.frame(matrix(ncol=3,nrow=0, 
                             dimnames=list(NULL, c("tcr", "fpeptide", 
                                                   "auc"))))

for (fpeptide in c(1,2,3,4,5,7,9,10)){
  auc_files <- list.files(path = "unpooled_batman_outputs/",
                          pattern = paste0("aucs_batman_unpooled_.*_",
                                           fpeptide,"_.*.csv"),
                          full.names = TRUE)
#  print(length(auc_files))
  if (length(auc_files)>0){
  auc_df <- data.frame()
  for (file in auc_files){
    # Read data
    file_data <- read.csv(file) 
    auc_df <- rbind(auc_df, file_data)
  }
  
  colnames(auc_df) <- c("tcr","auc")
  
  # Average over folds and random seed
  # Specify data frame
  auc_mean <- auc_df %>%
    # Specify group indicator, column, function
    group_by(tcr) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(auc),
                 list(auc = mean))
  
  auc_mean$fpeptide <- fpeptide/10
  
  auc_unpooled <- rbind(auc_unpooled,auc_mean)
  }
}

auc_mhcii <- auc_unpooled[!auc_unpooled$tcr %in% all_TCRs,]
auc_unpooled <- auc_unpooled[auc_unpooled$tcr %in% all_TCRs,] #discard MHCII TCRs

auc_unpooled$pooling <- 'BATMAN unpooled'
auc_mhcii$pooling <- 'BATMAN unpooled'

auc_all <- rbind(auc_all,auc_unpooled)

###############################
############ pTEAM AUCs #######
###############################
auc_pteam <- data.frame(matrix(ncol=3,nrow=0, 
                                  dimnames=list(NULL, c("tcr", "fpeptide", 
                                                        "auc"))))

for (fpeptide in c(1,2,3,4,5,7,9,10)){
  auc_files <- list.files(path = "pteam_outputs/",
                          pattern = paste0("aucs_pteam_.*_",
                                           fpeptide,"_.*.csv"),
                          full.names = TRUE)
  #  print(length(auc_files))
  if (length(auc_files)>0){
    auc_df <- data.frame()
    for (file in auc_files){
      # Read data
      file_data <- read.csv(file)[,2:3] 
      auc_df <- rbind(auc_df, file_data)
    }
    
    colnames(auc_df) <- c("tcr","auc")
    
    # Average over folds and random seed
    # Specify data frame
    auc_mean <- auc_df %>%
      # Specify group indicator, column, function
      group_by(tcr) %>%
      # Calculate the mean of the "Frequency" column for each group
      summarise_at(vars(auc),
                   list(auc = mean))
    
    auc_mean$fpeptide <- fpeptide/10
    
    auc_pteam <- rbind(auc_pteam,auc_mean)
  }
}

auc_pteam_mhcii <- auc_pteam[!auc_pteam$tcr %in% all_TCRs,
                                   ] #MHCII TCRs
auc_pteam <- auc_pteam[auc_pteam$tcr %in% all_TCRs,] #discard MHCII TCRs

auc_pteam$pooling <- 'pTEAM'
auc_pteam_mhcii$pooling <- 'pTEAM'

auc_all <- rbind(auc_all,auc_pteam)


# Load full 5-fold training aucs
within_aucs <- read.csv('../2d/pTEAM_AUCs.csv')
within_aucs <- within_aucs[within_aucs$type=="pTEAM_within",2:3]
colnames(within_aucs) <- c('tcr','auc')

#separate MHCI and II
within_aucs_mhcii <- within_aucs[!within_aucs$tcr %in% all_TCRs,]
within_aucs <- within_aucs[within_aucs$tcr %in% all_TCRs,]

within_aucs$fpeptide <- 1
within_aucs$pooling <- 'pTEAM'

within_aucs_mhcii$fpeptide <- 1
within_aucs_mhcii$pooling <- 'pTEAM'

auc_all <- rbind(auc_all,within_aucs)

auc_all_mhcii <- rbind(auc_mhcii,auc_pteam_mhcii,within_aucs_mhcii)


###################################
######## Plotting: MHCI ###########
###################################

p <- ggplot(auc_all, aes(y=auc,x=factor(fpeptide),
)) + 
  geom_boxplot(aes(color = pooling,
    fill = pooling),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1,
               position=position_dodge2(0.8,preserve = "total")) +
  scale_color_manual(values=c('#AA4A44'	,'black','#0066b2')) +
  scale_fill_manual(values=c("#FFD6D7","#C0C0C0","#F0FFFF")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels = all_TCRs), 
                 group=pooling), 
             alpha=0.85,
             size=2,
             position = position_jitterdodge(jitter.width=0.2,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)) +
  labs(y="TCR activation prediction AUC",
       x="Fraction of mutant peptides sampled",
       color="") +
  theme_cowplot() +
#  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
  ylim(0.6,0.95) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "methods_aucs_for_2e.pdf",
       width=9, height=4)


###################################
######## Plotting: MHCII ###########
###################################

p <- ggplot(auc_all_mhcii, aes(y=auc,x=factor(fpeptide),
)) + 
  geom_boxplot(aes(color = pooling,
                   fill = pooling),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1,
               position=position_dodge2(0.8,preserve = "total")) +
  scale_color_manual(values=c('black','#0066b2')) +
  scale_fill_manual(values=c("#C0C0C0","#F0FFFF")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels=TCR_list_mhcii), 
                 group=pooling), 
             alpha=0.85,
             size=2,
             position = position_jitterdodge(jitter.width=0.2,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)[19:22]) +
  labs(y="TCR activation prediction AUC",
       x="Fraction of mutant peptides sampled",
       color="") +
  theme_cowplot() +
  #  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
  ylim(0.6,0.95) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "methods_aucs_for_2e_mhcii.pdf",
       width=7, height=4)




# Save raw data
write.csv(auc_all, "raw_data_fig_2e.csv")
write.csv(auc_all_mhcii, "raw_data_fig_2e_mhcii.csv")










