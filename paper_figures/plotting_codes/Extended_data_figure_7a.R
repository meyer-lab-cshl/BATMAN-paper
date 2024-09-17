# # Code to plot scores from different TCR-pMHC methods

rm(list=ls())#Clear environment
# Load required libraries
library(readxl)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(forcats)
library(viridis)
library(loo)
library(pROC)
library(tidyr)

#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load peptide information
peptide_folds <- read.csv('fig_2_training_files/full_training_data.csv')
peptide_folds$activation <- peptide_folds$activation + 1
peptide_folds$fold <- peptide_folds$fold + 1 #Python to R

# List of 3 class 9-mer-binding TCRs
tcr_names <- c("18A2","NYE-S1","868Z11","TCR3","TCR6","TCR7","A6",
               "T1","FLT3DY","A23","TCR2-T")

# Load peptide data for selected TCRs

peptide_folds <- peptide_folds[peptide_folds$tcr %in% tcr_names,]


# peptide scores from BATMAN within-TCR and cross-TCR for full and symm matrices
peptide_score <- array(NA,dim = c(dim(peptide_folds)[1],4))

# BATMAN scores
output_data <- read.csv('fig_2_training_files/BATMAN_scores_full_matrix.csv')

peptide_score[,1] <- 
  as.double(output_data$within_tcr_score[match(peptide_folds$peptide,
                                               output_data$peptide)])
peptide_score[,2] <- 
  as.double(output_data$cross_tcr_score[match(peptide_folds$peptide,
                                              output_data$peptide)])

output_data <- read.csv('fig_2_training_files/BATMAN_scores_symm_matrix.csv')

peptide_score[,3] <- 
  as.double(output_data$within_tcr_score[match(peptide_folds$peptide,
                                               output_data$peptide)])
peptide_score[,4] <- 
  as.double(output_data$cross_tcr_score[match(peptide_folds$peptide,
                                              output_data$peptide)])

method_name <- c('BATMAN Full','BATMAN Symmetric')

# Calculate AUCs for within and cross-TCR modes
# Create and shape the full dataFrame for plotting
plot_data <- data.frame(matrix(nrow=0,ncol=4))
colnames(plot_data) <- c('TCR','method','mode','score')


for (tcr in tcr_names){
  training_folds_tcr <- peptide_folds[peptide_folds$tcr == tcr,]
  peptide_score_tcr <- peptide_score[peptide_folds$tcr == tcr,]
  
  # BATMAN scores
  # Get AUC for 5 folds
  AUCs <- array(NA,dim = c(3,5,2))
  for (fold in 1:5){
    # Full matrix
    ROC<-multiclass.roc(training_folds_tcr$activation[training_folds_tcr$fold==fold],
                        peptide_score_tcr[training_folds_tcr$fold==fold,
                                          1],direction=">")
    
    AUCs[,fold,1] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
    
    # Symm matrix
    ROC<-multiclass.roc(training_folds_tcr$activation[training_folds_tcr$fold==fold],
                        peptide_score_tcr[training_folds_tcr$fold==fold,
                                          3],direction=">")
    
    AUCs[,fold,2] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  }
  # Within-TCR scores
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN Full',
                         mode = 'Within TCR',
                         score = mean(AUCs[,,1]))
  plot_data <- rbind(plot_data,new_data)
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN Symmetric',
                         mode = 'Within TCR',
                         score = mean(AUCs[,,2]))
  plot_data <- rbind(plot_data,new_data)
  # Full TCR data
  ROC<-multiclass.roc(training_folds_tcr$activation,
                      peptide_score_tcr[,2],direction=">")
  AUCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN Full',
                         mode = 'Cross TCR',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
  
  ROC<-multiclass.roc(training_folds_tcr$activation,
                      peptide_score_tcr[,4],direction=">")
  AUCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN Symmetric',
                         mode = 'Cross TCR',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
}

###################
###### plotting ###
###################

# plot
p <- ggplot(plot_data, 
            aes(x=factor(method, levels = c("BATMAN Full", "BATMAN Symmetric")),
                y=score)) +
  geom_boxplot(aes(color=factor(mode, levels = c('Within TCR','Cross TCR'))),
               outlier.colour = NA,
               show.legend = FALSE, alpha=1,lwd=1
  ) +
  scale_color_manual(values=c('#DC143C','#007FFF')) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=TCR, 
                 group=factor(mode,levels = c('Within TCR','Cross TCR'))), 
             alpha=0.9,
             size=5, show.legend = FALSE,
             position = position_jitterdodge(jitter.width=0.1)) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
                                         labs(x="Method",
                                              y="AUC score",
                                              color="") +
  ylim(NA, 0.89) +
  theme_cowplot() +
  theme(legend.text = element_text(size=25)) +
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(size=25,angle=45,vjust=1),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=25, vjust = -1),
        axis.title.y = element_text(size=25))

ggsave(plot=p, "edfig7a.pdf",
       width=5, height=10)













