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

## Fig 2a       #############
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


# peptide scores from 7 methods + BATMAN within-TCR and cross-TCR
peptide_score <- array(NA,dim = c(dim(peptide_folds)[1],9))

# Load TCR-pMHC scores from 7 other methods
# We added a negative sign before ranks to use them as score
#

#ERGO II
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "ergo//ergo_output_mcpas.csv"))

peptide_score[,1] <- as.double(output_data$Score[match(peptide_folds$peptide,
                                                    output_data$Peptide)])



#epiTCR
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "epitcr//output.csv"))

peptide_score[,2] <- as.double(output_data$predict_proba[match(peptide_folds$peptide,
                                                       output_data$epitope)])

#ImRex
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "imrex//output-prediction.csv"))

peptide_score[,3] <- as.double(output_data$prediction_score[match(peptide_folds$peptide,
                                                  output_data$antigen.epitope)])

#pMTnet
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "pmtnet//pMTnet_prediction.csv"))

peptide_score[,4] <- -as.double(output_data$Rank[match(peptide_folds$peptide,
                                                          output_data$Antigen)])

#IEDB
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/iedb/",pattern = ".csv")
file_names <- paste0(rep("../codes/misc_methods/iedb/",
                         length(file_names)),file_names)
output_data <- read_csv(file_names)

peptide_score[,5] <- as.double(output_data$score[match(peptide_folds$peptide,
                                                      output_data$peptide)])

#NetTepi
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/nettepi/",pattern = "NetTepi_out_")
file_names <- paste0(rep("../codes/misc_methods/nettepi/",
                         length(file_names)),file_names)
output_data <- data.frame()
for (file in file_names){
  output_data <- rbind(output_data,read.table(file))
}

peptide_score[,6] <- as.double(output_data$V8[match(peptide_folds$peptide,
                                                       output_data$V2)])

#PRIME 2.0
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/prime/",pattern = "result_")
file_names <- paste0(rep("../codes/misc_methods/prime/",
                         length(file_names)),file_names)
output_data <- data.frame()
for (file in file_names){
  output_data <- rbind(output_data,read.table(file))
}

peptide_score[,7] <- as.double(output_data$V3[match(peptide_folds$peptide,
                                                    output_data$V1)])

# BATMAN scores
output_data <- read.csv('fig_2_training_files/BATMAN_scores_full_matrix.csv')

peptide_score[,8] <- 
  as.double(output_data$within_tcr_score[match(peptide_folds$peptide,
                                                    output_data$peptide)])
peptide_score[,9] <- 
  as.double(output_data$cross_tcr_score[match(peptide_folds$peptide,
                                                          output_data$peptide)])

method_name <- c('ERGO II','epiTCR','ImRex','pMTnet','IEDB','NetTepi',
                 'PRIME 2.0','BATMAN')

# Calculate AUCs for within and cross-TCR modes
# Create and shape the full dataFrame for plotting
plot_data <- data.frame(matrix(nrow=0,ncol=4))
colnames(plot_data) <- c('TCR','method','mode','score')


for (tcr in tcr_names){
  training_folds_tcr <- peptide_folds[peptide_folds$tcr == tcr,]
  peptide_score_tcr <- peptide_score[peptide_folds$tcr == tcr,]
  
  # BATMAN scores
  # Get AUC for 5 folds
  AUCs <- array(NA,dim = c(3,5))
  for (fold in 1:5){
    ROC<-multiclass.roc(training_folds_tcr$activation[training_folds_tcr$fold==fold],
                        peptide_score_tcr[training_folds_tcr$fold==fold,
                                          8],direction=">")
    
    AUCs[,fold] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  }
  # Within-TCR score
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN',
                         mode = 'Within TCR',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
  # Full TCR data
  ROC<-multiclass.roc(training_folds_tcr$activation,
                      peptide_score_tcr[,9],direction=">")
  AUCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  new_data <- data.frame(TCR = tcr,
                         method = 'BATMAN',
                         mode = 'Cross TCR',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
  
  for (method_index in 1:7){
    # Check if all scores are NA
    if (all(is.na(peptide_score_tcr[,method_index]))==FALSE){

    # Get AUC for 5 folds
    AUCs <- array(NA,dim = c(3,5))
    for (fold in 1:5){
    ROC<-multiclass.roc(training_folds_tcr$activation[training_folds_tcr$fold==fold],
                        peptide_score_tcr[training_folds_tcr$fold==fold,
                                      method_index],direction="<")

    AUCs[,fold] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
    }
    # Within-TCR score
    new_data <- data.frame(TCR = tcr,
                           method = method_name[method_index],
                           mode = 'Within TCR',
                           score = mean(AUCs))
    plot_data <- rbind(plot_data,new_data)
    # Full TCR data
    ROC<-multiclass.roc(training_folds_tcr$activation,
                        peptide_score_tcr[,method_index],direction="<")
    AUCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
    new_data <- data.frame(TCR = tcr,
                           method = method_name[method_index],
                           mode = 'Cross TCR',
                           score = mean(AUCs))
    plot_data <- rbind(plot_data,new_data)
    }
  }
}

# Read saved pTEAM scores
all_scores <- read.csv("../data/within_tcr1.csv")
plot_data <- rbind(plot_data,
                   data.frame(score=array(data.matrix(all_scores[all_scores$Row=='Atchley+RF',
                                                           2:12])),
                              TCR = tcr_names,
                              method = 'pTEAM',
                              mode = 'Within TCR'))

all_scores <- read.csv("../data/cross_tcr1.csv")
plot_data <- rbind(plot_data,
                   data.frame(score=array(data.matrix(all_scores[all_scores$Row=='Atchley+RF',
                                                           2:12])),
                              TCR = tcr_names,
                              method = 'pTEAM',
                              mode = 'Cross TCR'))





###################
###### plotting ###
###################

# plot
p <- ggplot(plot_data, 
            aes(x=factor(method, levels = c("BATMAN","PRIME 2.0","NetTepi",
                                            "IEDB","ImRex","pMTnet","epiTCR",
                                            "pTEAM","ERGO II")), y=score)) +
     geom_boxplot(aes(color=factor(mode, levels = c('Within TCR','Cross TCR'))),
                  outlier.colour = NA,
                  show.legend = FALSE, alpha=1,lwd=1
                  ) +
  scale_color_manual(values=c('#DC143C','#007FFF')) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=TCR, 
                 group=factor(mode,levels = c('Within TCR','Cross TCR'))), 
             alpha=0.9,
             size=5,
             position = position_jitterdodge(jitter.width=0.15)) +
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
  theme(axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=25, vjust = -1),
        axis.title.y = element_text(size=25))

ggsave(plot=p, "within_cross_by_method.pdf",
       width=20, height=10)

##        Fig 2b       #############
within_wider <- read_csv("../data/within_tcr1.csv")

# format data into long form
within <- within_wider %>%
  pivot_longer(-Row, names_to="tcr", values_to = "score") %>%
  rename(method=Row) %>%
  mutate(method=fct_inorder(method))

within <- within[within$method %in% c('Hamming','BLOSUM100','PAM10','Dayhoff',
                                                            'Gonnet','Atchley'),]

# Add BATMAN within-TCR data
batman_data <- plot_data[(plot_data$mode=='Within TCR') & (plot_data$method=='BATMAN'),
                         c('TCR','method','score')]
colnames(batman_data) <- c('tcr','method','score')

within <- rbind(within, batman_data)

# plot
p_within <- ggplot(within, 
                   aes(x=factor(method, levels = c('BATMAN','BLOSUM100','PAM10',
                                                   'Hamming','Atchley','Dayhoff',
                                                   'Gonnet')), y=score)) +
  geom_boxplot(outlier.colour = NA,
               show.legend = FALSE, alpha=1,lwd=0.6, color='#DC143C') +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=tcr), 
             alpha=0.9,
             size=10,
             position = position_jitter(width=0.1),show.legend = FALSE) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
  labs(x="AA distance matrix",y="AUC score",color="") +
  ylim(NA, 0.9) +
  theme_cowplot() +
  theme(legend.text = element_text(size=45)) +
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.position="right")+
  theme(axis.text.x = element_text(size=45,angle=45),
        axis.text.y = element_text(size=45),
        axis.title.x = element_text(size=45, vjust = 2),
        axis.title.y = element_text(size=45))
ggsave(plot=p_within, "within_by_distance_matrix.pdf",
       width=15, height=10)

# barplot for positional weights
weights <- data.frame(position=c(1:9),
                      w=c(0.038,0.466,0.333,0.265,0.531,0.811,0.395,0.084,0.444))
p_bar <- ggplot(weights, aes(x=position, y=w)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_cowplot() +
  labs(x="Position",
       y="Weight") +
  scale_x_reverse(breaks=c(1:9))

ggsave(plot=p_bar, "positional_weights.pdf",
       width=2, height=4)


  


