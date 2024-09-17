# Plots performance of difference pooling methods over TCRs grouped by index

rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)
library(readxl)
library(viridis)
library(pROC)
library(ggnewscale)



#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load functions
source("../codes/all_functions.R") 

# Arrays to store 3 types of AUCs for 66 TCRs
#(unpooled, and for both pooled within index and pooled within n-mer:
#weights only, symm, and full AA matrix - 7 in total)
# (initialize with NaN so that missing data is not plotted)
all_AUC_non_weak <- array(NaN,dim=c(66,7))
all_AUC_non_strong <- array(NaN,dim=c(66,7))
all_AUC_weak_strong <- array(NaN,dim=c(66,7))

# Unpooled results (take BLOSUM100 only, AA_index=16)
# Indices of 2-class TCRs
TCR_index_2_class <- c(10,11,22,24,26:32)
class_pair_2_class <- c(12,12,12,12,23,12,23,23,23,23,12)#1:non, 2: weak, 3: strong

#3-class TCRs only
for (tcr in setdiff(1:66,TCR_index_2_class)){ 
  # results for Bayesian inference
      load(paste("../data/within_tcr_unpooled_outputs//output_",
               16,"_",tcr,".rda",sep=""))
    all_AUC_non_weak[tcr,1]<-mean(classification_results$AUCs[1,])
    all_AUC_non_strong[tcr,1]<-mean(classification_results$AUCs[2,])
    all_AUC_weak_strong[tcr,1]<-mean(classification_results$AUCs[3,])
}

#2-class TCRs only
for (i in 1:length(TCR_index_2_class)){ 
  # results for Bayesian inference
    load(paste("../data/within_tcr_unpooled_outputs//output_",
               16,"_",TCR_index_2_class[i],".rda",sep=""))
    
    if (class_pair_2_class[i]==12){
      all_AUC_non_weak[TCR_index_2_class[i],1]<-
        mean(classification_results$AUCs)}
    
    if (class_pair_2_class[i]==13){
      all_AUC_non_strong[TCR_index_2_class[i],1]<-
        mean(classification_results$AUCs)}
    
    if (class_pair_2_class[i]==23){
      all_AUC_weak_strong[TCR_index_2_class[i],1]<-
        mean(classification_results$AUCs)}
    
}
# Load TCR names
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)
TCRs <- unique(TCR_data$tcr_name)

## Pooled results
# Load training fold data for peptides
peptide_data <- read_excel("../data/training_folds.xlsx")
#3-class TCRs only
for (tcr in setdiff(1:66,TCR_index_2_class)){
  print(tcr)
  # Find if the index peptide is 9-mer or 8-mer
  index_peptide <- unique(TCR_data$index_peptide[TCR_data$tcr_name==TCRs[tcr]])
  peptide_length <- length(unlist(strsplit(index_peptide,split="")))
  
  # Array to store fold-specific aucs
  aucs_folds <- array(NaN,dim=c(5,7,3)) #fold-by-method-by-class_pairs
  for (fold in 0:4){
    # Load test peptide data
    test_peptide <- peptide_data$peptide_list[peptide_data$TCR_name==TCRs[tcr] & 
                                peptide_data$training_folds==fold]
    test_peptide_binding_category <- 1 + peptide_data$peptide_binding_category[
                                              peptide_data$TCR_name==TCRs[tcr] & 
                                              peptide_data$training_folds==fold]
    
    ## Weights only ##  
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                 'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Pooling across all indices #
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/weights_only/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
    weights_all <- read.csv(paste0('../data/within_tcr_pooled/weights_only/',
                               'weights_all_',peptide_length,
                               '_mer_fold_',fold,'.csv'))
    weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                                         2:(1+peptide_length)]))
    
    # Calculate AUC
    ROC<-multiclass.roc(test_peptide_binding_category,
                        (peptide_feature%*%weights)[,1],direction="<")
    
    aucs_folds[fold+1,1,1:3] <- c(auc(ROC$rocs[[1]]),
                                  auc(ROC$rocs[[2]]),
                                  auc(ROC$rocs[[3]]))
    
    }
    # Pooling within an index #
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/weights_only/',
                           'weights_',index_peptide,
                           '_fold_',fold,'.csv'))){
    weights_all <- read.csv(paste0('../data/within_tcr_pooled/weights_only/',
                                   'weights_',index_peptide,
                                   '_fold_',fold,'.csv'))
    weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                        2:(1+peptide_length)]))
    
    # Calculate AUC
    ROC<-multiclass.roc(test_peptide_binding_category,
                        (peptide_feature%*%weights)[,1],direction="<")
    
    aucs_folds[fold+1,2,1:3] <- c(auc(ROC$rocs[[1]]),
                                  auc(ROC$rocs[[2]]),
                                  auc(ROC$rocs[[3]]))
    }
    
    ## Weights+ Symmetric AA matrix ##
    # Pooling across all indices #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                     'aa_matrix_all_',peptide_length,
                                     '_mer_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                     'weights_all_',peptide_length,
                                     '_mer_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,3,1:3] <- c(auc(ROC$rocs[[1]]),
                                    auc(ROC$rocs[[2]]),
                                    auc(ROC$rocs[[3]]))
      
    }
    
    # Pooling within an index #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                        'aa_matrix_',index_peptide,
                        '_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                     'weights_',index_peptide,
                                     '_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,4,1:3] <- c(auc(ROC$rocs[[1]]),
                                    auc(ROC$rocs[[2]]),
                                    auc(ROC$rocs[[3]]))
    }
    ## Weights+ Full AA matrix ##
    # Pooling across all indices #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                        'aa_matrix_all_',peptide_length,
                        '_mer_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                                    'weights_all_',peptide_length,
                                    '_mer_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,5,1:3] <- c(auc(ROC$rocs[[1]]),
                                    auc(ROC$rocs[[2]]),
                                    auc(ROC$rocs[[3]]))
      
    }
    
    # Pooling within an index #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                        'aa_matrix_',index_peptide,
                        '_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                                    'weights_',index_peptide,
                                    '_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,6,1:3] <- c(auc(ROC$rocs[[1]]),
                                    auc(ROC$rocs[[2]]),
                                    auc(ROC$rocs[[3]]))
    }    
  }
    # Average over folds and store
    all_AUC_non_weak[tcr,2:7]<-colMeans(aucs_folds[,1:6,1])
    all_AUC_non_strong[tcr,2:7]<-colMeans(aucs_folds[,1:6,2])
    all_AUC_weak_strong[tcr,2:7]<-colMeans(aucs_folds[,1:6,3])
    
    # for 8-mer-binder tcrs (all ova-specific), within-index pooling=across-8-mer pooling
    if (peptide_length==8){
      all_AUC_non_weak[tcr,c(3,5,7)] <- all_AUC_non_weak[tcr,c(2,4,6)]
      all_AUC_non_strong[tcr,c(3,5,7)] <- all_AUC_non_strong[tcr,c(2,4,6)]
      all_AUC_weak_strong[tcr,c(3,5,7)] <- all_AUC_weak_strong[tcr,c(2,4,6)]
    }
  }
##############################

######################################################
#2-class TCRs only
for (i in 1:length(TCR_index_2_class)){
  tcr <- TCR_index_2_class[i]
  print(tcr)
  # Find if the index peptide is 9-mer or 8-mer
  index_peptide <- unique(TCR_data$index_peptide[TCR_data$tcr_name==TCRs[tcr]])
  peptide_length <- length(unlist(strsplit(index_peptide,split="")))
  
  # Array to store fold-specific aucs
  aucs_folds <- array(NaN,dim=c(5,7)) #fold-by-method
  for (fold in 0:4){
    # Load test peptide data
    test_peptide <- peptide_data$peptide_list[peptide_data$TCR_name==TCRs[tcr] & 
                                                peptide_data$training_folds==fold]
    test_peptide_binding_category <- 1 + peptide_data$peptide_binding_category[
      peptide_data$TCR_name==TCRs[tcr] & 
        peptide_data$training_folds==fold]
    
    ## Weights only ##  
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Pooling across all indices #
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/weights_only/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      weights_all <- read.csv(paste0('../data/within_tcr_pooled/weights_only/',
                                     'weights_all_',peptide_length,
                                     '_mer_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,1] <- auc(ROC$rocs[[1]])
      
    }
    # Pooling within an index #
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/weights_only/',
                           'weights_',index_peptide,
                           '_fold_',fold,'.csv'))){
      weights_all <- read.csv(paste0('../data/within_tcr_pooled/weights_only/',
                                     'weights_',index_peptide,
                                     '_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,2] <- auc(ROC$rocs[[1]])
    }
    
    ## Weights+ Symmetric AA matrix ##
    # Pooling across all indices #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                        'aa_matrix_all_',peptide_length,
                        '_mer_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                    'weights_all_',peptide_length,
                                    '_mer_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,3] <- auc(ROC$rocs[[1]])
      
    }
    
    # Pooling within an index #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                        'aa_matrix_',index_peptide,
                        '_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                    'weights_',index_peptide,
                                    '_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,4] <- auc(ROC$rocs[[1]])
    }
    ## Weights+ Full AA matrix ##
    # Pooling across all indices #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                        'aa_matrix_all_',peptide_length,
                        '_mer_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'weights_all_',peptide_length,
                           '_mer_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                                    'weights_all_',peptide_length,
                                    '_mer_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,5] <- auc(ROC$rocs[[1]])
      
    }
    
    # Pooling within an index #
    # Load/create AA matrix
    AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
                'T','V','W','Y')
    load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
    AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
    
    # Load matrix multiplier
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      AA_multiplier <- 
        read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                        'aa_matrix_',index_peptide,
                        '_fold_',fold,'.csv'))
      AA_multiplier <- data.matrix(AA_multiplier[,2:21])
      AA_distance <- AA_distance*AA_multiplier}
    # Make peptide features
    peptide_feature <- generate_peptide_feature(test_peptide,index_peptide,
                                                AA_distance)
    
    # Load positional weight for selected TCR
    if (file.exists(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                           'aa_matrix_',index_peptide,
                           '_fold_',fold,'.csv'))){
      weights_all <-read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                                    'weights_',index_peptide,
                                    '_fold_',fold,'.csv'))
      weights <- array(unlist(weights_all[weights_all[,1]==TCRs[tcr],
                                          2:(1+peptide_length)]))
      
      # Calculate AUC
      ROC<-multiclass.roc(test_peptide_binding_category,
                          (peptide_feature%*%weights)[,1],direction="<")
      
      aucs_folds[fold+1,6] <- auc(ROC$rocs[[1]])
    }    
  }
  # Average over folds and store
  if (class_pair_2_class[i]==12){
  all_AUC_non_weak[tcr,2:7]<-colMeans(aucs_folds[,1:6])}
  
  if (class_pair_2_class[i]==13){
  all_AUC_non_strong[tcr,2:7]<-colMeans(aucs_folds[,1:6])}
  
  if (class_pair_2_class[i]==23){
  all_AUC_weak_strong[tcr,2:7]<-colMeans(aucs_folds[,1:6])}
  
  # for 8-mer-binder tcrs (all ova-specific), within-index pooling=across-8-mer pooling
  if (peptide_length==8){
    all_AUC_non_weak[tcr,c(3,5,7)] <- all_AUC_non_weak[tcr,c(2,4,6)]
    all_AUC_non_strong[tcr,c(3,5,7)] <- all_AUC_non_strong[tcr,c(2,4,6)]
    all_AUC_weak_strong[tcr,c(3,5,7)] <- all_AUC_weak_strong[tcr,c(2,4,6)]
  }
}

##############################
# order of methods: 
# unpooled, weights_only_all, weights_only_within_index,
# symm_only_all, symm_only_within_index
# full_only_all, full_only_within_index

methods <- c('Unpooled','BLOSUM_across','BLOSUM_within', 
             'Symmetric_across','Symmetric_within',
             'Full_across','Full_within')

# Add TCR and AA matrix names to array
colnames(all_AUC_non_weak) <- methods
colnames(all_AUC_non_strong) <- methods
colnames(all_AUC_weak_strong) <- methods
all_AUC_non_weak<- data.frame(TCRs,all_AUC_non_weak)
all_AUC_non_strong<- data.frame(TCRs,all_AUC_non_strong)
all_AUC_weak_strong<- data.frame(TCRs,all_AUC_weak_strong)

# Convert to long form for box plot
all_AUC_non_weak <- all_AUC_non_weak %>%
  pivot_longer(-TCRs, names_to="method", values_to = "score")
  
all_AUC_non_strong <- all_AUC_non_strong %>%
  pivot_longer(-TCRs, names_to="method", values_to = "score") 

all_AUC_weak_strong <- all_AUC_weak_strong %>%
  pivot_longer(-TCRs, names_to="method", values_to = "score") 

# Merge data and indicate class pair
all_AUC <- rbind.data.frame(all_AUC_non_weak,all_AUC_non_strong,all_AUC_weak_strong)
all_AUC$class_pair <- c(rep("12", each = dim(all_AUC_non_weak)[1]),
                        rep("13", each = dim(all_AUC_non_strong)[1]),
                        rep("23", each = dim(all_AUC_weak_strong)[1]))


# Indicate index peptides for each tcr in data point
# Array to store TCR names and index peptides
index_peptide <- array("X",dim=c(dim(all_AUC)[1]))
for (i in 1:dim(all_AUC)[1]){
  TCR_index <- match(all_AUC$TCRs[i],unique(all_AUC$TCRs)) #this is because R changed TCR names
  index_peptide[i] <- unique(TCR_data$index_peptide[TCR_data$tcr_name==TCRs[TCR_index]])
}

all_AUC$index_peptide <- index_peptide

# Discard data for one 10-mer-binding and one 11-mer-binding TCR
all_AUC <- all_AUC[!(all_AUC$index_peptide %in% c('ALWGPDPAAA','ALYDKTKRIFL')),]

all_AUC %>%
  mutate(method=fct_inorder(method))

level_order <- unique(all_AUC$method)

# For SIINFEKL TCRs, classify them into 2 categories
tcr_type <- array(NaN,dim=c(dim(all_AUC)[1]))

naive_tcrs <- c("B11","B15","B3", "F4","E8", "B13","H6", "G6", "F5", 
                "H5", "B2", "B6", "B5","E9", "E4", "G2", "B16","B14")
ed_tcrs <- c("OT1","Ed5","Ed8","Ed9","Ed10","Ed16-1","Ed16-30","Ed21","Ed23",
             "Ed28","Ed31","Ed33","Ed39","Ed40","Ed45","Ed46")

tcr_type[all_AUC$TCRs %in% naive_tcrs] <- "Naive"
tcr_type[all_AUC$TCRs %in% ed_tcrs] <- "Educated"

all_AUC$tcr_type <- tcr_type




###############
##  plot ######
###############

# Special plotting specs for SIINFEKL data
all_AUC <- all_AUC %>%
  mutate(color_by = TCRs,
         color_by = case_when(index_peptide == "SIINFEKL" ~ tcr_type,
                              TRUE ~ color_by))

auc_by_index <- split(all_AUC, f=all_AUC$index_peptide)

boxplot_by_index <- lapply(auc_by_index, function(x) {
  
  set_alpha <- 0.7 
  set_size <- 2.25
  set_jitter_width <- 0.15
  
  if (unique(x$index_peptide)=="SIINFEKL"){
    set_alpha <- 0.85
    set_size <- 1.6
    set_jitter_width <- 0.2
  }
  
  p <- ggplot(x, 
              aes(x=factor(method, levels = unique(all_AUC$method)), y=score)) +
    geom_boxplot(aes(color=factor(class_pair,levels = unique(class_pair))),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1,
                 position=position_dodge(0.75)) +
    scale_color_manual(values=c('#50C878','#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(color_by), 
                   group=factor(class_pair,levels = unique(class_pair))), 
               alpha=set_alpha,
               size=set_size,
               position = position_jitterdodge(jitter.width=set_jitter_width,
                                               dodge.width = 0.75)) +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(y="AUC",
         x="",
         color="") +
#    ylim(min(all_AUC$score,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=20))+
    facet_wrap(vars(index_peptide))+
    theme(strip.text = element_text(
      size = 25, color = "black")) +
    theme(legend.text = element_text(size=15)) +
    guides(color = guide_legend(ncol = 3)) +
    theme(legend.position="bottom")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=4,nrow=3,
                            align = "hv",axis="tb")

ggsave(plot=p_all, "../figures/extended_data_fig3/pooled_auc_by_index.pdf",
       width=20, height=20)