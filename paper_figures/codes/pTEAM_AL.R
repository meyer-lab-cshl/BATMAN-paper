# Uses AL with pTEAM for a particular TCR

rm(list=ls())#Clear environment

library(readxl)
library(pROC)
library(randomForest)
library(tidyverse)

setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load TCR and peptide information from the main TCR datafile
TCR_data_full <- read_csv("data/batman_train_data.csv")

TCR_data_full$activation <- TCR_data_full$activation + 1 

# filter based on List of 3 class 9-mer-binding TCRs
TCR_names <- c("18A2","NYE-S1","868Z11","TCR3","TCR6","TCR7","A6",
               "T1","FLT3DY","A23","TCR2-T")
TCR_data_full <- TCR_data_full[TCR_data_full$tcr %in% TCR_names,]

peptides_to_test <- 9 # number of peptides to test at an AL cycle
n_cycle <- 10 # number of AL cycle

auc_pTEAM <- array(0,dim=c(n_cycle,length(TCR_names)))# Array to store AUCs

# Loop over TCR
for (TCR_index in 1:length(TCR_names)) {
  
  # Extract data for a particular TCR
  TCR_data <- TCR_data_full[TCR_data_full$tcr==TCR_names[TCR_index],]
  
  # Identify the central(/index) peptide and take it out of data
  index_peptide <- unique(TCR_data$index)
  
  peptide_length <- 9
  n_peptide <- dim(TCR_data)[1] #number of peptides
                            
  #Add Atchley factor embeddings#
  # Read Atchley Factors
  AA_list_Atchley <- c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F',
                       'P','S','T','W','Y','V','-')
  Atchley_factors <- as.matrix(read.csv("data/AA_matrices/Atchley.csv")[,2:6])
  Atchley_factors <- rbind(Atchley_factors,matrix(0,nrow=1, ncol=5))# Gap is padded as zero
  
  # Generate embeddings for peptides with the Atchley factors
  peptide_embedding <- matrix(0,nrow=n_peptide,
                         ncol=(10*peptide_length))
  
  # peptide embeddings
  peptide_list <- TCR_data$peptide
  
  for (peptide in 1:length(peptide_list)) {
    for (position in 1:peptide_length) {
      # Follow recipe in pTEAM paper
      peptide_embedding[peptide,(position*5-4):(position*5)] <- 
        Atchley_factors[which(AA_list_Atchley==
                  unlist(strsplit(peptide_list[peptide],split=""))[position]),]
      
      peptide_embedding[peptide,(peptide_length*5+position*5-4):
                         (peptide_length*5+position*5)] <- 
      Atchley_factors[which(AA_list_Atchley==
                  unlist(strsplit(peptide_list[peptide],split=""))[position]),]-
      Atchley_factors[which(AA_list_Atchley==
                        unlist(strsplit(index_peptide,split=""))[position]),]
    }
  }
  
  # Separate train and test set for AL
  # Step 0
  train_peptide <- index_peptide #initialize train set with index peptide
  
  # Find mutants with largest BLOSUM62 distance from index
  blosum62 <- read.csv("data/AA_matrices/blosum62_original.csv",row.names =1)
  diag(blosum62) <- 'NA' #Take out diagonal elements
  
  for (position in 1:peptide_length){
    wt_aa <- unlist(strsplit(index_peptide,split=""))[position]
    chosen_aa <- names(which.max(blosum62[wt_aa,]))
    
    # construct and add mutant to training set
    chosen_mutant <- unlist(strsplit(index_peptide,split=""))
    chosen_mutant[position] <- chosen_aa
    train_peptide <- c(train_peptide,str_c(chosen_mutant,collapse=""))
  }
    # full train and test set and embeddings
    train_set <- TCR_data[TCR_data$peptide %in% train_peptide,]
    test_set <- TCR_data[!TCR_data$peptide %in% train_peptide,]
    
    train_embedding <- peptide_embedding[TCR_data$peptide %in% train_peptide,]
    test_embedding <- peptide_embedding[!TCR_data$peptide %in% train_peptide,]
    
    # Set up random forest for TCRs
    set.seed(300+TCR_index)
    peptide_random_forest_classifier <- randomForest(
      train_embedding,
      as.factor(train_set$activation)) #Train
    
    peptide_random_forest_prediction <- predict(
      peptide_random_forest_classifier, 
      test_embedding,type="prob") #Predict for left out TCR
    
    ROC <- multiclass.roc(test_set$activation,
                          peptide_random_forest_prediction,direction=">")
    
    if (length(ROC$rocs)==1){ #Missing a class in train or test
      auc_pTEAM[1,TCR_index] <- (auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2
      
    } else {
      auc_pTEAM[1,TCR_index] <- 
        mean(c((auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2,
                       (auc(ROC$rocs[[2]][[1]])+auc(ROC$rocs[[2]][[2]]))/2,
                       (auc(ROC$rocs[[3]][[1]])+auc(ROC$rocs[[3]][[2]]))/2))
    }
    
    for (AL_step in 2:n_cycle){ # run AL cycle
      
      print(c(TCR_index,AL_step))
      
    # Evaluate uncertainty to test new peptides
    peptide_uncertainty <- 
      rowMeans(1/abs(peptide_random_forest_prediction - 
            colMeans(peptide_random_forest_prediction)))
    
    # Append training data
    train_peptide <- 
      c(train_peptide, 
        test_set$peptide[order(peptide_uncertainty,
                               decreasing = TRUE)[1:peptides_to_test]])
    
    # Update train and test set
    # full train and test set and embeddings
    train_set <- TCR_data[TCR_data$peptide %in% train_peptide,]
    test_set <- TCR_data[!TCR_data$peptide %in% train_peptide,]
    
    train_embedding <- peptide_embedding[TCR_data$peptide %in% train_peptide,]
    test_embedding <- peptide_embedding[!TCR_data$peptide %in% train_peptide,]
    
    # Set up random forest for TCRs
    set.seed(300+TCR_index)
    peptide_random_forest_classifier <- randomForest(
      train_embedding,
      as.factor(train_set$activation)) #Train
    
    peptide_random_forest_prediction <- predict(
      peptide_random_forest_classifier, 
      test_embedding,type="prob") #Predict for left out TCR
    
    ROC <- multiclass.roc(test_set$activation,
                          peptide_random_forest_prediction,direction=">")
    
    if (length(ROC$rocs)==1){ #Missing a class in train or test
      auc_pTEAM[AL_step,TCR_index] <- (auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2
      
    } else {
      auc_pTEAM[AL_step,TCR_index] <- 
        mean(c((auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2,
               (auc(ROC$rocs[[2]][[1]])+auc(ROC$rocs[[2]][[2]]))/2,
               (auc(ROC$rocs[[3]][[1]])+auc(ROC$rocs[[3]][[2]]))/2))}
    
    }
  
}

# Save data
auc_pTEAM <- data.frame(auc_pTEAM)
colnames(auc_pTEAM) <- TCR_names
rownames(auc_pTEAM) <- as.character(0:9)

write.csv(auc_pTEAM,file="data/AL_results_pTEAM.csv")


