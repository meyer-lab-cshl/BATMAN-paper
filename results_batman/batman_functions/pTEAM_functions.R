# Align and featurize an array of TCR CDR3 sequences
generate_CDR3_feature <- function(cdr3_list) {
  
  library(muscle) #library for AA alignment

  # Extract unique sequences present in the data and store their indices
  cdr3_unique <- unique(cdr3_list)
  cdr3_unique_indices <- match(cdr3_list,cdr3_unique)
  
  # MUSCLE-align sequences
  cdr3_aligned <- as.character(muscle(stringset = AAStringSet(x=cdr3_unique, 
                                                            start=NA, 
                                                            end=NA, 
                                                            width=NA, 
                                                            use.names=FALSE)))
  cdr3_list_aligned <- cdr3_aligned[cdr3_unique_indices]
  
  # Featurize with Atchley factors
  # table 2 in https://www.pnas.org/content/102/18/6395
  atchley_factors_table <- t(array(c(
    c(-0.591,	-1.302,	-0.733,	1.570,	-0.146),
    c(-1.343,	0.465,	-0.862,	-1.020,	-0.255),
    c(1.050,	0.302,	-3.656,	-0.259,	-3.242),
    c(1.357,	-1.453,	1.477,	0.113,	-0.837),
    c(-1.006,	-0.590,	1.891,	-0.397,	0.412),
    c(-0.384,	1.652,	1.330,	1.045,	2.064),
    c(0.336,	-0.417,	-1.673,	-1.474,	-0.078),
    c(-1.239,	-0.547,	2.131,	0.393,	0.816),
    c(1.831,	-0.561,	0.533,	-0.277,	1.648),
    c(-1.019,	-0.987,	-1.505,	1.266,	-0.912),
    c(-0.663,	-1.524,	2.219,	-1.005,	1.212),
    c(0.945,	0.828,	1.299,	-0.169,	0.933),
    c(0.189,	2.081,	-1.628,	0.421,	-1.392),
    c(0.931,	-0.179,	-3.005,	-0.503,	-1.853),
    c(1.538,	-0.055,	1.502,	0.440,	2.897),
    c(-0.228,	1.399,	-4.760,	0.670,	-2.647),
    c(-0.032,	0.326,	2.213,	0.908,	1.313),
    c(-1.337,	-0.279,	-0.544,	1.242,	-1.262),
    c(-0.595,	0.009,	0.672,	-2.128,	-0.184),
    c(0.260,	0.830,	3.097,	-0.838,	1.512),
    c(0,0,0,0,0) #zeros for gap padding
    ), dim=c(5,21)))
  
  # Order of AAs
  aa_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
               'T','V','W','Y','-')
  
  atchley_factors_table <- data.frame(atchley_factors_table,row.names = aa_list)
  
  # Featurize aligned CDR3
  cdr3_feature <- data.frame(
    atchley_factors_table[unlist(strsplit(cdr3_list_aligned,split="")),]
    )
  
  cdr3_feature <- matrix(c(t(cdr3_feature)),nrow=length(cdr3_list),byrow=TRUE)
  
  return(cdr3_feature)
}
##################################################################

## Featurize an array pair of mutated and unmutated sequences
generate_peptide_feature <- function(mutant_list,index_list) {
  
  # Return error if unmutated peptide is present in mutant list
  if (sum(as.integer(mutant_list==index_list))>0){
    stop("Unmutated peptide cannot be in mutant list!")
  }
  
  #encode mutation position as one-hot vector
  mutations_position <- t(matrix(
    paste0(unlist(strsplit(mutant_list,split=""))) != 
      paste0(unlist(strsplit(index_list,split=""))),
    ncol=length(mutant_list)))*1
  
  # Get WT and mutated AA
  wt_aa <- paste0(unlist(strsplit(index_list,split="")))[
    paste0(unlist(strsplit(mutant_list,split=""))) != 
    paste0(unlist(strsplit(index_list,split="")))]
  
  mutant_aa <- paste0(unlist(strsplit(mutant_list,split="")))[
    paste0(unlist(strsplit(mutant_list,split=""))) != 
    paste0(unlist(strsplit(index_list,split="")))]
  
  # Get atchley embeddings
  # https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00238-6 for details
  # table 2 in https://www.pnas.org/content/102/18/6395
  atchley_factors_table <- t(array(c(
    c(-0.591,	-1.302,	-0.733,	1.570,	-0.146),
    c(-1.343,	0.465,	-0.862,	-1.020,	-0.255),
    c(1.050,	0.302,	-3.656,	-0.259,	-3.242),
    c(1.357,	-1.453,	1.477,	0.113,	-0.837),
    c(-1.006,	-0.590,	1.891,	-0.397,	0.412),
    c(-0.384,	1.652,	1.330,	1.045,	2.064),
    c(0.336,	-0.417,	-1.673,	-1.474,	-0.078),
    c(-1.239,	-0.547,	2.131,	0.393,	0.816),
    c(1.831,	-0.561,	0.533,	-0.277,	1.648),
    c(-1.019,	-0.987,	-1.505,	1.266,	-0.912),
    c(-0.663,	-1.524,	2.219,	-1.005,	1.212),
    c(0.945,	0.828,	1.299,	-0.169,	0.933),
    c(0.189,	2.081,	-1.628,	0.421,	-1.392),
    c(0.931,	-0.179,	-3.005,	-0.503,	-1.853),
    c(1.538,	-0.055,	1.502,	0.440,	2.897),
    c(-0.228,	1.399,	-4.760,	0.670,	-2.647),
    c(-0.032,	0.326,	2.213,	0.908,	1.313),
    c(-1.337,	-0.279,	-0.544,	1.242,	-1.262),
    c(-0.595,	0.009,	0.672,	-2.128,	-0.184),
    c(0.260,	0.830,	3.097,	-0.838,	1.512),
    c(0,0,0,0,0) #zeros for gap padding
  ), dim=c(5,21)))
  
  # Order of AAs
  aa_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
               'T','V','W','Y','-')
  
  atchley_factors_table <- data.frame(atchley_factors_table,row.names = aa_list)
  
  # Mutant features
  mutant_feature <- data.frame(
    atchley_factors_table[unlist(strsplit(mutant_list,split="")),])
  
  mutant_feature <- matrix(c(t(mutant_feature)),nrow=length(mutant_list),byrow=TRUE)
  
  # Difference between WT and mutant AA features
  wt_feature <- data.frame(
    atchley_factors_table[unlist(strsplit(index_list,split="")),])
  
  wt_feature <- matrix(c(t(wt_feature)),nrow=length(index_list),byrow=TRUE)
  
  mutant_feature_diff <- mutant_feature - wt_feature
  
  # WT AA feature
  wt_aa_feature <- data.frame(
    atchley_factors_table[wt_aa,])
  
  wt_aa_feature <- matrix(c(t(wt_aa_feature)),nrow=length(mutant_list),byrow=TRUE)
  
  # Mutant aa feature
  mutant_aa_feature <- data.frame(
    atchley_factors_table[mutant_aa,])
  
  mutant_aa_feature <- matrix(c(t(mutant_aa_feature)),
                              nrow=length(mutant_list),byrow=TRUE)
  
  # Concatenate all features together
  peptide_features <- cbind(mutant_feature,mutant_feature_diff,
                            mutations_position,wt_aa_feature,mutant_aa_feature)
  
  return(peptide_features)
}

############################################################################
## Generate peptide and (optionally) TCR features for train/test input
generate_tcr_pmhc_feature <- function(TCR_data,use_tcr=FALSE) {
  
  # Featurize inputs
  tcr_pmhc_features <- generate_peptide_feature(TCR_data$peptide,
                                                TCR_data$index_peptide)
  
  if (use_tcr==TRUE){ #Use embedded TCR sequences
    cdr3a_features <- generate_CDR3_feature(TCR_data$cdr3a)
    cdr3b_features <- generate_CDR3_feature(TCR_data$cdr3b)
    
    tcr_pmhc_features <- cbind(tcr_pmhc_features,cdr3a_features,cdr3b_features)
  }
  
  return(tcr_pmhc_features)
}

############################################################################


## Train the random forest in pTEAM and return test auc
train_pTEAM <- function(train_data,test_data,use_tcr=FALSE,seed=123) {
  
  set.seed(seed)
  library(randomForest)
  library(pROC)
  
  # Merge all data
  train_data$is_train <- 1
  test_data$is_train <- 0
  
  TCR_data <- rbind(train_data,test_data)
  
  # Featurize inputs
  tcr_pmhc_features <- generate_tcr_pmhc_feature(TCR_data,use_tcr)
  
  # Train Random Forest
  trained_random_forest <- randomForest(
    tcr_pmhc_features[TCR_data$is_train==1,],
    as.factor(TCR_data$activation[TCR_data$is_train==1]))
  
  # Get test AUC averaged over 3 classes
  # Test scores
  rf_prediction <- predict(
    trained_random_forest, 
    tcr_pmhc_features[TCR_data$is_train==0,],
    type="prob")
  
  # Calculate ROC and average AUC for 3 TCR activation classes
  ROC <- multiclass.roc(TCR_data$activation[TCR_data$is_train==0],
    rf_prediction,
    direction=">") #ROC curve
  
  if (length(ROC$rocs)==1){ #Missing a class
    auc_mean <- (auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2
    
  } else { # All 3 classes present
    auc_mean <- mean(
      c((auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2,
        (auc(ROC$rocs[[2]][[1]])+auc(ROC$rocs[[2]][[2]]))/2,
        (auc(ROC$rocs[[3]][[1]])+auc(ROC$rocs[[3]][[2]]))/2))
  }
  
  return(auc_mean)
}






