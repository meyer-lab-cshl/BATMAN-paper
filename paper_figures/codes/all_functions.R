## Load continuous activation data for mutants against selected TCR
load_tcr_activation <- function(TCR_index) {
  
  library(readxl)
  
  # set working directory to where this script is
  # setwd(getSrcDirectory(function(){})[1]) 
  
  # Load TCR and peptide information from the main TCR datafile
  suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
  )
  
  # Extract TCR names present in the data
  TCR_names <- unique(TCR_data$tcr_name)
  
  # Extract peptides and their activities for the particular TCR
  peptide_list <- TCR_data$peptide[TCR_data$tcr_name==TCR_names[TCR_index]]
  peptide_activity <- TCR_data$peptide_activity[
                                   TCR_data$tcr_name==TCR_names[TCR_index]]
  
  # Define the index peptide as the one with the highest activity
  index_peptide <- TCR_data$peptide[TCR_data$tcr_name==TCR_names[TCR_index]][
                                                  which.max(peptide_activity)]
  
  # Identify the index (/central/cognate) peptide
  index_peptide <- unique(TCR_data$index_peptide[
                                       TCR_data$tcr_name==TCR_names[TCR_index]])
  
  #Normalize peptide activities to the index activity
  peptide_activity <- peptide_activity/(peptide_activity[
                              which.max(peptide_activity)])
  
  #Discretize activation values to peptide binding categories
  #Initialize with all weak binders
  peptide_binding_category <- matrix(2,nrow=length(peptide_list),ncol=1)
  peptide_binding_category[peptide_activity<0.1] <- 1 #Non-binder
  peptide_binding_category[peptide_activity>=0.5] <- 3 #strong-binder
  peptide_binding_category <- as.integer(peptide_binding_category)
  
  #Discard binding class with <5 elements if 3 classes are present
  if (length(unique(peptide_binding_category))==3 & 
      min(table(peptide_binding_category)) <= 4){
    
    peptide_list<-peptide_list[peptide_binding_category!=which.min(
                                          table(peptide_binding_category))[[1]]]
    
    peptide_activity<-peptide_activity[peptide_binding_category!=which.min(
                                          table(peptide_binding_category))[[1]]]
    
    peptide_binding_category<-peptide_binding_category[
      peptide_binding_category!=which.min(table(peptide_binding_category))[[1]]]
     
}
  

  activation_data <- list("peptide_list" = peptide_list,
                          "peptide_activity" = peptide_activity,
                          "index_peptide" = index_peptide,
                          "peptide_binding_category"=peptide_binding_category)
}

##################################################################
# Load AA distance matrix
load_AA_distance <- function(AA_matrix_index) {

  library(readxl)
  
  # set working directory to where this script is
#  setwd(getSrcDirectory(function(){})[1])
  # List of AA distance functions
  AA_matrix_name_list <-c("hamming",paste(rep("blosum",13),seq(30,90,5),sep=""),
                           "blosum62","blosum100",
                           paste(rep("pam",50),seq(10,500,10),sep=""),
                           "dayhoff","gonnet","atchley_l2","atchley_cos");
  
  #Order of AAs
  AA_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
               'T','V','W','Y')
  
  #Name of the selected distance function
  AA_matrix_name <- AA_matrix_name_list[AA_matrix_index]
  
  #Load matrix
  load(paste("../data/distance_functions/",AA_matrix_name,".rda",sep=""))
  AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
}
###################################################################
# Generate position-dependent distance feature for mutant peptides
generate_peptide_feature <- function(peptide_list,index_peptide,AA_distance) {
 
  # Infer peptide length
  peptide_length <- length(unlist(strsplit(index_peptide,split="")))
  
  peptide_feature <- matrix(0,nrow=length(peptide_list),ncol=peptide_length)
  
  for (peptide in 1:length(peptide_list)) {
    for (position in 1:peptide_length) {
      peptide_feature[peptide,position] <- 
        (-1)*AA_distance[unlist(strsplit(index_peptide,split=""))[position],
                unlist(strsplit(peptide_list[peptide],split=""))[position]]
    }
  }
  return(peptide_feature)
}
####################################################################
# Perform regression and extract positional weights and correlation coefficients
# results from 5-fold cv stratified over TCR activation value deciles
regression_cv <- function(peptide_feature,peptide_activity,nfold,seed1,seed2) {

  library(rstan)
  library(brms)
  library(loo)
  library(dplyr) 
  
  # Indicate peptide deciles
  peptide_decile <- ntile(peptide_activity,10)
  peptide_data <- data.frame(peptide_feature,peptide_activity)
  
  #Infer peptide length
  peptide_length <- dim(peptide_feature)[2]
  
  # Do regression with brms package
  #specify initial parameter guesses for 4 chains
  initial_parameters1<-as.list(matrix(0.3,nrow=2+peptide_length,ncol=1))
  names(initial_parameters1) <- c("Intercept[1]","Sigma",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters2<-as.list(matrix(0.4,nrow=2+peptide_length,ncol=1))
  names(initial_parameters2) <- c("Intercept[1]","Sigma",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters3<-as.list(matrix(0.5,nrow=2+peptide_length,ncol=1))
  names(initial_parameters3) <- c("Intercept[1]","Sigma",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters4<-as.list(matrix(0.6,nrow=2+peptide_length,ncol=1))
  names(initial_parameters4) <- c("Intercept[1]","Sigma",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameter_list<-list(initial_parameters1,initial_parameters2,
                            initial_parameters3,initial_parameters4)
  
  # brms Regression object
  peptide_regression <- 
    brm(data = peptide_data, 
        family = gaussian,
        as.formula(paste("peptide_activity~1+",gsub(", ","+",
                  toString(colnames(peptide_data)[1:peptide_length])),sep="")),
        prior = c(prior(normal(0,1.5), class = Intercept),
                  prior(beta(2,2), lb=0.01, ub=0.99, class = b),
                  prior(exponential(1), class = sigma)),
        init = initial_parameter_list,
        iter = 2000, warmup = 1000, chains = 4,
        seed = seed1)
  
  #5-fold cross-validation sampling
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_decile)
  peptide_regression_cross_validated <- kfold(peptide_regression, K=nfold, 
                                          save_fits=TRUE, folds=training_folds)
  kfold_fits<-peptide_regression_cross_validated$fits
  
  #Extract position-dependent weights from the classifier for 5 folds
  position_dependent_weights<-array(0,dim=c(1+peptide_length,nfold))
  correlations<-array(0,dim=c(1,nfold))
  for (i in 1:(nfold)){
  position_dependent_weights[,i]<-fixef(kfold_fits[[i]])[1:(1+peptide_length),1]
  correlations[1,i]<-cor(predict(kfold_fits[[i]],
                                 peptide_data[training_folds==i,])[,1],
                          peptide_activity[training_folds==i],method="spearman")
  }
  
  regression_results <- list("correlations" = correlations,
                      "position_dependent_weights" = position_dependent_weights)
  
}
###############################################################################
####################################################################
# Perform regression and extract positional weights and correlation coefficients
# results from 5-fold cv stratified over TCR activation value deciles
classification_3_class_cv <- function(peptide_feature,
                                      peptide_binding_category,
                                      nfold,seed1,seed2) {
  
  library(rstan)
  library(brms)
  library(loo)
  library(dplyr)
  library(pROC)
  
  
  peptide_data <- data.frame(peptide_feature,peptide_binding_category)
  
  #Infer peptide length
  peptide_length <- dim(peptide_feature)[2]
  
  # Do regression with brms package
  #specify initial parameter guesses for 4 chains
  initial_parameters1<-as.list(matrix(0.3,nrow=2+peptide_length,ncol=1))
  names(initial_parameters1) <- c("Intercept[1]","Intercept[2]",
                                paste(rep.int('b[',peptide_length),
                        (1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters2<-as.list(matrix(0.4,nrow=2+peptide_length,ncol=1))
  names(initial_parameters2) <- c("Intercept[1]","Intercept[2]",
                                  paste(rep.int('b[',peptide_length),
                        (1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters3<-as.list(matrix(0.5,nrow=2+peptide_length,ncol=1))
  names(initial_parameters3) <- c("Intercept[1]","Intercept[2]",
                                  paste(rep.int('b[',peptide_length),
                        (1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters4<-as.list(matrix(0.6,nrow=2+peptide_length,ncol=1))
  names(initial_parameters4) <- c("Intercept[1]","Intercept[2]",
                                  paste(rep.int('b[',peptide_length),
                        (1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameter_list<-list(initial_parameters1,initial_parameters2,
                               initial_parameters3,initial_parameters4)
  
  # brms classifier object
  peptide_classifier <- 
    brm(data = peptide_data, 
        family = cumulative,
        as.formula(paste("peptide_binding_category~1+",gsub(", ","+",
                  toString(colnames(peptide_data)[1:peptide_length])),sep="")),
        prior = c(prior(normal(0,1.5), class = Intercept),
                  prior(beta(2,2), lb=0.01, ub=0.99, class = b)),
        init = initial_parameter_list,
        iter = 2000, warmup = 1000, chains = 4,
        seed = seed1)
  
  #5-fold cross-validation sampling
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_binding_category)
  peptide_classifier_cross_validated <- kfold(peptide_classifier, K=nfold, 
                                          save_fits=TRUE, folds=training_folds)
  kfold_fits<-peptide_classifier_cross_validated$fits
  
  #Extract position-dependent weights and AUCs from the classifier for 5 folds
  position_dependent_weights<-array(0,dim=c(peptide_length,nfold))
  AUCs <- array(0,dim=c(3,nfold))
  
  for (i in 1:(nfold)){
  position_dependent_weights[,i]<-fixef(kfold_fits[[i]])[3:(2+peptide_length),1]
  
  ROC<-multiclass.roc(peptide_binding_category[training_folds==i],
  (peptide_feature[training_folds==i,]%*%position_dependent_weights
                                                        [,i])[,1],direction="<")
    
  AUCs[,i] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  }
  
  classification_results <- list("AUCs" = AUCs,
                      "position_dependent_weights" = position_dependent_weights)
  
}
###############################################################################
####################################################################
# Perform regression and extract positional weights and correlation coefficients
# results from 5-fold cv stratified over TCR activation value deciles
classification_2_class_cv <- function(peptide_feature,
                                      peptide_binding_category,
                                      nfold,seed1,seed2) {
  
  library(rstan)
  library(brms)
  library(loo)
  library(dplyr) 
  library(pROC)
  
  
  peptide_data <- data.frame(peptide_feature,peptide_binding_category)
  
  #Infer peptide length
  peptide_length <- dim(peptide_feature)[2]
  
  # Do regression with brms package
  #specify initial parameter guesses for 4 chains
  initial_parameters1<-as.list(matrix(0.3,nrow=1+peptide_length,ncol=1))
  names(initial_parameters1) <- c("Intercept[1]",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters2<-as.list(matrix(0.4,nrow=1+peptide_length,ncol=1))
  names(initial_parameters2) <- c("Intercept[1]",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters3<-as.list(matrix(0.5,nrow=1+peptide_length,ncol=1))
  names(initial_parameters3) <- c("Intercept[1]",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameters4<-as.list(matrix(0.6,nrow=1+peptide_length,ncol=1))
  names(initial_parameters4) <- c("Intercept[1]",paste(rep.int('b[',
        peptide_length),(1:peptide_length),rep.int(']',peptide_length),sep=""))
  
  initial_parameter_list<-list(initial_parameters1,initial_parameters2,
                               initial_parameters3,initial_parameters4)
  
  # brms classification object
  peptide_classifier <- 
    brm(data = peptide_data, 
        family = bernoulli,
        as.formula(paste("peptide_binding_category~1+",gsub(", ","+",
                  toString(colnames(peptide_data)[1:peptide_length])),sep="")),
        prior = c(prior(normal(0,1.5), class = Intercept),
                  prior(beta(2,2), lb=0.01, ub=0.99, class = b)),
        init = initial_parameter_list,
        iter = 2000, warmup = 1000, chains = 4,
        seed = seed1)
  
  #5-fold cross-validation sampling
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_binding_category)
  peptide_classifier_cross_validated <- kfold(peptide_classifier, K=nfold, 
                                           save_fits=TRUE, folds=training_folds)
  kfold_fits<-peptide_classifier_cross_validated$fits
  
  #Extract position-dependent weights and AUCs from the classifier for 5 folds
  position_dependent_weights<-array(0,dim=c(peptide_length,nfold))
  AUCs <- array(0,dim=c(1,nfold))
  
  for (i in 1:(nfold)){
  position_dependent_weights[,i]<-fixef(kfold_fits[[i]])[2:(1+peptide_length),1]
    
  ROC<-multiclass.roc(peptide_binding_category[training_folds==i],
  (peptide_feature[training_folds==i,]%*%position_dependent_weights
                                                        [,i])[,1],direction="<")
    
  AUCs[,i] <- auc(ROC$rocs[[1]])
  }
  
  classification_results <- list("AUCs" = AUCs,
                      "position_dependent_weights" = position_dependent_weights)
  
}
###############################################################################
# Generate Atchley feature for mutant peptides

generate_peptide_feature_atchley <- function(peptide_list,index_peptide){
  
  library(readxl)
  
  # Read Atchley Factors
  AA_list <- c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S',
                                                                'T','W','Y','V')
  Atchley_factors<-as.matrix(read.csv("../data/atchley.csv")[,2:6])
  
  # Deduce peptide length
  peptide_length=length(unlist(strsplit(index_peptide,split="")))
  
  #Generate embeddings for peptides with the Atchley factors
  # according to https://www.biorxiv.org/content/10.1101/2023.05.10.540189v1
  peptide_feature<-matrix(0,nrow=length(peptide_list),ncol=10*peptide_length)
  
  for (peptide in 1:length(peptide_list)) {
    for (position in 1:peptide_length) {
      peptide_feature[peptide,(position*5-4):(position*5)]<-Atchley_factors[
    which(AA_list==unlist(strsplit(peptide_list[peptide],split=""))[position]),]
      
      peptide_feature[peptide,(peptide_length*5+position*5-4):
                        (peptide_length*5+position*5)]  <-  Atchley_factors[
   which(AA_list==unlist(strsplit(peptide_list[peptide],split=""))[position]),]-
                                                            Atchley_factors[
           which(AA_list==unlist(strsplit(index_peptide,split=""))[position]),]
    }
  }
  
  return(peptide_feature)
  
}
###########################################################################
# Perform regression with Atchley+RF and extract correlation coefficients
# Results from n-fold cv stratified over TCR activation value deciles
regression_cv_atchley <- function(peptide_feature,peptide_activity,
                                                 nfold,seed1,seed2){
  
  library(randomForest)
  library(loo) 
  library(pROC)
  library(dplyr)
  
  # Indicate peptide deciles
  peptide_decile<-ntile(peptide_activity,10)
  
  #Set up k folds for cross validation
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_decile)
  
  # Set up random forest for k folds
  correlations=array(0,dim=c(1,nfold))
  
  for (i in 1:nfold){
    #Regression Trees
    set.seed(seed1)
    peptide_random_forest_regression <-
      randomForest(peptide_feature[training_folds!=i,],
                   peptide_activity[training_folds!=i],ntrees=500)
    # Prediction
    peptide_activity_predicted <- 
      predict(peptide_random_forest_regression,
              peptide_feature[training_folds==i,])
    
    correlations[1,i]<-cor(peptide_activity_predicted,
                           peptide_activity[training_folds==i],
                           method="spearman")
  }
  return(correlations)
  
}
###########################################################################
# Perform classification with Atchley+RF and extract AUCs
# Results from n-fold cv stratified over TCR activation categories
# Seeds are for RF and creating folds for 3 class classification

classification_3_class_cv_atchley <- function(peptide_feature,
                                peptide_binding_category,nfold,seed1,seed2){
  
  #Set up k folds for cross validation
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_binding_category)
  
  # Set up random forest for k folds
  AUCs=array(0,dim=c(3,nfold))
  
  for (i in 1:nfold){
    #Classification trees
    set.seed(seed1)
    peptide_random_forest_classifier <- 
      randomForest(peptide_feature[training_folds!=i,],
                   as.factor(peptide_binding_category[training_folds!=i]),
                   ntrees=500)
    
    # Prediction scores
    peptide_random_forest_prediction <- 
      predict(peptide_random_forest_classifier, 
              peptide_feature[training_folds==i,],
              type="prob")
    
    # AUC calculation for 3 classes
    ROC<-multiclass.roc(peptide_binding_category[training_folds==i],
                        peptide_random_forest_prediction)
    AUCs[1:3,i]<-c((auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2,
                 (auc(ROC$rocs[[2]][[1]])+auc(ROC$rocs[[2]][[2]]))/2,
                 (auc(ROC$rocs[[3]][[1]])+auc(ROC$rocs[[3]][[2]]))/2)
  }
  return(AUCs)
}
##############################################################################
###########################################################################
# Perform classification with Atchley+RF and extract AUCs
# Results from n-fold cv stratified over TCR activation categories
# Seeds are for RF and creating folds for 3 class classification

classification_2_class_cv_atchley <- function(peptide_feature,
                                    peptide_binding_category,nfold,seed1,seed2){
  
  #Set up k folds for cross validation
  set.seed(seed2)
  training_folds=kfold_split_stratified(K = nfold, x = peptide_binding_category)
  
  # Set up random forest for k folds
  AUCs=array(0,dim=c(1,nfold))
  
  for (i in 1:nfold){
    #Classification trees
    set.seed(seed1)
    peptide_random_forest_classifier <- 
      randomForest(peptide_feature[training_folds!=i,],
                   as.factor(peptide_binding_category[training_folds!=i]),
                   ntrees=500)
    
    # Prediction scores
    peptide_random_forest_prediction <- 
      predict(peptide_random_forest_classifier, 
              peptide_feature[training_folds==i,],
              type="prob")
    
    # AUC calculation for 3 classes
    ROC<-multiclass.roc(peptide_binding_category[training_folds==i],
                        peptide_random_forest_prediction)
    AUCs[1,i]<-auc(ROC$rocs[[1]][[1]])
  }
  return(AUCs)
}
##############################################################################
# Function for LOO TCR performance calculation
auc_loo_tcr <- function(TCR,weight,aa_matrix){
  library(readxl)
  # Read TCR data
  suppressWarnings(
    TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
  )
  
  # Load peptide list and index peptide
  peptide_list <- TCR_data$peptide[TCR_data$tcr_name==TCR]
  index_peptide <- unique(TCR_data$index_peptide[TCR_data$tcr_name==TCR])
  
  # Assign peptide binding category
  peptide_activity <- TCR_data$peptide_activity[TCR_data$tcr_name==TCR]
  
  #Normalize peptide activities to the index activity
  peptide_activity <- peptide_activity/(peptide_activity[
    which.max(peptide_activity)]) 
  
  #Initialize with all weak binders
  peptide_binding_category <- matrix(2,nrow=length(peptide_list),ncol=1)
  peptide_binding_category[peptide_activity<0.1] <- 1 #Non-binder
  peptide_binding_category[peptide_activity>=0.5] <- 3 #strong-binder
  
  #Create peptide features
  peptide_features <- generate_peptide_feature(peptide_list,
                                               index_peptide,aa_matrix)
  
  # Calculate ROC and AUC
  ROC<-multiclass.roc(peptide_binding_category[,1],
                      (peptide_features%*%weight)[,1],
                      direction="<")
  
  if (dim(table(peptide_binding_category)) ==3){#3-class
  AUC <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  } else {#2-class
  AUC <- auc(ROC$rocs[[1]])  
  }
  return(AUC)
}












