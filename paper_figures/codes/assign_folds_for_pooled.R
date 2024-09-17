# Function to distribute TCR-pMHC pairs training folds identical to unpooled tasks
# The output spreadsheet is used by python scripts for pooled Bayesian inference

rm(list=ls())#Clear environment
library(loo)
library(writexl)
library(readxl)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

# set empty dataframe
output_data <- data.frame(matrix(ncol = 20, nrow = 0))
colnames(output_data) <- c("TCR_name","TCR_index","peptide_list",
                           "peptide_activity","index_peptide", 
                           "peptide_binding_category","training_folds",
                         paste(rep("feature.",11),seq(1,11,1),sep=""),
                         "aa_change_symm","aa_change_full")
# Extract TCR names
# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# Extract TCR names present in the data
TCR_names <- unique(TCR_data$tcr_name)

# Loop over all TCRs
for (TCR_index in 1:66){
  print(TCR_index)

# Load functions
source("all_functions.R")

# Load, normalize, and discretize activation of selected TCR by mutants
activation_data <- load_tcr_activation(TCR_index)
# This also fully discards a class with <5 elements

# Split data into 5folds
set.seed(200+TCR_index)
training_folds=kfold_split_stratified(K = 5, 
                                   x = activation_data$peptide_binding_category)
activation_data$training_folds <- training_folds-1 #-1 for python

# repeat index peptide and TCR name and index to make a column
activation_data$index_peptide <- rep(activation_data$index_peptide, 
                                    each = length(activation_data$peptide_list))
activation_data$TCR_index <- rep(TCR_index, 
                                 each = length(activation_data$peptide_list))

activation_data$TCR_name <- rep(TCR_names[TCR_index], 
                                 each = length(activation_data$peptide_list))

# Store BLOSUM100 peptide features
# Load AA distance matrix
AA_distance <- load_AA_distance(16) #BLOSUM100 has index 16

# Initialize all features with NaN
activation_data$feature <- array(NaN,
                                 dim=c(length(activation_data$peptide_list),11))
# Generate position-dependent distance feature for mutant peptides
peptide_feature <- generate_peptide_feature(activation_data$peptide_list,
                                            activation_data$index_peptide[1],
                                            AA_distance)
activation_data$feature[,1:dim(peptide_feature)[2]] <- (-1)*peptide_feature
#(For python, we do not need the minus sign)

activation_data$peptide_binding_category <- 
                                    activation_data$peptide_binding_category -1
#(For python)

# Indexes for AA change in a mutant (encodes which AA is substituted by which)
aa_change_symm <- matrix(NA,nrow=length(activation_data$peptide_list),ncol=1)
aa_change_full <- matrix(NA,nrow=length(activation_data$peptide_list),ncol=1)

#Order of AAs
AA_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
             'T','V','W','Y')

for (peptide in 1:length(activation_data$peptide_list)) {
  #Skip index peptide
  if (activation_data$index_peptide[1]!=activation_data$peptide_list[peptide]){
  # Location of mismatch
  mismatch_position <- which((unlist(strsplit(
    activation_data$index_peptide[1],split=""))!=
    unlist(strsplit(activation_data$peptide_list[peptide],split="")))==TRUE)
  
  initial_aa <-  unlist(strsplit(
                  activation_data$index_peptide[1],split=""))[mismatch_position]
  
  final_aa <-  unlist(strsplit(
    activation_data$peptide_list[peptide],split=""))[mismatch_position]
  
  
  aa_change_symm[peptide,1] <- 20*(min(which(AA_list==initial_aa),
                            which(AA_list==final_aa))-1) + 
                       (max(which(AA_list==initial_aa),
                            which(AA_list==final_aa))-1)
  
  aa_change_full[peptide,1] <- 20*(which(AA_list==initial_aa)-1) + 
                                  (which(AA_list==final_aa)-1)
  
  
  }
}
activation_data$aa_change_symm <- aa_change_symm
activation_data$aa_change_full <- aa_change_full

# add to dataframe
activation_data <- data.frame(activation_data)
output_data <- rbind.data.frame(output_data,activation_data)
}

# Write to spreadsheet
filename<-paste0("../data/training_folds.xlsx")
write_xlsx(output_data,filename,col_names=TRUE)