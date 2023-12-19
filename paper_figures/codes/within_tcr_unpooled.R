# Script taking a TCR index and a AA distance matrix index
#for within-TCR classification and regression and outputting correlation and AUC

rm(list=ls())#Clear environment
args = commandArgs(trailingOnly=TRUE)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

# Take AA distance function (1-70) and TCR index (1-66) as input (two integers)
AA_matrix_index <- as.integer(args[1])
TCR_index <- as.integer(args[2])

# Load functions
source("all_functions.R")

# Load, normalize, and discretize activation of selected TCR by mutants
activation_data <- load_tcr_activation(TCR_index)
# This also fully discards a class with <5 elements

# Load AA distance matrix
AA_distance <- load_AA_distance(AA_matrix_index)

# Generate position-dependent distance feature for mutant peptides
peptide_feature <- generate_peptide_feature(activation_data$peptide_list,
                                            activation_data$index_peptide,
                                            AA_distance)

# Perform regression and extract positional weights and correlation coefficients
# Results from n-fold cv stratified over TCR activation value deciles
# Seeds are for Bayes sampler and creating folds
regression_results <- regression_cv(peptide_feature,
                                    activation_data$peptide_activity,nfold=5,
                                    seed1=100+TCR_index, seed2=200+TCR_index)

# Perform classification and extract positional weights and AUCs
# Results from n-fold cv stratified over TCR activation categories
# Seeds are for Bayes sampler and creating folds

# Decide if the TCR has all 3 binding classes or 2 and classify accordingly
if(length(unique(activation_data$peptide_binding_category))==3){
  
classification_results <- classification_3_class_cv(peptide_feature,
                              activation_data$peptide_binding_category,nfold=5,
                              seed1=100+TCR_index, seed2=200+TCR_index)

} else {
  
classification_results <- classification_2_class_cv(peptide_feature,
                              activation_data$peptide_binding_category,nfold=5,
                              seed1=100+TCR_index, seed2=200+TCR_index)
}
##########################################################################
# Save outputs as RDA files
save(regression_results,classification_results,
     file=paste("../data/within_tcr_unpooled_outputs//output_",
                AA_matrix_index,"_",TCR_index,".rda",sep=""))
