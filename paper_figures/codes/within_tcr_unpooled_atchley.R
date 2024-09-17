# Script taking a TCR index and a AA distance matrix index
#for within-TCR classification and regression and outputting correlation and AUC

rm(list=ls())#Clear environment
args = commandArgs(trailingOnly=TRUE)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

# Take TCR index (1-66) as input (one integer)
TCR_index <- 1 #as.integer(args[1])

# Load functions
source("all_functions.R")

# Load, normalize, and discretize activation of selected TCR by mutants
activation_data <- load_tcr_activation(TCR_index)
# This also fully discards a class with <5 elements

# Generate Atchley feature for mutant peptides
peptide_feature <-generate_peptide_feature_atchley(activation_data$peptide_list,
                                            activation_data$index_peptide)

# Perform regression and extract correlation coefficients
# Results from n-fold cv stratified over TCR activation value deciles
# Seeds are for RF and creating folds
correlations <- regression_cv_atchley(peptide_feature,
                                    activation_data$peptide_activity,nfold=5,
                                    seed1=100+TCR_index,seed2=200+TCR_index)

# Perform classification and extract AUCs
# Results from n-fold cv stratified over TCR activation categories
# Seeds are for RF and creating folds

# Decide if the TCR has all 3 binding classes or 2 and classify accordingly
if(length(unique(activation_data$peptide_binding_category))==3){
  
AUCs <- classification_3_class_cv_atchley(peptide_feature,
                              activation_data$peptide_binding_category,nfold=5,
                              seed1=300+TCR_index, seed2=200+TCR_index)

} else {
  
AUCs <- classification_2_class_cv_atchley(peptide_feature,
                              activation_data$peptide_binding_category,nfold=5,
                              seed1=300+TCR_index, seed2=200+TCR_index)
}
##########################################################################
# Save outputs as RDA files
save(correlations,AUCs,
     file=paste("../data/within_tcr_unpooled_outputs_atchley//output_",
                TCR_index,".rda",sep=""))
