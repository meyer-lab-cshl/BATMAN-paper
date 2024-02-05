# # Code to assign folds and create training data

rm(list=ls())#Clear environment
# Load required libraries
library(readxl)
library(writexl)
library(loo)
library(pROC)
library(dplyr)


#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# Extract TCR names present in the data
TCR_names <- unique(TCR_data$tcr_name)

# set empty dataframe
peptide_folds <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(peptide_folds) <- c("peptide",
                             "peptide_activity",
                             "index",
                             "activation",
                             "training_folds",
                             "tcr")


# Loop over all TCRs
for (TCR_index in 1:66){
  print(TCR_index)
  
  # Load functions
  source("../codes/all_functions.R")
  
  # Load, normalize, and discretize activation of selected TCR by mutants
  activation_data <- load_tcr_activation(TCR_index)
  # This also fully discards a class with <5 elements
  
  # Split data into 5folds
  set.seed(200+TCR_index)
  training_folds=kfold_split_stratified(K = 5,
                                        x = activation_data$peptide_binding_category)
  activation_data$training_folds <- training_folds-1 #for python
  
  
  # repeat TCR name to make a column
  activation_data$tcr <- rep(TCR_names[TCR_index],
                                  each = length(activation_data$peptide_list))
  
  activation_data$peptide_binding_category <- 
    activation_data$peptide_binding_category-1 #For python
  
  # Add to full data
  
  peptide_folds <- rbind(peptide_folds,data.frame(activation_data))
}

colnames(peptide_folds) <- c('peptide','activity','index_peptide','activation','fold','tcr')

# List of 3 class 9-mer-binding TCRs
tcr_names <- c("18A2","NYE-S1","868Z11","TCR3","TCR6","TCR7","A6",
               "T1","FLT3DY","A23","TCR2-T")

# Load peptide data for selected TCRs

peptide_folds <- select(peptide_folds[peptide_folds$tcr %in% tcr_names,],
                        c('tcr',
                          'index_peptide',
                          'peptide',
                          'activation',
                          'fold'))

# Write to csv
write.csv(peptide_folds, "fig_2_training_files\\full_training_data.csv",
          row.names = FALSE)
