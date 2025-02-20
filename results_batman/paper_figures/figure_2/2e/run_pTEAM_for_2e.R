rm(list=ls())#Clear environment
args = commandArgs(trailingOnly=TRUE)

library(readxl)
library(dplyr)


# Inputs
rand_seed <- as.integer(args[1]) #random seed
fpeptide <- as.integer(args[2]) #10*fraction of peptides per mutation site per TCR
fold <- as.integer(args[3]) #test fold number

# source pTEAM functions
source("../../../batman_functions/pTEAM_functions.R")

# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel(paste0("../../../tcr_epitope_datasets/",
                        "mutational_scan_datasets/train_test_data_folds.xlsx"))
)

# Subset selected TCRs
TCR_list_9mer <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11') 

TCR_list_10mer <- c('TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va')

TCR_list_mhcii <- c('TCR-F5','TCR-3598-2','MBP-TCR','B3K508')

TCR_list <- c(TCR_list_9mer,TCR_list_10mer,TCR_list_mhcii)

TCR_data <- TCR_data[TCR_data$tcr %in% TCR_list,]

# Remove unmutated peptides
TCR_data <- TCR_data[TCR_data$peptide!=TCR_data$index_peptide,]


##########################################################
# Within TCR fold
##########################################################
AUC_within_TCR <- array(0,dim=c(length(TCR_list),1))# Array for AUCs for all TCRs

for (TCR_index in 1:length(TCR_list)){
  
  # Subsample train data, stratify on mutation position
  train_data <- TCR_data[(TCR_data$tcr==TCR_list[TCR_index]) &
                           (TCR_data$training_fold!=fold),]
  set.seed(rand_seed)
  train_data <- train_data %>%
    group_by(mutation_position) %>%
    sample_frac(fpeptide/10)
  
  
  # Test data
  test_data <- TCR_data[(TCR_data$tcr==TCR_list[TCR_index]) &
                          (TCR_data$training_fold==fold),]

  # Train RF for every TCR without using TCR sequence
  AUC_within_TCR[TCR_index] <- train_pTEAM(train_data,test_data,seed = 111)
}

# Make a DataFrame
AUC_within_TCR <- data.frame(tcr = TCR_list,
                             AUC=AUC_within_TCR)

# export
write.csv(AUC_within_TCR,
          file=paste0("pteam_outputs/aucs_pteam_",rand_seed,"_",
                        fpeptide,"_",fold,".csv"))



