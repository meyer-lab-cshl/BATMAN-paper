rm(list=ls())#Clear environment

library(readxl)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

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
# Within TCR 5-fold CV
##########################################################
AUC_within_TCR <- array(0,dim=c(length(TCR_list),5))# Array for AUCs for 5 folds

for (TCR_index in 1:length(TCR_list)){
  for (fold in 1:5){
    
    print(c(TCR_index,fold)) #Print progress
    
    # Train RF for every TCR and fold without using TCR sequence
    AUC_within_TCR[TCR_index,fold] <- train_pTEAM(
      TCR_data[(TCR_data$tcr==TCR_list[TCR_index]) & # Train set
                 (TCR_data$training_fold!=fold-1),],
      TCR_data[(TCR_data$tcr==TCR_list[TCR_index]) & # Test set
                 (TCR_data$training_fold==fold-1),],
      seed = 111)
  }
}
# Average over 5 folds
AUC_within_TCR <- rowMeans(AUC_within_TCR)

# Make a DataFrame
AUC_within_TCR <- data.frame(tcr = TCR_list,
                             AUC=AUC_within_TCR,
                             type='pTEAM_within')

###################################################################
# Leave-one-TCR-out tests without using TCR seq
###################################################################

# 9-mers
AUC_LOO_TCR_9mers <- array(0,dim=c(length(TCR_list_9mer)))# Array for AUCs
TCR_data_9mer <- TCR_data[TCR_data$tcr %in% TCR_list_9mer,]

# Loop over TCRs
for (TCR_index in 1:length(TCR_list_9mer)){
  
    print(c(TCR_index)) #Print progress
    
    # Train RF for LOO TCR without using TCR sequence
    AUC_LOO_TCR_9mers[TCR_index] <- train_pTEAM(
      TCR_data_9mer[TCR_data_9mer$tcr!=TCR_list_9mer[TCR_index],], #Train set
      TCR_data_9mer[TCR_data_9mer$tcr==TCR_list_9mer[TCR_index],], #Test set
      seed = 111)
    
}

# 10-mers
AUC_LOO_TCR_10mers <- array(0,dim=c(length(TCR_list_10mer)))# Array for AUCs
TCR_data_10mer <- TCR_data[TCR_data$tcr %in% TCR_list_10mer,]

# Loop over TCRs
for (TCR_index in 1:length(TCR_list_10mer)){
  
  print(c(TCR_index)) #Print progress
  
  # Train RF for LOO TCR without using TCR sequence
  AUC_LOO_TCR_10mers[TCR_index] <- train_pTEAM(
    TCR_data_10mer[TCR_data_10mer$tcr!=TCR_list_10mer[TCR_index],], #Train set
    TCR_data_10mer[TCR_data_10mer$tcr==TCR_list_10mer[TCR_index],], #Test set
    seed = 111)
}

# Store data
AUC_LOO_TCR_without_seq <- data.frame(tcr = c(TCR_list_9mer,TCR_list_10mer),
 AUC=c(AUC_LOO_TCR_9mers,AUC_LOO_TCR_10mers),
 type='pTEAM_LOO_without_seq')

###################################################################
# Leave-one-TCR-out tests with using TCR seq
###################################################################

# Discard TCRs without CDR3a/b seqs
TCR_data <- TCR_data[!is.na(TCR_data$cdr3b),]
TCR_data <- TCR_data[TCR_data$cdr3b!='NA',]

# 9-mers
TCR_data_9mer <- TCR_data[TCR_data$tcr %in% TCR_list_9mer,]
TCR_list_9mer <- unique(TCR_data_9mer$tcr)
AUC_LOO_TCR_9mers <- array(0,dim=c(length(TCR_list_9mer)))# Array for AUCs

# Loop over TCRs
for (TCR_index in 1:length(TCR_list_9mer)){
  
  print(c(TCR_index)) #Print progress
  
  # Train RF for LOO TCR without using TCR sequence
  AUC_LOO_TCR_9mers[TCR_index] <- train_pTEAM(
    TCR_data_9mer[TCR_data_9mer$tcr!=TCR_list_9mer[TCR_index],], #Train set
    TCR_data_9mer[TCR_data_9mer$tcr==TCR_list_9mer[TCR_index],], #Test set
    use_tcr = TRUE,
    seed = 111)
  
}

# 9-mers
TCR_data_10mer <- TCR_data[TCR_data$tcr %in% TCR_list_10mer,]
TCR_list_10mer <- unique(TCR_data_10mer$tcr)
AUC_LOO_TCR_10mers <- array(0,dim=c(length(TCR_list_10mer)))# Array for AUCs

# Loop over TCRs
for (TCR_index in 1:length(TCR_list_10mer)){
  
  print(c(TCR_index)) #Print progress
  
  # Train RF for LOO TCR without using TCR sequence
  AUC_LOO_TCR_10mers[TCR_index] <- train_pTEAM(
    TCR_data_10mer[TCR_data_10mer$tcr!=TCR_list_10mer[TCR_index],], #Train set
    TCR_data_10mer[TCR_data_10mer$tcr==TCR_list_10mer[TCR_index],], #Test set
    use_tcr = TRUE,
    seed = 111)
  
}

# Store data
AUC_LOO_TCR_with_seq <- data.frame(tcr = c(TCR_list_9mer,TCR_list_10mer),
  AUC=c(AUC_LOO_TCR_9mers,AUC_LOO_TCR_10mers),
  type='pTEAM_LOO_with_seq')

# Combine all data
pTEAM_AUCs <-rbind(AUC_within_TCR,AUC_LOO_TCR_without_seq,AUC_LOO_TCR_with_seq)

# export
write.csv(pTEAM_AUCs,file='pTEAM_AUCs.csv')



