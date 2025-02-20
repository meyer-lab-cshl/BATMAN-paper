rm(list=ls())#Clear environment

library(readxl)
library(dplyr)

# source pTEAM functions
source("../../../batman_functions/pTEAM_functions.R")

# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  tcr_pmhc_data <- read_excel(paste0("../../../tcr_epitope_datasets/",
                        "mutational_scan_datasets/train_test_data_folds.xlsx"))
)

# Subset selected index peptides
index_peptide_list <- c('FRDYVDRFYKTLRAEQASQE', 'FEAQKAKANKAVD', 'SIINFEKL',
                        'TPQDLNTML', 'VPSVWRSSL', 'FMNKFIYEI', 'SLLMWITQC',
                        'IMDQVPFSV', 'NLVPMVATV', 'VVVGAVGVGK')

tcr_pmhc_data <- tcr_pmhc_data %>% 
  filter(index_peptide %in% index_peptide_list)

# Subset for TCRs with CDR3a/b sequences
tcr_pmhc_data <- tcr_pmhc_data %>% 
  filter(!is.na(cdr3b))

# For Ova peptide, subset for educated TCRs and OT1
tcr_pmhc_data <- tcr_pmhc_data %>% 
  filter(index_peptide != 'SIINFEKL' | 
           grepl('Ed', tcr) | 
           tcr == 'OT1')

# Empty dataframe to store aucs
auc_loo_no_tcr_seq <- data.frame(tcr = character(), 
                                        index_peptide = character(), 
                                        auc = numeric(), 
                                        stringsAsFactors = FALSE)

auc_loo_with_tcr_seq <- data.frame(tcr = character(), 
                                 index_peptide = character(), 
                                 auc = numeric(), 
                                 stringsAsFactors = FALSE)
# Remove unmutated peptides
tcr_pmhc_data <- tcr_pmhc_data[tcr_pmhc_data$peptide!=tcr_pmhc_data$index_peptide,]


##########################################################
# Run LOO-TCR predictions for BATMAN
##########################################################
for (index_peptide in index_peptide_list){
  
  # Subset to antigen-specific TCR repertoire
  tcr_data <- tcr_pmhc_data[tcr_pmhc_data$index_peptide==index_peptide,]
  
  # LOO TCR loop
  for (tcr in unique(tcr_data$tcr)){
    
    # LOO train and test data
    train_data <- tcr_data[tcr_data$tcr!=tcr,]
    test_data <- tcr_data[tcr_data$tcr==tcr,]
    
    # Train RF for every TCR with and without using TCR sequence
    auc_no_tcr_seq <- train_pTEAM(train_data,test_data,seed = 111)
    
    auc_with_tcr_seq <- train_pTEAM(train_data,test_data,
                                    use_tcr=TRUE,seed = 111)
    
    # Add to df
    new_data_no_tcr_seq <- data.frame(tcr=tcr,index_peptide=index_peptide,
                                      auc=auc_no_tcr_seq)
    
    new_data_with_tcr_seq <- data.frame(tcr=tcr,index_peptide=index_peptide,
                                      auc=auc_with_tcr_seq)
    
    auc_loo_no_tcr_seq <- rbind(auc_loo_no_tcr_seq,new_data_no_tcr_seq)
    auc_loo_with_tcr_seq <- rbind(auc_loo_with_tcr_seq,new_data_with_tcr_seq)
    
  }
}

# Merge and save data
auc_loo_no_tcr_seq$method <- 'pTEAM'
auc_loo_with_tcr_seq$method <- 'pTEAM+TCR'

auc_all <- rbind(auc_loo_no_tcr_seq,auc_loo_with_tcr_seq)

# export
write.csv(auc_all,file='pTEAM_within_repertoire_LOO_AUCs.csv')



