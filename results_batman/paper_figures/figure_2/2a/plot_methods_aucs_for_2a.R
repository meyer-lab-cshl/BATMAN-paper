rm(list=ls())#Clear environment

# Packages
library(readxl)
library(pROC)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)



setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_pMHC_data <- read_excel(paste0("../../../tcr_epitope_datasets/",
  "mutational_scan_datasets/train_test_data_folds.xlsx"))
)

# Subset selected TCRs and columns
TCR_list_9mer <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11')

TCR_list_10mer <- c('TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va')

TCR_list_mhcii<- c('TCR-F5','TCR-3598-2','MBP-TCR','B3K508')

all_TCRs <- c(TCR_list_9mer,TCR_list_10mer,TCR_list_mhcii)

TCR_pMHC_data <- TCR_pMHC_data[TCR_pMHC_data$tcr %in% all_TCRs,]
TCR_pMHC_data <- select(TCR_pMHC_data,c('tcr','index_peptide','peptide','activation'))

# open empty df for storing AUC
AUC_all_methods <- data.frame(matrix(ncol=3,nrow=0, 
                                     dimnames=list(NULL, c("tcr", "method", 
                                                           "auc"))))

# calculate AUCs for different methods

############
# ATM-TCR 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/atmtcr/',
                              'pred_original_atmtcr_input.csv'), 
                   header=FALSE,
                   sep="\t")
# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$V1

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$V1,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$V5

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "ATM-TCR", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# AttnTAP 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/attntap/',
                             'attntap_output_mcpas.csv')) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$antigen

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$antigen,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "AttnTAP", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# DeepTR 
############

# UnZip files first
# get names of all the zip files
zipfile_list <- list.files(path = "tcr_pmhc_methods_predictions/deeptr/",
                   pattern = "*.zip", full.names = TRUE)

# unzip all files and extract data
scores_df <- data.frame()

for (zipfile in zipfile_list){
  
# Unzip file
unzip(zipfile)

# Read data
file_data <- read.csv('DeepTRP_predictions.tsv',sep=",") 
scores_df <- rbind(scores_df, file_data)
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Antigens

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Antigens,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Score_pMHC.TCR

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "DeepTR", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}


############
# epiTCR 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/epitcr/',
                             'epitcr_output_mhci.csv'),skip=14) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$epitope

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$epitope,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$predict_proba

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "epiTCR", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# ERGO-II 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/ergo/',
                             'ergo_mcpas_output.csv')) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Score

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "ERGO-II", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# HeteroTCR 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/heterotcr/',
                             'heterotcr_output.tsv'),sep='\t') 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$probability

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "HeteroTCR", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
#   IEDB 
############

# MHCI data
#load data from all CSV files
output_files <- list.files(path = "tcr_pmhc_methods_predictions/iedb/",
           pattern = "iedb_output_peptides_.*.csv", full.names = TRUE)

scores_df <- data.frame()

for (file in output_files){
  # Read data
  file_data <- read.csv(file) 
  scores_df <- rbind(scores_df, file_data)
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$immunogenicity.score

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "IEDB", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

# MHC II
#load data from all txt files
output_files <- list.files(path = "tcr_pmhc_methods_predictions/iedb/",
                           pattern = "iedb_output_peptides_.*.txt",
                           full.names = TRUE)

scores_df <- data.frame()

for (file in output_files){
  # Read data
  file_data <- read.csv(file,sep="\t",skip=3) 
  scores_df <- rbind(scores_df,file_data)
  scores_df <- unique(scores_df) #drop duplicate
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$score

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "IEDB", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# ImRex 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/imrex/',
                             'imrex_output_predictions.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$antigen.epitope

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$antigen.epitope,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction_score

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "ImRex", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# iTCep 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/itcep/',
                             'itcep_output_with_tcr.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Probability

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "iTCep", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

############
# NetTCR-2.2 
############

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/nettcr/',
                             'nettcr_predictions.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "NetTCR-2.2", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}


################
#   NetTepi-1.0 
###############

#load data from all CSV files
output_files <- list.files(path = "tcr_pmhc_methods_predictions/nettepi/",
                           pattern = "nettepi_output_.*.xls",
                           full.names = TRUE)

scores_df <- data.frame()

for (file in output_files){
  # Read data
  file_data <- read.csv(file,sep="\t") 
  scores_df <- rbind(scores_df, file_data)
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Tcell

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "NetTepi-1.0", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

################
# pMTnet_Omni 
###############

# With full TCR seq data, using pMTnet_Omni
#load data from all CSV files
output_files <- list.files(path = "tcr_pmhc_methods_predictions/pmtnet/",
                           pattern = "pMTnet_2025.*.csv",
                           full.names = TRUE)

scores_df <- data.frame()

for (file in output_files){
  # Read data
  file_data <- read.csv(file) 
  scores_df <- rbind(scores_df, file_data)
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$logit

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "pMTnet-Omni", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

# With only CDR3 seq data, using pMTnet v1

scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/pmtnet/',
                             'pMTnet_v1_20250105190604__prediction.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Antigen

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Antigen,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Rank

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction=">") #direction flipped since it is rank
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "pMTnet-Omni", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

################
## PRIME-2.1 ####
################

#load data from all txt files
output_files <- list.files(path = "tcr_pmhc_methods_predictions/prime/",
                           pattern = "prime_output_peptides_.*.txt",
                           full.names = TRUE)

scores_df <- data.frame()

for (file in output_files){
  # Read data
  file_data <- read.csv(file,sep="\t",skip=11)
  
  # Rename cols to harmonize data
  colnames(file_data) <- c("peptide","Rank_bestAllele","Score_bestAllele",
                           "RankBinding_bestAllele","BestAllele","Rank","Score",
                           "RankBinding")
    
  scores_df <- rbind(scores_df,file_data)
}

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Score

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "PRIME-2.1", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

################
# TCRPrediction 
################

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/tcrprediction/',
                             'tcrprediction_input_peptides_entire_cross_newemb.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$pred1

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "TCRPrediction", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

################
# TITAN 
################

# Load peptide order
peptide_ID <- read.csv(paste0('tcr_pmhc_methods_predictions/titan/',
                              'epitope_seq_ID.csv'),
                       sep = '\t',header = FALSE)
peptide_list <- peptide_ID$V1

#load data
scores_df <- read.csv(paste0('tcr_pmhc_methods_predictions/titan/',
                             'titan_output_converted.csv')) 

# set rownames to peptide
rownames(scores_df) <- peptide_list

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% peptide_list,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$X0

# Calculate AUC for each TCR
for (tcr in unique(TCR_data_subset$tcr)) {
  
  # Get ROC
  ROC <- multiclass.roc(as.factor(TCR_data_subset$activation[TCR_data_subset$tcr==tcr]),
                        TCR_data_subset$score[TCR_data_subset$tcr==tcr],
                        direction="<")
  # Mean AUC for 3 class-pairs
  auc <- mean(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  # Store data
  newdata <- data.frame(tcr = tcr, method = "TITAN", auc = auc)
  AUC_all_methods <- rbind(AUC_all_methods,newdata)
}

# Add MHC info

AUC_all_methods$MHC_type <- 'MHCI'
AUC_all_methods$MHC_type[AUC_all_methods$tcr %in% TCR_list_mhcii] <- 'MHCII'

# Sort methods by AUC
methods_mean <- AUC_all_methods %>%
  # Specify group indicator, column, function
  group_by(method) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(auc),
               list(mean_auc = mean))
methods_mean <- methods_mean[order(methods_mean$mean_auc,decreasing=TRUE),]



###################################
######## Plotting #################
###################################

p <- ggplot(AUC_all_methods, aes(y=auc,
                                 x=factor(method,levels = methods_mean$method),
                                 )) + 
  geom_boxplot(aes(color = MHC_type,fill = MHC_type),
               outlier.colour = NA,
               show.legend = FALSE, alpha=1,
               position=position_dodge2(0.8,preserve = "single")) +
  scale_color_manual(values=c('#5D3FD3','#66a61e')) +
  scale_fill_manual(values=c("#beacdf", "#bed8c0")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels = all_TCRs), 
                 group=MHC_type), 
             alpha=0.85,
             size=3,
             position = position_jitterdodge(jitter.width=0.08,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)) +
  labs(y="TCR activation prediction AUC",
       x="Pretrained TCR-pMHC method",
       color="") +
  theme_cowplot() +
  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "methods_aucs_for_2a.pdf",
       width=10, height=5)

# Save raw data
write.csv(AUC_all_methods, "raw_data_fig_2a.csv")


















