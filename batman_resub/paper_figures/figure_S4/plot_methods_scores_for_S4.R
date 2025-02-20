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
  TCR_pMHC_data <- read_excel(paste0("../../tcr_epitope_datasets/",
  "mutational_scan_datasets/train_test_data_folds.xlsx"))
)

# Subset selected TCRs and columns
TCR_list_9mer <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                 'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11')

TCR_list_10mer <- c('TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va')

TCR_list_mhcii<- c('TCR-F5','TCR-3598-2','MBP-TCR','B3K508')

all_TCRs <- c(TCR_list_9mer,TCR_list_10mer,TCR_list_mhcii)

TCR_pMHC_data <- TCR_pMHC_data[TCR_pMHC_data$tcr %in% all_TCRs,]
TCR_pMHC_data <- select(TCR_pMHC_data,c('tcr','index_peptide','peptide','peptide_activity'))

# open empty df for storing AUC
score_all_methods <- data.frame(matrix(ncol=6,nrow=0, 
                                     dimnames=list(NULL, c('tcr','index_peptide',
                                                           'peptide','peptide_activity',
                                                           'score','method'))))
# stores scores for different methods

############
# ATM-TCR 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/atmtcr/',
                              'pred_original_atmtcr_input.csv'), 
                   header=FALSE,
                   sep="\t")
# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$V1

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$V1,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$V5

TCR_data_subset$method <- "ATM-TCR"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)
  
  
############
# AttnTAP 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/attntap/',
                             'attntap_output_mcpas.csv')) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$antigen

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$antigen,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction

TCR_data_subset$method <- "AttnTAP"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# DeepTR 
############

# UnZip files first
# get names of all the zip files
zipfile_list <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/deeptr/",
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

TCR_data_subset$method <- "DeepTR"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# epiTCR 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/epitcr/',
                             'epitcr_output_mhci.csv'),skip=14) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$epitope

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$epitope,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$predict_proba

TCR_data_subset$method <- "epiTCR"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)
############
# ERGO-II 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/ergo/',
                             'ergo_mcpas_output.csv')) 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Score

TCR_data_subset$method <- "ERGO-II"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# HeteroTCR 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/heterotcr/',
                             'heterotcr_output.tsv'),sep='\t') 
#mcpas performs better than vdjdb

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$probability

TCR_data_subset$method <- "HeteroTCR"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
#   IEDB 
############

# MHCI data
#load data from all CSV files
output_files <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/iedb/",
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

TCR_data_subset$method <- "IEDB"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

# MHC II
#load data from all txt files
output_files <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/iedb/",
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

TCR_data_subset$method <- "IEDB"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# ImRex 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/imrex/',
                             'imrex_output_predictions.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$antigen.epitope

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$antigen.epitope,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction_score

TCR_data_subset$method <- "ImRex"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# iTCep 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/itcep/',
                             'itcep_output_with_tcr.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Probability

TCR_data_subset$method <- "iTCep"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

############
# NetTCR-2.2 
############

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/nettcr/',
                             'nettcr_predictions.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$prediction

TCR_data_subset$method <- "NetTCR-2.2"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

################
#   NetTepi-1.0 
###############

#load data from all CSV files
output_files <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/nettepi/",
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

TCR_data_subset$method <- "NetTepi-1.0"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

################
# pMTnet_Omni 
###############

# With full TCR seq data, using pMTnet_Omni
#load data from all CSV files
output_files <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/pmtnet/",
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

TCR_data_subset$method <- "pMTnet-Omni"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

# With only CDR3 seq data, using pMTnet v1

scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/pmtnet/',
                             'pMTnet_v1_20250105190604__prediction.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$Antigen

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$Antigen,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$Rank

TCR_data_subset$method <- "pMTnet-Omni"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

################
## PRIME-2.1 ####
################

#load data from all txt files
output_files <- list.files(path = "../figure_2/2a/tcr_pmhc_methods_predictions/prime/",
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

TCR_data_subset$method <- "PRIME-2.1"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

################
# TCRPrediction 
################

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/tcrprediction/',
                             'tcrprediction_input_peptides_entire_cross_newemb.csv')) 

# set rownames to peptide (unique in data)
rownames(scores_df) <- scores_df$peptide

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% scores_df$peptide,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$pred1

TCR_data_subset$method <- "TCRPrediction"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)


################
# TITAN 
################

# Load peptide order
peptide_ID <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/titan/',
                              'epitope_seq_ID.csv'),
                       sep = '\t',header = FALSE)
peptide_list <- peptide_ID$V1

#load data
scores_df <- read.csv(paste0('../figure_2/2a/tcr_pmhc_methods_predictions/titan/',
                             'titan_output_converted.csv')) 

# set rownames to peptide
rownames(scores_df) <- peptide_list

# Subset of TCRs with scores available
TCR_data_subset <- TCR_pMHC_data[TCR_pMHC_data$peptide %in% peptide_list,]

#extract score
TCR_data_subset$score <- scores_df[TCR_data_subset$peptide,]$X0

TCR_data_subset$method <- "TITAN"

# Add to data
score_all_methods <- rbind(score_all_methods,TCR_data_subset)

# Add MHC info
score_all_methods$MHC_type <- 'MHCI'
score_all_methods$MHC_type[score_all_methods$tcr %in% TCR_list_mhcii] <- 'MHCII'



###################################
######## Plotting #################
###################################

p <- ggplot(score_all_methods, aes(y=score,
                                 x=peptide_activity,
                                 )) + 
  geom_point(aes(color=factor(tcr,levels = all_TCRs)), 
             alpha=0.85,
             size=0.5) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)) +
  labs(y="TCR-pMHC score",
       x="Normalized TCR activation",
       color="") +
  theme_cowplot() +
  facet_wrap(vars(method),scales = "free")+
  theme(strip.text = element_text(
    size = 20, color = "black")) +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=20)) +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(size=20))
#  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "figS4.pdf",
       width=15, height=15)

# Save raw data
write.csv(score_all_methods, "raw_data_fig_S4.csv")


















