rm(list=ls())#Clear environment
# Load required libraries
library(readxl)
library(R.matlab)
library(pROC)

# List of 3 class 9mer-binding TCRs
tcr_list <- c(6,12,52,55,58,59,60,62,64,65,68)

aucs_pmtnet <- array(0,dim=c(1,length(tcr_list)))


for (TCR in setdiff(1:11,c(8,9))){
  
  TCR_index <- tcr_list[TCR]
  
  # Load TCR and peptide information from the main TCR datafile
  TCR_data <- read_excel("C:/Users/amita/Downloads/epitope_analysis/antigen-availability-master/antigen-availability-master/TCR-crossreactivity/TCR_epitope_database.xlsx")
  
  # Extract TCR names present in the data
  TCR_names <- unique(TCR_data$tcr_name)
  
  # Extract peptides and their activities for a particular TCR
  peptide_list <- TCR_data$peptide[TCR_data$tcr_name==TCR_names[TCR_index]]
  peptide_activity <- TCR_data$peptide_activity[TCR_data$tcr_name==TCR_names[TCR_index]]
  
  # Define the index peptide as the one with the highest activity
  index_peptide <- TCR_data$peptide[TCR_data$tcr_name==TCR_names[TCR_index]][which.max(peptide_activity)]
  
  #Normalize peptide activities to the index activity
  peptide_activity <- peptide_activity/(peptide_activity[which.max(peptide_activity)])
  
  # Indicate peptide binding categories
  peptide_binding_category <- matrix(1,nrow=length(peptide_list),ncol=1)#Initialize with all weak binders
  peptide_binding_category[peptide_activity<0.1]<- 0 #Non-binder
  peptide_binding_category[peptide_activity>=0.5]<- 2 #strong-binder
  peptide_binding_category<- as.integer(peptide_binding_category)
  
  set.seed(200+TCR_index)
  nfold<- 5#number of training folds
  training_folds<- kfold_split_stratified(K = nfold, x = peptide_binding_category)#generate test indices for folds
  
  #Load file
  output_data <- read.csv(paste0("C:/Users/amita/Downloads/epitope_analysis/",
             "antigen-availability-master/antigen-availability-master/",
             "TCR-crossreactivity/final_codes/local_outputs/",
             "misc_methods/pmtnet//pMTnet_prediction.csv"))
  
  peptide_rank <- as.double(output_data$Rank[match(peptide_list,output_data$Antigen)]) #find peptide ranks
  
  #Calculate AUC
  
    ROC<-multiclass.roc(peptide_binding_category,peptide_rank,direction=">")
    AUROCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  
  aucs_pmtnet[1,TCR] <- mean(AUROCs)
}

save(aucs_pmtnet,file=paste0("C:/Users/amita/Downloads/epitope_analysis/",
                      "antigen-availability-master/antigen-availability-master/",
                      "TCR-crossreactivity/final_codes/local_outputs/",
                      "misc_methods/pmtnet/pmtnet_outputs.rda"))