rm(list=ls())#Clear environment

library(readxl)
library(pROC)
library(randomForest)
library(muscle)

setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is


# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# Discard TCRs without CDR3a/b seqs
TCR_data <- TCR_data[TCR_data$cdr3a!="NA" & TCR_data$cdr3b!="NA",]

# Filter data
# Use only all 9-mer-binding TCRs
peptide_length <- 9
TCR_data <- TCR_data[nchar(TCR_data$index_peptide)==peptide_length,]

# Filter based on index peptide
index_peptide_selected <- 'TPQDLNTML'
TCR_data <- TCR_data[TCR_data$index_peptide==index_peptide_selected,]

# Extract TCR names present in the data
TCR_names <- unique(TCR_data$tcr_name)

# Extract and align TCR CDR3a/b seqs for selected TCRs
cdr3a_list <- rep('X',length(TCR_names))
cdr3b_list <- rep('X',length(TCR_names))

for (tcr in 1:length(TCR_names)){
  cdr3a_list[tcr] <- unique(TCR_data$cdr3a[TCR_data$tcr_name==TCR_names[tcr]])
  cdr3b_list[tcr] <- unique(TCR_data$cdr3b[TCR_data$tcr_name==TCR_names[tcr]])
}

cdr3a_list <- as.character(muscle(stringset = AAStringSet(x=cdr3a_list, 
                                                          start=NA, 
                                                          end=NA, 
                                                          width=NA, 
                                                          use.names=FALSE)))

cdr3b_list <- as.character(muscle(stringset = AAStringSet(x=cdr3b_list, 
                                                          start=NA, 
                                                          end=NA, 
                                                          width=NA, 
                                                          use.names=FALSE)))

# Extract all peptides for all selected TCRs to get the total number of peptides
n_peptide <- dim(TCR_data)[1]

#aligned cdr3a, cdr3b length
cdr3a_length <- length(unlist(strsplit(cdr3a_list[1],split="")))
cdr3b_length <- length(unlist(strsplit(cdr3b_list[1],split="")))
cdr3_length <- cdr3a_length + cdr3b_length

# Create arrays to store tcr-peptide info for training and testing

#array to store CRD3a/b and peptide Atchley embeddings
tcr_peptide_embedding <- matrix(0,nrow=n_peptide,
                                ncol=(10*peptide_length+5*cdr3_length));

peptide_binding_category <- as.integer(matrix(0,nrow=n_peptide,ncol=1));

#which TCR the peptide belongs to
tcr_number <- as.integer(matrix(0,nrow=n_peptide,ncol=1));

# generate train/test data by looping over TCRs
peptide_count <- 0
for (TCR_index in 1:length(TCR_names)) {
  
  # Extract peptides and their activities for a particular TCR
  peptide_list <- TCR_data$peptide[TCR_data$tcr_name==TCR_names[TCR_index]]
  peptide_activity <- TCR_data$peptide_activity[
                                        TCR_data$tcr_name==TCR_names[TCR_index]]
  
  # Identify the central(/index) peptide
  central_peptide <- unique(TCR_data$index_peptide
                            [TCR_data$tcr_name==TCR_names[TCR_index]])
  
  #Normalize peptide activities to the maximum activity
  peptide_activity <- peptide_activity/
                                 (peptide_activity[which.max(peptide_activity)])
  
  # Indicate peptide binding categories
  #Initialize with all weak binders
  binding_category <- matrix(2,nrow=length(peptide_list),ncol=1)
  binding_category[peptide_activity<0.1] <- 1 #Non-binder
  binding_category[peptide_activity>=0.5] <- 3 #strong-binder
  binding_category <- as.integer(binding_category)
  
  ## Atchley factor classifier ##
  
  # Read Atchley Factors
  AA_list_Atchley <- c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F',
                       'P','S','T','W','Y','V','-')
  Atchley_factors <- as.matrix(read.csv("../data/Atchley.csv")[,2:6])
  Atchley_factors <- rbind(Atchley_factors,matrix(0,nrow=1, ncol=5))# Gap is padded as zero
  
  # Generate embeddings for peptides with the Atchley factors
  embeddings <- matrix(0,nrow=length(peptide_list),
                         ncol=(10*peptide_length+5*cdr3_length))
  
  # peptide embeddings
  for (peptide in 1:length(peptide_list)) {
    for (position in 1:peptide_length) {
      embeddings[peptide,(position*5-4):(position*5)] <- 
        Atchley_factors[which(AA_list_Atchley==
                  unlist(strsplit(peptide_list[peptide],split=""))[position]),]
      
      embeddings[peptide,(peptide_length*5+position*5-4):
                         (peptide_length*5+position*5)] <- 
        
      Atchley_factors[which(AA_list_Atchley==
                  unlist(strsplit(peptide_list[peptide],split=""))[position]),]-
        
      Atchley_factors[which(AA_list_Atchley==
                        unlist(strsplit(central_peptide,split=""))[position]),]
    }
  }
  
  # CDR3a embeddings
  for (position in 1:cdr3a_length) {
    embeddings[,(peptide_length*10+position*5-4):
                 (peptide_length*10+position*5)] <-
     Atchley_factors[which(AA_list_Atchley==
                   unlist(strsplit(cdr3a_list[TCR_index],split=""))[position]),]
  }
  
  # CDR3b embeddings
  for (position in 1:cdr3b_length) {
    embeddings[,(peptide_length*10+cdr3a_length*5+position*5-4):
                 (peptide_length*10+cdr3a_length*5+position*5)] <- 
      Atchley_factors[which(AA_list_Atchley==
                   unlist(strsplit(cdr3b_list[TCR_index],split=""))[position]),]
  }
  
  # Store full data for the TCR
  tcr_peptide_embedding[(peptide_count+1):
                        (peptide_count+length(peptide_list)),]<- embeddings
  
  peptide_binding_category[(peptide_count+1):
                      (peptide_count+length(peptide_list))] <- binding_category
  
  tcr_number[(peptide_count+1):(peptide_count+length(peptide_list))]<- TCR_index
  
  peptide_count <- peptide_count + length(peptide_list)
}

# Set up random forest for LOO TCRs

loo_tcr <- array(0,dim=c(3,length(TCR_names)))# Array to store AUCs

for (i in 1:length(TCR_names)){
  print(i)
  #Classification trees
  set.seed(300+TCR_index)
  peptide_random_forest_classifier <- randomForest(
    tcr_peptide_embedding[tcr_number!=i,],
    as.factor(peptide_binding_category[tcr_number!=i])) #Train
  
  peptide_random_forest_prediction <- predict(
    peptide_random_forest_classifier, 
    tcr_peptide_embedding[tcr_number==i,],type="prob") #Predict for left out TCR
  
  ROC <- multiclass.roc(peptide_binding_category[tcr_number==i],
                      peptide_random_forest_prediction,direction=">") #ROC
  
  if (length(ROC$rocs)==1){ #Missing a class
    loo_tcr[1,i] <- (auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2
    
  } else {
  loo_tcr[,i] <- c((auc(ROC$rocs[[1]][[1]])+auc(ROC$rocs[[1]][[2]]))/2,
                   (auc(ROC$rocs[[2]][[1]])+auc(ROC$rocs[[2]][[2]]))/2,
                   (auc(ROC$rocs[[3]][[1]])+auc(ROC$rocs[[3]][[2]]))/2)
  }
}

save(TCR_names,loo_tcr,
     file=paste0("../data/loo_tcr/atchley/auc_loo_tcr_atchley_",
                 index_peptide_selected,".rda"))


