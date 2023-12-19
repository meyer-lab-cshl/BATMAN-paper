rm(list=ls())#Clear environment
library(readxl)
library(pROC)


tcr_list <- c(6,12,52,55,58,59,60,62,64,65,68)
aucs_msk <- array(0,dim=c(1,length(tcr_list)))

for (tcr in 1:length(tcr_list)){
  
TCR_index <- tcr_list[tcr]

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
peptide_activity<-peptide_activity/(peptide_activity[which.max(peptide_activity)])

# Indicate peptide binding categories
peptide_binding_category <- matrix(2,nrow=length(peptide_list),ncol=1)#Initialize with all weak binders
peptide_binding_category[peptide_activity<0.1]=1 #Non-binder
peptide_binding_category[peptide_activity>=0.5]=3 #strong-binder
peptide_binding_category<- as.integer(peptide_binding_category)

# load distances from the index peptide
epitope_distance_data <- read.csv("C:/Users/amita/Downloads/epitope_analysis/antigen-availability-master/antigen-availability-master/TCR-crossreactivity/final_codes/local_outputs/misc_methods/msk/NeoantigenEditing-main/NeoantigenEditing-main/distances.csv")
a<-epitope_distance_data$seq
index_peptide_distances <- matrix(0,nrow=length(peptide_list),ncol=1)
for (peptide in 1:length(peptide_list)) {
  index_peptide_distances[peptide,1]<-(-1)*epitope_distance_data[which(epitope_distance_data$seq%in%index_peptide),1+which(epitope_distance_data$seq%in%peptide_list[peptide])]
}
  
#Calculate AUROCs on the test sets of each fold

  ROC<-multiclass.roc(peptide_binding_category,index_peptide_distances[,1],direction="<")
  AUROCs<-c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))

aucs_msk[1,tcr] <- mean(AUROCs)
}

save(aucs_msk,tcr_list,
     file=paste0("C:/Users/amita/Downloads/epitope_analysis/",
                 "antigen-availability-master/antigen-availability-master/",
                 "TCR-crossreactivity/final_codes/msk_outputs.rda"))