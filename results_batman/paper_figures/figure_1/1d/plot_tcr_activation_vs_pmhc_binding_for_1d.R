rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)
library(viridis)
library(loo)
library(pROC)

############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to current

# List of 9-mer HLA-A*02:01 restricted index peptides with all positions mutagenized
index_peptide_list <- c('NLVPMVATV','SLYNTVATL','SLLMWITQC','IMDQVPFSV','LLFGYPVYV',
                        'ALYDKTKRI','YIMSDSNYV','SLFLGILSV','FMNKFIYEI')

# Get TCR activation data
tcr_data <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                            '/TCR_pMHCI_mutational_scan_database.xlsx'))

# Filter by index peptide
tcr_data <- tcr_data[tcr_data$index_peptide %in% index_peptide_list,]

# Assign TCR activation category
tcr_data$tcr_activation <- "Weak" #Initialize
tcr_names <- unique(tcr_data$tcr) # List all TCRs

# Loop over TCRs
for (tcr in tcr_names) {
  
  # Extract peptides and their activities for a particular TCR
  peptide_list <- tcr_data$peptide[tcr_data$tcr==tcr]
  peptide_activity <- tcr_data$peptide_activity[tcr_data$tcr==tcr]
  
  # Normalize peptide activities to the maximum activity
  peptide_activity <- peptide_activity/max(peptide_activity)
  
  # Indicate TCR activation categories
  #Initialize with all weak binders
  activation <- matrix("Weak",nrow=length(peptide_list),ncol=1)
  activation[peptide_activity<0.1] <- "No" #Non-binder
  activation[peptide_activity>=0.5] <- "Strong" #strong-binder
  
  # Add to TCR data
  tcr_data$tcr_activation[tcr_data$tcr==tcr] <- activation
}


# Get mutation location
# First flatten WT and mutated seqs, then do elementwise comparison
mutation_location <- unlist(strsplit(paste0(tcr_data$index_peptide,collapse=''),'')) !=
                     unlist(strsplit(paste0(tcr_data$peptide,collapse=''),''))
# Multiply by digits array to assign location and reshape
mutation_location <- matrix(mutation_location*rep(c(1:9),
                                                  length(mutation_location)/9),
                            ncol = 9, byrow = TRUE)
# Sum to get the nonzero mutation location
mutation_location <- rowSums(mutation_location)
tcr_data$mutation_location <- mutation_location

# Load NetMHCPan predictions
pmhc_data <- read.csv(paste0('../../../tcr_epitope_datasets/netmhcpan_predictions',
                              '/netMHCpan_formatted_predictions.csv'))

# subset by HLA-A*02:01
pmhc_data <- pmhc_data[pmhc_data$MHC=='HLA-A*0201',]

# Set Peptides as row names
rownames(pmhc_data) <- pmhc_data$Peptide

# Change missing BindLevel to "NB"
pmhc_data$BindLevel[pmhc_data$BindLevel==""] <- "NB"

# Get MHC binding for desired mutant peptides
tcr_data$BindLevel <- pmhc_data[tcr_data$peptide,'BindLevel']

# Remove WT peptides
tcr_data <- tcr_data[tcr_data$mutation_location != 0,]

# Count number of data entries for different pMHC binding and TCR activation
counts <- tcr_data %>% dplyr::count(BindLevel, tcr_activation)

############
## Plot ####
############

p <- ggplot(tcr_data,aes(y=factor(tcr_activation,levels=c('No','Weak','Strong')),
                   x=factor(BindLevel,levels=c('NB','WB','SB')))) + 
  geom_jitter(aes(color=factor(mutation_location)),size=0.5,alpha=0.5)+
  
   scale_color_manual(values=c('#b30000','#EC5800','#0096FF','#5D3FD3',
                             '#B2BEB5','#66a61e','#50C878','#0047AB','#a6761d')) +
                               
                                         labs(x="Predicted MHC binding",
                                              y="TCR activation",
                                              color="")+
  theme_cowplot() +
  theme(axis.text.x = element_text(size=35),
        axis.title.y = element_text(size=35),
        axis.title.x = element_text(size=35),
        axis.text.y = element_text(size=35)) +
  
  annotate("text", x=counts$BindLevel, y=counts$tcr_activation, 
           label= counts$n,size=7)

ggsave(plot=p, "tcr_activation_mhc_binding_0201.pdf",
       width=7.5, height=6)

write.csv(tcr_data,'raw_data_fig_1d.csv')




