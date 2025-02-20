rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)
library(readxl)
library(viridis)
library(ggseqlogo)
library(muscle)

#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Get TCR activation data
TCR_data <- read_excel(paste0('../../tcr_epitope_datasets/mutational_scan_datasets',
                              '/train_test_data_folds.xlsx'))


# Discard TCRs without CDR3b sequence
TCR_data_full <- TCR_data[!is.na(TCR_data$cdr3b),]

# List index_peptides
index_peptides <- unique(TCR_data_full$index_peptide)

# List for CDR sequences
cdr3a_list_full <- list()
cdr3b_list_full <- list()

# Filter based on index peptide
for (index_peptide_selected in index_peptides){
  TCR_data <- TCR_data_full[TCR_data_full$index_peptide==index_peptide_selected,]
  
  # Extract TCR names present in the data
  TCR_names <- unique(TCR_data$tcr)
  
  # Extract and align TCR CDR3a/b seqs for selected TCRs
  cdr3a_list <- rep('X',length(TCR_names))
  cdr3b_list <- rep('X',length(TCR_names))
  
  for (tcr in 1:length(TCR_names)){
    cdr3a_list[tcr] <- unique(TCR_data$cdr3a[TCR_data$tcr==TCR_names[tcr]])
    cdr3b_list[tcr] <- unique(TCR_data$cdr3b[TCR_data$tcr==TCR_names[tcr]])
  }
  
  title <- paste0(index_peptide_selected,' (n=',length(TCR_names),')')
  
  cdr3a_list_full[[title]] <- 
    as.character(muscle(stringset = AAStringSet(x=cdr3a_list, 
                                                start=NA, 
                                                end=NA, 
                                                width=NA, 
                                                use.names=FALSE)))
  
  cdr3b_list_full[[title]] <- 
    as.character(muscle(stringset = AAStringSet(x=cdr3b_list, 
                                                start=NA, 
                                                end=NA, 
                                                width=NA, 
                                                use.names=FALSE)))
}

#######################
#### Plot #############
#######################

# Create custom colour scheme
color_scheme = make_col_scheme(chars=c('R','H','K','D','E','W','F','Y','L','I','M','V','A',
                              'S','T','N','Q','C','G','P'),
                      
                      groups=c('charged','charged','charged','charged','charged',
                               'hydrophobic aromatic','hydrophobic aromatic','hydrophobic aromatic',
                               'hydrophobic non-aromatic','hydrophobic non-aromatic',
                               'hydrophobic non-aromatic','hydrophobic non-aromatic',
                               'hydrophobic non-aromatic',
                               'polar uncharged','polar uncharged',
                               'polar uncharged','polar uncharged',
                               'special case','special case','special case'), 
                      
                      cols=c('#00AEEF','#00AEEF','#00AEEF','#00AEEF','#00AEEF',
                              '#7F3F98','#7F3F98','#7F3F98',
                              '#8DC63F','#8DC63F','#8DC63F','#8DC63F','#8DC63F',
                              '#F78F3D','#F78F3D','#F78F3D','#F78F3D',
                              'black','black','black'))

p1<- ggseqlogo(cdr3a_list_full, seq_type='aa',method='prob',ncol=2,
               col_scheme=color_scheme) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme(strip.text = element_text(
    size = 15, color = "black")) +
  theme(legend.text=element_text(size=15))

p2<- ggseqlogo(cdr3b_list_full, seq_type='aa',method='prob',ncol=2,
               col_scheme=color_scheme) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  theme(strip.text = element_text(
    size = 15, color = "black")) +
  theme(legend.text=element_blank())

p_all <- plot_grid(p1, p2, ncol=2, labels=c('CDR3a','CDR3b')) 


ggsave(plot=p_all, "figS1.pdf",
       width=15, height=15)
