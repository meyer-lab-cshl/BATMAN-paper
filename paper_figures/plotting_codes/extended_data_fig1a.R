# Plots performance of difference pooling methods over TCRs grouped by index

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
# Load data
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# Discard patented TCRs
TCR_data_full <- TCR_data[!(TCR_data$tcr_name %in% c('T1','T3','FLT3DY')),]

# List index_peptides
index_peptides <- unique(TCR_data_full$index_peptide)

# List for CDR sequences
cdr3a_list_full <- list()
cdr3b_list_full <- list()

# Filter based on index peptide
for (index_peptide_selected in index_peptides){
TCR_data <- TCR_data_full[TCR_data_full$index_peptide==index_peptide_selected,]

# Extract TCR names present in the data
TCR_names <- unique(TCR_data$tcr_name)

# Extract and align TCR CDR3a/b seqs for selected TCRs
cdr3a_list <- rep('X',length(TCR_names))
cdr3b_list <- rep('X',length(TCR_names))

for (tcr in 1:length(TCR_names)){
  cdr3a_list[tcr] <- unique(TCR_data$cdr3a[TCR_data$tcr_name==TCR_names[tcr]])
  cdr3b_list[tcr] <- unique(TCR_data$cdr3b[TCR_data$tcr_name==TCR_names[tcr]])
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
p1<- ggseqlogo(cdr3a_list_full, seq_type='aa',method='prob',ncol=1) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  theme(strip.text = element_text(
          size = 20, color = "black"))

p2<- ggseqlogo(cdr3b_list_full, seq_type='aa',method='prob',ncol=1) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  theme(strip.text = element_text(
    size = 20, color = "black"))

p_all <- plot_grid(p1, p2, ncol=2) 


ggsave(plot=p_all, "../figures/extended_data_fig1/seq_logo.pdf",
       width=10, height=20)

