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

############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to current

# Get TCR activation data for MHC I and II
tcr_activation_data_mhcI <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                            '/TCR_pMHCI_mutational_scan_database.xlsx'))

tcr_activation_data_mhcII <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                                   '/TCR_pMHCII_mutational_scan_database.xlsx'))

tcr_data <- rbind(tcr_activation_data_mhcI,tcr_activation_data_mhcII)

# Only use TCRs with available CDR3ab sequence
tcr_data <- tcr_data[!is.na(tcr_data$cdr3b) & !is.na(tcr_data$cdr3a),]


# Assign TCR activation category,mutation location and WT and mutant AA
tcr_data$tcr_activation <- "Weak" #Initialize
tcr_data$wt_aa <- NA
tcr_data$mutant_aa <- NA
tcr_data$mutation_location <- NA
tcr_names <- unique(tcr_data$tcr) # List all TCRs

# Loop over TCRs for assigning activation category
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

# Remove index peptides from data
tcr_data <- tcr_data[tcr_data$index_peptide!=tcr_data$peptide,]

# Loop over TCRs for recording mutation position and mutated aa
for (tcr in tcr_names) {  
  #Subset full TCR data to the current tcr
  tcr_data_subset <- tcr_data[tcr_data$tcr==tcr,]
  
  # Get mutation location
  # First flatten WT and mutated seqs, then do elementwise comparison
  mutation_location <- unlist(strsplit(paste0(tcr_data_subset$index_peptide,
                                              collapse=''),'')) !=
    unlist(strsplit(paste0(tcr_data_subset$peptide,collapse=''),''))
  
  # Record WT and mutant AA at mutation location
  wt_aa <- unlist(strsplit(paste0(tcr_data_subset$index_peptide,
                                  collapse=''),''))[mutation_location]
  
  mutant_aa <- unlist(strsplit(paste0(tcr_data_subset$peptide,
                                  collapse=''),''))[mutation_location]
  
  # Multiply by digits array to assign mutation location and reshape
  peptide_length <- nchar(unique(tcr_data_subset$index_peptide))
  mutation_location <- matrix(mutation_location*rep(c(1:peptide_length),
                                                    length(mutation_location)/peptide_length),
                              ncol = peptide_length, byrow = TRUE)
  
  # Sum to get the nonzero mutation location
  mutation_location <- rowSums(mutation_location)
  
  # Record mutation_position, WT and mutant amino acid
  tcr_data$mutation_location[tcr_data$tcr==tcr] <- mutation_location
  tcr_data$wt_aa[tcr_data$tcr==tcr] <- wt_aa
  tcr_data$mutant_aa[tcr_data$tcr==tcr] <- mutant_aa
}

# Assign WT and mutant AA hydrophobicity
# Read hydrophobicity data
aa_hydrophobicity <- read_excel(paste0('../../../tcr_epitope_datasets/amino_acid_features',
                                              '/interfacial_hydrophobicity.xlsx'))
aa_hydrophobicity <- data.frame(aa_hydrophobicity, row.names = 1)

# Record WT, Mutant hydrophobicity
tcr_data$wt_aa_hydrophobicity <- aa_hydrophobicity[tcr_data$wt_aa,1]
tcr_data$mutant_aa_hydrophobicity <- aa_hydrophobicity[tcr_data$mutant_aa,1]


# Calculate average hydrophobicity of Strong, Weak, and Non-activators and CDR3ab for each tcr
tcr_hydrophobicity <- data.frame(unique(tcr_data$tcr))
colnames(tcr_hydrophobicity)<- 'tcr'

tcr_hydrophobicity$mean_strong <- NA
tcr_hydrophobicity$mean_non <- NA
tcr_hydrophobicity$index_peptide <- NA
tcr_hydrophobicity$sum_cdr3a <- NA
tcr_hydrophobicity$sum_cdr3b <- NA
tcr_hydrophobicity$cdr3a_length <- NA
tcr_hydrophobicity$cdr3b_length <- NA
tcr_hydrophobicity$percent_strong_activator <- NA



for (tcr in tcr_hydrophobicity$tcr) {
  
  # Record Mean hydrophibicity of mutant AA residues for SB and NB
  tcr_hydrophobicity$mean_strong[tcr_hydrophobicity$tcr==tcr] <-
    mean(tcr_data$mutant_aa_hydrophobicity[tcr_data$tcr==tcr &
                                          tcr_data$tcr_activation=="Strong"])
  
  tcr_hydrophobicity$mean_non[tcr_hydrophobicity$tcr==tcr] <-
    mean(tcr_data$mutant_aa_hydrophobicity[tcr_data$tcr==tcr &
                                             tcr_data$tcr_activation=="No"])
  
  tcr_hydrophobicity$index_peptide[tcr_hydrophobicity$tcr==tcr] <- unique(tcr_data$index_peptide[tcr_data$tcr==tcr])
  
  # Calculate total CDR3ab hydrophobicity
  cdr3a <- unique(tcr_data$cdr3a[tcr_data$tcr==tcr])
  
  tcr_hydrophobicity$sum_cdr3a[tcr_hydrophobicity$tcr==tcr] <-
    sum(aa_hydrophobicity[unlist(strsplit(cdr3a,split="")),1])
  
  cdr3b <- unique(tcr_data$cdr3b[tcr_data$tcr==tcr])
  
  tcr_hydrophobicity$sum_cdr3b[tcr_hydrophobicity$tcr==tcr] <-
    sum(aa_hydrophobicity[unlist(strsplit(cdr3b,split="")),1])
  
  # Get CDR3ab lengths
  tcr_hydrophobicity$cdr3a_length[tcr_hydrophobicity$tcr==tcr] <-  nchar(cdr3a)
  tcr_hydrophobicity$cdr3b_length[tcr_hydrophobicity$tcr==tcr] <-  nchar(cdr3b)
  
  # Get % strong activator
  tcr_hydrophobicity$percent_strong_activator[tcr_hydrophobicity$tcr==tcr] <-
    100*nrow(tcr_data[tcr_data$tcr==tcr & 
                        tcr_data$tcr_activation=="Strong",])/nrow(tcr_data[tcr_data$tcr==tcr,])
  
  }

# Hydrophobicity difference of non and strong binders
tcr_hydrophobicity$mean_diff <- tcr_hydrophobicity$mean_non - 
                                          tcr_hydrophobicity$mean_strong 
# Peptide type
tcr_hydrophobicity$mhc_type <- "MHCI"
tcr_hydrophobicity$mhc_type[tcr_hydrophobicity$index_peptide %in% unique(
  tcr_activation_data_mhcII$index_peptide)] <- "MHCII"

############
## Plot ####
############

p <- ggplot(tcr_hydrophobicity,aes(y=mean_diff,
                   x=sum_cdr3b,color=mhc_type)) + 
  geom_point(size=3,alpha=0.9) +
  geom_smooth(data = tcr_hydrophobicity, aes(y=mean_diff, x=sum_cdr3b),
              colour="gray",
              method = "lm") +
  scale_color_manual(values=c('#5D3FD3','#66a61e')) +
  labs(x="CDR3b hydrophobicity",y="Mutant AA Hydrophobicity difference (non-strong)",color="")+
  theme_cowplot() +
  geom_hline(yintercept=0, linetype='dotted', col = 'black') +
  geom_vline(xintercept=0, linetype='dotted', col = 'black') +
  theme(axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20)) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5))

ggsave(plot=p, "tcr_hydrophobicity_allowed_subs.pdf",
       width=5, height=5)

write.csv(tcr_hydrophobicity,'raw_data_fig_1g.csv')


