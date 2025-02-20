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

# Record WT, Mutant, and change in hydrophobicity
tcr_data$wt_aa_hydrophobicity <- aa_hydrophobicity[tcr_data$wt_aa,1]
tcr_data$mutant_aa_hydrophobicity <- aa_hydrophobicity[tcr_data$mutant_aa,1]

tcr_data$hydrophobicity_change <- tcr_data$mutant_aa_hydrophobicity -
                                                  tcr_data$wt_aa_hydrophobicity

# Calculate average and sd hydrophobicity of Strong, Weak, and Strong-activators for each WT AA
mutant_hydrophobicity <- data.frame(rownames(aa_hydrophobicity))
colnames(mutant_hydrophobicity)<- 'wt_aa'

mutant_hydrophobicity$mean_strong <- NA
mutant_hydrophobicity$mean_non <- NA

mutant_hydrophobicity$sd_strong <- NA
mutant_hydrophobicity$sd_non <- NA

mutant_hydrophobicity$percent_strong_activator <- NA



for (aa in rownames(aa_hydrophobicity)) {
  mutant_hydrophobicity$mean_strong[mutant_hydrophobicity$wt_aa==aa] <-
    mean(tcr_data$mutant_aa_hydrophobicity[tcr_data$wt_aa==aa &
                                          tcr_data$tcr_activation=="Strong"])
  
  mutant_hydrophobicity$sd_strong[mutant_hydrophobicity$wt_aa==aa] <-
    sd(tcr_data$mutant_aa_hydrophobicity[tcr_data$wt_aa==aa &
                                             tcr_data$tcr_activation=="Strong"])
  
  mutant_hydrophobicity$mean_non[mutant_hydrophobicity$wt_aa==aa] <-
    mean(tcr_data$mutant_aa_hydrophobicity[tcr_data$wt_aa==aa &
                                             tcr_data$tcr_activation=="No"])
  
  mutant_hydrophobicity$sd_non[mutant_hydrophobicity$wt_aa==aa] <-
    sd(tcr_data$mutant_aa_hydrophobicity[tcr_data$wt_aa==aa &
                                             tcr_data$tcr_activation=="No"])
  
  mutant_hydrophobicity$percent_strong_activator[mutant_hydrophobicity$wt_aa==aa] <-
    100*nrow(tcr_data[tcr_data$wt_aa==aa & 
                          tcr_data$tcr_activation=="Strong",])/nrow(tcr_data[tcr_data$wt_aa==aa,])
}

mutant_hydrophobicity$mean_diff <- mutant_hydrophobicity$mean_non - 
                                          mutant_hydrophobicity$mean_strong
mutant_hydrophobicity$sd_diff <- sqrt((mutant_hydrophobicity$sd_non)^2 +
                                        (mutant_hydrophobicity$sd_strong)^2) 

############
## Plot ####
############
# Value used to transform the data
coeff <- 75

p <- ggplot(mutant_hydrophobicity,aes(x=factor(wt_aa,
                                    levels = rownames(aa_hydrophobicity)))) + 
  geom_point(aes(y=mean_diff),size=5,color='#D22B2B',show.legend=FALSE) +
  # geom_point(aes(y=percent_strong_activator/coeff),size=2,color='#4682B4',
  #            show.legend=FALSE) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "hp(AAsB-AAnon)[kcal/mol]",
    
    # Add a second axis and specify its features
    #sec.axis = sec_axis(~.*coeff, name="% Strong-activating mutants")
  ) +
  theme(
    axis.title.y = element_text(color = '#B47846', size=20)
    #axis.title.y.right = element_text(color = '#4682B4', size=20)
  ) +
  
  labs(x="WT AA",color="") +
  theme_cowplot() +
  geom_hline(yintercept=0, linetype='dashed', col = '#966919') +
  geom_vline(xintercept=9.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

ggsave(plot=p, "peptide_hydrophobicity_allowed_subs.pdf",
       width=5, height=5)

#Save raw data
write.csv(mutant_hydrophobicity,'raw_data_fig_1f.csv')


