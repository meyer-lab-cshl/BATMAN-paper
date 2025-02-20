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
library(scales)

############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to current

# Get TCR activation data for MHC I and II
tcr_activation_data_mhcI <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                            '/TCR_pMHCI_mutational_scan_database.xlsx'))

tcr_activation_data_mhcII <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                                   '/TCR_pMHCII_mutational_scan_database.xlsx'))

tcr_activation_data <- rbind(tcr_activation_data_mhcI,tcr_activation_data_mhcII)

# Separately get c259 TCR activation data for the index peptide for which Kd is available
tcr_activation_data <- tcr_activation_data[tcr_activation_data$tcr!='c259',]
tcr_activation_data_c259 <- read_excel('tcr_pmhc_binding_datasets/c259/activation_data_c259.xlsx')
tcr_activation_data <- rbind(tcr_activation_data,tcr_activation_data_c259)

# Read TCR binding data
tcr_binding_data <- read_excel('tcr_pmhc_binding_datasets/tcr_pmhc_kd.xlsx')

# Subset to common TCRs for which both activation and Kd are available
tcr_activation_data <- tcr_activation_data[
                     tcr_activation_data$tcr %in% unique(tcr_binding_data$tcr),]

# Normalize TCR activation
tcr_names <- unique(tcr_activation_data$tcr) # List all TCRs

# Loop over TCRs
for (tcr in tcr_names) {
  
  # Extract peptides activities for a particular TCR
  peptide_activity <- tcr_activation_data$peptide_activity[
                                                    tcr_activation_data$tcr==tcr]
  
  # Normalize peptide activities to the maximum activity
  peptide_activity <- peptide_activity/max(peptide_activity)
  
  # Add to TCR data
  tcr_activation_data$peptide_activity[tcr_activation_data$tcr==tcr] <- peptide_activity
}

# Merge binding and activation data
tcr_data <- merge(tcr_activation_data,tcr_binding_data,
                  by=c("tcr","index_peptide","peptide"))

# Indicate index peptide
tcr_data$is_index <- 0
tcr_data$is_index[tcr_data$index_peptide==tcr_data$peptide] <- 1

############
## Plot ####
############

p <- ggplot(tcr_data,aes(y=peptide_activity,
                   x=Kd)) + 
  geom_point(aes(color=factor(tcr,levels=c('868Z11','c259','TIL1383I','B3K506',
                                           'B3K508','TCR2W1S12-20.4','TCR75-1',
                                           'YAe5-62.8')),shape=as.factor(is_index)),
             size=3*tcr_data$is_index+3,alpha=0.8)+
  
   scale_color_manual(values=c('#b30000','#808080','#FAA0A9','#0096FF','#5D3FD3',
                                        '#0047AB','#bc80bd','#50C878')) +
                               
                                         labs(x="TCR-pMHC binding affinity (Kd in M)",
                                              y="Normalized TCR activation",
                                              color="")+
  theme_cowplot() +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
  theme(axis.text.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        axis.text.y = element_text(size=30)) +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(size=25))

ggsave(plot=p, "tcr_activation_binding.pdf",
       width=7, height=7)

write.csv(tcr_data,'raw_data_fig_1e.csv')


