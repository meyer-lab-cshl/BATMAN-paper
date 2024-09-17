# Code to plot scores from different TCR-pMHC methods against actual activation

rm(list=ls())#Clear environment
# Load required libraries
library(readxl)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(forcats)
library(viridis)


#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load data
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# List of 3 class 9-mer-binding TCRs
tcr_names <- c("18A2","NYE-S1","868Z11","TCR3","TCR6","TCR7","A6",
               "T1","FLT3DY","A23","TCR2-T")

# Load peptide data for selected TCRs
TCR_data <- select(TCR_data[TCR_data$tcr_name %in% tcr_names,],
                  c('tcr_name', 'peptide', 'peptide_activity'))

# Normalize peptide activity
for (TCR in tcr_names){
  TCR_data$peptide_activity[TCR_data$tcr_name==TCR] <-
    TCR_data$peptide_activity[TCR_data$tcr_name==TCR]/
    max(TCR_data$peptide_activity[TCR_data$tcr_name==TCR])
}

# peptide scores from BATMAN+7 methods
peptide_score <- array(NA,dim = c(dim(TCR_data)[1],8))
peptide_score[,1] <- TCR_data$peptide_activity
  
# Load TCR-pMHC scores from 7 other methods
# We added a negative sign before ranks to use them as score
#



#ERGO II
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "ergo//ergo_output_mcpas.csv"))
  
peptide_score[,2] <- as.double(output_data$Score[match(TCR_data$peptide,
                                                    output_data$Peptide)])


#epiTCR
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "epitcr//output.csv"))

peptide_score[,3] <- as.double(output_data$predict_proba[match(TCR_data$peptide,
                                                       output_data$epitope)])

#ImRex
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "imrex//output-prediction.csv"))

peptide_score[,4] <- as.double(output_data$prediction_score[match(TCR_data$peptide,
                                                  output_data$antigen.epitope)])

#pMTnet
output_data <- read.csv(paste0("../codes/misc_methods/",
                               "pmtnet//pMTnet_prediction.csv"))

peptide_score[,5] <- -as.double(output_data$Rank[match(TCR_data$peptide,
                                                          output_data$Antigen)])

#IEDB
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/iedb/",pattern = ".csv")
file_names <- paste0(rep("../codes/misc_methods/iedb/",
                         length(file_names)),file_names)
output_data <- read_csv(file_names)

peptide_score[,6] <- as.double(output_data$score[match(TCR_data$peptide,
                                                      output_data$peptide)])

#NetTepi
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/nettepi/",pattern = "NetTepi_out_")
file_names <- paste0(rep("../codes/misc_methods/nettepi/",
                         length(file_names)),file_names)
output_data <- data.frame()
for (file in file_names){
  output_data <- rbind(output_data,read.table(file))
}

peptide_score[,7] <- as.double(output_data$V8[match(TCR_data$peptide,
                                                       output_data$V2)])

#PRIME 2.0
#Load outputs for all MHCs
file_names <- list.files("../codes/misc_methods/prime/",pattern = "result_")
file_names <- paste0(rep("../codes/misc_methods/prime/",
                         length(file_names)),file_names)
output_data <- data.frame()
for (file in file_names){
  output_data <- rbind(output_data,read.table(file))
}

peptide_score[,8] <- as.double(output_data$V3[match(TCR_data$peptide,
                                                    output_data$V1)])
# Create and shape the full dataFrame for plotting
plot_data <- data.frame(matrix(nrow=0,ncol=5))
colnames(plot_data) <- c('TCR','peptide','activation','method','score')

method_name <- c('ERGO II','epiTCR','ImRex','pMTnet','IEDB','NetTepi',
                 'PRIME 2.0')

corrs = array(NA,dim=7) # Array to store correlation coefficients

for (method_index in 1:7){
  TCR <- TCR_data$tcr_name
  peptide <- TCR_data$peptide
  activation <- peptide_score[,1]
  method <- rep(method_name[method_index],length(TCR_data$tcr_name))
  score <- peptide_score[,method_index+1]
  new_data <- data.frame(TCR, peptide, activation, method, score)
  colnames(new_data) <- c('TCR','peptide','activation','method','score')
  
  plot_data <- rbind(plot_data,new_data)
  corrs[method_index] <- cor(peptide_score[,method_index+1],peptide_score[,1],
                             use="complete.obs")^2
}

###################
###### plotting ###
###################

p <- ggplot(plot_data, aes(y=score, x=activation))
p <- p +
  geom_point(aes(color=TCR), alpha=0.5) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
  labs(y="TCR-pMHC score",
       x="Normalized TCR activation",
       color="") +
  xlim(NA, 0.9) +
  theme_cowplot() +
  facet_wrap(vars(method),scales = "free")+
  theme(strip.text = element_text(
    size = 25, color = "black")) +
  theme(legend.text = element_text(size=20)) +
  guides(color = guide_legend(ncol = 6)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30)) +
  scale_x_continuous(limits=c(0, 1), breaks = c(0,0.25,0.5,0.75,1))


ggsave(plot=p, "../figures/extended_data_fig6/score_by_method.pdf",
       width=20, height=15)


  


