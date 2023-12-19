# Plots performance of difference AA distance functions over all TCRs

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
library(ggrepel)

## Extended Data Figure 2a: average AUC (over 3 or 1 class-pairs and 5 folds) for TCRs
#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is
all_AUC <- array(0,dim=c(71,66))

# results for Bayesian inference
for (AA_distance in 1:70) {
  print(AA_distance)
  for (tcr in 1:66){
    load(paste("../data/within_tcr_unpooled_outputs//output_",
                                            AA_distance,"_",tcr,".rda",sep=""))
    all_AUC[AA_distance,tcr]<-mean(classification_results$AUCs)
  }
}

# results for Atchley
for (tcr in 1:66){
  load(paste("../data/within_tcr_unpooled_outputs_atchley//output_",
                                                  tcr,".rda",sep=""))
  all_AUC[71,tcr]<-mean(AUCs)
}

# Names of TCRs and AA distance matrices
# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)
# Extract TCR names present in the data
TCRs <- unique(TCR_data$tcr_name)

AA_distances <-
c("Hamming",paste(rep("BLOSUM",13),seq(30,90,5),sep=""),
                        "BLOSUM62","BLOSUM100",
                        paste(rep("PAM",50),seq(10,500,10),sep=""),
                        "Dayhoff","Gonnet","Atchley_l2","Atchley_cos","Atchley+RF")

# Add TCR and AA matrix names to array
colnames(all_AUC) <- TCRs

all_AUC <- data.frame(AA_distances,all_AUC)
###############
## Analysis ###
###############
# Convert to long form for box plot
AUC <- all_AUC %>%
  pivot_longer(-AA_distances, names_to="tcr", values_to = "score") %>%
  mutate(AA_distances=fct_inorder(AA_distances))                 


###############
##  plot ######
###############
p <- ggplot(AUC, aes(x=AA_distances, y=score))
p <- p + geom_boxplot(outlier.colour = NA,show.legend = FALSE) +
  geom_jitter(aes(color=as.factor(tcr)), alpha=0.5, width=0.1, height=0, show.legend = TRUE) +
  scale_colour_manual(values = palette(viridis(66))) +
  labs(y="Average AUC",
       x="",
       color="") +
#  xlim(NA, 0.9) +
  theme_cowplot() +
  theme_minimal_grid(14) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=15),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  theme(legend.position = "bottom")+ 
  guides(color = guide_legend(nrow = 9)) 

ggsave(plot=p, "../figures/extended_data_fig2/all_methods_with_legend.pdf",
       width=20, height=6.5)
