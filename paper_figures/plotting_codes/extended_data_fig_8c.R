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
library(viridis)
library(loo)
library(pROC)

############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# List of 3 class 9-mer-binding TCRs
tcr_names <- c("18A2","868Z11","FLT3DY","NYE-S1","A23","A6","T1",
               "TCR2-T","TCR3","TCR6","TCR7")

plot_data <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(plot_data) <- c('npeptide','tcr','auc')

for (tcr_name in tcr_names){

# Load TCR-specific OPT AUC data
auc_opt <- read_csv(paste0("../data/AL_OPT/",tcr_name,"_OPT.csv"))
colnames(auc_opt)[1] <- 'npeptide'
auc_opt$tcr <- tcr_name
auc_opt <-  na.omit(auc_opt) #remove NA
auc_opt <- auc_opt[4:(dim(auc_opt)[1]-1),] #first 2 and Last points are reducdant
auc_opt <- auc_opt[,c(1,2,4)] #remove peptide seq

plot_data <- rbind(plot_data,auc_opt)
}

################
### Plotting  ##
################

p <- ggplot(data=plot_data, aes(x=npeptide, y=auc, group=tcr)) +
  geom_line(aes(color=tcr)) +
  geom_point(aes(color=tcr,size=1,alpha=0.9)) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
                                         labs(y="AUC",
                                              x="Number of peptides sampled",
                                              color="") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=25),
        axis.title.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        axis.text.y = element_text(size=25))

ggsave(plot=p, "../figures/extended_data_fig8/auc_opt.pdf",width=15, height=5)







