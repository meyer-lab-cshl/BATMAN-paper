# Plots performance of LOO TCR tasks for BATMAN and Atchley

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
library(pROC)
library(ggnewscale)



#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load functions
source("../codes/all_functions.R") 

# Arrays to store 3 types of AUCs for 66 TCRs
#(Atchley+RF, weights_only, full and symmetric AA matrix,
#all for both within index peptide and across n-mers)
# (initialize with NaN so that missing data is not plotted)
all_AUC_non_weak <- array(NaN,dim=c(66,8))
all_AUC_non_strong <- array(NaN,dim=c(66,8))
all_AUC_weak_strong <- array(NaN,dim=c(66,8))

# Indices of 2-class TCRs
TCR_index_2_class <- c(26,29,30,31)
class_pair_2_class <- c(23,23,23,23)#1:non, 2: weak, 3: strong

# Load TCR names
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)
TCR <- unique(TCR_data$tcr_name)

# Loop over TCR index
for (TCR_index in 1:66){
  print(TCR_index)
  
  # Index peptide and its length
  index_peptide <- 
               unique(TCR_data$index_peptide[TCR_data$tcr_name==TCR[TCR_index]])
  
  index_peptide_length <- nchar(index_peptide)
    
  ## Atchley+RF results for across n-mer #######################################
  if (file.exists(paste0("../data/loo_tcr/atchley/auc_loo_tcr_atchley_all_",
                        index_peptide_length,"_mer.rda"))){
  load(paste0("../data/loo_tcr/atchley/auc_loo_tcr_atchley_all_",
              index_peptide_length,"_mer.rda"))
  if (TCR[TCR_index] %in% TCR_names){
  if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
    all_AUC_weak_strong[TCR_index,1] <- 
                                     loo_tcr[1,which(TCR_names==TCR[TCR_index])]
  } else {
    all_AUC_non_weak[TCR_index,1] <- loo_tcr[1,which(TCR_names==TCR[TCR_index])]
    all_AUC_non_strong[TCR_index,1]<-loo_tcr[2,which(TCR_names==TCR[TCR_index])]
    all_AUC_weak_strong[TCR_index,1]<-loo_tcr[3,which(TCR_names==TCR[TCR_index])]
  }
  }
  }
  
  ## Atchley+RF results for within index #######################################
  if (file.exists(paste0("../data/loo_tcr/atchley/auc_loo_tcr_atchley_",
                         index_peptide,".rda"))){
    load(paste0("../data/loo_tcr/atchley/auc_loo_tcr_atchley_",
                index_peptide,".rda"))
    if (TCR[TCR_index] %in% TCR_names){
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,2] <- 
        loo_tcr[1,which(TCR_names==TCR[TCR_index])]
    } else {
      all_AUC_non_weak[TCR_index,2] <- loo_tcr[1,which(TCR_names==TCR[TCR_index])]
      all_AUC_non_strong[TCR_index,2]<-loo_tcr[2,which(TCR_names==TCR[TCR_index])]
      all_AUC_weak_strong[TCR_index,2]<-loo_tcr[3,which(TCR_names==TCR[TCR_index])]
    }
    }
  }
  ##### BATMAN results #########################################################
  # Load BLOSUM100 matrix
  AA_list <-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',
              'T','V','W','Y')
  load("../data/distance_functions/BLOSUM100.rda")# Load BLOSUM100
  AA_distance <- matrix(aa_mat,nrow=20,ncol=20,dimnames=list(AA_list,AA_list))
  
  # weights-only, across index ##########################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/weights_only/weights_all_",
                         index_peptide_length,"_mer_",TCR[TCR_index],".csv"))){
    suppressMessages(
    weights<-unlist(read_csv(paste0("../data/loo_tcr/weights_only/weights_all_",
                index_peptide_length,"_mer_",TCR[TCR_index],".csv"),
                col_types = cols())[2:(1+index_peptide_length)])
    )
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance)
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,3] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,3]<-AUC[1]
      all_AUC_non_strong[TCR_index,3]<-AUC[2]
      all_AUC_weak_strong[TCR_index,3]<-AUC[3]
    }
  }
  # weights-only, within index ##########################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/weights_only/weights_",
                         index_peptide,"_",TCR[TCR_index],".csv"))){
    suppressMessages(
      weights<-unlist(read_csv(paste0("../data/loo_tcr/weights_only/weights_",
                                      index_peptide,"_",TCR[TCR_index],".csv"),
                               col_types = cols())[2:(1+index_peptide_length)])
    )
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance)
    
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,4] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,4]<-AUC[1]
      all_AUC_non_strong[TCR_index,4]<-AUC[2]
      all_AUC_weak_strong[TCR_index,4]<-AUC[3]
    }
  }
  
  # weights+symm AA matrix, across index #######################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/symm_aa_matrix/weights_all_",
                         index_peptide_length,"_mer_",TCR[TCR_index],".csv"))){
    suppressMessages(
      weights<-unlist(read_csv(paste0("../data/loo_tcr/symm_aa_matrix/weights_all_",
                            index_peptide_length,"_mer_",TCR[TCR_index],".csv"),
                               col_types = cols())[2:(1+index_peptide_length)])
    )
      
    AA_multiplier<-read.csv(paste0("../data/loo_tcr/symm_aa_matrix/aa_matrix_all_",
                            index_peptide_length,"_mer_",TCR[TCR_index],".csv"))
    AA_multiplier <- data.matrix(AA_multiplier[,2:21])
    
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance*AA_multiplier)
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,5] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,5]<-AUC[1]
      all_AUC_non_strong[TCR_index,5]<-AUC[2]
      all_AUC_weak_strong[TCR_index,5]<-AUC[3]
    }
  }
  # weights+symm AA matrix, within index ##########################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/symm_aa_matrix/weights_",
                         index_peptide,"_",TCR[TCR_index],".csv"))){
    suppressMessages(
      weights<-unlist(read_csv(paste0("../data/loo_tcr/symm_aa_matrix/weights_",
                                      index_peptide,"_",TCR[TCR_index],".csv"),
                               col_types = cols())[2:(1+index_peptide_length)])
    )
    AA_multiplier<-read.csv(paste0("../data/loo_tcr/symm_aa_matrix/aa_matrix_",
                                    index_peptide,"_",TCR[TCR_index],".csv"))
    AA_multiplier <- data.matrix(AA_multiplier[,2:21])
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance*AA_multiplier)
    
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,6] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,6]<-AUC[1]
      all_AUC_non_strong[TCR_index,6]<-AUC[2]
      all_AUC_weak_strong[TCR_index,6]<-AUC[3]
    }
  }
  
  # weights+full AA matrix, across index #######################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/full_aa_matrix/weights_all_",
                         index_peptide_length,"_mer_",TCR[TCR_index],".csv"))){
    suppressMessages(
      weights<-unlist(read_csv(paste0("../data/loo_tcr/full_aa_matrix/weights_all_",
                                      index_peptide_length,"_mer_",TCR[TCR_index],".csv"),
                               col_types = cols())[2:(1+index_peptide_length)])
    )
    
    AA_multiplier<-read.csv(paste0("../data/loo_tcr/full_aa_matrix/aa_matrix_all_",
                                   index_peptide_length,"_mer_",TCR[TCR_index],".csv"))
    AA_multiplier <- data.matrix(AA_multiplier[,2:21])
    
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance*AA_multiplier)
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,7] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,7]<-AUC[1]
      all_AUC_non_strong[TCR_index,7]<-AUC[2]
      all_AUC_weak_strong[TCR_index,7]<-AUC[3]
    }
  }
  # weights+full AA matrix, within index ##########################################
  #load weights
  if (file.exists(paste0("../data/loo_tcr/full_aa_matrix/weights_",
                         index_peptide,"_",TCR[TCR_index],".csv"))){
    suppressMessages(
      weights<-unlist(read_csv(paste0("../data/loo_tcr/full_aa_matrix/weights_",
                                      index_peptide,"_",TCR[TCR_index],".csv"),
                               col_types = cols())[2:(1+index_peptide_length)])
    )
    AA_multiplier<-read.csv(paste0("../data/loo_tcr/full_aa_matrix/aa_matrix_",
                                   index_peptide,"_",TCR[TCR_index],".csv"))
    AA_multiplier <- data.matrix(AA_multiplier[,2:21])
    # Get AUC
    AUC <- auc_loo_tcr(TCR[TCR_index],weights,AA_distance*AA_multiplier)
    
    # Save  
    if (TCR_index %in% TCR_index_2_class){#2-class TCR, all in weak-strong
      all_AUC_weak_strong[TCR_index,8] <-AUC 
      
    } else {
      all_AUC_non_weak[TCR_index,8]<-AUC[1]
      all_AUC_non_strong[TCR_index,8]<-AUC[2]
      all_AUC_weak_strong[TCR_index,8]<-AUC[3]
    }
  }
  
  
  # for 8-mer-binder tcrs (all ova-specific), within-index pooling=across-8-mer pooling
  if (index_peptide_length==8){
    all_AUC_non_weak[TCR_index,c(2,4,6,8)] <- all_AUC_non_weak[TCR_index,c(1,3,5,7)]
    all_AUC_non_strong[TCR_index,c(2,4,6,8)] <- all_AUC_non_strong[TCR_index,c(1,3,5,7)]
    all_AUC_weak_strong[TCR_index,c(2,4,6,8)] <- all_AUC_weak_strong[TCR_index,c(1,3,5,7)]
  }  
  
  
  
  
  
} #TCR loop

# order of methods: 

methods <- c('pTEAM_across','pTEAM_within','BLOSUM_across','BLOSUM_within', 
             'Symmetric_across','Symmetric_within',
             'Full_across','Full_within')

# Add TCR and AA matrix names to array
colnames(all_AUC_non_weak) <- methods
colnames(all_AUC_non_strong) <- methods
colnames(all_AUC_weak_strong) <- methods
all_AUC_non_weak<- data.frame(TCR,all_AUC_non_weak)
all_AUC_non_strong<- data.frame(TCR,all_AUC_non_strong)
all_AUC_weak_strong<- data.frame(TCR,all_AUC_weak_strong)

# Convert to long form for box plot
all_AUC_non_weak <- all_AUC_non_weak %>%
  pivot_longer(-TCR, names_to="method", values_to = "score")

all_AUC_non_strong <- all_AUC_non_strong %>%
  pivot_longer(-TCR, names_to="method", values_to = "score") 

all_AUC_weak_strong <- all_AUC_weak_strong %>%
  pivot_longer(-TCR, names_to="method", values_to = "score") 

# Merge data and indicate class pair
all_AUC <- rbind.data.frame(all_AUC_non_weak,all_AUC_non_strong,all_AUC_weak_strong)
all_AUC$class_pair <- c(rep("12", each = dim(all_AUC_non_weak)[1]),
                        rep("13", each = dim(all_AUC_non_strong)[1]),
                        rep("23", each = dim(all_AUC_weak_strong)[1]))


# Indicate index peptides for each tcr in data point
# Array to store TCR names and index peptides
index_peptide <- array("X",dim=c(dim(all_AUC)[1]))
for (i in 1:dim(all_AUC)[1]){
  TCR_index <- match(all_AUC$TCR[i],unique(all_AUC$TCR)) 
  #this is because R changed TCR names
  index_peptide[i] <- 
               unique(TCR_data$index_peptide[TCR_data$tcr_name==TCR[TCR_index]])
}

all_AUC$index_peptide <- index_peptide

# Discard data for one 10-mer-binding and one 11-mer-binding TCR
all_AUC <- all_AUC[!(all_AUC$index_peptide %in% c('ALWGPDPAAA','ALYDKTKRIFL')),]

all_AUC %>%
  mutate(method=fct_inorder(method))

level_order <- unique(all_AUC$method)
# For SIINFEKL TCRs, classify them into 2 categories
tcr_type <- array(NaN,dim=c(dim(all_AUC)[1]))

naive_tcrs <- c("B11","B15","B3", "F4","E8", "B13","H6", "G6", "F5", 
                "H5", "B2", "B6", "B5","E9", "E4", "G2", "B16","B14")
ed_tcrs <- c("OT1","Ed5","Ed8","Ed9","Ed10","Ed16-1","Ed16-30","Ed21","Ed23",
             "Ed28","Ed31","Ed33","Ed39","Ed40","Ed45","Ed46")

tcr_type[all_AUC$TCR %in% naive_tcrs] <- "Naive"
tcr_type[all_AUC$TCR %in% ed_tcrs] <- "Educated"

all_AUC$tcr_type <- tcr_type

###############
##  plot ######
###############

# Special plotting specs for SIINFEKL data
all_AUC <- all_AUC %>%
  mutate(color_by = TCR,
         color_by = case_when(index_peptide == "SIINFEKL" ~ tcr_type,
                              TRUE ~ color_by))

auc_by_index <- split(all_AUC, f=all_AUC$index_peptide)

boxplot_by_index <- lapply(auc_by_index, function(x) {
  
  set_alpha <- 0.7 
  set_size <- 2.25
  set_jitter_width <- 0.15
  
  if (unique(x$index_peptide)=="SIINFEKL"){
    set_alpha <- 0.85
    set_size <- 1.6
    set_jitter_width <- 0.2
  }
  
  p <- ggplot(x, 
              aes(x=factor(method, levels = unique(all_AUC$method)), y=score)) +
    geom_boxplot(aes(color=factor(class_pair,levels = unique(class_pair))),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1,
                 position=position_dodge(0.75)) +
    scale_color_manual(values=c('#50C878','#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(color_by), 
                   group=factor(class_pair,levels = unique(class_pair))), 
               alpha=set_alpha,
               size=set_size,
               position = position_jitterdodge(jitter.width=set_jitter_width,
                                               dodge.width = 0.75)) +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(y="AUC",
         x="",
         color="") +
    #    ylim(min(all_AUC$score,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=20))+
    facet_wrap(vars(index_peptide))+
    theme(strip.text = element_text(
      size = 25, color = "black")) +
    theme(legend.text = element_text(size=15)) +
    guides(color = guide_legend(ncol = 3)) +
    theme(legend.position="bottom")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=4,nrow=3,
                            align = "hv",axis="tb")

ggsave(plot=p_all, "../figures/extended_data_fig4/loo_tcr_by_index.pdf",
       width=20, height=20)

