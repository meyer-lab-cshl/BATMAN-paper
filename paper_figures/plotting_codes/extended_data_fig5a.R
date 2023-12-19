# Positional weights of TCRs
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

#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is
# Load data
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)

# List of TCRs
TCR <- unique(TCR_data$tcr_name)

#Selected TCRs to plot
TCR_names <- c('18A2','868Z11','A23','A6','FLT3DY','NYE-S1','T1',
               'TCR2-T','TCR3','TCR6','TCR7') 
TCR_index <- match(TCR_names,TCR)

# 11 TCRs of Fig 2a, 9 positional weights for each of the 3 types (3 bars/position): 
#unpooled classification with hamming,blosum100,pam10,... 
#unpooled regression with hamming,blosum100,pam10,... 
#pooled classification with blosum100, full and symmetric matrix (pooled within 9-mer)

names <- c('position','weight','is_pooled','aa_method','TCR')
weights_all <- vector("list",length=length(names))
names(weights_all) <- names

#Unpooled results
methods_index <- c(1,16,17,67,68,69) #Hamming,blosum100,pam10,Dayhoff,Gonnet, Atchley
methods_name <- c('Hamming','BLOSUM100','PAM10','Dayhoff','Gonnet','Atchley')

peptide_length <- 9
for (method in 1:length(methods_index)){
  for (tcr in 1:length(TCR_index)){
load(paste0('../data/within_tcr_unpooled_outputs/output_',methods_index[method],'_',
            TCR_index[tcr],'.rda'))
# classification weights
    weights <- rowMeans(classification_results$position_dependent_weights)
    weights <- weights/max(weights) #normalize
    
    weights_all[['position']] <- append(weights_all[['position']],1:peptide_length)
    weights_all[['weight']] <- append(weights_all[['weight']],weights)
    weights_all[['is_pooled']] <- append(weights_all[['is_pooled']],
                                         rep('unpooled_classification',peptide_length))
    weights_all[['aa_method']] <- append(weights_all[['aa_method']],
                                         rep(methods_name[method],peptide_length))
    weights_all[['TCR']] <- append(weights_all[['TCR']],
                                         rep(TCR_names[tcr],peptide_length))
    
# regression weights
    weights <- rowMeans(regression_results$position_dependent_weights)[2:(peptide_length+1)]
    weights <- weights/max(weights) #normalize
    
    weights_all[['position']] <- append(weights_all[['position']],1:peptide_length)
    weights_all[['weight']] <- append(weights_all[['weight']],weights)
    weights_all[['is_pooled']] <- append(weights_all[['is_pooled']],
                                         rep('unpooled_regression',peptide_length))
    weights_all[['aa_method']] <- append(weights_all[['aa_method']],
                                         rep(methods_name[method],peptide_length))
    weights_all[['TCR']] <- append(weights_all[['TCR']],
                                   rep(TCR_names[tcr],peptide_length))
    
  }
}

# Pooled results
for (tcr in TCR_names){
  index_peptide <- unique(TCR_data$index_peptide[TCR_data$tcr_name==tcr])
  
  # Weights only
  weight_folds <- array(0,dim = c(peptide_length,5))
  for (fold in 1:5){
    weights_all_tcr <-read.csv(paste0('../data/within_tcr_pooled/weights_only/',
                                  'weights_all_',peptide_length,
                                  '_mer_fold_',fold-1,'.csv'))
    weight_folds[,fold] <- array(unlist(weights_all_tcr[weights_all_tcr[,1]==tcr,
                                        2:(1+peptide_length)]))
  }
  weights <- rowMeans(weight_folds)
  weights <- weights/max(weights) #normalize
  
  weights_all[['position']] <- append(weights_all[['position']],1:peptide_length)
  weights_all[['weight']] <- append(weights_all[['weight']],weights)
  weights_all[['is_pooled']] <- append(weights_all[['is_pooled']],
                                       rep('pooled_classification',peptide_length))
  weights_all[['aa_method']] <- append(weights_all[['aa_method']],
                                       rep('BLOSUM100',peptide_length))
  weights_all[['TCR']] <- append(weights_all[['TCR']],
                                 rep(tcr,peptide_length))
  
  # Weights+ Symm AA matrix
  weight_folds <- array(0,dim = c(peptide_length,5))
  for (fold in 1:5){
    weights_all_tcr <-read.csv(paste0('../data/within_tcr_pooled/symm_aa_matrix/',
                                      'weights_all_',peptide_length,
                                      '_mer_fold_',fold-1,'.csv'))
    weight_folds[,fold] <- array(unlist(weights_all_tcr[weights_all_tcr[,1]==tcr,
                                                        2:(1+peptide_length)]))
  }
  weights <- rowMeans(weight_folds)
  weights <- weights/max(weights) #normalize
  
  weights_all[['position']] <- append(weights_all[['position']],1:peptide_length)
  weights_all[['weight']] <- append(weights_all[['weight']],weights)
  weights_all[['is_pooled']] <- append(weights_all[['is_pooled']],
                                       rep('pooled_classification',peptide_length))
  weights_all[['aa_method']] <- append(weights_all[['aa_method']],
                                       rep('Inferred_symmetric',peptide_length))
  weights_all[['TCR']] <- append(weights_all[['TCR']],
                                 rep(tcr,peptide_length))
  
  # Weights+ Full AA matrix
  weight_folds <- array(0,dim = c(peptide_length,5))
  for (fold in 1:5){
    weights_all_tcr <-read.csv(paste0('../data/within_tcr_pooled/full_aa_matrix/',
                                      'weights_all_',peptide_length,
                                      '_mer_fold_',fold-1,'.csv'))
    weight_folds[,fold] <- array(unlist(weights_all_tcr[weights_all_tcr[,1]==tcr,
                                                        2:(1+peptide_length)]))
  }
  weights <- rowMeans(weight_folds)
  weights <- weights/max(weights) #normalize
  
  weights_all[['position']] <- append(weights_all[['position']],1:peptide_length)
  weights_all[['weight']] <- append(weights_all[['weight']],weights)
  weights_all[['is_pooled']] <- append(weights_all[['is_pooled']],
                                       rep('pooled_classification',peptide_length))
  weights_all[['aa_method']] <- append(weights_all[['aa_method']],
                                       rep('Inferred_full',peptide_length))
  weights_all[['TCR']] <- append(weights_all[['TCR']],
                                 rep(tcr,peptide_length))
  
}

weights_all <- data.frame(weights_all)

weights_all %>%
mutate(is_pooled=fct_inorder(is_pooled))

levels <- unique(weights_all$is_pooled)

method_levels <- unique(weights_all$aa_method)
##########################
####    Plotting #########
##########################
p<- ggplot(weights_all, aes(x=factor(position), y=weight, 
                            fill=factor(is_pooled,levels=levels))) +
  geom_bar(fun = "median",stat='summary',show.legend = FALSE, width = 0.4, 
           position = position_dodge(width = 0.5),alpha=1) +
  geom_point(aes(color=factor(aa_method, levels = method_levels),
                 group=factor(is_pooled,levels=levels)),
             alpha=0.85, size=2,
             show.legend = TRUE,
             position = position_jitterdodge(dodge.width=0.45)) +
    scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000',
                                '#0047AB','#E9C44F','#808080','#bc80bd')) + 
      stat_summary(fun = "median",fun.min=min, fun.max=max,
               geom = "errorbar",
               position = position_dodge(0.5),
               width = 0.3) +
      labs(y="Positional weight",
       x="Position",
       color="") +
  #  xlim(NA, 0.9) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=1, size=25),
        axis.text.y = element_text(size=25),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=30),
        legend.text =element_text(size=20)) +
  facet_wrap(vars(TCR),ncol = 2)+
  theme(strip.text = element_text(
    size = 25, color = "black")) +
  theme(legend.position = "bottom")+ 
  guides(color = guide_legend(nrow = 2), size=25)

ggsave(plot=p, "../figures/extended_data_fig5/positional_weights_by_tcr.pdf",
       width=20, height=15)




