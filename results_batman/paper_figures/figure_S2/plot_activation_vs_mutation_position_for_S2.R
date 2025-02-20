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
library(RColorBrewer)


#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Get TCR activation data
TCR_data <- read_excel(paste0('../../tcr_epitope_datasets/mutational_scan_datasets',
                              '/train_test_data_folds.xlsx'))

# Discard data from index peptides with not all positions mutated
TCR_data <- TCR_data[!TCR_data$index_peptide=='FEAQKAKANKAVD',]

#Discard data from index peptides themselves
TCR_data <- TCR_data[TCR_data$mutation_position!=0,] 

#  Indicate mutation in anchor residues
TCR_data$is_anchor <- array('non-anchor',dim=dim(TCR_data)[1]) #initialize
TCR_data$is_anchor[nchar(TCR_data$peptide)==10 & 
                     TCR_data$mutation_position %in% c(1,2,10)]='anchor'

TCR_data$is_anchor[nchar(TCR_data$index_peptide)=='ALYDKTKRIFL' & 
                     TCR_data$mutation_position %in% c(1,2,11)] = 'anchor'

TCR_data$is_anchor[nchar(TCR_data$peptide)==9 & 
                     TCR_data$mutation_position %in% c(1,2,9)] = 'anchor'

TCR_data$is_anchor[nchar(TCR_data$peptide)==8 & 
                     TCR_data$mutation_position %in% c(2,3,5,8)] = 'anchor'

TCR_data$is_anchor[TCR_data$index_peptide=='FRDYVDRFYKTLRAEQASQE' & 
                     TCR_data$mutation_position %in% c(9,12)] = 'anchor'

TCR_data$is_anchor[TCR_data$index_peptide=='LPVPGVLLKEFTVSGNILTI' & 
                     TCR_data$mutation_position %in% c(7,10,12,15)] = 'anchor'

TCR_data$is_anchor[TCR_data$index_peptide=='ASQKRPSQRSK' & 
                     TCR_data$mutation_position %in% c(4,5)] = 'anchor'

TCR_data %>%
  mutate(index_peptide=fct_inorder(index_peptide))

# For SIINFEKL TCRs, classify them into 2 categories
tcr_type <- array(NaN,dim=c(dim(TCR_data)[1]))

naive_tcrs <- c("B11","B15","B3", "F4","E8", "B13","H6", "G6", "F5", 
                "H5", "B2", "B6", "B5","E9", "E4", "G2", "B16","B14")
ed_tcrs <- c("OT1","Ed5","Ed8","Ed9","Ed10","Ed16-1","Ed16-30","Ed21","Ed23",
             "Ed28","Ed31","Ed33","Ed39","Ed40","Ed45","Ed46")

tcr_type[TCR_data$tcr %in% naive_tcrs] <- "Naive"
tcr_type[TCR_data$tcr %in% ed_tcrs] <- "Educated"

TCR_data$tcr_type <- tcr_type

#######################
#### Plot #############
#######################

# Special plotting specs for SIINFEKL data
TCR_data <- TCR_data %>%
  mutate(color_by = tcr,
         color_by = case_when(index_peptide == "SIINFEKL" ~ tcr_type,
                              TRUE ~ color_by))

TCR_data_by_index <- split(TCR_data, 
                           factor(TCR_data$index_peptide,
                                  levels = unique(TCR_data$index_peptide)))

boxplot_by_index <- lapply(TCR_data_by_index, function(x) {
  
  set_alpha <- 0.75 
  set_size <- 1
  set_jitter_width <- 0.15
  
  if (unique(x$index_peptide)=="SIINFEKL"){
    set_alpha <- 0.25
    set_size <- 0.6
    set_jitter_width <- 0.3
  }
  
  p <- ggplot(x, 
              aes(x=factor(mutation_position), y=peptide_activity)) +
    geom_boxplot(aes(color=factor(is_anchor)),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1) +
    scale_color_manual(values=c('#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(tcr)), 
               alpha=set_alpha,
               size=set_size,
               position = position_jitter(width=set_jitter_width)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(15, "Dark2"))(42)) +
    labs(y="TCR activation",
         x="Peptide mutation position",
         color="") +
    ylim(0,1)+
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(size=15),
          axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size=20))+
    facet_wrap(vars(factor(index_peptide,levels=unique(index_peptide))))+
    theme(strip.text = element_text(
      size = 10, color = "black")) +
    theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
    theme(legend.text = element_text(size=15))+
      guides(color = guide_legend(ncol = 3))
  #    theme(legend.position="bottom")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=3)

ggsave(plot=p_all, "figS2.pdf",
       width=20, height=20)





