# TCR activation by mutated peptides
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

# Discard patented TCRs
TCR_data <- TCR_data[!(TCR_data$tcr_name %in% c('T1','T3','FLT3DY')),]

# List of TCRs
TCR <- unique(TCR_data$tcr_name)

# Column to store mutation position
TCR_data$position <- 0*TCR_data$index_peptide_activity

# Normalize peptide activity data
for (TCR_name in TCR){
  TCR_data$peptide_activity[TCR_data$tcr_name==TCR_name] <-
    (TCR_data$peptide_activity[TCR_data$tcr_name==TCR_name])/
    max(TCR_data$peptide_activity[TCR_data$tcr_name==TCR_name])


# Record mutation position
index_peptide <- TCR_data$index_peptide[TCR_data$tcr_name==TCR_name]
mutant_peptide <- TCR_data$peptide[TCR_data$tcr_name==TCR_name]

peptide_length <- nchar(unique(index_peptide))
n_peptide <- length(mutant_peptide) #number of peptides

TCR_data$position[TCR_data$tcr_name==TCR_name] <- 
  rowSums((matrix(unlist(strsplit(index_peptide,split='')),ncol=peptide_length,byrow=TRUE)!=
  matrix(unlist(strsplit(mutant_peptide,split='')),ncol=peptide_length,byrow=TRUE))*
  t(array(1:peptide_length,dim=c(peptide_length,n_peptide))))

}

TCR_data <- TCR_data[TCR_data$position!=0,] #Discard index peptides

#  Indicate anchor residues
TCR_data$is_anchor <- array('non-anchor',dim=dim(TCR_data)[1])
TCR_data$is_anchor[nchar(TCR_data$peptide)==10 & 
                     TCR_data$position %in% c(1,2,10)]='anchor'
TCR_data$is_anchor[nchar(TCR_data$peptide)==9 & 
                     TCR_data$position %in% c(1,2,9)]='anchor'
TCR_data$is_anchor[nchar(TCR_data$peptide)==8 & 
                     TCR_data$position %in% c(2,3,5,8)]='anchor'

TCR_data %>%
  mutate(index_peptide=fct_inorder(index_peptide))

# For SIINFEKL TCRs, classify them into 2 categories
tcr_type <- array(NaN,dim=c(dim(TCR_data)[1]))

naive_tcrs <- c("B11","B15","B3", "F4","E8", "B13","H6", "G6", "F5", 
                "H5", "B2", "B6", "B5","E9", "E4", "G2", "B16","B14")
ed_tcrs <- c("OT1","Ed5","Ed8","Ed9","Ed10","Ed16-1","Ed16-30","Ed21","Ed23",
             "Ed28","Ed31","Ed33","Ed39","Ed40","Ed45","Ed46")

tcr_type[TCR_data$tcr_name %in% naive_tcrs] <- "Naive"
tcr_type[TCR_data$tcr_name %in% ed_tcrs] <- "Educated"

TCR_data$tcr_type <- tcr_type

#######################
#### Plot #############
#######################

# Special plotting specs for SIINFEKL data
TCR_data <- TCR_data %>%
  mutate(color_by = tcr_name,
         color_by = case_when(index_peptide == "SIINFEKL" ~ tcr_type,
                              TRUE ~ color_by))

TCR_data_by_index <- split(TCR_data, 
                           factor(TCR_data$index_peptide,
                                  levels = unique(TCR_data$index_peptide)))

boxplot_by_index <- lapply(TCR_data_by_index, function(x) {
  
  set_alpha <- 0.5 
  set_size <- 1
  set_jitter_width <- 0.15
  
  if (unique(x$index_peptide)=="SIINFEKL"){
    set_alpha <- 0.25
    set_size <- 0.6
    set_jitter_width <- 0.3
  }
  
  p <- ggplot(x, 
              aes(x=factor(position), y=normalized_peptide_activity)) +
    geom_boxplot(aes(color=factor(is_anchor)),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1) +
    scale_color_manual(values=c('#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(color_by)), 
               alpha=set_alpha,
               size=set_size,
               position = position_jitter(width=set_jitter_width)) +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(y="Normalized peptide activity",
         x="",
         color="") +
    ylim(0,1)+
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(size=10),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=10))+
    facet_wrap(vars(factor(index_peptide,levels=unique(index_peptide))))+
    theme(strip.text = element_text(
      size = 25, color = "black")) +
    theme(legend.text = element_text(size=15))
#    guides(color = guide_legend(ncol = 3)) +
#    theme(legend.position="bottom")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=1,
                            align = "hv",axis="tb")

ggsave(plot=p_all, "../figures/extended_data_fig1/positional_mutation_activity.pdf",
       width=6, height=20)

