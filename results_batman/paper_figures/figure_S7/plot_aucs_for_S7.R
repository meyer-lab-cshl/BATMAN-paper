rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(viridis)
library(forcats)
library(tidyr)
library(RColorBrewer)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# Load unpooled AUCs
auc_all <- read.csv('aucs_Full_unpooled.csv')
colnames(auc_all)[1] <- 'tcr'

# Load auc values for different BATMAN modes
auc_new <- read.csv('aucs_Symmetric_across.csv')
colnames(auc_new)[1] <- 'tcr'
auc_all <- rbind(auc_all, auc_new)

auc_new <- read.csv('aucs_Symmetric_within.csv')
colnames(auc_new)[1] <- 'tcr'
auc_all <- rbind(auc_all, auc_new)

auc_new <- read.csv('aucs_Full_across.csv')
colnames(auc_new)[1] <- 'tcr'
auc_all <- rbind(auc_all, auc_new)

auc_new <- read.csv('aucs_Full_within.csv')
colnames(auc_new)[1] <- 'tcr'
auc_all <- rbind(auc_all, auc_new)

###################################################################
### Load index peptides corresponding to TCR names      ###########
TCR_data <- read_excel(paste0('../../tcr_epitope_datasets/mutational_scan_datasets',
                              '/train_test_data_folds.xlsx'))
TCR_data <- unique(TCR_data[c('tcr','index_peptide')])
TCR_data <- data.frame(TCR_data, row.names = 1)

# Add index peptides
auc_all$index_peptide <- TCR_data[auc_all$tcr,]

auc_all <- pivot_longer(auc_all,
                              cols = c('WS','NS','NW'),
                              names_to = 'pair',
                              values_to = 'auc')

# Split df by index peptide to plot them separately
auc_by_index <- split(auc_all, f=auc_all$index_peptide)

boxplot_by_index <- lapply(auc_by_index, function(x) {
  
  set_alpha <- 0.85 
  set_size <- 1
  set_jitter_width <- 0.15
  
  p <- ggplot(x, aes(x=fct_inorder(factor(type)), y=auc)) +
    geom_boxplot(aes(color=factor(pair,levels = unique(pair))),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1,
                 position=position_dodge(0.75)) +
    scale_color_manual(values=c('#50C878','#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(tcr), 
                   group=factor(pair,levels = unique(pair))), 
               alpha=set_alpha,
               size=set_size,show.legend = FALSE,
               position = position_jitterdodge(jitter.width=set_jitter_width,
                                               dodge.width = 0.75)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(15, "Dark2"))(42)) +
    labs(y="Classification AUC",
         x="",
         color="") +
    #    ylim(min(auc_df_subset,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=15),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=20))+
    facet_wrap(~index_peptide,scales='free_y')+
    theme(strip.text = element_text(
      size = 20, color = "black")) +
    theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
    theme(legend.text = element_text(size=10)) +
    guides(color = guide_legend(ncol = 2)) +
    theme(legend.position="right")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index,ncol=5,align = "hv",axis="tb")

# Save plot
ggsave(plot=p_all, "aucs_for_figS7.pdf",
       width=15, height=12)

# Save raw data
write.csv(auc_all, "raw_data_fig_S7.csv")











