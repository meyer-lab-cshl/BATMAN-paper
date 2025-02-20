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

###############################
## BATMAN with AA matrices#####
###############################
# List all files
auc_files <- list.files(path = "aucs_TCR_AA_matrix/",
                           pattern = "auc_.*.csv",
                           full.names = TRUE)
# open empty df for storing
auc_df <- data.frame()
for (file in auc_files){
  # Read data
  file_data <- read.csv(file) 
  auc_df <- rbind(auc_df, file_data)
}

auc_df$mean_auc <- (auc_df$WS + auc_df$NS +auc_df$NW)/3

# Make AA matrix names in all caps
auc_df$aa_matrix <- toupper(auc_df$aa_matrix)

###################################
######## Plotting #################
###################################

p <- ggplot(auc_df,aes(y=mean_auc,x=fct_inorder(factor(aa_matrix)))) + 
  geom_boxplot(outlier.colour = NA,
               show.legend = FALSE) +
  geom_jitter(aes(color=tcr,group=aa_matrix),
              alpha=0.5, width=0.1, height=0,size=0.1, show.legend = FALSE) +
  #  scale_colour_manual(values = palette(viridis(151))) +
  labs(y="Average classification AUC",
       x="",
       color="") +
  #  xlim(NA, 0.9) +
  theme_cowplot() +
  theme_minimal_grid(14) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=5),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank())
#  theme(legend.position = "bottom")+ 
#  guides(color = guide_legend(nrow = 9))

# Save plot
ggsave(plot=p, "auc_for_S5b.pdf",
       width=7, height=3)


###################################################################
### Load index peptides corresponding to TCR names      ###########
TCR_data <- read_excel(paste0('../../../tcr_epitope_datasets/mutational_scan_datasets',
                              '/train_test_data_folds.xlsx'))
TCR_data <- unique(TCR_data[c('tcr','index_peptide')])
TCR_data <- data.frame(TCR_data, row.names = 1)

# Plot for selected AA matrices
auc_df_subset <- auc_df[auc_df$aa_matrix %in% c('HAMMING','BLOSUM100','PAM10',
                                                'DAYHOFF','GONNET'),]

# Add index peptides
auc_df_subset$index_peptide <- TCR_data[auc_df_subset$tcr,]

auc_df_subset <- pivot_longer(auc_df_subset,
                              cols = c('WS','NS','NW'),
                              names_to = 'pair',
                              values_to = 'auc')

# Split df by index peptide to plot them separately
auc_by_index <- split(auc_df_subset, f=auc_df_subset$index_peptide)

boxplot_by_index <- lapply(auc_by_index, function(x) {
  
  set_alpha <- 0.85 
  set_size <- 1
  set_jitter_width <- 0.15
  
  p <- ggplot(x, aes(x=fct_inorder(factor(aa_matrix)), y=auc)) +
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
    labs(y="AUC",
         x="",
         color="") +
#    ylim(min(auc_df_subset,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=15))+
    facet_wrap(vars(index_peptide))+
    theme(strip.text = element_text(
      size = 20, color = "black")) +
    theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
    theme(legend.text = element_text(size=10)) +
    guides(color = guide_legend(ncol = 2)) +
    theme(legend.position="right")
  
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index,ncol=5,align = "hv",axis="tb")

# Save plot
ggsave(plot=p_all, "auc_for_S5c.pdf",
       width=15, height=15)

# Save raw data
write.csv(auc_df, "raw_data_fig_S5b-c.csv")











