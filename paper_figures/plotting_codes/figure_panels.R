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
library(ggnewscale)

###############
##  plot ######
###############
# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])
load("all_AUC.rda")

# Special plotting specs for SIINFEKL data
all_AUC <- all_AUC %>%
    mutate(color_by = tcr,
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
  
    p <- ggplot(x, aes(x=AA_distances, y=score)) +
        geom_boxplot(aes(color=factor(class_pair,levels = unique(class_pair))),
                     outlier.colour = NA,
                     show.legend = FALSE, alpha=1,lwd=1,
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
        ylim(min(all_AUC$score,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
        theme_cowplot() +
        # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
        #       axis.text.y = element_text(size=20),
        #       axis.title.y = element_text(size=20)) +
        theme(axis.text.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=20))+
        facet_wrap(vars(index_peptide))+
        theme(strip.text = element_text(
            size = 25, color = "black")) +
        theme(legend.text = element_text(size=15)) +
       guides(color = guide_legend(ncol = 3)) +
       theme(legend.position="bottom",
             legend.box.background = element_rect(colour = "black"))
      # 
})



# boxplot_by_index <- lapply(seq_along(boxplot_by_index), function(x) {
#     if (x < (length(boxplot_by_index)+1)) {
#         boxplot_by_index[[x]] +
#             theme(axis.text.x = element_blank())
#     } else {
#         boxplot_by_index[[x]]
#     }
#     })



p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=5,nrow=3,
                   align = "hv",axis="tb")

ggsave(plot=p_all, "auc_by_index.pdf",
       width=20, height=12)
