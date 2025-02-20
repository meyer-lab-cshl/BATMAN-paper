rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

##################
### Load Data ####
##################
batman_aucs <- read.csv('BATMAN_within_repertoire_LOO_AUCs.csv')
batman_aucs$method <- 'BATMAN'
  
pteam_aucs <- read.csv('pTEAM_within_repertoire_LOO_AUCs.csv')

auc_all <- rbind(batman_aucs,pteam_aucs)

###################################
######## Plotting: MHCI ###########
###################################

p <- ggplot(auc_all, aes(y=auc,x=factor(index_peptide,
                                        levels=c("FRDYVDRFYKTLRAEQASQE","FEAQKAKANKAVD",
                                                 "VVVGAVGVGK","TPQDLNTML","VPSVWRSSL","FMNKFIYEI",           
                                                 "SLLMWITQC","IMDQVPFSV","NLVPMVATV",
                                                 "SIINFEKL")),
)) + 
  geom_boxplot(aes(color = method,
    fill = method),
               outlier.colour = NA,
               show.legend = TRUE, alpha=1,
               position=position_dodge(0.8)) +
  scale_color_manual(values=c('#AA4A44'	,'#5072A7','#0066b2')) +
  scale_fill_manual(values=c("#FFD6D7","#F0F8FF","#B9D9EB")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(index_peptide,
                              levels=c("FRDYVDRFYKTLRAEQASQE","FEAQKAKANKAVD",
                                       "VVVGAVGVGK","TPQDLNTML","VPSVWRSSL","FMNKFIYEI",           
                                       "SLLMWITQC","IMDQVPFSV","NLVPMVATV",
                                       "SIINFEKL")), 
                 group=method), 
             show.legend = TRUE,
             alpha=0.8,
             size=2,
             position = position_jitterdodge(jitter.width=0.2,
                                             dodge.width = 0.8)) +
   labs(y="TCR activation prediction AUC",
       x="Index peptide",
       color="") +
  theme_cowplot() +
#  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_blank(),#text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
#  ylim(0.6,0.95) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 1))

# Save plot
ggsave(plot=p, "methods_aucs_for_2f.pdf",
       width=9, height=2.7)



# Save raw data
write.csv(auc_all, "raw_data_fig_2f.csv")










