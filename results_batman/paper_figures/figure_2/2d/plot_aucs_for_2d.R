rm(list=ls())#Clear environment

# Packages
library(readxl)
library(pROC)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# selected TCRs
TCR_list_9mer <- c('a3a','TIL1383I','TCR81-14','18A2','NYE-S1','TCR6','TCR7',
                   'A6','T1','FLT3DY','A23','TCR2-T','R24','868Z11')

TCR_list_10mer <- c('TCR-1E6', 'APN-TCR', 'EWW-TCR', 'A11Va')

TCR_list_mhcii<- c('TCR-F5','TCR-3598-2','MBP-TCR','B3K508')

all_TCRs <- c(TCR_list_9mer,TCR_list_10mer,TCR_list_mhcii)

# Read pTEAM AUCs
all_aucs <- read.csv('pTEAM_AUCs.csv')[,2:4]
all_aucs$test <- 'LOO TCR'
all_aucs$test[all_aucs$type=='pTEAM_within'] <- 'Within TCR'

all_aucs$type[all_aucs$type %in% c('pTEAM_within','pTEAM_LOO_without_seq')] <- 'pTEAM'
all_aucs$type[all_aucs$type=='pTEAM_LOO_with_seq'] <- 'pTEAM+TCR'

# Read BATMAN AUCs
batman_aucs_1 <- read.csv('BATMAN_AUCs.csv')[,1:2]
colnames(batman_aucs_1) <- c('tcr','AUC')
batman_aucs_1$type <- 'BATMAN'
batman_aucs_1$test <- 'Within TCR'

batman_aucs_2 <- read.csv('BATMAN_AUCs.csv')[,c(1,3)]
colnames(batman_aucs_2) <- c('tcr','AUC')
batman_aucs_2$type <- 'BATMAN'
batman_aucs_2$test <- 'LOO TCR'

all_aucs <- rbind(all_aucs,batman_aucs_1,batman_aucs_2)

# Read BATMAN (no MHC) AUCs
batman_aucs_1 <- read.csv('BATMAN_no_mhc_AUCs.csv')[,1:2]
colnames(batman_aucs_1) <- c('tcr','AUC')
batman_aucs_1$type <- 'BATMAN (no MHC)'
batman_aucs_1$test <- 'Within TCR'

batman_aucs_2 <- read.csv('BATMAN_no_mhc_AUCs.csv')[,c(1,3)]
colnames(batman_aucs_2) <- c('tcr','AUC')
batman_aucs_2$type <- 'BATMAN (no MHC)'
batman_aucs_2$test <- 'LOO TCR'

all_aucs <- rbind(all_aucs,batman_aucs_1,batman_aucs_2)

# Read BATMAN (BLOSUM100) AUCs
batman_aucs_1 <- read.csv('BATMAN_BLOSUM100_AUCs.csv')[,1:2]
colnames(batman_aucs_1) <- c('tcr','AUC')
batman_aucs_1$type <- 'BATMAN (BLOSUM100)'
batman_aucs_1$test <- 'Within TCR'

batman_aucs_2 <- read.csv('BATMAN_BLOSUM100_AUCs.csv')[,c(1,3)]
colnames(batman_aucs_2) <- c('tcr','AUC')
batman_aucs_2$type <- 'BATMAN (BLOSUM100)'
batman_aucs_2$test <- 'LOO TCR'

all_aucs <- rbind(all_aucs,batman_aucs_1,batman_aucs_2)

###################################
######## Plotting #################
###################################

p <- ggplot(all_aucs, aes(y=AUC,x=factor(type,levels=c("BATMAN","BATMAN (no MHC)",
                                                       "BATMAN (BLOSUM100)","pTEAM",
                                                       "pTEAM+TCR")),
)) + 
  geom_boxplot(aes(color = factor(test,levels=c('Within TCR','LOO TCR')),
                   fill = factor(type, levels=c("BATMAN","BATMAN (no MHC)",
                                                "BATMAN (BLOSUM100)","pTEAM",
                                                "pTEAM+TCR"))),
               outlier.colour = NA,
               show.legend = FALSE, alpha=1,
               position=position_dodge2(0.8,preserve = "single")) +
  scale_color_manual(values=c('#007FFF','black')) +
  scale_fill_manual(values=c("#FFD6D7","#FFD6D7","#FFD6D7", "#F0FFFF","#F0FFFF")) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color=factor(tcr,levels = all_TCRs), 
                 group=factor(test,levels=c('Within TCR','LOO TCR'))), 
             alpha=0.85,
             size=2,
             position = position_jitterdodge(jitter.width=1.4,
                                             dodge.width = 0.8)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set1"))(22)) +
  labs(y="TCR activation prediction AUC",
       x="",
       color="") +
  theme_cowplot() +
  geom_hline(yintercept=0.5, linetype='dashed', col = '#966919') +
  theme(axis.text.x = element_text(size=15,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "methods_aucs_for_2d.pdf",
       width=6.5, height=5)

# Save raw data
write.csv(all_aucs, "raw_data_fig_2d.csv")





  







