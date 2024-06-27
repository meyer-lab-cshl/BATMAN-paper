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


auc_all <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(auc_all) <- c('method','tcr','score','sd')


# within and pan TCR results
# Load peptide information
peptide_folds <- read.csv('fig_2_training_files/full_training_data.csv')
peptide_folds$activation <- peptide_folds$activation + 1
peptide_folds$fold <- peptide_folds$fold + 1 #Python to R

# List of 3 class 9-mer-binding TCRs
tcr_names <- c("18A2","NYE-S1","868Z11","TCR3","TCR6","TCR7","A6",
               "T1","FLT3DY","A23","TCR2-T")

# Load peptide data for selected TCRs
peptide_folds <- peptide_folds[peptide_folds$tcr %in% tcr_names,]

# peptide scores from BATMAN within-TCR and cross-TCR for full and symm matrices
peptide_score <- array(NA,dim = c(dim(peptide_folds)[1],2))

# BATMAN scores
output_data <- read.csv('fig_2_training_files/BATMAN_scores_full_matrix.csv')

peptide_score[,1] <- 
  as.double(output_data$within_tcr_score[match(peptide_folds$peptide,
                                               output_data$peptide)])
peptide_score[,2] <- 
  as.double(output_data$cross_tcr_score[match(peptide_folds$peptide,
                                              output_data$peptide)])

method_name <- c('BATMAN Full')

# Calculate AUCs for within and cross-TCR modes
# Create and shape the full dataFrame for plotting
plot_data <- data.frame(matrix(nrow=0,ncol=4))
colnames(plot_data) <- c('tcr','method','score','type')

for (tcr in tcr_names){
  training_folds_tcr <- peptide_folds[peptide_folds$tcr == tcr,]
  peptide_score_tcr <- peptide_score[peptide_folds$tcr == tcr,]
  
  # BATMAN scores
  # Get AUC for 5 folds
  AUCs <- array(NA,dim = c(3,5))
  for (fold in 1:5){
    # Full matrix
    ROC<-multiclass.roc(training_folds_tcr$activation[training_folds_tcr$fold==fold],
                        peptide_score_tcr[training_folds_tcr$fold==fold,
                                          1],direction=">")
    
    AUCs[,fold] <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
    
  }
  # Within-TCR scores
  new_data <- data.frame(tcr = tcr,
                         method = 'Within TCR',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
  
  # Full TCR data
  ROC<-multiclass.roc(training_folds_tcr$activation,
                      peptide_score_tcr[,2],direction=">")
  AUCs <- c(auc(ROC$rocs[[1]]),auc(ROC$rocs[[2]]),auc(ROC$rocs[[3]]))
  new_data <- data.frame(tcr = tcr,
                         method = 'Leave one TCR out',
                         score = mean(AUCs))
  plot_data <- rbind(plot_data,new_data)
}

plot_data$type <- 'No AL'


#################################
### Load and average AL aucs  ###
#################################

##### only random ##############

count <- 0 #counts number of random runs
AL_auc_mean <- data.frame(array(0,dim = c(10,12)))
colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                           "TCR6","TCR7","T1","FLT3DY")

for (AL_seed in 1:100){
  filename <- paste0('../data/AL_cluster/data/auc_random_',AL_seed,'.csv')
  if (file.exists(filename)){
    count <- count +1
    AL_data_wider <- read.csv(filename,check.names = FALSE)
    AL_auc_mean <- AL_auc_mean + AL_data_wider
  }
}
AL_auc_mean <- AL_auc_mean/count
colnames(AL_auc_mean)[1] <- "X"

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$method <- (AL_data$method+1)*9
AL_data$type <- 'All random'

# add to plotting data
plot_data <- rbind(plot_data,AL_data)

##### AL and random alternating ##############

count <- 0 #counts number of random runs
AL_auc_mean <- data.frame(array(0,dim = c(10,12)))
colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                           "TCR6","TCR7","T1","FLT3DY")

for (AL_seed in 1:100){
  filename <- paste0('../data/AL_cluster/data/auc_active_random_alternate_',AL_seed,'.csv')
  if (file.exists(filename)){
    count <- count +1
    AL_data_wider <- read.csv(filename,check.names = FALSE)
    AL_auc_mean <- AL_auc_mean + AL_data_wider
  }
}
AL_auc_mean <- AL_auc_mean/count
colnames(AL_auc_mean)[1] <- "X"

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$method <- (AL_data$method+1)*9
AL_data$type <- 'Even steps random'

# add to plotting data
plot_data <- rbind(plot_data,AL_data)

##### AL and then all random##############

count <- 0 #counts number of random runs
AL_auc_mean <- data.frame(array(0,dim = c(10,12)))
colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                           "TCR6","TCR7","T1","FLT3DY")

for (AL_seed in 1:100){
  filename <- paste0('../data/AL_cluster/data/auc_active_then_random_',AL_seed,'.csv')
  if (file.exists(filename)){
    count <- count +1
    AL_data_wider <- read.csv(filename,check.names = FALSE)
    AL_auc_mean <- AL_auc_mean + AL_data_wider
  }
}
AL_auc_mean <- AL_auc_mean/count
colnames(AL_auc_mean)[1] <- "X"

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$method <- (AL_data$method+1)*9
AL_data$type <- 'All random from 2nd step'

# add to plotting data
plot_data <- rbind(plot_data,AL_data)

##### AL and random alternating ##############

count <- 0 #counts number of random runs
AL_auc_mean <- data.frame(array(0,dim = c(10,12)))
colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                           "TCR6","TCR7","T1","FLT3DY")

for (AL_seed in 1:100){
  filename <- paste0('../data/AL_cluster/data/auc_active_random_mixed_',AL_seed,'.csv')
  if (file.exists(filename)){
    count <- count +1
    AL_data_wider <- read.csv(filename,check.names = FALSE)
    AL_auc_mean <- AL_auc_mean + AL_data_wider
  }
}
AL_auc_mean <- AL_auc_mean/count
colnames(AL_auc_mean)[1] <- "X"

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$method <- (AL_data$method+1)*9
AL_data$type <- '3 random each step'

# add to plotting data
plot_data <- rbind(plot_data,AL_data)

##### AL and random alternating ##############

count <- 0 #counts number of random runs
AL_auc_mean <- data.frame(array(0,dim = c(10,12)))
colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                           "TCR6","TCR7","T1","FLT3DY")

for (AL_seed in 1:100){
  filename <- paste0('../data/AL_cluster/data/auc_active_',AL_seed,'.csv')
  if (file.exists(filename)){
    count <- count +1
    AL_data_wider <- read.csv(filename,check.names = FALSE)
    AL_auc_mean <- AL_auc_mean + AL_data_wider
  }
}
AL_auc_mean <- AL_auc_mean/count
colnames(AL_auc_mean)[1] <- "X"

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$method <- (AL_data$method+1)*9
AL_data$type <- 'No random'

# add to plotting data
plot_data <- rbind(plot_data,AL_data)




################
### Plotting  ##
################

# Plotting Extended Data Figure
p <- ggplot(plot_data, aes(y=score, x=factor(method, levels = c("Leave one TCR out",
                                                                "9","18","27","36","45",
                                                                "54","63","72","81","90",
                                                                "Within TCR"))))

p <- p + geom_boxplot(aes(color=factor(type, levels = c("No AL","No random",
                                                        "All random","Even steps random",
                                                        "All random from 2nd step",
                                                        "3 random each step"))),
                      outlier.colour = NA,
                      show.legend = FALSE, alpha=1,lwd=0.4,
                      position=position_dodge(0.9)) +
  scale_color_manual(values=c('#000000','#984ea3','#7570b3','#e7298a','#e6ab02','#a6761d')) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color=tcr,
                  group = factor(type, levels = c("No AL","No random",
                                                  "All random","Even steps random",
                                                  "All random from 2nd step",
                                                  "3 random each step"))), 
              alpha=0.9, 
              position = position_jitterdodge(jitter.width=0.05,
                                              dodge.width = 0.9),size=2) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
                                         labs(y="Average AUC",
                                              x="Number of peptides sampled",
                                              color="") +
  ylim(NA, 0.9) +
  theme_cowplot() +
  #  theme(legend.direction="vertical", legend.position="right") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=40),
        axis.title.y = element_text(size=45),
        axis.title.x = element_text(size=45),
        axis.text.y = element_text(size=45))


ggsave(plot=p, "../figures/extended_data_fig8/auc_al_strategies.pdf",width=15, height=10)

