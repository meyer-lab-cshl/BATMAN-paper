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

plot_data <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(plot_data) <- c('method','tcr','score','type')

#################################
### Load and average AL aucs  ###
#################################

### AL all steps where only >1 mutant/site is taken #######

for (npeptide in c(2,3,4,6,8,12)){
  
  count <- 0 #counts number of random runs
  
  AL_auc_mean <- data.frame(array(0,dim = c(as.integer(12/npeptide),12)))
  
  colnames(AL_auc_mean) <- c("","18A2","868Z11","A23","A6","NYE-S1","TCR2-T","TCR3",
                             "TCR6","TCR7","T1","FLT3DY")
  
  for (AL_seed in 1:100){
    filename <- paste0('../data/AL_cluster/data/auc_npeptide_',npeptide,'_nstep_',as.integer(12/npeptide),'_',AL_seed,'.csv')
    if (file.exists(filename)){
      count <- count +1
      AL_data_wider <- read.csv(filename,check.names = FALSE)
      AL_auc_mean <- AL_auc_mean + AL_data_wider
    }
  }
  AL_auc_mean <- AL_auc_mean/count
  colnames(AL_auc_mean)[1] <- "X"
  
  # Convert to steps
  AL_auc_mean$X <- (AL_auc_mean$X+1)
  
  # convert data to long form
  AL_data <- AL_auc_mean %>%
    pivot_longer(-X, names_to="tcr", values_to = "score") %>%
    rename(method=X)
  
  total_npeptide <-as.integer(AL_data$method)*npeptide*9
  
  AL_data$type <- 'Total #peptides = '
  AL_data$type <- paste0(AL_data$type,as.character(total_npeptide))
  
  # add to plotting data
  if (count>0){
    plot_data <- rbind(plot_data,AL_data)
  }
  
}

##### AL ##############

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

AL_auc_mean <- AL_auc_mean[4,]

# convert data to long form
AL_data <- AL_auc_mean %>%
  pivot_longer(-X, names_to="tcr", values_to = "score") %>%
  rename(method=X)
AL_data$type <- 'Total #peptides = 36'
AL_data$method <- 4

# add to plotting data
plot_data <- rbind(plot_data,AL_data)

# filter for only certain total #peptide
plot_data <- plot_data[plot_data$type %in% c("No AL","Total #peptides = 36",
                                             "Total #peptides = 72",
                                             "Total #peptides = 108"),]

plot_data <- plot_data[plot_data$method %in% c("Leave one TCR out",
                                               "1","2","4",
                                               "Within TCR"),]


################
### Plotting  ##
################

# Plotting Extended Data Figure
p <- ggplot(plot_data, aes(y=score, x=factor(method, levels = c("1","2","4"))))

p <- p + geom_boxplot(aes(color=factor(type, levels = c("Total #peptides = 36",
                                                        "Total #peptides = 72",
                                                        "Total #peptides = 108"))),
                      outlier.colour = NA,
                      show.legend = FALSE, alpha=1,lwd=0.6,
                      position=position_dodge(0.8)) +
  scale_color_manual(values=c('#fec44f','#fe9929','#ec7014')) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color=tcr,
                  group = factor(type)), 
              alpha=0.9, 
              position = position_jitterdodge(jitter.width=0.05,
                                              dodge.width = 0.8),size=4) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB',
                                       '#EC5800','#E4D00A','#FAA0A9','#808080',
                                       '#bc80bd','#50C878')) + 
                                         labs(y="Average AUC",
                                              x="Number of experimental rounds",
                                              color="") +
  ylim(0.475, 0.9) +
  theme_cowplot() +
  #  theme(legend.direction="vertical", legend.position="right") +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=1, size=45),
        axis.title.y = element_text(size=45),
        axis.title.x = element_text(size=45),
        axis.text.y = element_text(size=45))


ggsave(plot=p, "../figures/fig2/auc_al_npeptide.pdf",width=15, height=10)

