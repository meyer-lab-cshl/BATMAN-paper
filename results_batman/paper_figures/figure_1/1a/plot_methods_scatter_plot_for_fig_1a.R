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
library(ggrepel)


############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is
stats <- read_csv("methods_stats.csv")
stats <- stats %>% select(method,pairs_per_TCR,data_balance)
################
## analysis ####
################
# plot
p <- ggplot(stats, aes(y=pairs_per_TCR, x=data_balance, label=method)) + 
  geom_point(color='#412e8e',size=9) +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma = 0.1,base = 10),breaks=c(0,1,10,50,100,175)) +
  scale_x_continuous(breaks=c(1,6,11),limits=c(0.8,11)) +
#  scale_x_continuous(trans=scales::pseudo_log_trans(sigma = 0.1,base = 10),limits=c(0,175),breaks=c(0,1,2,5,10,50,150))+
  theme_cowplot() + geom_text_repel(max.overlaps=27,size = 9) +
  xlab("#Total/|#Positive-#Negative|") + ylab("#pMHC-TCR pairs per TCR") +
  theme(axis.text.x = element_text(size=40),
        axis.text.y = element_text(size=40),
        axis.title.x = element_text(size=40),
        axis.title.y = element_text(size=40))
  
ggsave(plot=p, "methods_scatter_plot_fig_1a.pdf",
       width=8, height=10)

write.csv(stats,'raw_data_fig_1a.csv')
