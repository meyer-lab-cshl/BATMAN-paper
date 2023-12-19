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
stats <- read_csv("../data/methods_stats.csv")
stats <- stats %>% select(method,pairs_per_TCR,negative_by_total,positive_by_total)
################
## analysis ####
################
# plot
p <- ggplot(stats, aes(y=pairs_per_TCR, x=1/abs(positive_by_total-negative_by_total), label=method)) + 
  geom_point(color='#412e8e',size=9) +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma = 0.1,base = 10),breaks=c(0,1,10,50,100,175)) +
  scale_x_continuous(breaks=c(1,6,11),limits=c(0.8,11)) +
#  scale_x_continuous(trans=scales::pseudo_log_trans(sigma = 0.1,base = 10),limits=c(0,175),breaks=c(0,1,2,5,10,50,150))+
  theme_cowplot(font_size = 30) + geom_text_repel(max.overlaps=27,size = 9) +
  xlab("Full data/|Positive-Negative|") + ylab("#pMHC-TCR pairs per TCR")
  
ggsave(plot=p, "../figures/fig1/methods.pdf",
       width=8, height=10)
