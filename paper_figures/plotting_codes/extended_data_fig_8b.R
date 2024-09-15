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
colnames(plot_data) <- c('method','tcr','score','sd')

####### pTEAM results original ########
# Add pTEAM results
within_wider <- read_csv("../data/within_tcr1.csv")
cross_wider <- read_csv("../data/cross_tcr1.csv")

# format data into long form
within <- within_wider %>%
  pivot_longer(-Row, names_to="tcr", values_to = "score") %>%
  rename(method=Row) %>%
  mutate(method=fct_inorder(method))

within <- within[within$method=="Atchley+RF",]
within$method <- "Within TCR"
within$type <- "pTEAM (original)"
plot_data <- rbind(plot_data,within)

cross <- cross_wider %>%
  pivot_longer(-Row, names_to="tcr", values_to = "score") %>%
  rename(method=Row) %>%
  mutate(method=fct_inorder(method))

cross <- cross[cross$method=="Atchley+RF",]
cross$method <- "Leave one TCR out"
cross$type <- "pTEAM (original)"
plot_data <- rbind(plot_data,cross)

pTEAM_AL <- read_csv("../data/AL_results_pTEAM.csv")[1:10,]
pTEAM_AL <- pTEAM_AL %>%
  pivot_longer(-...1, names_to="tcr", values_to = "score") %>%
  rename(method=...1)
pTEAM_AL$method <- (pTEAM_AL$method+1)*9
pTEAM_AL$type <- 'pTEAM (original)'

plot_data <- rbind(plot_data,pTEAM_AL)

pTEAM_AL <- read_csv("../data/AL_results_pTEAM_v2.csv")[1:10,]
pTEAM_AL <- pTEAM_AL %>%
  pivot_longer(-...1, names_to="tcr", values_to = "score") %>%
  rename(method=...1)
pTEAM_AL$method <- (pTEAM_AL$method+1)*9
pTEAM_AL$type <- 'pTEAM (modified)'

plot_data <- rbind(plot_data,pTEAM_AL)



################
### Plotting  ##
################

# Plotting Extended Data Figure
p <- ggplot(plot_data, aes(y=score, x=factor(method, levels = c("Leave one TCR out",
                                                                "9","18","27","36","45",
                                                                "54","63","72","81","90",
                                                                "Within TCR"))))

p <- p + geom_boxplot(aes(color=factor(type, levels = c("pTEAM (original)",
                                                        "pTEAM (modified)"))),
                      outlier.colour = NA,
                      show.legend = TRUE, alpha=1,lwd=0.6,
                      position=position_dodge(0.8)) +
  scale_color_manual(values=c('#1b9e77','#b2df8a')) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color=tcr,
                  group = factor(type, levels = c("pTEAM (original)",
                                                  "pTEAM (modified)"))), 
              alpha=0.9, 
              position = position_jitterdodge(jitter.width=0.05,
                                              dodge.width = 0.8),size=4) +
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


ggsave(plot=p, "../figures/extended_data_fig8/auc_al_pTEAM.pdf",width=15, height=10)

