#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(forcats)

############
## data ####
############
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is
within_wider <- read_csv("../data/within_tcr1.csv")
cross_wider <- read_csv("../data/cross_tcr1.csv")
################
## analysis ####
################

# format data into long form
within <- within_wider %>%
  pivot_longer(-Row, names_to="tcr", values_to = "score") %>%
  rename(method=Row) %>%
  mutate(method=fct_inorder(method))

cross <- cross_wider %>%
  pivot_longer(-Row, names_to="tcr", values_to = "score") %>%
  rename(method=Row) %>%
  mutate(method=fct_inorder(method))

# plot
p_within <- ggplot(within, aes(y=method, x=score))
p_within <- p_within + geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(color=tcr), alpha=0.9, width=0, height=0.25) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB','#EC5800','#E4D00A','#FAA0A9','#808080','#bc80bd','#50C878')) + 
  labs(y="Distance function",
       x="AUC score",
       color="") +
  xlim(NA, 0.9) +
  theme_cowplot() +
  theme(legend.direction="horizontal", legend.position="bottom")
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  
p_cross <- ggplot(cross, aes(y=method, x=score))
p_cross <- p_cross + geom_boxplot() +
  geom_jitter(aes(color=tcr), alpha=0.9, width=0, height=0.25) +
  scale_color_manual(values=c('#088F8F','#0096FF','#5D3FD3','#b30000','#0047AB','#EC5800','#E4D00A','#FAA0A9','#808080','#bc80bd','#50C878')) + 
  labs(y="Prediction method",
       x="AUC score",
       color="TCR clone") +
  xlim(NA, 0.9) +
  theme_cowplot()+
  theme(legend.direction="horizontal", legend.position="bottom")
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

 p_all <- plot_grid(p_within, p_cross, ncol=1,
           labels=c("a","b")) 

 ggsave(plot=p_all, "../figures/fig2/auc_tcrs.pdf",
        width=5, height=7)
 
 # barplot for positional weights
 weights <- data.frame(position=c(1:9),
            w=c(0.038,0.466,0.333,0.265,0.531,0.811,0.395,0.084,0.444))
 p_bar <- ggplot(weights, aes(x=position, y=w)) +
   geom_bar(stat="identity") +
   coord_flip() +
   theme_cowplot() +
   labs(x="Position",
        y="Weight") +
   scale_x_reverse(breaks=c(1:9))
 
 ggsave(plot=p_bar, "../figures/fig2/positional_weights.pdf",
        width=2, height=4)