rm(list=ls())#Clear environment

# Packages
library(readxl)
library(tidyverse)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)


# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# Get data
nlv_data <- read.csv('nlv_activity_distance.csv')

# Pivot longer
nlv_data <- pivot_longer(nlv_data,
                          cols=c("pan.TCR","NLV2","NLV3","TCR2","TCR52.10","TCR82.14"),
                          names_to = "tcr")

# Select TCRs to plot
nlv_data <- nlv_data[nlv_data$tcr %in% c("NLV3","TCR2","TCR52.10",
                                         "pan.TCR"),]

###################################
######## Plotting: MHCI ###########
###################################

p <- ggplot(nlv_data, aes(y=peptide_activity,x=value)) + 
  geom_smooth(data = nlv_data, aes(y=peptide_activity,x=value),
              colour="gray",
              method = "lm") +
  geom_point(color='#3A85A8',alpha=0.7,show.legend = FALSE,
             size=1) +
  labs(y="Peptide enrichment against NLV-expanded T cells",
       x="Peptide-to-index distance",
       color="") +
  theme_cowplot() +
  facet_wrap(~factor(tcr,levels=c("NLV3","TCR2","TCR52.10",
                                  "pan.TCR")),
             nrow=1,scales = "free_x") +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm"))) +
    theme(axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12)) +
    ylim(0,1.2)
#  theme(legend.position="right")
#  guides(color = guide_legend(ncol = 2))

# Save plot
ggsave(plot=p, "activity_vs_distance_for_5c.pdf",
       width=5, height=3)


# Save raw data
write.csv(nlv_data, "raw_data_fig_5c.csv")