# # Code to plot AA distances 

rm(list=ls())#Clear environment

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(viridis)
library(tidyr)

# Read Data
peptide_folds <- read.csv('aa_within.csv',row.names=1)

# AA lists
all_aa <- row.names(peptide_folds)
#hydrophobic_aa <- c("A","V","I","L","M","F","C","W")
hydrophobic_aa <- c("I","L","V")
non_hydrophobic_aa <- c("R","N","D","Q","E","K")
non_hydrophobic_aa <- setdiff(all_aa,hydrophobic_aa)

hydrophobic2non <- 
  data.frame(distance = as.vector(as.matrix(peptide_folds[hydrophobic_aa,non_hydrophobic_aa])),
             type = 'Hydrophobic to hydrophilic')

non2hydrophobic <- 
  data.frame(distance = as.vector(as.matrix(peptide_folds[non_hydrophobic_aa,hydrophobic_aa])),
             type = 'Hydrophilic to hydrophobic')

plot_data <- rbind(hydrophobic2non,non2hydrophobic)


# Plot
p <- ggplot(plot_data, 
aes(x=type, y=distance)) +
  geom_boxplot(aes(color=type),outlier.colour = NA,
               show.legend = FALSE, alpha=1,lwd=1
  ) +
  geom_point(alpha=0.6,
             size=5, position = position_jitter(width=0.1)) + 
                                         labs(x="",
                                              y="Distance",
                                              color="") +
#  ylim(NA, 0.89) +
  theme_cowplot() +
#  theme(legend.text = element_text(size=25)) +
#  guides(color = guide_legend(ncol = 1)) +
#  theme(legend.position="right")+
  theme(axis.text.x = element_text(size=20,angle=45,hjust=1),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25))

ggsave(plot=p, "AA_distances_hydro.pdf",
       width=8, height=10)
