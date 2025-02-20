rm(list=ls())#Clear environment

# Packages
library(readxl)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(viridis)
library(forcats)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

###############################
## BATMAN with AA matrices#####
###############################
# List all files
corr_files <- list.files(path = "Corrs_TCR_AA_matrix/",
                           pattern = "corr_.*.csv",
                           full.names = TRUE)
# open empty df for storing
corr_df <- data.frame()
for (file in corr_files){
  # Read data
  file_data <- read.csv(file) 
  corr_df <- rbind(corr_df, file_data)
}

colnames(corr_df) <- c("X","corr",'aa','tcr')

# Make AA matrix names in all caps
corr_df$aa <- toupper(corr_df$aa)


###################################
######## Plotting #################
###################################

p <- ggplot(corr_df,aes(y=corr,x=fct_inorder((factor(aa))))) + 
  geom_boxplot(outlier.colour = NA,
               show.legend = FALSE) +
  geom_jitter(aes(color=tcr,group=aa),
              alpha=0.5, width=0.1, height=0,size=0.1, show.legend = FALSE) +
#  scale_colour_manual(values = palette(viridis(151))) +
  labs(y="Average Spearman r",
       x="",
       color="") +
  #  xlim(NA, 0.9) +
  theme_cowplot() +
  theme_minimal_grid(14) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=5),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank())
#  theme(legend.position = "bottom")+ 
#  guides(color = guide_legend(nrow = 9))

# Save plot
ggsave(plot=p, "corr_for_S5a.pdf",
       width=7, height=3)


# Save raw data
write.csv(corr_df, "raw_data_fig_S5a.csv")

