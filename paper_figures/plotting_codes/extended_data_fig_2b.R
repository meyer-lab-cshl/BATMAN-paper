# Plots performance of difference AA distance functions over TCRs grouped by index

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
library(readxl)
library(viridis)
library(ggnewscale)



## Extended Data Figure 2b-o: AUC (average of 5 folds) for TCRs
#############################
## Gather and store data ####
#############################
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is
selected_AA_distance <- c(1,16,17,67,68,69) #Plot only some AA distance functions

# Arrays to store 3 types of AUCs 
# (initialize with NaN so that missing data is not plotted)
all_AUC_non_weak <- array(NaN,dim=c(length(selected_AA_distance)+1,66))
all_AUC_non_strong <- array(NaN,dim=c(length(selected_AA_distance)+1,66))
all_AUC_weak_strong <- array(NaN,dim=c(length(selected_AA_distance)+1,66))

# Indices of 2-class TCRs
TCR_index_2_class <- c(10,11,22,24,26:32)
class_pair_2_class <- c(12,12,12,12,23,12,23,23,23,23,12)#1:non, 2: weak, 3: strong

#3-class TCRs only
for (tcr in setdiff(1:66,TCR_index_2_class)){ 
  # results for Bayesian inference
  for (AA_distance in 1:length(selected_AA_distance)) {
    load(paste("../data/within_tcr_unpooled_outputs//output_",
               selected_AA_distance[AA_distance],"_",tcr,".rda",sep=""))
    all_AUC_non_weak[AA_distance,tcr]<-mean(classification_results$AUCs[1,])
    all_AUC_non_strong[AA_distance,tcr]<-mean(classification_results$AUCs[2,])
    all_AUC_weak_strong[AA_distance,tcr]<-mean(classification_results$AUCs[3,])
  }
  # results for Atchley
  load(paste("../data/within_tcr_unpooled_outputs_atchley//output_",
             tcr,".rda",sep=""))
  all_AUC_non_weak[length(selected_AA_distance)+1,tcr]<-mean(AUCs[1,])
  all_AUC_non_strong[length(selected_AA_distance)+1,tcr]<-mean(AUCs[2,])
  all_AUC_weak_strong[length(selected_AA_distance)+1,tcr]<-mean(AUCs[3,])
}

#2-class TCRs only
for (i in 1:length(TCR_index_2_class)){ 
  # results for Bayesian inference
  for (AA_distance in 1:length(selected_AA_distance)) {
    load(paste("../data/within_tcr_unpooled_outputs//output_",
               selected_AA_distance[AA_distance],"_",TCR_index_2_class[i],
               ".rda",sep=""))
    
    if (class_pair_2_class[i]==12){
      all_AUC_non_weak[AA_distance,TCR_index_2_class[i]]<-
        mean(classification_results$AUCs)}
    
    if (class_pair_2_class[i]==13){
      all_AUC_non_strong[AA_distance,TCR_index_2_class[i]]<-
        mean(classification_results$AUCs)}
    
    if (class_pair_2_class[i]==23){
      all_AUC_weak_strong[AA_distance,TCR_index_2_class[i]]<-
        mean(classification_results$AUCs)}
    
    }
  # results for Atchley
  load(paste("../data/within_tcr_unpooled_outputs_atchley//output_",
             TCR_index_2_class[i],".rda",sep=""))
  
  if (class_pair_2_class[i]==12){
  all_AUC_non_weak[length(selected_AA_distance)+1,TCR_index_2_class[i]]<-mean(AUCs)
  }
  if (class_pair_2_class[i]==13){
  all_AUC_non_strong[length(selected_AA_distance)+1,TCR_index_2_class[i]]<-mean(AUCs)
  }
  if (class_pair_2_class[i]==23){
  all_AUC_weak_strong[length(selected_AA_distance)+1,TCR_index_2_class[i]]<-mean(AUCs)
  }
}

# Names of TCRs and AA distance matrices
# Load TCR and peptide information from the main TCR datafile
suppressWarnings(
  TCR_data <- read_excel("../data/TCR_epitope_database.xlsx")
)
# Extract TCR names present in the data
TCRs <- unique(TCR_data$tcr_name)

AA_distances <-
c("Hamming",paste(rep("BLOSUM",13),seq(30,90,5),sep=""),
                        "BLOSUM62","BLOSUM100",
                        paste(rep("PAM",50),seq(10,500,10),sep=""),
                        "Dayhoff","Gonnet","Atchley","Atchley_cos","pTEAM")
AA_distances <- AA_distances[c(1,16,17,67,68,69,71)] #add Atchley
# Add TCR and AA matrix names to array
colnames(all_AUC_non_weak) <- TCRs
colnames(all_AUC_non_strong) <- TCRs
colnames(all_AUC_weak_strong) <- TCRs
all_AUC_non_weak<- data.frame(AA_distances,all_AUC_non_weak)
all_AUC_non_strong<- data.frame(AA_distances,all_AUC_non_strong)
all_AUC_weak_strong<- data.frame(AA_distances,all_AUC_weak_strong)
###############
## Analysis ###
###############
# Convert to long form for box plot
all_AUC_non_weak <- all_AUC_non_weak %>%
  pivot_longer(-AA_distances, names_to="tcr", values_to = "score") %>%
  mutate(AA_distances=fct_inorder(AA_distances))

all_AUC_non_strong <- all_AUC_non_strong %>%
  pivot_longer(-AA_distances, names_to="tcr", values_to = "score") %>%
  mutate(AA_distances=fct_inorder(AA_distances))

all_AUC_weak_strong <- all_AUC_weak_strong %>%
  pivot_longer(-AA_distances, names_to="tcr", values_to = "score") %>%
  mutate(AA_distances=fct_inorder(AA_distances))

# Merge data and indicate class pair
all_AUC <- rbind.data.frame(all_AUC_non_weak,all_AUC_non_strong,all_AUC_weak_strong)
all_AUC$class_pair <- c(rep("12", each = dim(all_AUC_non_weak)[1]),
                        rep("13", each = dim(all_AUC_non_strong)[1]),
                        rep("23", each = dim(all_AUC_weak_strong)[1]))

# Indicate index peptides for each tcr in data point
# Array to store TCR names and index peptides
index_peptide <- array("X",dim=c(dim(all_AUC)[1]))
for (i in 1:dim(all_AUC)[1]){
  TCR_index <- match(all_AUC$tcr[i],unique(all_AUC$tcr)) #this is because R changed TCR names
  index_peptide[i] <- unique(TCR_data$index_peptide[TCR_data$tcr_name==TCRs[TCR_index]])
}

all_AUC$index_peptide <- index_peptide

# For SIINFEKL TCRs, classify them into 2 categories
tcr_type <- array(NaN,dim=c(dim(all_AUC)[1]))

naive_tcrs <- c("B11","B15","B3", "F4","E8", "B13","H6", "G6", "F5", 
                "H5", "B2", "B6", "B5","E9", "E4", "G2", "B16","B14")
ed_tcrs <- c("OT1","Ed5","Ed8","Ed9","Ed10","Ed16.1","Ed16.30","Ed21","Ed23",
             "Ed28","Ed31","Ed33","Ed39","Ed40","Ed45","Ed46")

tcr_type[all_AUC$tcr %in% naive_tcrs] <- "Naive"
tcr_type[all_AUC$tcr %in% ed_tcrs] <- "Educated"

all_AUC$tcr_type <- tcr_type




###############
##  plot ######
###############

# Special plotting specs for SIINFEKL data
all_AUC <- all_AUC %>%
  mutate(color_by = tcr,
         color_by = case_when(index_peptide == "SIINFEKL" ~ tcr_type,
                              TRUE ~ color_by))

auc_by_index <- split(all_AUC, f=all_AUC$index_peptide)

boxplot_by_index <- lapply(auc_by_index, function(x) {
  
  set_alpha <- 0.7 
  set_size <- 2.25
  set_jitter_width <- 0.15
  
  if (unique(x$index_peptide)=="SIINFEKL"){
    set_alpha <- 0.85
    set_size <- 1.6
    set_jitter_width <- 0.2
  }
  
  p <- ggplot(x, aes(x=AA_distances, y=score)) +
    geom_boxplot(aes(color=factor(class_pair,levels = unique(class_pair))),
                 outlier.colour = NA,
                 show.legend = FALSE, alpha=1,
                 position=position_dodge(0.75)) +
    scale_color_manual(values=c('#50C878','#DC143C','#007FFF')) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color=as.factor(color_by), 
                   group=factor(class_pair,levels = unique(class_pair))), 
               alpha=set_alpha,
               size=set_size,
               position = position_jitterdodge(jitter.width=set_jitter_width,
                                               dodge.width = 0.75)) +
    scale_color_brewer(type="qual", palette = "Dark2") +
    labs(y="AUC",
         x="",
         color="") +
    ylim(min(all_AUC$score,na.rm = TRUE), max(all_AUC$score,na.rm = TRUE)) +
    theme_cowplot() +
    # theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
    #       axis.text.y = element_text(size=20),
    #       axis.title.y = element_text(size=20)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=20),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=20))+
    facet_wrap(vars(index_peptide))+
    theme(strip.text = element_text(
      size = 25, color = "black")) +
    theme(legend.text = element_text(size=15)) +
    guides(color = guide_legend(ncol = 3)) +
    theme(legend.position="bottom")
   
})

p_all <- cowplot::plot_grid(plotlist = boxplot_by_index, ncol=5,nrow=3,
                            align = "hv",axis="tb")

ggsave(plot=p_all, "../figures/extended_data_fig2/auc_by_index.pdf",
       width=20, height=17)





