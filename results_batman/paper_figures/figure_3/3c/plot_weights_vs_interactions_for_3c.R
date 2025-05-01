#################
## libraries ####
#################
library(tidyr)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(ggrepel)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

############
## data  ###
############

# peptide contacts with MHC and TCR measured from pdb structures
contacts <- readxl::read_xlsx("TCR-peptide-contacts.xlsx",
                              sheet="tidy-summary")

# mapping file between TCR name and pdb id
mapping <- readxl::read_xlsx("PDB_ID_list.xlsx")

################
## analyses ####
################

contacts[is.na(contacts)] <- 0
contacts <- contacts %>%
    select(-total) %>%
    mutate(pdb_id=tolower(pdb_id)) %>%
    pivot_longer(MHC:Peptide, names_to="interaction",
                 values_to="contacts")

# get BATMAN inferred positional weights for TCRs with available structure
pepweights <- lapply(mapping$tcr[mapping$pdb_id %in% contacts$pdb_id], function(tcr) {
  tmp  <- read_csv(paste0("inferred_weights/inferred_weights_", tcr, ".csv"))
  tibble(tcr=as.vector(tmp[[1]]),
         position=as.numeric(names(tmp)[-1]) + 1,
         weights=unlist(tmp[-1]))
}) %>%
  bind_rows

# subset weights for positions in structure and normalise to max weight per TCR
pepweights <- pepweights  %>%
  filter(!(tcr == "TCR-F5" & position %in% c(1:6, 20)),
         !(tcr == "F24" & position %in% c(1:6, 20))) %>%
  mutate(position = case_when(tcr %in% c("TCR-F5", "F24") ~ position - 6,
                              TRUE ~ position)) %>%
  group_by(tcr) %>%
  mutate(weights_norm = weights/max(weights))


contacts <- mapping %>%
    select(tcr, pdb_id) %>%
    inner_join(contacts) %>%
    left_join(pepweights) %>%
    select(-pdb_id)

# subset to interactions (H-bonds and others) with TCR
tcr_interactions <- contacts %>%
  filter(measure=="all", interaction == "TCR")

# Add residue AAs
tcr_epitope <- mapping[c('tcr','peptide')]
tcr_epitope$peptide[tcr_epitope$peptide == 'FRDYVDRFYKTLRAEQASQE'] <- 'RFYKTLRAEQASQ'
rownames(tcr_epitope) <- tcr_epitope$tcr

peptides <- tcr_epitope[tcr_interactions$tcr,]
tcr_interactions$pepresidue <- substr(peptides$peptide,tcr_interactions$position,
                                      tcr_interactions$position)
tcr_interactions$pepresidue <- paste0(tcr_interactions$position,
                                      tcr_interactions$pepresidue)

# Save raw data
write.csv(tcr_interactions, "raw_data_fig_3c.csv")


# Show label only for 1E6
tcr_interactions$pepresidue[tcr_interactions$tcr!='TCR-1E6'] <- NA

# Show label only for certain residues of 1E6
tcr_interactions$pepresidue[!tcr_interactions$pepresidue %in% c('4G','5P',
                                                                '6D','7P')] <- NA

#### Plotting ########

p <- ggplot(tcr_interactions,aes(x=contacts, y=weights_norm,label=pepresidue))
p <- p + geom_point(aes(x=contacts, y=weights_norm,
                        color=factor(tcr,levels=c("TCR-1E6", "47BE7","TIL1383I",
                                                  "NYE-S1", "NYE-S2", "NYE-S3",
                                                  "TCR-F5", "F24"))),
                    size=4,alpha=0.95,show.legend = FALSE) +
  scale_color_manual(values = c("#FF7F00","#002244","#B23648","#3A85A8","#00BFFF",
                                         "#007FFF","#F2E631","#5A4FCF")) +
  labs(color="",
       y="Normalized positional weight",
       x="Number of TCR interactions") +
  geom_smooth(aes(x=contacts, y=weights_norm,
                  group=factor(tcr,levels=c("TCR-1E6", "47BE7","TIL1383I","NYE-S1",
                                            "NYE-S2", "NYE-S3","TCR-F5", "F24")),
                  color=factor(tcr,levels=c("TCR-1E6", "47BE7","TIL1383I",
                                            "NYE-S1", "NYE-S2", "NYE-S3",
                                            "TCR-F5", "F24"))),
              formula= y ~ x,
              method ="lm",
              se=FALSE,show.legend = FALSE,size=0.5) +
  cowplot::theme_cowplot() +
  facet_wrap(vars(factor(tcr,levels=c("TCR-1E6", "47BE7","TIL1383I",
                                      "NYE-S1", "NYE-S2", "NYE-S3",
                                      "TCR-F5", "F24"))),ncol=2) +
  theme(strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")))+
  geom_text_repel(size = 7) +
  theme(text = element_text(size=35),
        axis.text = element_text(size=25)) +
  theme(legend.position = "bottom")


ggsave(plot=p, "structure_weight_relationship_3c.pdf",
       width=6, height=11)

