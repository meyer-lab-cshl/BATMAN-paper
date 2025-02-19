#################
## libraries ####
#################
library(tidyr)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)

############
## data  ###
############

# peptide contacts with MHC and TCR measured from pdb structures
contacts <- readxl::read_xlsx("data/TCR-peptide-contacts.xlsx",
                              sheet="tidy-summary")

# mapping file between TCR name and pdb id
mapping <- readxl::read_xlsx("data/PDB_ID_list.xlsx")

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
  tmp  <- read_csv(paste0("data/inferred_weights/inferred_weights_", tcr, ".csv"))
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

p <- ggplot(tcr_interactions)
p <- p + geom_point(aes(x=contacts, y=weights_norm,
                   color=tcr)) +
    facet_wrap(~tcr, ncol=4) +
  scale_color_brewer(type="qual", palette = "Dark2") +
  labs(color="TCR",
       y="Normalised weights",
       x="Number of interactions") +
  geom_smooth(aes(x=contacts, y=weights_norm, group=tcr, color=tcr),
              formula= y ~ x,
              method ="lm",
              se=FALSE) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20))

ggsave(plot=p, "../figures/extended_data_fig8/structure-weight-relationship.pdf",
       width=15, height=4)
