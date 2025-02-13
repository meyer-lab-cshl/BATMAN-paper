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
contacts <- readxl::read_xlsx("../data/TCR-peptide-contacts.xlsx",
                              sheet="tidy-summary")

# mapping file between TCR name and pdb id
mapping <- readxl::read_xlsx("../data/PDB_ID_list.xlsx")

################
## analyses ####
################

# summarise interactions across chains for TCR and MHC
contacts <- contacts %>%
  mutate(primary_interaction = gsub("-.*", "", primary_interaction)) %>%
  group_by(position, primary_interaction, measure, pdb_id) %>%
  summarise(interactions=sum(interactions, na.rm=TRUE))

# get BATMAN inferred positional weights for TCRs with available structure
pepweights <- lapply(mapping$tcr[mapping$pdb_id %in% contacts$pdb_id], function(tcr) {
  tmp  <- read_csv(paste0("../data/inferred_weights/inferred_weights_", tcr, ".csv"))
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

# create full grid of possible interactions and merge all dataframes
allmhci <- expand.grid(c(1:9),
                        c("MHC", "TCR"),
                        c("H-bonds", "all"),
                        mapping$tcr[!grepl("DRB", mapping$mhc)])
allmhcii <- expand.grid(c(1:13),
                       c("MHC", "TCR"),
                       c("H-bonds", "all"),
                       mapping$tcr[grepl("DRB", mapping$mhc)])
allcombs <- bind_rows(allmhci, allmhcii) %>%
  as_tibble
colnames(allcombs) <- c("position", "primary_interaction", "measure", "tcr")


contacts <- mapping %>%
  select(tcr, pdb_id) %>%
  inner_join(contacts) %>%
  select(-pdb_id)

contacts <- allcombs %>%
  left_join(contacts) %>%
  left_join(pepweights) %>%
  mutate(interactions = case_when(is.na(interactions) ~ 0,
                                  TRUE ~ interactions))

# subset to interactions (H-bonds and others) with TCR 
tcr_interactions <- contacts %>%
  filter(measure=="all", primary_interaction == "TCR")

p <- ggplot(tcr_interactions)
p <- p + geom_point(aes(x=interactions, y=weights_norm,
                   color=tcr)) +
  scale_color_brewer(type="qual", palette = "Dark2") +
  labs(color="TCR",
       y="Normalised weights",
       x="Number of interactions") +
  geom_smooth(aes(x=interactions, y=weights_norm, group=tcr, color=tcr),
              formula= y ~ x,
              method ="lm",
              se=FALSE) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20))
           
ggsave(plot=p, "../figures/extended_data_fig8/structure-weight-relationship.pdf",
       width=15, height=4)
