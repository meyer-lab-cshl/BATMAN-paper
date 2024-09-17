library(tidyr)
library(readr)
library(dplyr)
library(seqinr)

peptides <- read_csv("data/netmhcpan_inputs_mhc_i.csv")[, -1]

# format MHC names for netMHCpan
peptides <- peptides %>%
  mutate(mhc = case_when(mhc == "H2-Db" ~ "H-2-Db",
                         mhc == "H2-Kb" ~ "H-2-Kb", 
        # remove * and : : cluster submission can't deal with these in filenames
                         TRUE ~ gsub("[\\*:]","", mhc)),
         length=nchar(peptide),
         combined = paste0(mhc, "_", length))

peptides_list <- split(peptides, f=peptides$combined)


tt <- lapply(seq_along(peptides_list), function(x) {
  combined <- names(peptides_list)[x]
  tmp <- peptides_list[[x]] %>%
    mutate(id=paste0(combined, "_", 1:n()))
  write.fasta(sequences=as.list(tmp$peptide),
              names=tmp$id,
              file.out = file.path("data", paste0("peptides_", combined, ".fasta")),
              as.string=TRUE)
})
