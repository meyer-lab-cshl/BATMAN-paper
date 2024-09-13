library(tidyr)
library(readr)
library(dplyr)
library(seqinr)

peptides <- read_csv("data/netmhcpan_inputs_mhc_i.csv")[, -1]
peptides_list <- split(peptides, f=peptides$mhc)

tt <- lapply(seq_along(peptides_list), function(x) {
  mhc <- names(peptides_list)[x]
  tmp <- peptides_list[[x]] %>%
    mutate(id=paste0(mhc, "_", 1:n()))
  write.fasta(sequences=as.list(tmp$peptide),
              names=tmp$id,
              file.out = file.path("data", paste0("peptides_", mhc, ".fasta")),
              as.string=TRUE)
})
