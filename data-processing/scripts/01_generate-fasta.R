library(tidyr)
library(readr)
library(dplyr)
library(seqinr)

#fn"data/netmhcpan_inputs_mhci_1dhamming.csv"
fn <- snakemake@input[['peptides']])
distance <- gsub(".*_(.*).csv", "\\1")
peptides <- read_csv(fn)[, -1]
peptides_list <- split(peptides, f=peptides$mhc)


tt <- lapply(seq_along(peptides_list), function(x) {
  mhc <- names(peptides_list)[x]
  tmp <- peptides_list[[x]] %>%
    mutate(id=paste0(mhc, "_", 1:n()))
  write.fasta(sequences=as.list(tmp$peptide),
              names=tmp$id,
              file.out = file.path("data",
                                   paste0("peptides_", mhc, ".fasta")
                                   ),
              as.string=TRUE)
})
