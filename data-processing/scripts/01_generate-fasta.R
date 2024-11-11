library(tidyr)
library(readr)
library(dplyr)
library(seqinr)

#fn"data/netmhcpan_inputs_mhci_1dhamming.csv"
fn <- snakemake@input[['peptides']]
distance <- gsub(".*_(.*).csv", "\\1", fn)
peptides <- read_csv(fn)

fn.out <- file.path("data", "fastafiles", distance)
if (!dir.exists(fn.out)) dir.create(fn.out, recursive=TRUE)

# format MHC names for netMHCpan
peptides <- peptides %>%
  mutate(mhc = case_when(mhc == "H2-Db" ~ "H-2-Db",
                         mhc == "H2-Kb" ~ "H-2-Kb",
                         mhc ==  "I-Ab" ~ "H-2-IAb",
                         mhc ==  "I-Au" ~"H-2-IAu",
        # remove * and : : cluster submission can't deal with these in filenames
                         grepl("DR", mhc) ~ gsub("\\*","_", mhc),
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
              file.out = file.path(fn.out, paste0("peptides_", combined, ".fasta")),
              as.string=TRUE)
})
