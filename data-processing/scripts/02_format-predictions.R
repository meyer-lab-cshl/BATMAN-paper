#################
## libraries ####
#################
library(tidyverse)
library(stringr)
library(optparse)

############
## data  ###
############

## command line arguments ####
option_list <- list(
    make_option(c("-p", "--dirpred"), action="store",
               dest="dirpred",
               type="character", help="Path the directory with netMHCpan results
                [default: %default].", default=NULL),
    make_option(c("-i", "--sampleid"), action="store",
               dest="id",
               type="character", help="sample id as indicated in filename
                [default: %default].", default=NULL),
    make_option(c("-c", "--class"), action="store",
                dest="class",
                type="character", help="MHC class
                [default: %default].", default=NULL),
    make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out ",
               "[default: %default]."),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script",
                "[default: %default].")
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$dirpred <- "data-processing/data/predictions/mhcii-1dhamming"
    args$id <- "H-2-IAb_13"
    args$verbose <- TRUE
    args$class <- "classII"
}

################
## analysis ####
################

format_predictions <- function(args) {
    if (args$verbose){
        message("Read binding predictions for ", args$id)
    }
    filepath <- file.path(args$dirpred,
                          str_c("netMHCpan_", args$id, ".txt"))
    dat <- read_lines(filepath, skip_empty_rows=TRUE)

    if (args$verbose) message("Extract information")
    ## Remove comment and separator lines ####
    dat <- dat[!grepl("^#", dat)]
    dat <- dat[!grepl("^-", dat)]

    ## Extract info about training data ####
    dist_training <- dat[grepl("Distance to training data", dat)]
    dat <- dat[!grepl("Distance to training data", dat)]

    ## Check for error messages ####
    error <- dat[grepl("Error", dat)]
    dat <- dat[!grepl("Error", dat)]

    ## Extract predictions summary lines ####
    index_proteins <- which(grepl("Number of weak binders", dat))
    proteins <- dat[index_proteins]
    dat <- dat[-index_proteins]

    if (args$verbose) message("Format data")
    ## Format vector of entries into tibble ####
    header <- dat[1]
    dat <- dat[!grepl("Pos", dat)]

    columncount <- ifelse(args$class == "classI", 17, 15)
    header <- header %>%
        str_squish() %>%
        str_split_fixed(pattern=" ", n=columncount)
    dat <- dat %>%
        str_squish() %>%
        str_split_fixed(pattern=" ", n=columncount)
    colnames(dat) <- header

    # Rank_EL and Score_EL names new in netMHCpan4.1; this formating won't work
    # for earlier netMHCpan versions
    # aff column abbreviation different for netMHCpan I and II: make consistent
    aff_column <- header[grepl("Aff", header)]
    
    dat <- dat[1:nrow(dat),] %>%
        as_tibble %>%
        dplyr::rename(Rank_EL=`%Rank_EL`) %>%
        dplyr::rename(Rank_BA=`%Rank_BA`) %>%
        dplyr::rename(Aff_nM={aff_column}) %>%
        mutate(BindLevel = gsub("<= ", "", BindLevel),
               BindLevel = case_when(BindLevel =="" ~ "NB",
                                     TRUE ~ BindLevel))

}

dat <- format_predictions(args)

if (args$verbose) message("Write results")
write_csv(dat,
          file.path(args$dirpred,
                    str_c("netMHCpan_", args$id, "_formated.csv")),
          na="")


