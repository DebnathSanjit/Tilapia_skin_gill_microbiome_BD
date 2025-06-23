#!/usr/bin/env Rscript --vanilla
#mamba activate R

library(dada2, quietly = T)
library(ggplot2)
library(svglite)
library(dplyr) # Required or error not enough reads appears
source("scripts/novaseq_error_correction_models.R")
set.seed(1234)

args <- commandArgs(trailingOnly = T)
path <- args[1]

filtFs <- file.path(sort(list.files(path, pattern = ".*filt.fwd.fq.gz", 
                                    full.names = T, recursive = T)))
filtFs <- filtFs[sapply(filtFs, file.size) > 20]


errF_3 <- try(learnErrors(filtFs, multithread = T, nbases = 1e10,
                      errorEstimationFunction = loessErrfun_mod3, 
                      verbose = T))
saveRDS(errF_3, paste0(path,"/errF_3.rds"))

errF_3.p <- plotErrors(errF_3, nominalQ = T)
ggsave(paste0(path, "/errF_3.svg"), plot = errF_3.p, 
              units = "cm", width = 25, height = 25)
