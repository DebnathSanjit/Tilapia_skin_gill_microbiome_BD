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

filtRs <- file.path(sort(list.files(path, pattern = ".*filt.rev.fq.gz", 
                                    full.names = T, recursive = T)))
filtRs <- filtRs[sapply(filtRs, file.size) > 20]


errR_1 <- try(learnErrors(filtRs, multithread = T, nbases = 1e10,
                      errorEstimationFunction = loessErrfun_mod1, 
                      verbose = T))
saveRDS(errR_1, paste0(path,"/errR_1.rds"))

errR_1.p <- plotErrors(errR_1, nominalQ = T)
ggsave(paste0(path, "/errR_1.svg"), plot = errR_1.p, 
              units = "cm", width = 25, height = 25)
