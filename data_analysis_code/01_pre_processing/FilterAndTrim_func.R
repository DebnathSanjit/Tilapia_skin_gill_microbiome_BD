#!/usr/bin/env Rscript --vanilla

#mamba activate R

library(dada2, quietly = T)
set.seed(1234)

args <- commandArgs(trailingOnly = T)

input.fwd <- args[1]
input.rev <- args[2]
output.fwd <- args[3]
output.rev <- args[4]

try(filterAndTrim(fwd=input.fwd, filt=output.fwd, 
              rev=input.rev, filt.rev=output.rev,
              maxN=0, maxEE=c(2,2), truncQ=2, 
              rm.phix=T, compress=T, multithread=F, verbose = T))
