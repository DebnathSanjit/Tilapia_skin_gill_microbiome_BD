# Set seed
set.seed(1234)

#mamba activate R
library(dada2, quietly = T)
library(phyloseq, quietly = T)

# Define args
args <- commandArgs(trailingOnly = T)
path <- args[1]
ps <- readRDS(args[2])
silva <- args[3]
silva_species <- args[4]

taxa <- assignTaxonomy(as(refseq(ps), "character"), silva, multithread=T, 
                       tryRC = T)
taxa <- addSpecies(taxa, silva_species, tryRC = T)

row.names(taxa) <- names(rownames(taxa))
tax_table(ps) <- taxa

saveRDS(ps, paste0(path, "/phyloseq_taxa.rds"))
