#!/usr/bin/env Rscript --vanilla
#mamba activate R

set.seed(1234)

library(dada2, quietly = T)
library(dplyr, verbose = F)
library(phyloseq)
library(Biostrings)
library(stringr)

# We have reads in multiple orientations that need to be reverse complimented 
# Before being combined. See issue below
# https://github.com/benjjneb/dada2/issues/938

# Define args
args <- commandArgs(trailingOnly = T)
path <- args[1]
errF <- readRDS(args[2]) # pass in which error model to use
errR <- readRDS(args[3])
primer.combos <- read.csv(args[4], col.names = c("fwd", "rev"))
lookup_table <- setNames(primer.combos$fwd, primer.combos$rev)

# This needs to be redone as not all samples pass filterAndTrim
# Define filtFs and filtRs
filtFs <- file.path(sort(list.files(path, pattern = ".*filt.fwd.fq.gz", 
                                    full.names = T)))
filtRs <- file.path(sort(list.files(path, pattern = ".*filt.rev.fq.gz", 
                                    full.names = T)))

filtFs <- filtFs[sapply(filtFs, file.size) > 20]
filtRs <- filtRs[sapply(filtRs, file.size) > 20]

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose=T)
derepRs <- derepFastq(filtRs, verbose=T)

dadaFs <- dada(derepFs, err=errF, multithread=T, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=T, pool="pseudo")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)

seqtab <- makeSequenceTable(mergers)
#print(table(nchar(colnames(seqtab))))

seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(252,254)]

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=T)

reorient <- seqtab.nochim %>%
  as.data.frame() %>%
  filter(grepl(".*806rcbc", rownames(.))) %>%
  `colnames<-`(reverseComplement(DNAStringSet(colnames(.)))) %>%
  `rownames<-`(str_replace_all(rownames(.), lookup_table)) %>%
  as.matrix()

without.rev.first <- seqtab.nochim %>%
  as.data.frame() %>%
  filter(!grepl(".*806rcbc", rownames(.))) %>%
  as.matrix()

combined.seqtab <- mergeSequenceTables(reorient, without.rev.first, 
                                       repeats = "sum", tryRC = T)

ps <- otu_table(combined.seqtab, taxa_are_rows = F)

dna <- DNAStringSet(taxa_names(ps))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(refseq(dna), ps)

saveRDS(ps, paste0(path, "/phyloseq.rds"))
writeXStringSet(refseq(ps), paste0(path, "/ASVs.fna"))
