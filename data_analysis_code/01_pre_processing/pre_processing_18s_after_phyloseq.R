# Date 20240202
# Sanjit Debnath

library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")

# set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/Pre_processing_18S")

# Load phyloseq obect
ps.18s <- readRDS("Sanjit_18S_PS.rds")
ps.18s # 3052 taxa and 65 samples

#to write the samples names into csv file
seq.names <- sample_names(ps.18s) %>% data.frame()
#write.csv(seq.names, "sequence_names.csv", row.names = F)

# add metadata
metadata <- read.csv("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/Metadata_FSS_BD/FSS_BD_metadata_v5_18S.csv")
row.names(metadata) <- metadata[,"Seq_name"]
sample_data(ps.18s) <- metadata

# Check phyloseq object
ps.18s

# save phyloseq object with added metadata
#saveRDS(ps.18s, "phyloseq_with_metadata_20240202.rds")

# read phyloseq object
ps.18s <- readRDS("phyloseq_with_metadata_20240202.rds")
ps.18s

# now need to remove the contaminants and do further analysis
# Check contaminants----

##Inspect Library Sizes----
df <- as.data.frame(sample_data(ps.18s)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.18s)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Plot library size
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + 
  geom_point()

##Transpose the otu table
# as the taxa are in column in my previous otu table, need to transpose the otu_table using following code
# otherwise, will not work at get vector of contaminants stage
otu_table(ps.18s) <- t(otu_table(ps.18s, taxa_are_rows = "T"))

#do some filteering like on this https://f1000research.com/articles/5-1492
### prevalence filter----
#This function will remove taxa (ASVs) with low prevalence, where prevalence is the fraction of total samples in which an ASV is observed.
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps.18s), MARGIN = ifelse(taxa_are_rows(ps.18s), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps.18s), tax_table(ps.18s))

# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps1 <- prune_taxa(keepTaxa, ps.18s)
ps1 #2984 taxa and 65 samples

# Save file
#saveRDS(ps1, "ps_no_contam1.rds")

# Identify Contaminants - Prevalence
# for checking and removing contaminants, as prevalence seems works better, 
# I will only use prevalence filter (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
sample_data(ps1)$is.neg <- sample_data(ps1)$Sample_or_control == "Negative"
ps_water_contamdf.prev <- isContaminant(ps1, method="prevalence", neg="is.neg") #default prevalence 0.1
table(ps_water_contamdf.prev$contaminant)
#FALSE  TRUE 
#2982     2 

# check the contaminats
which(ps_water_contamdf.prev$contaminant)
# 387 547

# Increase the prevalence to 0.5
ps_water_contamdf.prev05 <- isContaminant(ps1, method="prevalence", neg="is.neg", threshold=0.5)
table(ps_water_contamdf.prev05$contaminant)
#FALSE  TRUE 
# 2977     7  

which(ps_water_contamdf.prev05$contaminant)
# 384  387  444  547 1954 2757 2893

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps1, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Negative", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=ps_water_contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Make data.frame of prevalence in positive and negative samples (threshold=0.5)
df.pa05 <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=ps_water_contamdf.prev05$contaminant)
ggplot(data=df.pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#get vector of contaminants
wat.contam <- ps1@otu_table[ps_water_contamdf.prev05$contaminant,] %>% row.names() 
wat.contam

#"ASV911"  "ASV2008" "ASV588"  "ASV25"   "ASV2821" "ASV1705" "ASV130" 

view(ps1@tax_table[wat.contam,])

# We can manually check these asv if they are actual contaminat or not

# If want to check individual likely contaminats
# Check possible Contaminants in controls (mock, negs) as then often appear in control generally with higher counts compared to samples.
# check consistency of possible Contaminants in controls and samples as contaminates may show varying counts in different samples, while genuine biological signals are often more consistent
# check the presence of possible contaminants in negs, as they may be introduced during the DNA extraction or PCR steps.

### But we can also just remove these asvs without checking, as they don't make much difference 
# Remove contaminants----
water.nocontam <- prune_taxa(!(taxa_names(ps1) %in% wat.contam), ps1)
water.nocontam # 2977 taxa and 65 samples

# Save this as RDS file
#saveRDS(water.nocontam, "phyloseq_18s_no_contam_20240309.rds")

# load latest phyloseq obect
ps <- readRDS("phyloseq_18s_no_contam_20240309.rds")
ps #2977 taxa and 65 samples

## subset only mock to check the accuracy later
ps_mock <- ps %>% 
  ps_filter(Sample_type == "mock")
ps_mock #5 taxa and 2 samples
saveRDS(ps_mock, "phyloseq_18s_mock_20240309.rds")

# as previously we have removed contaminants using negs, and also have separated 
#mock for checking accuracy later, we do not need them any more, so remove them
ps.no.control <- ps %>% 
  ps_filter(Sample_type == "Pond_water")
ps.no.control # 2977 taxa and 59 samples


## Remove anything Class == craniata which removes fish and humans----
### Removing Craniata----
ps_no_Craniata <- ps.no.control %>% subset_taxa(Class !="Craniata"  | is.na(Class))
ps_no_Craniata #2961 taxa and 59 samples

### Removing Teleostei----
ps_no_Teleostei <- ps_no_Craniata %>% subset_taxa(Family !="Teleostei"  | is.na(Family))
ps_no_Teleostei #2961 taxa and 59 samples

### Removing NA----
ps_no_NA <- ps_no_Teleostei %>% subset_taxa(Domain !='NA') 
ps_no_NA #2961 taxa and 59 samples

# check the count number per sample and add in the metadata
# Calculate the count number per sample
count_per_sample <- sample_sums(ps_no_NA)

# Add the count number to the sample metadata
sample_data(ps_no_NA)$count_per_sample <- count_per_sample

# Check the updated sample metadata
sample_data(ps_no_NA)

# Remove low read count samples----
ps_high_counts <- prune_samples(sample_sums(ps_no_NA)>=1000, ps_no_NA) %>% tidyPS() #remove low read count (less than 1000 reads) samples that have failed
ps_high_counts # 2961 taxa and 58 samples 

# Rename this and saved as all_filtered as everything has been filtered from this 
ps_all_filtered <- ps_high_counts
ps_all_filtered # 2961 taxa and 58 samples 

# Phylogenetic tree has already been added, so this is the final phyloseq object 
#saveRDS(ps_all_filtered, "phyloseq_18s_filtered_final_20240309.rds")


