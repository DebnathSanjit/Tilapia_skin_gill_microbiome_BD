# date: 27/11/2023
# Author: Sanjit Debnath

#Load libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")
library(decontam); packageVersion("decontam")
library(microViz); packageVersion("microViz")

#set working dictionary
setwd("C:/Users/Scd_Research_Data/Data/Field_study_samples_BD/R_scripts")

#small function to tidy ps object after subsetting to remove any zero taxa or samples, 
#optionally transforms counts to relative abundance (JMcM)
tidyPS <- function(ps, RA = FALSE){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  if(RA == TRUE){
    ps <- transform_sample_counts(ps, function(x){x / sum(x)})
  }
  return(ps)}

#set theme
theme_set(theme_bw())

#Load phyloseq object
ps <- readRDS("C:/Users/Scd_Research_Data/Data/Field_study_samples_BD/R_scripts/phyloseq_taxa.rds")
ps # 127266 taxa and 384 samples

# Write the samples names into csv file
#metadata <- sample_names(ps) %>% data.frame()
#write.csv(metadata, "metadata_old.csv", row.names = F)

# Add metadata/sample data to the phyloseq object
metadata <- read.csv("FSS_BD_metadata_v5.csv")
#metadata.seq <- read.csv("")

#combined_df <- metadata.new %>% #it is the whole metadata with sample names and other details
#  left_join(., metadata.seq, by = "Plate_ID") #it is the name with the sequence, need to import it first

row.names(metadata) <- metadata[,"Seq_name"]
sample_data(ps) <- metadata


#rownames(metadata) <- metadata$Seq_name 
#sample_data(ps) <- metadata

#add reads number into the metadata
sample_data(ps)$Read_depth <- sample_sums(ps)
ps # 127266 taxa and 352 samples 

#saveRDS(ps, "C:/Users/Scd_Research_Data/Data/Field_study_samples_BD/R_scripts/phyloseq_taxa2.rds")

# 1. Remove Eukaryotes, chlorophyl, mitocondria----
# Check the tax_table for any Eukaryotes, chloroplast, mitochondria
view(ps@tax_table)

##removing eukaryotes
ps_no_euk <- ps %>% subset_taxa(Kingdom !="Eukaryota"  | is.na(Kingdom))
ps_no_euk # 127266 taxa and 352 samples

##removing chloroplast
ps_no_chl <- ps_no_euk %>% subset_taxa(Order !="Chloroplast"  | is.na(Order))
ps_no_chl #  125519 taxa and 352 samples

##removing mitochondria
ps_no_mit <- ps_no_chl %>% subset_taxa(Family !='Mitochondria'  | is.na(Family)) 
ps_no_mit # 124507 taxa and 352 samples

#removing NA
ps_no_NA <- ps_no_mit %>% subset_taxa(Kingdom !='NA') 
ps_no_NA #124490 taxa and 352 samples 

# 2. prevalence filtering----
ps0 <- readRDS("ps_no_NA.rds")
ps0 #124490 taxa and 352 samples

##Inspect Library Sizes----
df <- as.data.frame(sample_data(ps0)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps0)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

##Transpose the otu table----
# as the taxa are in column in my previous otu table, need to transpose the otu_table using following code
otu_table(ps0) <- t(otu_table(ps0, taxa_are_rows = "T"))

#do some filteering like on this https://f1000research.com/articles/5-1492
### prevalence filter----
#This function will remove taxa (ASVs) with low prevalence, where prevalence is the fraction of total samples in which an ASV is observed.
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps0), MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps0), tax_table(ps0))
# Define prevalence threshold as present ASV in 2 samples or more
prevalenceThreshold <- 2
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps1 <- prune_taxa(keepTaxa, ps0)
ps1 #51153 taxa and 352 samples

# Save file
#saveRDS(ps1, "ps_no_contam1.rds") 

# 3. Phylogenetic tree----
#function to root the tree
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

# Load tree
tree <- read_tree ("C:/Users/Scd_Research_Data/Data/Field_study_samples_BD/R_scripts/ASVs.msa.treefile")

out_group <- pick_new_outgroup(tree)
rooted_tree <- ape::root(tree, outgroup=out_group, resolve.root=TRUE)

phy_tree(ps) <- rooted_tree
#saveRDS(ps, "ps_filtered_with_tree.rds")

# 4. Removing contaminants----
# date: 16/10/2024 (updated from 27/11/2023)

#set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/Pre_processing_16S")

# Load phyloseq object with tree
ps.16s <- readRDS("ps_filtered_with_tree.rds")
ps.16s # 51153 taxa and 352 samples 

# There were some problem with the previous metadata, some information was missing
# so I updated metadata. Later I have updated the metadata again to add stocking density column, also changed 
## "Non.diseased" into "Non-diseased". Also changed PCR_batch to PCR_library. As i am redoing contaminant step again, it's better to add the metadata here
# Load the updated metadata/sample data to the phyloseq object
# load the new metadata
metadata.16s <- read.csv("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/Metadata_FSS_BD/FSS_BD_metadata_v7_16S.csv")
# add the new metadata into the phyloseq object
row.names(metadata.16s) <- metadata.16s[,"Seq_name"]
sample_data(ps.16s) <- metadata.16s

# check the new ps
ps.16s #51153 taxa and 352 samples 

## Rename the column according to sample name----
# Extract the sample data
sample_data_df <- as.data.frame(sample_data(ps.16s))

# Set the row names to be the values in the Sample column
row.names(sample_data_df) <- sample_data_df$Sample  # Use the Sample column for row names

# Update sample names in the phyloseq object to match sample_data_df
sample_names(ps.16s) <- row.names(sample_data_df)

# Update the sample data in the phyloseq object
sample_data(ps.16s) <- sample_data(sample_data_df)

#Identify Contaminants - Prevalence----
# I checked using frequency method, but that's not very helpful. So I will only use prevalence

#check the contamination using blanks and remove them (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
sample_data(ps.16s)$is.neg <- sample_data(ps.16s)$Sample_or_control == "Negative"    #the $ operator is used to extract or subset a specific part of a data object in R
contamdf.prev <- isContaminant(ps.16s, method="prevalence", neg="is.neg") # default threshold 0.1
table(contamdf.prev$contaminant)
#FALSE  TRUE 
#51126    27
head(which(contamdf.prev$contaminant)) # check the contaminates
#  478  725 3025 5134 5584 7371

# check by increasing threshold=0.5
contamdf.prev05 <- isContaminant(ps.16s, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
#FALSE  TRUE 
#51019   134 
head(which(contamdf.prev05$contaminant)) #to check which ASVs are contaminated
# 27 317 322 342 392 478

# Make phyloseq object of presence-absence in negative controls and true samples
ps.16s.pa <- transform_sample_counts(ps.16s, function(abund) 1*(abund>0))
ps.16s.pa.neg <- prune_samples(sample_data(ps.16s.pa)$Sample_or_control == "Negative", ps.16s.pa)
ps.16s.pa.pos <- prune_samples(sample_data(ps.16s.pa)$Sample_or_control == "Sample", ps.16s.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.16s.pa.pos), pa.neg=taxa_sums(ps.16s.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Make data.frame of prevalence in positive and negative samples (threshold=0.5)
df.pa05 <- data.frame(pa.pos=taxa_sums(ps.16s.pa.pos), pa.neg=taxa_sums(ps.16s.pa.neg),
                      contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa05, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove contaminants from final phyloseq----
# Default prevalence
contam_asvs <- ps.16s@otu_table[contamdf.prev$contaminant,] %>% row.names() #get vector of contaminants
# Previously in this step an error appeard because of the otu_table not correct,
# So I have transpose the otu_table and then it worked.
contam_asvs
#[1] "ASV13737"  "ASV50732"  "ASV19142"  "ASV101438" "ASV18942"  "ASV60716"  "ASV29737" 
#[8] "ASV12947"  "ASV91995"  "ASV58559"  "ASV38109"  "ASV76979"  "ASV59411"  "ASV36357" 
#[15] "ASV39620"  "ASV65884"  "ASV46366"  "ASV5149"   "ASV39123"  "ASV41913"  "ASV25158" 
#[22] "ASV64643"  "ASV45529"  "ASV76787"  "ASV92426"  "ASV30449"  "ASV31453" 

#ps.16s@otu_table[1:5,1:5]
#contamdf.prev$contaminant %>% head(100)#to view which ASVs are possible contaminants (TRUE)

#Keep only possible contaminants in the phyloseq object
prune_taxa(contam_asvs, ps.16s) #prune_taxa function keeps the taxa from phyloseq object that we are interested in

#Checl taxonomy of these possilbe contaminants
#ps.16s@tax_table[contam_asvs,] #check taxonomy, expect or dodgy. Ralstonia always comes up - bad
view(ps.16s@tax_table[contam_asvs,])
view(t(ps.16s@otu_table[contam_asvs,])) # check the otu table, contaminant: ASV consistently shows high abundance in negative controls, low abundance in true samples, and is absent or minimal in the mock sample

# Threshold = 0.5 
#0.5 threshold does more strong clean-up than default (0.1) 
#get vector of contaminants
contam_asvs_05 <- ps.16s@otu_table[contamdf.prev05$contaminant,] %>% row.names() #get vector of contaminants
#contamdf.prev05$contaminant %>% head(100)#to view which ASVs are possible contaminants (TRUE)
contam_asvs_05 # 134 ASVs
view(ps.16s@tax_table[contam_asvs_05,])

# Find out which ASVs have higher reads in negative samples compared to true samples
# Extract contaminant ASVs
contam_asvs05 <- ps.16s@otu_table[contamdf.prev05$contaminant, ] 

# Convert to a data frame for easier manipulation
contam_df05 <- as.data.frame(contam_asvs05)

# Extract sample names for negative and true samples
neg_samples <- sample_names(ps.16s)[sample_data(ps.16s)$Sample_or_control == "Negative"]
true_samples <- sample_names(ps.16s)[sample_data(ps.16s)$Sample_or_control == "Sample"]
mock_samples <- sample_names(ps.16s)[sample_data(ps.16s)$Sample_or_control == "Positive"]

## Calculate total reads in negative and true samples
#contam_df05$neg_reads <- rowSums(contam_df05[, neg_samples], na.rm = TRUE)
#contam_df05$true_reads <- rowSums(contam_df05[, true_samples], na.rm = TRUE)
#contam_df05$mock_reads <- rowSums(contam_df05[, mock_samples], na.rm = TRUE)

# selects only ASVs with at least one negative sample count higher than all true sample counts. 
# Identify contaminants based on the condition
contaminants05 <- rownames(contam_df05)[apply(contam_df05, 1, function(asv_counts) {
  max(asv_counts[neg_samples]) > max(asv_counts[true_samples])
})]

# Filter the original dataframe to retain only identified contaminants
contam_df05_filtered <- contam_df05[contaminants05, ]

# Display the identified contaminant ASVs
contam_df05_transpose <- t(contam_df05_filtered) # 13 ASVs has higher reads in negative samples compared to true and mock samples. these will be removed as contaminant
# for these 13 ASVs, Negatives have higher higher read than true samples and mocks. Remove these

# Besides check the tax table for contam_asvs05 and and if 
# Acinetobacter, Burkholderia, Bacillus, Corynebacterium, Escherichia/Shigella, 
# Pseudomonas, Propionibacterium/Cutibacterium, Ralstonia, Sphingomonas, Streptococcus,  are present with high read in negative. If so, they will be contaminant, 
view(ps.16s@tax_table[contam_asvs_05,])

# ASV2371 Acinetobacter
ASV2371 <- prune_taxa("ASV2371", ps.16s) %>% tidyPS() %>% otu_table() 
ASV2371 # relatively high is true samples, no contaminant

# ASV9736 Corynebacterium
ASV9736 <- prune_taxa("ASV9736", ps.16s) %>% tidyPS() %>% otu_table() 
ASV9736 # relatively high is true samples, no contaminant

# ASV11106 Corynebacterium
ASV11106 <- prune_taxa("ASV11106", ps.16s) %>% tidyPS() %>% otu_table() 
ASV11106 # relatively high is true samples, no contaminant

# ASV11318 Corynebacterium
ASV11318 <- prune_taxa("ASV11318", ps.16s) %>% tidyPS() %>% otu_table() 
ASV11318 # relatively high is true samples, no contaminant

# ASV5003 Corynebacterium
ASV5003 <- prune_taxa("ASV5003", ps.16s) %>% tidyPS() %>% otu_table() 
ASV5003 # relatively high is true samples, no contaminant

# ASV8890 Corynebacterium
ASV8890 <- prune_taxa("ASV8890", ps.16s) %>% tidyPS() %>% otu_table() 
ASV8890 # relatively high is true samples, no contaminant

# ASV8071 Corynebacterium
ASV8071 <- prune_taxa("ASV8071", ps.16s) %>% tidyPS() %>% otu_table() 
ASV8071 # relatively high is true samples, no contaminant

# ASV5149 Escherichia-Shigella
ASV5149 <- prune_taxa("ASV5149", ps.16s) %>% tidyPS() %>% otu_table() 
ASV5149 # Relatively high in one samples than the negative, no contaminant

# ASV26507 Pseudomonas
ASV26507 <- prune_taxa("ASV26507", ps.16s) %>% tidyPS() %>% otu_table() 
ASV26507 # Relatively high in true samples, no contaminant

# ASV9175 Pseudomonas
ASV9175 <- prune_taxa("ASV9175", ps.16s) %>% tidyPS() %>% otu_table() 
ASV9175 # Relatively high in true samples, no contaminant

# ASV5150 Pseudomonas 
ASV5150 <- prune_taxa("ASV5150", ps.16s) %>% tidyPS() %>% otu_table() 
ASV5150 # relatively high is true samples, no contaminant

# ASV6873 Pseudomonas
ASV6873 <- prune_taxa("ASV6873", ps.16s) %>% tidyPS() %>% otu_table() 
ASV6873 # relatively high is true samples, no contaminant

# ASV11744 Ralstonia
ASV11744 <- prune_taxa("ASV11744", ps.16s) %>% tidyPS() %>% otu_table() 
ASV11744 # relatively high is true samples, no contaminant

# so 13 ASVs in contam_df05_filtered are possible contaminants. Remove these

# Remove contaminants----
#ps.nocontam05 <- prune_taxa(!(taxa_names(ps.16s) %in% contaminants05), ps.16s)
#ps.nocontam05 # 51153 taxa and 352 samples

# As ASVs are in datafreame now, identify them to remove
asvs_to_remove <- rownames(contam_df05_filtered) 

# Remove ASVs from the phyloseq object
ps.16s_nocontam05 <- prune_taxa(!(taxa_names(ps.16s) %in% asvs_to_remove), ps.16s)
ps.16s_nocontam05 # 51140 taxa and 352 samples

#add reads number into the metadata----
sample_data(ps.16s_nocontam05)$Read_depth <- sample_sums(ps.16s_nocontam05)
ps.16s_nocontam05 # 51140 taxa and 352 samples

# Save this as RDS file
#saveRDS(ps.16s_nocontam05, "phyloseq_FSS_BD_no_contam_16s_20241016.rds")

###################
# load latest phyloseq obect----
ps.16s1 <- readRDS("phyloseq_FSS_BD_no_contam_16s_20241016.rds")
ps.16s1 # 51140 taxa and 352 samples

#Remove low read count samples----
ps.16s_high_counts <- prune_samples(sample_sums(ps.16s1)>=2000, ps.16s1) %>% tidyPS() #remove low read count (less than 2000 reads) samples that have failed
ps.16s_high_counts # 51140 taxa and 343 samples 

# subset only mock to check the accuracy later----
ps.16s_mock <- ps.16s_high_counts %>% 
  ps_filter(Sample_type %in% c("Mock_DNA", "Mock_extract"))
ps.16s_mock # 67 taxa and 4 samples

# save as rds file
#saveRDS(ps.16s_mock, "phyloseq_FSS_BD_mock_16s_20241016.rds")

# as previously we have removed contaminants using negs, and also I have separated 
# mock for checking accuracy later, we do not need them any more, so remove them
ps.16s.no.control <- ps.16s_high_counts %>% 
  ps_filter(Sample_type %in% c("Gill_swab", "Skin_swab", "Pond_water"))
ps.16s.no.control # 51140 taxa and 334 samples

# save this final phyloseq object for downstream analysis
#saveRDS(ps.16s.no.control, "phyloseq_16S_all_filtered_with_tree_16s_20241016.rds")

### As the 16s data has over 50k taxa, I will filter the taxa
#### with prevalence >= 2 (same as before) and with a TotalAbundance >= 30
# Load the updated phyloseq object
ps.16s <- readRDS("phyloseq_16S_all_filtered_with_tree_16s_20241016.rds")

#do some filteering like on this https://f1000research.com/articles/5-1492
### prevalence filter----
#This function will remove taxa (ASVs) with low prevalence, where prevalence is the fraction of total samples in which an ASV is observed.
# Compute prevalence of each feature, store as data.frame
# Filter to remove the rare taxa
prevdf <- apply(X = otu_table(ps.16s), MARGIN = ifelse(taxa_are_rows(ps.16s), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps.16s))

# keep taxa those are present in at least in 2 sample and with a minimum abundace of 30
keepTaxa30a <- rownames(prevdf)[(prevdf$Prevalence >= 2 & prevdf$TotalAbundance >= 30)]
ps1 <- prune_taxa(keepTaxa30a, ps.16s)
ps1 # 33312 taxa and 334 samples

# save the new phyloseq object
#saveRDS(ps1, "phyloseq_FSS_BD_metadata_v7_16S_abd30_20241016.rds")



