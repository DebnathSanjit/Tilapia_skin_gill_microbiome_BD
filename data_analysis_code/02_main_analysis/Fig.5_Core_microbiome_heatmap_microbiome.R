# Phyloseq -> Subset -> core (microbiome package)----
# 29/10/2024
# Author: Sanjit Debnath

# This script is from Leo Lahti, Sudarshan Shetty et al. for core microbiome. 
# core microbiome calculation (https://microbiome.github.io/tutorials/Core.html)

# I'm calculating gill and core microbiota, not for water. Combine both and then subset these taxa from the original phyloseq object and plot these

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(microbiome); packageVersion("microbiome")
library(viridis); packageVersion("viridis")
library(RColorBrewer); packageVersion("RColorBrewer")
#library(openxlsx); packageVersion("openxlsx") # to save table as excel file
library(writexl); packageVersion("writexl") # to save table as excel file

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

#to make the variable using my desired color #set variables desired color
sample.colors <- c("Gill_swab"= "#4363d8", "Skin_swab" = "#8F7700FF", "Pond_water"= "#5EC747",
                   "Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5",
                   "Gill_Diseased" = "#F35E5A", "Gill_Non-diseased" = "#9DCC00",
                   "Skin_Diseased" = "#B68A06", "Skin_Non-diseased" = "#7cb9cb",
                   "Water_Non-diseased" = "#619CFF", "Water_Diseased" = "#FDBF6F")

# Prokaryotes
#load phyloseq object
ps <- readRDS("phyloseq_FSS_BD_metadata_v7_16S_abd30_20241016.rds")
ps # 33312 taxa and 334 samples

# Remove others except tilapia, combine sample_type and disease state, rename sample type 
pseq <- ps %>%
  subset_taxa(Kingdom !="Archaea"  | is.na(Kingdom)) %>% #remove archaea, over 1000 taxa
  ps_filter(Sampled_species == "Tilapia" | Sample_type == "Pond_water") %>%
  ps_mutate(
    CombinedType = paste(Sample_type, Reported_disease, sep = "_"),
    Sample_Disease = ifelse(Sample_type == "Gill_swab", "Gill", 
                            ifelse(Sample_type == "Skin_swab", "Skin", 
                                   ifelse(Sample_type == "Pond_water", "Water", Sample_type))) %>% 
      paste(., Reported_disease, sep = "_")
  ) %>% tax_fix() # use tax_fix to fix the unclassified rank
pseq #  32127 taxa and 298 samples

# there were some problem with the tax_table. so need to fix it first (https://github.com/joey711/phyloseq/issues/1551)
#Problematic Genus values detected in tax_table:
# Incertae Sedis / Unknown Family Family
# Fix taxonomic labels in your phyloseq object
pseq_fixed <- tax_fix(pseq, unknowns = c("Incertae Sedis", "Unknown Family Family", "Unknown Family"))

## 1. Calculate relative abundance----
pseq.rel <- microbiome::transform(pseq_fixed, "compositional")

## Subset gill----
gill_swab <- pseq.rel %>% 
  ps_filter(Sample_type == "Gill_swab")
gill_swab # 31264 taxa and 118 samples

## Calculate core----
# With A full phyloseq object of the core microbiota is obtained as follow
gill_core <- core(gill_swab, detection = 0.0001, prevalence = 0.90)
gill_core  # 21 taxa and 118 samples

### Calculate Phylum and genus
# Extract the taxonomy table and convert it to a data frame
#tax_table_gill <- as.data.frame(tax_table(gill_core))
tax_table_gill <- as.data.frame(phyloseq::tax_table(gill_core))

# Convert row names to a column named "ASV"
tax_table_gill <- tax_table_gill %>% 
  rownames_to_column(var = "ASV")

# Count the number of unique Phylum
gill_phylum <- length(unique(tax_table_gill$Phylum))
gill_phylum # 7

# Count the number of unique Genus
gill_genus <- length(unique(tax_table_gill$Genus))
gill_genus # 14

## Subset skin----
skin_swab <- pseq.rel %>% 
  ps_filter(Sample_type == "Skin_swab")
skin_swab # 29594 taxa and 120 samples

## Calculate core----
# With A full phyloseq object of the core microbiota is obtained as follow
skin_core <- core(skin_swab, detection = 0.0001, prevalence = 0.90)
skin_core  # 17 taxa and 118 samples

# Extract the taxonomy table and convert it to a data frame
#tax_table_skin <- as.data.frame(tax_table(skin_core))
tax_table_skin <- as.data.frame(phyloseq::tax_table(skin_core))

# Convert row names to a column named "ASV"
tax_table_skin <- tax_table_skin %>% 
  rownames_to_column(var = "ASV")

# Count the number of unique Phylum
skin_phylum <- length(unique(tax_table_skin$Phylum))
skin_phylum # 7

# Count the number of unique Genus
skin_genus <- length(unique(tax_table_skin$Genus))
skin_genus # 11

##### Create a list of core tax tables----
core_tax_table <- list(Gill = tax_table_gill, 
                    Skin = tax_table_skin)

##### Save as Excel file----
#write_xlsx(core_tax_table, "Gill and skin core microbiota tax tables_20241029.xlsx")

## Combine skin and gill core
# Remove the phylogenetic trees as with different number of tree nodes, i can't combine phyloseq object
gill_core_no_tree <- phyloseq(otu_table(gill_core), sample_data(gill_core), tax_table(gill_core))
skin_core_no_tree <- phyloseq(otu_table(skin_core), sample_data(skin_core), tax_table(skin_core))

# Merge the phyloseq objects without trees
combined_swab_core <- merge_phyloseq(gill_core_no_tree, skin_core_no_tree)
combined_swab_core # 23 taxa and 238 samples

## Subset these taxa from the original phyloseq object----
# Extract taxa names from merged_phyloseq
combined_swab_taxa_names <- taxa_names(combined_swab_core)

# Use these names in subset_taxa to filter s16_ps_gen
swab_core <- subset_taxa(pseq.rel, rownames(tax_table(pseq.rel)) %in% combined_swab_taxa_names)
swab_core # 23 taxa and 298 samples

### Plot heatmap----
# Now plot heatmap using microViz package. It's much easy and I can modify as i want the plot
# If i make the variable factor, the heatmap does work. So I added a new column with sample names with number so that they appear as i want
swab_core <- swab_core %>%
  ps_mutate(
    Sample_type2 = ifelse(Sample_type == "Gill_swab", "01.Gill", 
                          ifelse(Sample_type == "Skin_swab", "02.Skin", 
                                 ifelse(Sample_type == "Pond_water", "03.Water", Sample_type)))
  ) %>% 
  tax_fix()  # use tax_fix to fix the unclassified rank
swab_core #  23 taxa and 298 samples 

#### Genus----
heatmap.gen <- swab_core %>%
  # sort all samples by similarity
  ps_seriate(rank = "Genus", tax_transform = "compositional", dist = "bray") %>%
  # arrange the samples into Disease State groups
  ps_arrange(Sample_type2, Reported_disease) %>%
  #tax_filter(#min_prevalence = 0.90, # taxa must be present in at least 10% of the samples
  #  min_total_abundance = 0.1) %>% # RA should be > 1%
  tax_transform("compositional", rank = "Genus") %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(2, "cm")) # to adjust prevalence bar
    ),
    row_names_gp = grid::gpar(fontsize = 15, fontface = "italic"),
   # colors = heat_palette(palette = "YlGnBu", rev = TRUE), # If I want to change the color. try colorspace::hcl_palettes() and choose a color palette. "Rocket" is default
    #colors = heat_palette(palette = "Rocket", rev = TRUE),
    colors = heat_palette(palette = "GnBu", rev = TRUE),
   #colors = heat_palette(palette = "YlGn", rev = TRUE),
    grid_col = NA,
    sample_anno = sampleAnnotation(
      "Sample type" = anno_sample("Sample_type2"),
      "Reported disease" = anno_sample("Reported_disease"),
      col = cols,
      border = FALSE
    ),
    sample_seriation = "Identity" # suppress sample reordering
  )
heatmap.gen

# Save as 1500*900

