# Date: 16/10/2024
# This script is to plot the expected mock community vs mock community in sequenced samples

# Load Libraries----
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

# Expected mock composition
# https://files.zymoresearch.com/protocols/_d6305_d6306_zymobiomics_microbial_community_dna_standard.pdf
# https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-dna-standard

# Prokaryotes----
# Create the data frame with species names and their proportions mentioned in the zymo website
zymo.com <- data.frame(
  Species = c("Pseudomonas aeruginosa","Escherichia coli","Salmonella enterica",
              "Lactobacillus fermentum","Enterococcus faecalis",
              "Staphylococcus aureus","Listeria monocytogenes","Bacillus subtilis"), 
  #"Saccharomyces cerevisiae", "Cryptococcus neoformans"), # since these are eukaryotes, they will not be present in 16s, so remove them and their proportion
  Proportion = c(4.2,10.1,10.4,18.4,9.9,15.5,14.1,17.4),
  Genus = c("Pseudomonas","Escherichia","Salmonella","Lactobacillus","Enterococcus","Staphylococcus","Listeria","Bacillus"),
  Genus2 = c("Bacillus", "Enterococcus","Escherichia-Shigella","Limosilactobacillus","Listeria","Pseudomonas","Salmonella","Staphylococcus"))

# Sort the genus as they were present in the zymo bouchure
zymo.com <- zymo.com[order(zymo.com$Genus2), ]

# set colours
#species.color <- c("Listeria monocytogenes" = "#FDBF6F", "Pseudomonas aeruginosa" = "#b249d5", 
#                   "Bacillus subtilis" = "#1F78B4", "Escherichia coli" = "#33A02C", 
#                   "Salmonella enterica" = "#B2DF8A", "Lactobacillus fermentum" = "#A6CEE3", 
#                   "Enterococcus faecalis" = "#FF7F00", "Staphylococcus aureus" = "#FB9A99") 
genus.color <- c("Listeria" = "#FDBF6F", "Pseudomonas" = "#b249d5", 
                   "Bacillus" = "#1F78B4", "Escherichia" = "#33A02C", 
                   "Salmonella" = "#B2DF8A", "Lactobacillus" = "#A6CEE3", 
                   "Enterococcus" = "#FF7F00", "Staphylococcus" = "#FB9A99") 

#"Saccharomyces cerevisiae" = "#A6761D", "Cryptococcus neoformans" = "#5c47b8")

# Plot the stacked bar chart
zymo.exp <- ggplot(zymo.com, aes(x = "", y = Proportion/100, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus.color) +
  labs(title = "Expected Zymo mock community", 
       x = "Expected \nmock community", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", margin = margin(t = 5)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.key.height = unit(0.8, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill = guide_legend(ncol = 1))


# Load the phyloseq object
ps_mock <- readRDS("phyloseq_FSS_BD_mock_16s_20241016.rds")
ps_mock <- ps_mock %>% tax_fix()
ps_mock # 67 taxa and 4 samples

### bar plot (microviz)----
# my palette
myPal.gen <- tax_palette(
  data = ps_mock, rank = "Genus", n = 25, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
myPal.gen["Pseudomonas"] <- "#b249d5" # Override existing color if any color is not good
tax_palette_plot(myPal.gen)

# Barplot
zymo.extra <- ps_mock %>%
  comp_barplot(tax_order = sort(zymo.com$Genus2, decreasing = T),
               tax_level = "Genus",
               n_taxa = 8,
               tax_sort = "name",
               sample_order = "asis", 
               bar_width = 0.9,
               bar_outline_colour = NA,
               merge_other = FALSE,
               label = "Sample",
               palette = myPal.gen) +
  labs(title = "Relative abundance of mock community",
       y = "Relative Abundance") +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", margin = margin(t = 5)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.key.height = unit(0.8, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(reverse=T)) # Reverse the legend
zymo.extra

# Combine plots
zymo.16s <- cowplot::plot_grid(zymo.exp,
                   zymo.extra,
                   ncol = 2,
                   labels = "AUTO",
                   rel_widths = c(1.3,2),
                   align = "h")
# Save as 1500*800

# Microeukaryotes----
# The expected genus are not present in mock samples

# Create the data frame with species names and their proportions mentioned in the zymo website
zymo.com.18s <- data.frame(
  Species = c("Saccharomyces cerevisiae", "Cryptococcus neoformans"), # since these are eukaryotes, they will not be present in 16s, so remove them and their proportion
  Proportion = c(73.81, 26.19), # the website says 9.3 and 3.3%. Since I'm comparing these two only, I transform them as relative abundance (9.3*100)/12.6
  Genus = c("Saccharomyces","Cryptococcus"),
  Genus2 = c("Saccharomyces","Cryptococcus"),
  Family = c("Saccharomycetaceae","Cryptococcaceae"))

# Sort the genus as they were present in the zymo bouchure
zymo.com.18s <- zymo.com.18s[order(zymo.com.18s$Family), ]

# set colours
#species.color.18s <- c("Saccharomyces cerevisiae" = "#A6761D", "Cryptococcus neoformans" = "#5c47b8")
#genus.color.18s <- c("Saccharomyces" = "#A6761D", "Cryptococcus" = "#5c47b8")
fam.color.18s <- c("Saccharomycetaceae" = "#A6761D", "Cryptococcaceae" = "#5c47b8")

# Plot the stacked bar chart
zymo.exp.18s <- ggplot(zymo.com.18s, aes(x = "", y = Proportion/100, fill = Family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fam.color.18s) +
  labs(title = "Expected Zymo mock community", 
       x = "Expected \nmock community", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", margin = margin(t = 5)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.key.height = unit(0.8, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill = guide_legend(ncol = 1))
zymo.exp.18s

# Load the phyloseq object
ps.18s_mock <- readRDS("phyloseq_18s_mock_20240309.rds")
ps.18s_mock <- ps.18s_mock %>% tax_fix()
ps.18s_mock # 5 taxa and 2 samples

# #removing Craniata----
ps_no_Craniata <- ps.18s_mock %>% subset_taxa(Class !="Craniata"  | is.na(Class))
ps_no_Craniata # 4 taxa and 2 sample

#removing Teleostei
ps_no_Teleostei <- ps_no_Craniata %>% subset_taxa(Family !="Teleostei"  | is.na(Family))
ps_no_Teleostei #  4 taxa and 2 samples 

mock.18s <- ps_no_Teleostei

### bar plot (microviz)----
## my palette
#myPal.fam.18s <- tax_palette(
#  data = mock.18s, rank = "Family", n = 25, pal = "brewerPlus", # At genus level is difficult to get the community for 16S. so do it for family level
#  add = c(Other = "lightgrey"))
##myPal.gen["Pseudomonas"] <- "#b249d5" # Override existing color if any color is not good
#tax_palette_plot(myPal.fam.18s)
#genus.color.18s <- c("Saccharomyces" = "#A6761D", "Cryptococcus" = "#5c47b8")
#fam.color.18s <- c("Saccharomycetaceae" = "#A6761D", "Cryptococcaceae" = "#5c47b8")
fam.color.18s2 <- c("Saccharomycetales" = "#A6761D", "Tremellomycetes" = "#5c47b8", "Other" = "grey")

# Barplot
zymo.extra.18s <- mock.18s %>%
  comp_barplot(#tax_order = sort(zymo.com.18s$Genus2, decreasing = T),
    tax_level = "Family",
    n_taxa = 2,
    tax_sort = "name",
    sample_order = "asis", 
    bar_width = 0.9,
    bar_outline_colour = NA,
    merge_other = FALSE,
    label = "Sample",
    palette = fam.color.18s2) +
  labs(title = "Relative abundance of mock community",
       y = "Relative Abundance") +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold", margin = margin(t = 5)),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.key.height = unit(0.8, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(reverse=T)) # Reverse the legend
zymo.extra.18s

# Combine plots
zymo.18s <- cowplot::plot_grid(zymo.exp.18s,
                   zymo.extra.18s,
                   ncol = 2,
                   labels = c("C", "D"),
                   rel_widths = c(1.3,2),
                   align = "h")
# Save as 1500*800

# Combine 16s and 18s together
cowplot::plot_grid(zymo.16s,
                   zymo.18s,
                   ncol = 1)
# Save as 2000*1200
