# Date: 14/11/2024 
# author: Sanjit Debnath

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(patchwork); packageVersion("patchwork")
library(colorblindr); packageVersion("colorblindr")
library(paletteer); packageVersion("paletteer")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

# Define custom themes function
my_custom_theme1 <- function() {
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)), 
    legend.position = "none"
  ) 
} # For first plot

my_custom_theme2 <- function() {
  theme(plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_blank(),#element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),#element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.8, "cm"))  # Adjust the height of legend keys
} # for pond wise plot with legend

my_custom_theme3 <- function() {
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), 
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_blank(), #text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),#text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "none"
  )
} # for with mean plot without legend

# Genus
my_custom_theme4 <- function() {
  theme(plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_blank(),#text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),#text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.key.height = unit(0.8, "cm"))
} # for pond wise plot with legend

#Prokaryotes----
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
  )
pseq # 32127 taxa and 298 samples

# Define the correct order of Pond_name
pond_order <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20")

# Modify Pond_name factor levels in the sample data
sample_data(pseq)$Pond_name <- factor(sample_data(pseq)$Pond_name, levels = pond_order)

# Define the correct order of Sample_Disease
sample_order <- c("Gill_Diseased", "Gill_Non-diseased", "Skin_Diseased", "Skin_Non-diseased", "Water_Diseased", "Water_Non-diseased")

# Modify Sample_Disease factor levels in the sample data
sample_data(pseq)$Sample_Disease <- factor(sample_data(pseq)$Sample_Disease, levels = sample_order)

# there were some problem with the tax_table that needs to be fixed for microViz package. so need to fix it first (https://github.com/joey711/phyloseq/issues/1551)
# Fix taxonomic labels in your phyloseq object
pseq_fixed <- tax_fix(pseq, unknowns = c("Incertae Sedis", "Unknown Family Family", "Unknown Family"))

# Rank = Phylum----
# Define a named vector of colors for different phyla
phyla.colors <- c(
  "Proteobacteria" = "#00D9FF", "Bacteroidota" = "#875692",  "Planctomycetota" = "#f38400",
  "Actinobacteriota" = "#a1caf1", "Verrucomicrobiota" = "#f3c300","Firmicutes" = "#c2b280",
  "Fusobacteriota" = "#ABE496",  "Cyanobacteria" = "#008856",  "Chloroflexi" = "#e68fac",
  "Acidobacteriota" = "#0067a5",  "Patescibacteria" = "#0A47FFFF",  "Bdellovibrionota" = "#604e97",
  "Deinococcota" = "#f6a600",  "Bacteria Kingdom" = "#b3446c",  "Myxococcota" = "#dcd300",
  "Desulfobacterota" = "#882d17",  "Gemmatimonadota" = "#f99379",  "SAR324 clade(Marine group B)" = "#654522",
  "Armatimonadota" = "#e25822",  "Halobacterota" = "#8db600",  "Other" = "lightgrey") # Default color for any other phyla

## Non-diseased gill----
ndgill <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Gill_Non-diseased")
#ndgill
topPhy.ndgill <- ndgill %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

### Average###
p1.avg <- ndgill %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.ndgill, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Gill_swab" = "G\nND")) + # this is necessary to make the plots with same height
  labs(title = "GND",
       y = "Relative abundance")+
  #x = "Gill\n(non-diseased)") +
  guides(fill = guide_legend(ncol = 1)) +
  my_custom_theme1()

### Pond-wise###
p1.pond <- ndgill %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.ndgill, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "GND",
       y = "Relative abundance",
       x = "Pond name") +
  guides(fill = guide_legend(ncol = 1)) +
  my_custom_theme2()

## Diseased gill----
dgill <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Gill_Diseased")
#dgill
topPhy.dgill <- dgill %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

# Average
p2.avg <- dgill %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.dgill, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Gill_swab" = "G\nD"))+
  labs(title = "GD") +
  my_custom_theme3()

# Pondwise
p2.pond <- dgill %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.dgill, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "GD") +
  my_custom_theme2()

## Non-diseased skin----
ndskin <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Skin_Non-diseased")
#ndskin
topPhy.ndskin <- ndskin %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

# Average
p3.avg <- ndskin %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.ndskin, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Skin_swab" = "S\nND"))+
  labs(title = "SND")+
  #y = "Relative abundance") +
  my_custom_theme3()

# Pond-wise
p3.pond <- ndskin %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.ndskin, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "SND")+
  #y = "Relative abundance") +
  my_custom_theme2()

## Diseased skin----
dskin <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Skin_Diseased")
#dskin
topPhy.dskin <- dskin %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

# Average
p4.avg <- dskin %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.dskin, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Skin_swab" = "S\nD"))+
  labs(title = "SD")+
  #y = "Relative abundance") +
  my_custom_theme3()

# Pond-wise
p4.pond <- dskin %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.dskin, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "SD")+
  #y = "Relative abundance") +
  my_custom_theme2()

## Non-diseased water----
ndwater <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Water_Non-diseased")
#ndwater
topPhy.ndwater <- ndwater %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

# Average
p5.avg <- ndwater %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.ndwater, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nND"))+
  labs(title = "WND") +
  my_custom_theme3()

# Pond-wise
p5.pond <- ndwater %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.ndwater, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "WND") +
  my_custom_theme2()

## Diseased water----
dwater <- pseq_fixed %>% 
  ps_filter(Sample_Disease == "Water_Diseased")
#dwater
topPhy.dwater <- dwater %>%
  tax_top(n = 10, rank = "Phylum") %>%
  sort() # this makes them alphabetical

# Average
p6.avg <- dwater %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topPhy.dwater, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nD"))+
  labs(title = "WD") +
  my_custom_theme3()

# Pond-wise
p6.pond <- dwater %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topPhy.dwater, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "WD") +
  my_custom_theme2()

#### Combine Phylum plot----
comb.phy <- cowplot::plot_grid(p1.avg,
                               p2.avg,
                               p3.avg,
                               p4.avg,
                               p5.avg,
                               p6.avg,
                               p1.pond + theme(legend.position = "none"),
                               p2.pond + theme(legend.position = "none"),
                               p3.pond + theme(legend.position = "none"),
                               p4.pond + theme(legend.position = "none"),
                               p5.pond + theme(legend.position = "none"),
                               p6.pond + theme(legend.position = "none"),
                               rel_widths = c(3.6,2,2,2,2,2,4,4,4,4,4,4),
                               nrow = 1)
comb.phy

# Rank = Genus----
# add specific color for the interested genera
gen.colors <- c(
  "Acinetobacter" = "#66A61E",
  "Aeromonas" = "#0A47FFFF",
  "Alcaligenaceae Family" = "#654522",
  "Bacteria Kingdom" = "lightgreen",
  "Brevibacterium" = "#e68fac",
  "Cetobacterium" = "#B2DF8A",
  "Chryseobacterium" = "#A6761D",
  "Clostridiaceae Family" = "#83d5af",
  "Clostridium sensu stricto 1" = "#1B9E77",
  "CL500-3" = "#A6CEE3",
  "CL500-29 marine group" = "#5e79b2",
  "Comamonadaceae Family" = "#FFFF99",
  "Cyanobium PCC-6307" = "#6A3D9A",
  "Enhydrobacter" = "#e25822",
  "env.OPS 17 Family" = "#5c47b8",
  "Flavobacterium" = "#00D9FF",
  "Fluviicola" = "#b249d5",
  "hgcI clade" ="#1F78B4", 
  "Isosphaeraceae Family" = "#CCCCCCFF",
  "Kapabacteriales Order" = "#cfd251",
  "Lactobacillus" = "#f99379",
  "LD29" = "#FDBF6F",
  "Methylacidiphilaceae Family" = "#B15928",
  "MWH-UniP1 aquatic group Family" = "#5f7b35",
  "Mycobacterium" = "#CAB2D6",
  "PeM15 Order" = "#666666",
  "Pirellulaceae Family" = "#33A02C",
  "Phycisphaeraceae Family" = "lightyellow",
  "Polynucleobacter" = "#E31A1C",
  "Prevotella_7" = "#ff69b4",
  "Prevotellaceae Family" = "#773a27",
  "Proteobacteria Phylum" = "#7cb9cb",
  "Rhodocyclaceae Family" = "#E7298A",
  "Rubinisphaeraceae Family" = "#E6AB02",
  "SH3-11" = "#7570B3",
  "Sphingomonas" = "#4C005C",
  "Verrucomicrobiae Class" = "#7edc45",
  "Vibrio" = "#8db600",
  "Vogesella" = "#532a5a",
  "2013Ark19i" = "#d3c4a8",
  "Other" = "lightgrey"
)

# To calculate the top 20 genera, if I consider the whole phyloseq object 
## and also if I split the disease and non-diseased, several genera differ.
### So, I need to subset disease and non-diseased, calculate top 20 genera for 
#### both diseased and non-diseased samples. This is different than phyla

## Non-diseased gill----
topGen.ndgill <- ndgill %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p1.avg.g <- ndgill %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.ndgill, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Gill_swab" = "G\nND")) +
  labs(title = "GND",
       y = "Relative abundance")+
  #x = "Gill\n(non-diseased)") +
  my_custom_theme1()
p1.avg.g

# Pond-wise
p1.pond.g <- ndgill %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.ndgill, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "GND",
       y = "Relative abundance") +
  my_custom_theme4() +
  guides(fill = guide_legend(ncol = 1)) #+ # Adjust the number of rows in the legend


## Diseased gill----
topGen.dgill <- dgill %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p2.avg.g <- dgill %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.dgill, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Gill_swab" = "G\nD"))+
  labs(title = "GD") +
  my_custom_theme3()

# Pond-wise
p2.pond.g <- dgill %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.dgill, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "GD") +
  my_custom_theme4()


## Non-diseased skin----
topGen.ndskin <- ndskin %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p3.avg.g <- ndskin %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.ndskin, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Skin_swab" = "S\nND"))+
  labs(title = "SND")+
  #y = "Relative abundance") +
  my_custom_theme3()

# Pond-wise
p3.pond.g <-  ndskin %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.ndskin, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "SND")+
  #y = "Relative abundance") +
  my_custom_theme4()

## Diseased skin----
topGen.dskin <- dskin %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p4.avg.g <- dskin %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.dskin, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Skin_swab" = "S\nD"))+
  labs(title = "SD")+
  #y = "Relative abundance") +
  my_custom_theme3()

# Pond-wise
p4.pond.g <-  dskin %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Skin_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.dskin, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "SD")+
  #y = "Relative abundance") +
  my_custom_theme4()

## Non-diseased water----
topGen.ndwater <- ndwater %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p5.avg.g <- ndwater %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.ndwater, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nND"))+
  labs(title = "WND") +
  my_custom_theme3()

# Pond-wise
p5.pond.g <-  ndwater %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.ndwater, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "WND") +
  my_custom_theme4()

## Diseased water----
topGen.dwater <- dwater %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
p6.avg.g <- dwater %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.dwater, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nD"))+
  labs(title = "WD") +
  my_custom_theme3()

# Pond-wise
p6.pond.g <-  dwater %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Sample_type == "Pond_water") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.dwater, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = gen.colors) +
  #coord_flip() +
  labs(title = "WD") +
  my_custom_theme4()

#### Combine genus----
comb.gen <- cowplot::plot_grid(p1.avg.g, 
                               p2.avg.g,
                               p3.avg.g,
                               p4.avg.g,
                               p5.avg.g,
                               p6.avg.g,
                               p1.pond.g + theme(legend.position = "none"),
                               p2.pond.g + theme(legend.position = "none"),
                               p3.pond.g + theme(legend.position = "none"),
                               p4.pond.g + theme(legend.position = "none"),
                               p5.pond.g + theme(legend.position = "none"),
                               p6.pond.g + theme(legend.position = "none"),
                               rel_widths = c(3.6,2,2,2,2,2,4,4,4,4,4,4),
                               nrow = 1)
comb.gen

# Microeukaryotes----
ps.18s <- readRDS("phyloseq_FSS_BD_metadata_v7_18S_20240703.rds")
ps.18s <- ps.18s %>% tax_fix() 
ps.18s # 2961 taxa and 58 samples

# Define the correct order of Pond_name
pond_order <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20")

# Modify Pond_name factor levels in the sample data
sample_data(ps.18s)$Pond_name <- factor(sample_data(ps.18s)$Pond_name, levels = pond_order)

# As 18s has only water sample, we can directly plot diseased vs non-diseased
# my palette
myPal.div.18s <- tax_palette(
  data = ps.18s, rank = "Division", n = 20, pal = "brewerPlus",
  add = c(Other = "lightgrey"))
myPal.div.18s["Opisthokonta"] <- "#00D9FF" # Override existing color if any color is not good
#myPal.phy.16s["Fusobacteriota"] <- "#ABE496"
tax_palette_plot(myPal.div.18s)

## Non-diseased----
topTaxa.ndis.18s <- ps.18s %>%
  ps_filter(Reported_disease == "Non-diseased") %>%
  tax_top(n = 10, rank = "Division") %>%
  sort() # this makes them alphabetical

## Non-diseased----
# Average
water_div.ndis <- ps.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Sample_type) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Non-diseased") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topTaxa.ndis.18s, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nND"))+
  labs(title = "WND",
       y = "Relative abundance") +
  my_custom_theme1()

# Pond-wise
water_div.ndis2 <- ps.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Non-diseased") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topTaxa.ndis.18s, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  labs(title = "WND",
       y = "Relative abundance") +
  my_custom_theme2()
#water_div.ndis2

## Diseased----
topTaxa.dis.18s <- ps.18s %>%
  ps_filter(Reported_disease == "Diseased") %>%
  tax_top(n = 10, rank = "Division") %>%
  sort() # this makes them alphabetical

# Average
water_div.dis <- ps.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Sample_type) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Diseased") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topTaxa.dis.18s, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 10, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nD"))+
  labs(title = "WD",
       y = "Relative abundance") +
  my_custom_theme3()
#water_div.dis

# Pond-wise
water_div.dis2 <-  ps.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Diseased") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topTaxa.dis.18s, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 10, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  labs(title = "WD",
       y = "Relative abundance") +
  my_custom_theme2()
#water_div.dis2

### Combine----
comb.div.18s <- cowplot::plot_grid(water_div.ndis,
                                   water_div.dis,
                                   water_div.ndis2 + theme(legend.position = "none"),
                                   water_div.dis2 + theme(legend.position = "none"),
                                   ncol = 4,
                                   rel_widths = c(2.5,2,4,4))
comb.div.18s

# Rank = Genus----
# my palette
myPal.gen.18s <- tax_palette(
  data = ps.18s, rank = "Genus", n = 20, pal = "kelly",
  add = c(Other = "lightgrey"))
myPal.gen.18s["Cryptomonadales_X Family"] <- "#5e79b2" # Override existing color if any color is not good
myPal.gen.18s["Chrysophyceae Class"] <- "#426600" 
myPal.gen.18s[ "Mediophyceae Class"] <- "#FFA8BB"
myPal.gen.18s["Tovellia"] <- "#ABE496"
#myPal.gen.18s["Dinophyceae Class"] <- 
tax_palette_plot(myPal.gen.18s)

## Non-diseased----
topGen.ndis.18s <- ps.18s %>%
  ps_filter(Reported_disease == "Non-diseased") %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
water_gen.ndis <- ps.18s %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Sample_type) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Non-diseased") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.ndis.18s, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.gen.18s) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nND"))+
  labs(title = "WND",
       y = "Relative abundance") +
  my_custom_theme1()
#water_gen.ndis

# Pond-wise
water_gen.ndis2 <- ps.18s %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Non-diseased") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.ndis.18s, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.gen.18s) +
  #coord_flip() +
  labs(title = "WND",
       y = "Relative abundance") +
  my_custom_theme4()
#water_gen.ndis2

## Diseased----
topGen.dis.18s <- ps.18s %>%
  ps_filter(Reported_disease == "Diseased") %>%
  tax_top(n = 15, rank = "Genus") %>%
  sort() # this makes them alphabetical

# Average
water_gen.dis <- ps.18s %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Sample_type) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Diseased") %>%
  phyloseq::merge_samples(group = "Sample_type") %>%
  comp_barplot(tax_order = topGen.dis.18s, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.gen.18s) +
  #coord_flip() +
  scale_x_discrete(labels = c("Pond_water" = "W\nD"))+
  labs(title = "WD") +
  my_custom_theme3()
#water_gen.dis

# Pond-wise
water_gen.dis2 <-  ps.18s %>%
  tax_sort(by = sum, at = "Genus") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(Reported_disease == "Diseased") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topGen.dis.18s, # this brings the named taxa to the front
               tax_level = "Genus", n_taxa = 15, 
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.gen.18s) +
  #coord_flip() +
  labs(title = "WD",
       y = "Relative abundance") +
  my_custom_theme4()
#water_gen.dis2

### Combine----
comb.gen.18s <- cowplot::plot_grid(water_gen.ndis,
                                   water_gen.dis,
                                   water_gen.ndis2 + theme(legend.position = "none"),
                                   water_gen.dis2 + theme(legend.position = "none"),
                                   ncol = 4,
                                   rel_widths = c(2.5,2,4,4))
comb.gen.18s

## legend are not common is diseased and non-diseased. So I added one now and then later will add the remaining taxa using inkscape
#legend.gen.18s <- get_legend(water_gen.ndis)
#
#comb.gen.18s2 <- cowplot::plot_grid(comb.gen.18s,
#                                    legend.gen.18s,
#                                    rel_widths = c(2, 0.5),
#                                    nrow = 1)
#comb.gen.18s2

# Combine all----
comb.barplot <- cowplot::plot_grid(comb.phy,
                                   comb.gen,
                                   comb.div.18s,
                                   comb.gen.18s,
                                   labels = "AUTO",
                                   align = "v",
                                   ncol = 1)
comb.barplot

# Add one legend with the biggest name to add the space. Then in inkscape, finalise the legends and add them
legend6 <- get_legend(p6.pond.g)

# Final plot----
comb.barplot2 <- cowplot::plot_grid(comb.barplot,
                                    legend6,
                                    ncol = 2,
                                    rel_widths = c(2, 0.35))
comb.barplot2
# saved as 2100*2970, then combine all the legend for genus following code below
## keep only unique genus and add this legend to the combine plot
# This one is mostly similar as the previous one. As it will take some time to finalise in inkscape, just leave as the current figure. later add this one



#### Prokaryotes genus legends----
legend1 <- get_legend(p1.pond.g)
legend2 <- get_legend(p2.pond.g)
legend3 <- get_legend(p3.pond.g)
legend4 <- get_legend(p4.pond.g)
legend5 <- get_legend(p5.pond.g)
legend6 <- get_legend(p6.pond.g)

combined_legend <- plot_grid(legend1, 
                             legend2, 
                             legend3, 
                             legend4, 
                             legend5, 
                             legend6, ncol = 6)
combined_legend

#### Microeukaryotes genus legends----
legend1.18s <- get_legend(water_gen.ndis2)
legend2.18s <- get_legend(water_gen.dis2)


combined_legend.18s <- plot_grid(legend1.18s, 
                                 legend2.18s, 
                                 ncol = 2)
combined_legend.18s

## edit in inkscape. to remove separate border and put a single border for 
# all groups, use 1.33 stroke width with 50% transparency

# Figures for legend----
# to use the legend, plot with taxa that are present in all sample types
# Alphabetical sorting of top 10 phyla
topTaxa.phy <- pseq_fixed %>%
  tax_top(n = 11, rank = "Phylum") %>%
  sort() # this makes them alphabetical
# Plot
top.phy <- pseq_fixed %>%
  tax_sort(by = sum, at = "Phylum") %>% # this orders all genera by abundance
  ps_select(Sample_type, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Sample_type == "Gill_swab") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topTaxa.phy, # this brings the named taxa to the front
               tax_level = "Phylum", n_taxa = 11,
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = phyla.colors) +
  #coord_flip() +
  labs(title = "Gill (non-diseased)",
       y = "Relative abundance",
       x = "Pond name") +
  theme(plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.8, "cm")) + # Adjust the height of legend keys
  guides(fill = guide_legend(ncol = 1)) #+ # Adjust the number of rows in the legend
top.phy # This one is necessary for the legend

# Rank = Division----
# for legend
# set up for alphabetical sorting
topTaxa.div.18s <- ps.18s %>%
  tax_top(n = 11, rank = "Division") %>%
  sort() # this makes them alphabetical

water_div.all <-  ps.18s %>%
  tax_sort(by = sum, at = "Division") %>% # this orders all genera by abundance
  ps_select(Reported_disease, Pond_name) %>% # avoids lots of phyloseq::merge_samples warnings
  #ps_filter(Reported_disease == "Diseased") %>%
  phyloseq::merge_samples(group = "Pond_name") %>%
  comp_barplot(tax_order = topTaxa.div.18s, # this brings the named taxa to the front
               tax_level = "Division", n_taxa = 11, # RA > 2%
               tax_sort = "name",
               sample_order = "asis", 
               merge_other = TRUE, other_name = "Other",
               bar_width = 0.9,
               bar_outline_colour = NA,
               palette = myPal.div.18s) +
  #coord_flip() +
  labs(title = "WD",
       y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_blank(), #element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_text(vjust = 0.5, hjust = 0.5, size = 14, face = "bold", margin = margin(t = 5)),  
        axis.title.y = element_blank(), #element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
        legend.position = "right",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_legend(ncol = 1))  # Adjust the number of rows in the legend
