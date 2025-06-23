# 05/11/2024
# Author: Sanjit Debnath

# Script for Co_occurrence network
# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(file2meco); packageVersion("file2meco") # to convent into microtable
library(microeco); packageVersion("microeco")
library(WGCNA); packageVersion("WGCNA")
library(magrittr); packageVersion("magrittr")
library(cowplot); packageVersion("cowplot")

## Setup working dictionary first
#set working dictionary
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

# Prokaryotes
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

# Filter out low-abundnace taxa----
# Step 1: Transform counts to relative abundances
pseq_relab <- transform_sample_counts(pseq, function(x) x / sum(x))

# Step 2: Identify taxa with relative abundance > 1% in any sample
taxa_above_1pct <- taxa_names(filter_taxa(pseq_relab, function(x) max(x) > 0.01, TRUE))

# Step 3: Identify taxa present in more than 30% of samples
taxa_above_30pct <- taxa_names(filter_taxa(pseq_relab, function(x) sum(x > 0) > 0.3 * length(x), TRUE))

# Step 4: Find the intersection of both conditions
taxa_to_keep <- intersect(taxa_above_1pct, taxa_above_30pct)

# Subset the phyloseq object to retain only the desired taxa
pseq_filtered <- prune_taxa(taxa_to_keep, pseq)

# Check the filtered phyloseq object
pseq_filtered # 300 taxa and 298 samples

# Make a barplot so that the legend can be used later
phyla.colors <- c(
  "Proteobacteria" = "#00e0ffff", "Bacteroidota" = "#ff7200ff", "Actinobacteriota" = "#ff82ffff", 
  "Firmicutes" = "#68d500ff","Planctomycetota" = "#025d43ff", "Verrucomicrobiota" = "#70446eff",
  "Chloroflexi" = "#e9a400ff", "Cyanobacteria" = "#78abffff", "Acidobacteriota" = "#85481aff",
  "Fusobacteriota" = "#00ceb9ff", "Patescibacteria" = "#00cb5dff", "Latescibacterota" = "#ff66b7ff",
  "Desulfobacterota" = "#8aa434ff","Margulisbacteria" = "#e1bdb2ff","Deinococcota" = "#009cc5ff", "Other" = "lightgrey") 
#comp_barplot(pseq_filtered, tax_level = "Phylum", n_taxa = 15, palette = phyla.colors)


# To plot network if we aggregate at Genus level and then use Phylum as classification, we will have both Phylum and genus on the plot
## Aggregate at Genus rank----
pseq_filtered.gen <- tax_glom(pseq_filtered, taxrank = "Genus")
pseq_filtered.gen # 156 taxa and 298 samples

# Subset each sample disease type
# Gill non-disease----
gnd.glom <- pseq_filtered.gen %>% ps_filter(Sample_Disease == "Gill_Non-diseased")
gnd.glom # 156 taxa and 58 samples

# convert to microtable
dataset.gnd.glom <- phyloseq2meco(gnd.glom) 

## Make a network ----
# using WGCNA package is faster and optimized for larger dataset and often provides a more robust and efficient network structure, especially when analyzing large microbiome datasets
net.gnd.glom <- trans_network$new(dataset = dataset.gnd.glom, 
                                  cor_method = "spearman", # non-parametric, doesn't need check normality
                                  use_WGCNA_pearson_spearman = TRUE) #
                                  #filter_thres = 0.01) # Relative abundance threshold, mean RA 0.1% across all samples
# After filtering, 156 features are remained

### Construct network; require igraph package----
net.gnd.glom$cal_network(COR_p_thres = 0.001, # the p value threshold
                         COR_p_adjust = "BH", # p value adjustment method
                         COR_cut = 0.6, # correlation coefficient threshold (rho)
                         add_taxa_name = c("Phylum", "Genus"), #one or more taxonomic rank name
                         COR_optimization = TRUE)
net.gnd.glom$res_network # IGRAPH c42f0fd UNW- 149 1269

### Partition network in modules----
# Modules are clusters of ASVs/OTUs that co-occur more frequently together across samples than with others.
## These groups might represent functionally related taxa, shared ecological niches, or taxa that interact with one another within the microbial ecosystem.
# invoke igraph cluster_fast_greedy function for this undirected network 
net.gnd.glom$cal_module(method = "cluster_fast_greedy",
                        module_name_prefix = "Module") # Totally, 7 modules are identified
#net.gnd.glom$plot_network(method = "igraph", layout = layout_with_kk)

### Save calculated network----
#net.gnd.glom$save_network(filepath = "network_gnd_glom2_with_WGCNA_p001_20241105.gexf")
# now go to Gephi software and plot the network for module/cluster

#### Network attributes####
net.gnd.glom$cal_network_attr()
net.gnd.glom.attr <- net.gnd.glom$res_network_attr
net.gnd.glom.attr

# To save this, add the rownames as a new column
net.gnd.glom.attr$Features <- rownames(net.gnd.glom.attr)
# save as excel file
#writexl::write_xlsx(net.gnd.glom.attr, "network_gnd_glom2_attributes.xlsx")

#### get node properties
net.gnd.glom$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
net.gnd.glom$get_edge_table()
# return t1$res_edge_table 
net.gnd.glom$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# Number of edges
gnd.edges <- table(net.gnd.glom$res_edge_table$label)
print(gnd.edges)

# eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module
net.gnd.glom$cal_eigen()
# return t1$res_eigen

##### Save this attribute table####
gnd.attibute.table <- net.gnd.glom$res_node_table
#writexl::write_xlsx(gnd.attibute.table, "net_gnd_glom_attribute_table_20241106.xlsx")

# Gill disease----
gd.glom <- pseq_filtered.gen %>% ps_filter(Sample_Disease == "Gill_Diseased")
gd.glom # 156 taxa and 60 samples

# convert to microtable
dataset.gd.glom <- phyloseq2meco(gd.glom) 

## Make a network ----
# using WGCNA package is faster and optimized for larger dataset and often provides a more robust and efficient network structure, especially when analyzing large microbiome datasets
net.gd.glom <- trans_network$new(dataset = dataset.gd.glom, 
                                 cor_method = "spearman", # non-parametric, doesn't need check normality
                                 use_WGCNA_pearson_spearman = TRUE) #
                                 #filter_thres = 0.001) # Relative abundance threshold, mean RA 0.1% across all samples
# After filtering, 156 features are remained

### Construct network; require igraph package----
net.gd.glom$cal_network(COR_p_thres = 0.001, # the p value threshold
                        COR_p_adjust = "BH", # p value adjustment method
                        COR_cut = 0.6, # correlation coefficient threshold (rho)
                        add_taxa_name = c("Phylum", "Genus"), #one or more taxonomic rank name
                        COR_optimization = TRUE)
net.gd.glom$res_network # IGRAPH 24513c8 UNW- 151 803

### Partition network in modules----
# Modules are clusters of ASVs/OTUs that co-occur more frequently together across samples than with others.
## These groups might represent functionally related taxa, shared ecological niches, or taxa that interact with one another within the microbial ecosystem.
# invoke igraph cluster_fast_greedy function for this undirected network 
net.gd.glom$cal_module(method = "cluster_fast_greedy",
                       module_name_prefix = "Module") # Totally, 9 modules are identified
# Save the calculated network to plot using gephi
#net.gd.glom$save_network(filepath = "network_gd_glom2_with_WGCNA_p001_20241105.gexf")
# now go to Gephi software and plot the network for module/cluster

#### Network attributes####
net.gd.glom$cal_network_attr()
net.gd.glom.attr <- net.gd.glom$res_network_attr
net.gd.glom.attr

# Add the rownames as a new column
net.gd.glom.attr$Features <- rownames(net.gd.glom.attr)
# save as excel file
#writexl::write_xlsx(net.gd.glom.attr, "network_gd_glom2_attributes.xlsx")

#### Node properties####
net.gd.glom$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# Access and save the result from the object's res_node_table attribute
node_table <- net.gd.glom$res_node_table # 127 nodes

# View the head of the node table
head(node_table)

# get edge properties
net.gd.glom$get_edge_table()
# return t1$res_edge_table 
net.gd.glom$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# Number of edges
gd.edges <- table(net.gd.glom$res_edge_table$label)
print(gd.edges)

# eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module
net.gd.glom$cal_eigen()
# return t1$res_eigen

##### Save this attribute table####
gd.attibute.table <- net.gd.glom$res_node_table
#writexl::write_xlsx(gd.attibute.table, "net_gd_glom_attribute_table_20241106.xlsx")


# Skin non-disease----
snd.glom <- pseq_filtered.gen %>% ps_filter(Sample_Disease == "Skin_Non-diseased")
snd.glom # 156 taxa and 58 samples

# convert to microtable
dataset.snd.glom <- phyloseq2meco(snd.glom) 

## Make a network ----
# using WGCNA package is faster and optimized for larger dataset and often provides a more robust and efficient network structure, especially when analyzing large microbiome datasets
net.snd.glom <- trans_network$new(dataset = dataset.snd.glom, 
                                  cor_method = "spearman", # non-parametric, doesn't need check normality
                                  use_WGCNA_pearson_spearman = TRUE) #
                                  #filter_thres = 0.001) # Relative abundance threshold, mean RA 0.1% across all samples
# After filtering, 156 features are remained

### Construct network; require igraph package----
net.snd.glom$cal_network(COR_p_thres = 0.001, # the p value threshold
                         COR_p_adjust = "BH", # p value adjustment method
                         COR_cut = 0.6, # correlation coefficient threshold (rho)
                         add_taxa_name = c("Phylum", "Genus"), #one or more taxonomic rank name
                         COR_optimization = TRUE)
net.snd.glom$res_network # IGRAPH 3f3d99c UNW- 140 835

### Partition network in modules----
# Modules are clusters of ASVs/OTUs that co-occur more frequently together across samples than with others.
## These groups might represent functionally related taxa, shared ecological niches, or taxa that interact with one another within the microbial ecosystem.
# invoke igraph cluster_fast_greedy function for this undirected network 
net.snd.glom$cal_module(method = "cluster_fast_greedy",
                        module_name_prefix = "Module") # Totally, 5 modules are identified
#net.snd.glom$plot_network(method = "igraph", layout = layout_with_kk)

# Save the calculated network to plot using gephi
# require rgexf package to be installed
#net.snd.glom$save_network(filepath = "network_snd_glom2_with_WGCNA_p001_20241105.gexf")
# now go to Gephi software and plot the network for module/cluster

#### Network attributes####
net.snd.glom$cal_network_attr()
net.snd.glom.attr <- net.snd.glom$res_network_attr
net.snd.glom.attr

# Add the rownames as a new column
net.snd.glom.attr$Features <- rownames(net.snd.glom.attr)
# save as excel file
#writexl::write_xlsx(net.snd.glom.attr, "network_snd_glom2_attributes.xlsx")

# get node properties
net.snd.glom$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
net.snd.glom$get_edge_table()
# return t1$res_edge_table 
net.snd.glom$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# Number of edges
snd.edges <- table(net.snd.glom$res_edge_table$label)
print(snd.edges)

# eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module
net.snd.glom$cal_eigen()
# return t1$res_eigen

##### Save this attribute table####
snd.attibute.table <- net.snd.glom$res_node_table
#writexl::write_xlsx(snd.attibute.table, "net_snd_glom_attribute_table_20241106.xlsx")

# Skin disease----
sd.glom <- pseq_filtered.gen %>% ps_filter(Sample_Disease == "Skin_Diseased")
sd.glom # 156 taxa and 62 samples

# convert to microtable
dataset.sd.glom <- phyloseq2meco(sd.glom) 

## Make a network ----
# using WGCNA package is faster and optimized for larger dataset and often provides a more robust and efficient network structure, especially when analyzing large microbiome datasets
net.sd.glom <- trans_network$new(dataset = dataset.sd.glom, 
                                 cor_method = "spearman", # non-parametric, doesn't need check normality
                                 use_WGCNA_pearson_spearman = TRUE) #
                                 #filter_thres = 0.001) # Relative abundance threshold, mean RA 0.1% across all samples
# After filtering, 156 features are remained

### Construct network; require igraph package----
net.sd.glom$cal_network(COR_p_thres = 0.001, # the p value threshold
                        COR_p_adjust = "BH", # p value adjustment method
                        COR_cut = 0.6, # correlation coefficient threshold (rho)
                        add_taxa_name = c("Phylum", "Genus"), #one or more taxonomic rank name
                        COR_optimization = TRUE)
net.sd.glom$res_network # IGRAPH 82985ed UNW- 115 498

### Partition network in modules----
# Modules are clusters of ASVs/OTUs that co-occur more frequently together across samples than with others.
## These groups might represent functionally related taxa, shared ecological niches, or taxa that interact with one another within the microbial ecosystem.
# invoke igraph cluster_fast_greedy function for this undirected network 
net.sd.glom$cal_module(method = "cluster_fast_greedy",
                       module_name_prefix = "Module") # Totally, 16 modules are identified
#net.sd.glom$plot_network(method = "igraph", layout = layout_with_kk)

# Save the calculated network to plot using gephi
# require rgexf package to be installed
#net.sd.glom$save_network(filepath = "network_sd_glom2_with_WGCNA_p001_20241105.gexf")
# now go to Gephi software and plot the network for module/cluster

#### Network attributes####
net.sd.glom$cal_network_attr()
net.sd.glom.attr <- net.sd.glom$res_network_attr
net.sd.glom.attr

# Add the rownames as a new column
net.sd.glom.attr$Features <- rownames(net.sd.glom.attr)
# save as excel file
#writexl::write_xlsx(net.sd.glom.attr, "network_sd_glom2_attributes.xlsx")

#### get node properties####
net.sd.glom$get_node_table(node_roles = TRUE)
# return t1$res_node_table

# get edge properties
net.sd.glom$get_edge_table()
# return t1$res_edge_table 
net.sd.glom$get_adjacency_matrix()
# return t1$res_adjacency_matrix

# Number of edges
sd.edges <- table(net.sd.glom$res_edge_table$label)
print(sd.edges)

# eigengene of a module, i.e. the first principal component of PCA, represents the main variance of the abundance in the species of the module
net.sd.glom$cal_eigen()
# return t1$res_eigen

##### Save this attribute table####
sd.attibute.table <- net.sd.glom$res_node_table
#writexl::write_xlsx(sd.attibute.table, "net_sd_glom_attribute_table_20241106.xlsx")

# combine diseased and non-diseased gill and skin network (Fig. 6)

# Fig.S5----
# Set phyla colours
phyla.colors2 <- c(
  "p__Proteobacteria" = "#00e0ffff", "p__Bacteroidota" = "#ff7200ff", "p__Actinobacteriota" = "#ff82ffff", 
  "p__Firmicutes" = "#68d500ff","p__Planctomycetota" = "#025d43ff", "p__Verrucomicrobiota" = "#70446eff",
  "p__Chloroflexi" = "#e9a400ff", "p__Cyanobacteria" = "#78abffff", "p__Acidobacteriota" = "#85481aff",
  "p__Fusobacteriota" = "#00ceb9ff", "p__Patescibacteria" = "#00cb5dff", "p__Latescibacterota" = "#ff66b7ff",
  "p__Desulfobacterota" = "#8aa434ff","p__Margulisbacteria" = "#e1bdb2ff","p__Deinococcota" = "#009cc5ff", 
  "p__Bacteria Kingdom" = "lightgrey") 

# Keystone species----
# Plot betweenness centrality vs. degree to determine keystone species----
## Non-diseased gill----
gnd.keystone <- ggplot(gnd.attibute.table, aes(x = degree, y = betweenness_centrality, color = Phylum)) +
  geom_point(size = 4) +  # Adjust size as needed
  scale_color_manual(values = phyla.colors2) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x = "Node degree", y = "Betweenness Centrality", title = "Non-diseased gill")

# For interactive plot
library(plotly); packageVersion("plotly")
#ggplotly(gnd.keystone)

## Diseased gill----
gd.keystone <- ggplot(gd.attibute.table, aes(x = degree, y = betweenness_centrality, color = Phylum)) +
  geom_point(size = 4) +  # Adjust size as needed
  scale_color_manual(values = phyla.colors2) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x = "Node degree", y = "Betweenness Centrality", title = "Diseased gill")

# For interactive plot
#ggplotly(gd.keystone)

## Non-diseased skin----
snd.keystone <- ggplot(snd.attibute.table, aes(x = degree, y = betweenness_centrality, color = Phylum)) +
  geom_point(size = 4) +  # Adjust size as needed
  scale_color_manual(values = phyla.colors2) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x = "Node degree", y = "Betweenness Centrality", title = "Non-diseased skin")

# For interactive plot
#ggplotly(snd.keystone)

## Diseased skin----
sd.keystone <- ggplot(sd.attibute.table, aes(x = degree, y = betweenness_centrality, color = Phylum)) +
  geom_point(size = 4) +  # Adjust size as needed
  scale_color_manual(values = phyla.colors2) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x = "Node degree", y = "Betweenness Centrality", title = "Diseased skin")

# For interactive plot
ggplotly(sd.keystone)

### Combine plots####
keystone.taxa <- cowplot::plot_grid(gnd.keystone + theme(legend.position = "none"),
                                    gd.keystone + theme(legend.position = "none"),
                                    snd.keystone + theme(legend.position = "none"),
                                    sd.keystone + theme(legend.position = "none"),
                                    ncol = 2,
                                    labels = "AUTO",
                                    align = "v")
phyla.legend <- get_legend(gnd.keystone)

keystone.taxa <- cowplot::plot_grid(keystone.taxa,
                                    phyla.legend,
                                    rel_widths = c(2, 0.2))
keystone.taxa
# save as 2000*1200, then add the species name Prioritize Nodes with High Betweenness Centrality  as these taxa serve as crucial connectors across different parts of the network
#After selecting high-BC nodes as primary keystone indicators, also include high-degree nodes as secondary keystone indicators

