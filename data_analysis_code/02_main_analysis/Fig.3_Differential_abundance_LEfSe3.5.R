# 25/05/2024 (updated from 08/07/2024)
# Author: Sanjit Debnath

# Bacterial signatures in different samples
#  Differentially abundance taxa
## most commonly used DA methods can be divided into three main categories: 
# a) simple statistical tests; b) RNA-seq based methods; c) metagenomic based methods
# we are not trying simple statistical tests

# deseq2 and edger shows the figure for two variables as I want. for more than 
# two variable, all looks same. But lefse is the most popular methods recently,
# Also for biomarker identification, LEfSe better
## based on all of the above method, i will follow lefse

# If i use all rank which is necessary for the cladogram, I have lots of taxa.
# So i decided to do the lefse at genus rank and use a lda cutoof of 3,
# so any taxa less tha 3 lda score will be not shown
# instead of showing all together, i think if i compare disease vs non-diseased
# across samples, that would be better

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(microbiomeMarker); packageVersion("microbiomeMarker")
# https://yiluheihei.github.io/microbiomeMarker/articles/microbiomeMarker-vignette.html#introduction
library(tibble); packageVersion("tibble")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

# Prokaryotes
#load phyloseq object
ps <- readRDS("phyloseq_FSS_BD_metadata_v7_16S_abd30_20241016.rds")
ps # 33214 taxa and 334 samples

# Remove others except tilapia, combine sample_type and disease state, rename sample type 
# All Samples----
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

# Gill----
gill <- pseq %>% 
  ps_filter(Sample_type == "Gill_swab")
gill # 31264 taxa and 118 samples

## LEfSe----
### Genus rank----
gill_lefse3.5 <- run_lefse(
  gill,
  group = "Reported_disease",
  taxa_rank = "Genus", # to plot cladogram, we need to avoid this
  wilcoxon_cutoff = 0.05, 
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 3.5)
gill_lefse3.5 # # 31 microbiome markers with 5 variables ( with p < 0.01, 23 microbiome markers with 5 variables) 2216 taxa and  118 samples

# It takes some time, so I saved it as rds, so that next time doesn't need to run again
#saveRDS(gill_lefse3.5, "lefse3.5_microbiomeMarker_gill_p0.05.rds")
gill_lefse3.5 <- readRDS("lefse3.5_microbiomeMarker_gill_p0.05.rds")
gill_lefse3.5 # 31 microbiome markers with 5 variables 2216 taxa and  118 samples

#### Bar plot for effect size----
plot_ef_bar(gill_lefse3.5)

#### Make it better with ggplot----
# Extract the LEfSe results and convert to a tibble
gill_lefse3.5_results <- as_tibble(marker_table(gill_lefse3.5))

# Check the structure of the LEfSe results
head(gill_lefse3.5_results)

# Add a column to indicate the direction (positive for Diseased, negative for Non-diseased)
gill_lefse3.5_results <- gill_lefse3.5_results %>%
  mutate(direction = ifelse(enrich_group == "Diseased", "Positive", "Negative"),
         lda_score = ifelse(enrich_group == "Diseased", ef_lda, -ef_lda))

# Custom plot using ggplot2
gill_lefse3.5.plot <- ggplot(gill_lefse3.5_results, aes(x = reorder(feature, lda_score), y = lda_score, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, face = "italic", margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_blank(),#text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14)) +
  labs(fill = "Enriched group",
       title = "Gill",
       x = "Enriched taxa",
       y = "LDA Score (log10)") +
  guides(color = guide_legend(nrow = 1))+
  ylim(-5,5)
gill_lefse3.5.plot

#### Abundance box plot----
gill_lefse3.5_abd <- plot_abundance(gill_lefse3.5, group = "Reported_disease")
gill_lefse3.5_abd + scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

#### Dot plot----
gill_lefse3.5_dot <- plot_ef_dot(gill_lefse3.5)
gill_lefse3.5_dot + scale_color_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

# Skin----
skin <- pseq %>% 
  ps_filter(Sample_type == "Skin_swab")
skin # 29496 taxa and 120 samples

## LEfSe----
### Genus rank----
skin_lefse3.5 <- run_lefse(
  skin,
  group = "Reported_disease",
  taxa_rank = "Genus", # to plot cladogram, we need to avoid this
  wilcoxon_cutoff = 0.05, # 30 microbiome markers with 5 variable
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 3.5)
skin_lefse3.5 # 30 microbiome markers with 5 variable (with p<0.01, 24 microbiome markers with 5 variables) 2188 taxa and  120 samples

#saveRDS(skin_lefse3.5, "lefse3.5_microbiomeMarker_skin_p0.05.rds")
skin_lefse3.5 <- readRDS("lefse3.5_microbiomeMarker_skin_p0.05.rds")
skin_lefse3.5 # 30 microbiome markers with 5 variables 2188 taxa and  120 samples

#### Bar plot for effect size----
plot_ef_bar(skin_lefse3.5)

#### Make it better with ggplot----
# Extract the LEfSe results and convert to a tibble
skin_lefse3.5_results <- as_tibble(marker_table(skin_lefse3.5))

# Check the structure of the LEfSe results
head(skin_lefse3.5_results)

# Add a column to indicate the direction (positive for Diseased, negative for Non-diseased)
skin_lefse3.5_results <- skin_lefse3.5_results %>%
  mutate(direction = ifelse(enrich_group == "Diseased", "Positive", "Negative"),
         lda_score = ifelse(enrich_group == "Diseased", ef_lda, -ef_lda))

# Custom plot using ggplot2
skin_lefse3.5.plot <- ggplot(skin_lefse3.5_results, aes(x = reorder(feature, lda_score), y = lda_score, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, face = "italic", margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_blank(),#text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14)) +
  labs(fill = "Enriched group",
       title = "Skin",
       x = "Enriched taxa",
       y = "LDA Score (log10)") +
  guides(color = guide_legend(nrow = 1))+
  ylim(-5,5)
skin_lefse3.5.plot

#### Abundance box plot----
skin_lefse3.5_abd <- plot_abundance(skin_lefse3.5, group = "Reported_disease")
skin_lefse3.5_abd + scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

#### Dot plot----
skin_lefse3.5_dot <- plot_ef_dot(skin_lefse3.5)
skin_lefse3.5_dot + scale_color_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

# Water----
water <- pseq %>% 
  ps_filter(Sample_type == "Pond_water")
water # 17973 taxa and 60 samples

## LEfSe----
### Genus rank----
water_lefse3.5 <- run_lefse(
  water,
  group = "Reported_disease",
  taxa_rank = "Genus", # to plot cladogram, we need to avoid this
  wilcoxon_cutoff = 0.05, 
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 3.5)
water_lefse3.5 # 28 microbiome markers with 5 variables (with p<0.01, 22 microbiome markers with 5 variables) 1712 taxa and 60 samples 
#saveRDS(water_lefse3.5, "lefse3.5_microbiomeMarker_water_p0.05.rds")
water_lefse3.5 <- readRDS("lefse3.5_microbiomeMarker_water_p0.05.rds")
water_lefse3.5 # 28 microbiome markers with 5 variables 1712 taxa and 60 samples 

#### Bar plot for effect size----
plot_ef_bar(water_lefse3.5)

#### Make it better with ggplot----
# Extract the LEfSe results and convert to a tibble
water_lefse3.5_results <- as_tibble(marker_table(water_lefse3.5))

# Check the structure of the LEfSe results
head(water_lefse3.5_results)

# Add a column to indicate the direction (positive for Diseased, negative for Non-diseased)
water_lefse3.5_results <- water_lefse3.5_results %>%
  mutate(direction = ifelse(enrich_group == "Diseased", "Positive", "Negative"),
         lda_score = ifelse(enrich_group == "Diseased", ef_lda, -ef_lda))

# Custom plot using ggplot2
water_lefse3.5.plot <- ggplot(water_lefse3.5_results, aes(x = reorder(feature, lda_score), y = lda_score, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5")) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, face = "italic", margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_blank(),#text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14))+
  labs(fill = "Enriched group",
       title = "Water",
       x = "Enriched taxa",
       y = "LDA Score (log10)") +
  guides(color = guide_legend(nrow = 1)) +
  ylim(-5,5)
water_lefse3.5.plot

#### Abundance box plot----
water_lefse3.5_abd <- plot_abundance(water_lefse3.5, group = "Reported_disease")
water_lefse3.5_abd + scale_fill_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

#### Dot plot----
water_lefse3.5_dot <- plot_ef_dot(water_lefse3.5)
water_lefse3.5_dot + scale_color_manual(values = c("Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5"))
# but I'm not using this plot

# Combine plots----
comb.bar.3.5 <- cowplot::plot_grid(#lefse.st.plot + theme(legend.position = "none"),
  #lefse.st3.5.plot,
  gill_lefse3.5.plot,# + theme(legend.position = "none"),
  skin_lefse3.5.plot,# + theme(legend.position = "none"),
  water_lefse3.5.plot,#+ theme(legend.position = "none"),
  labels = c("A", "B", "C"),
  ncol = 3)
comb.bar.3.5
# save as 2000*1200


