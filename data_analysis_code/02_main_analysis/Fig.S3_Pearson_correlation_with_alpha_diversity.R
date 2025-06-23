# date: 20240824
# Author: Sanjit Debnath

# With these code, I'm plotting the correlation between different sample types
# These all woks, just need to run these codes, save the final plot,
# just a little modification in inkscape. 
# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(microbiome); packageVersion("microbiome") # data analysis and visualisation
library(DT); packageVersion("DT") # interactive tables in html and markdown
library(data.table); packageVersion("data.table") # alternative to data.frame
library(picante); packageVersion("picante") # for PD

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

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

## A. Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy <- rarefy_even_depth(pseq, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#8187OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates <- estimate_richness(ps_rarefy, measures = c("Chao1", "Shannon"))
alpha_estimates <- cbind(alpha_estimates, sample_data(pseq))

### Phylogenetic diversity----
# we need to separately calculte the phylogenetic diversity
# follow: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html#diversities

#### Calculate phylogenetic diversity----
rarefy.asvtab <- as.data.frame(ps_rarefy@otu_table)
rarefy.tree <- ps_rarefy@phy_tree

# We first need to check if the tree is rooted or not 
ps_rarefy@phy_tree
# Rooted; includes branch lengths

df.pd <- pd(t(rarefy.asvtab), rarefy.tree, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
#datatable(df.pd)

# to plot the PD, we need to add the PD result with the previous alpha divesity dataframe
# get the metadata out 
alpha_estimate_withPD <- meta(alpha_estimates)

# Add the results of PD to this dataframe and then plot.
alpha_estimate_withPD$Phylogenetic_Diversity <- df.pd$PD 

## Subset the dataframe based on sample type----
gill_df <- subset(alpha_estimate_withPD, Sample_type == "Gill_swab")
skin_df <- subset(alpha_estimate_withPD, Sample_type == "Skin_swab")
water_df <- subset(alpha_estimate_withPD, Sample_type == "Pond_water")

#use dplyr summarise and group_by to group by the pond site, and take a mean value for each group
gill_df_mean <- gill_df %>% group_by(Pond_name) %>%
  summarise(Shannon_mean=mean(Shannon), Shannon_sd=sd(Shannon), 
            Chao1_mean=mean(Chao1), Chao1_sd=sd(Chao1),
            PD_mean=mean(Phylogenetic_Diversity), PD_sd=sd(Phylogenetic_Diversity))

skin_df_mean <- skin_df %>% group_by(Pond_name) %>%
  summarise(Shannon_mean=mean(Shannon), Shannon_sd=sd(Shannon), 
            Chao1_mean=mean(Chao1), Chao1_sd=sd(Chao1),
            PD_mean=mean(Phylogenetic_Diversity), PD_sd=sd(Phylogenetic_Diversity))

water_df_mean <- water_df %>% group_by(Pond_name) %>%
  summarise(Shannon_mean=mean(Shannon), Shannon_sd=sd(Shannon), 
            Chao1_mean=mean(Chao1), Chao1_sd=sd(Chao1),
            PD_mean=mean(Phylogenetic_Diversity), PD_sd=sd(Phylogenetic_Diversity))

#create a df for each indices with a single value per pond site, comparing the two sample types
richness_tab <- data.frame("Gill" = gill_df_mean$Chao1_mean,
                           "Skin" = skin_df_mean$Chao1_mean,
                           "Water" = water_df_mean$Chao1_mean)

shannon_tab <- data.frame("Gill" = gill_df_mean$Shannon_mean,
                          "Skin" = skin_df_mean$Shannon_mean,
                          "Water" = water_df_mean$Shannon_mean)

pd_tab <- data.frame("Gill" = gill_df_mean$PD_mean,
                      "Skin" = skin_df_mean$PD_mean,
                      "Water" = water_df_mean$PD_mean)


# Test to decide preason or spearman
shapiro.test(richness_tab$Gill) # p = 0.5122, normal distribution,
shapiro.test(richness_tab$Skin) # p = 0.1376, normal distribution,
shapiro.test(richness_tab$Water) # p = 0.6493, normal distribution. Okay to use pearson

shapiro.test(shannon_tab$Gill) # p = 0.3717, normal distribution
shapiro.test(shannon_tab$Skin) # p = 0.5661, normal distribution
shapiro.test(shannon_tab$Water) # p = 0.8548, normal distribution. Okay to use pearson

shapiro.test(pd_tab$Gill) # p = 0.9913, normal distribution,
shapiro.test(pd_tab$Skin) # p = 0.8378, normal distribution,
shapiro.test(pd_tab$Water) # p = 0.8034, normal distribution. Okay to use pearson

# Plot----
## Chao1 Richness----
p1_rich <- ggscatter(richness_tab, x = "Gill", y = "Skin",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "Tilapia gill", ylab = "Tilapia skin",
                     title = "Chao1 richness")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p1_rich

p2_rich <- ggscatter(richness_tab, y = "Gill", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia gill", xlab = "Pond water",
                     title = "Chao1 richness")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p2_rich

p3_rich <- ggscatter(richness_tab, y = "Skin", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia skin", xlab = "Pond water",
                     title = "Chao1 richness")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p3_rich

## Shannon divesity----
p1_shan <- ggscatter(shannon_tab, x = "Gill", y = "Skin",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "Tilapia gill", ylab = "Tilapia skin",
                     title = "Shannon diversity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p1_shan

p2_shan <- ggscatter(shannon_tab, y = "Gill", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia gill", xlab = "Pond water",
                     title = "Shannon diversity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p2_shan

p3_shan <- ggscatter(shannon_tab, y = "Skin", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia skin", xlab = "Pond water",
                     title = "Shannon diversity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p3_shan

## PD diversity----
p1_pd <- ggscatter(pd_tab, x = "Gill", y = "Skin",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     xlab = "Tilapia gill", ylab = "Tilapia skin",
                     title = "Phylosegenetic diversity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p1_pd

p2_pd <- ggscatter(pd_tab, y = "Gill", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia gill", xlab = "Pond water",
                     title = "Phylosegenetic diversity")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p2_pd

p3_pd <- ggscatter(pd_tab, y = "Skin", x = "Water",
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     ylab = "Tilapia skin", xlab = "Pond water",
                     title = "Phylosegenetic diversity") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
  )
#p3_pd

# Combine plots----
fig.s3 <- cowplot::plot_grid(p1_rich,
                   p1_shan,
                   p1_pd,
                   p2_rich,# + theme(plot.title = element_blank()),
                   p2_shan,# + theme(plot.title = element_blank()),
                   p2_pd,
                   p3_rich,# + theme(plot.title = element_blank()),
                   p3_shan,# + theme(plot.title = element_blank()),
                   p3_pd,
                   align = "v",
                   labels = "AUTO",
                   ncol = 3)
fig.s3
# save as 1500*835, do a little modification to finalise this


