# date: 20250616 (updated from 20241101)
# Author: Sanjit Debnath

# with these code, I plotted alpha-beta diversity, shared taxa (asv) and differentially abundant genera by lefse
# using venn diagram for prokaryotes across different sample types
# These all woks, just need to run these codes, save the final plot,
# just a little modification in inkscape. not using the chao1 result but plooting in case needed

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
#install.packages("ggpubr")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(stringr); packageVersion("stringr") # to wrap text
library(gridExtra); packageVersion("gridExtra")
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
library(dunn.test); packageVersion("dunn.test")

#to make the variable using my desired color, set variables desired color
sample.colors <- c("Gill_swab"= "#4363d8", "Skin_swab" = "#8F7700FF", "Pond_water"= "#5EC747",
                   "Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5")

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

# Prokaryotes----
## A. Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy <- rarefy_even_depth(pseq, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#8187OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates <- estimate_richness(ps_rarefy, measures = c("Chao1", "Shannon"))
alpha_estimates <- cbind(alpha_estimates, sample_data(pseq))

# Reorder the factor levels
alpha_estimates$Sample_type <- factor(alpha_estimates$Sample_type, 
                                      levels = c("Gill_swab", "Skin_swab", "Pond_water"))

### Richness----
# Chao1: An estimator of total species richness that accounts for the presence of rare species.
# Create the plot
st.c <- ggplot(data = alpha_estimates, aes(y = Chao1, x = Sample_type)) +
  geom_boxplot(aes(color = Sample_type), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_type), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = sample.colors, labels = c("Gill", "Skin", "Water")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_blank(),#text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    #axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    #axis.line.y = element_line(color = "black")   # Optional: Add border to panel
    axis.line = element_line(color = "black")
    ) +
  scale_x_discrete(labels = c("Gill_swab" = "Gill", 
                              "Skin_swab" = "Skin", 
                              "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Sample type",  # Customize legend title
       y = "Chao1",  # Add y-axis title to the plot
       x = "Sample type",
       title = "Chao1") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend

#### Perform Kw test----
kruskal.test(Chao1 ~ Sample_type, data = alpha_estimates)

##### perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.st.c <- dunn_test(Chao1 ~ Sample_type, data = alpha_estimates, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.st.c2 <- stat.test.st.c %>% add_xy_position(x = "Sample_type")
st.c2 <- st.c + stat_pvalue_manual(stat.test.st.c2, label = "p.adj.signif", 
                                   step.increase = 0.05, tip.length = 0.005, 
                                   size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
st.c2

### Diversity----
# Shannon Diversity Index: Incorporates both richness and evenness of species. It considers the number of species and their relative abundances.
st.s <- ggplot(data = alpha_estimates, aes(y = Shannon, x = Sample_type)) +
  geom_boxplot(aes(color = Sample_type), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_type), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = sample.colors, labels = c("Gill", "Skin", "Water")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_blank(),#text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    #axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    #axis.line.y = element_line(color = "black")   # Optional: Add border to panel
    axis.line = element_line(color = "black")
    ) +
  scale_x_discrete(labels = c("Gill_swab" = "Gill", 
                              "Skin_swab" = "Skin", 
                              "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Sample type",  # Customize legend title
       y = "Shannon",  # Add y-axis title to the plot
       x = "Sample type",
       title = "Shannon") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend

#### Perform Kw test----
kruskal.test(Shannon ~ Sample_type, data = alpha_estimates)

##### perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.st.s <- dunn_test(Shannon ~ Sample_type, data = alpha_estimates, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.st.s2 <- stat.test.st.s %>% add_xy_position(x = "Sample_type")
st.s2 <- st.s + stat_pvalue_manual(stat.test.st.s2, label = "p.adj.signif", 
                                   step.increase = 0.05, tip.length = 0.005, 
                                   size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
st.s2

### Phylogenetic diversity----
# we need to separately calculte the phylogenetic diversity
# follow: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html#diversities
library(microbiome); packageVersion("microbiome") # data analysis and visualisation
library(DT); packageVersion("DT") # interactive tables in html and markdown
library(data.table); packageVersion("data.table") # alternative to data.frame
library(picante); packageVersion("picante") # for PD

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

# Reorder the factor levels
alpha_estimate_withPD$Sample_type <- factor(alpha_estimate_withPD$Sample_type, 
                                            levels = c("Gill_swab", "Skin_swab", "Pond_water"))


#### Plot Phylogenetic diversity----
st.pd <- ggplot(data = alpha_estimate_withPD, aes(y = Phylogenetic_Diversity, x = Sample_type)) +
  geom_boxplot(aes(color = Sample_type), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_type), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = sample.colors, labels = c("Gill", "Skin", "Water")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_blank(),#text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_blank(),#text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    #axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    #axis.line.y = element_line(color = "black")   # Optional: Add border to panel
    axis.line = element_line(color = "black")
    ) +
  scale_x_discrete(labels = c("Gill_swab" = "Gill", 
                              "Skin_swab" = "Skin", 
                              "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Sample type",  # Customize legend title
       y = "PD",  # Add y-axis title to the plot
       x = "Sample type",
       title = "Faith's PD") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
st.pd

##### Perform Kw test----
kruskal.test(Phylogenetic_Diversity ~ Sample_type, data = alpha_estimate_withPD)

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.st.pd <- dunn_test(Phylogenetic_Diversity ~ Sample_type, data = alpha_estimate_withPD, p.adjust.method = "BH")

# add the value on the plotBox plot
stat.test.st.pd <- stat.test.st.pd %>% add_xy_position(x = "Sample_type")
st.pd2 <- st.pd + stat_pvalue_manual(stat.test.st.pd, label = "p.adj.signif", 
                                     step.increase = 0.05, tip.length = 0.005, 
                                     size = 8)  # Adjust the size of the asterisk here
#hide.ns = TRUE)
st.pd2

## B. Beta diversity (Bray-Curtis)----
st_bray <- pseq %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") 

# Sampling_month
st.beta <- st_bray %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Sample_type", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Sample_type)
  ) +
  scale_color_manual(values = sample.colors, labels = c("Gill", "Water", "Skin")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
    legend.position = "bottom",
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
    axis.line = element_line(color = "black")  # Correctly specify axis line color for both axes
  ) +
  #scale_x_discrete(labels = c("Gill_swab" = "Gill", 
  #                            "Skin_swab" = "Skin", 
  #                            "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Sample type",  # Customize legend title
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
st.beta

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_st <- sample_data(pseq) %>%
  data.frame() %>%
  tibble()
bray_dist.st = phyloseq::distance(pseq, method="bray")
PERM_st_bray <- adonis2(bray_dist.st ~ Sample_type, data = metadata_st)

# to add test result on the plot
st.beta2 <- st.beta + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_st_bray["Model","R2"], 3), 
                                             ", p = ", PERM_st_bray["Model", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
st.beta2

# pairwise adonis
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

#### Pairwise adonis----
pairwise_adonis_st <- pairwise.adonis(bray_dist.st, metadata_st$Sample_type, p.adjust.m = "BH")
pairwise_adonis_st

# Add both PERMANOVA and pairwise adonis results on the plot
st.beta3 <-st.beta + 
  annotate(
    geom = "label",
    label = paste(
      "PERMANOVA: R² = ", round(PERM_st_bray["Model", "R2"], 3), 
      ", p =", PERM_st_bray["Model", "Pr(>F)"], "\n",
      "Skin vs Gill: R² =", round(pairwise_adonis_st[1, "R2"], 3),
      ", FDR =", pairwise_adonis_st[1, "p.adjusted"], "\n",
      "Skin vs Water: R² =", round(pairwise_adonis_st[2, "R2"], 3),
      ", FDR =", pairwise_adonis_st[2, "p.adjusted"], "\n",
      "Gill vs Water: R² =", round(pairwise_adonis_st[3, "R2"], 3),
      ", FDR =", pairwise_adonis_st[3, "p.adjusted"]
    ),
    x = Inf, y = Inf, hjust = 1.05, vjust = 0.8
  )
st.beta3

## C. Shared ASVs----
# load libraries
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess") # needed to make the vennlist
library(VennDiagram); packageVersion("VennDiagram") # for venn diagram
library(gridExtra); packageVersion("gridExtra")# need to visualise venn diagram properly

### Sample types----
vennlist.16s.st <- get_vennlist(obj = pseq, factorNames = "Sample_type")

# to rename the name of the lists
names(vennlist.16s.st) <- c("Gill", "Water", "Skin")
#vennlist.16s.st

# Plot
venn.16s.st <- venn.diagram(vennlist.16s.st,
                            filename = NULL, # can't make it work
                            disable.logging = TRUE, # stop log file output in dictionary
                            main = "Shared bacterial ASVs", 
                            main.pos = c(0.5, 1.0), 
                            main.fontface = "bold",
                            main.fontfamily = "serif", main.col = "black",
                            main.cex = 2, main.just = c(0.5, 1),
                            fill = c("#4363d8", "#5EC747", "#8F7700FF"),
                            alpha = 0.5, 
                            fontfamily = "serif",
                            fontface = "bold",
                            cex = 1.6,
                            cat.col = "black",
                            cat.cex = 1.6,
                            cat.default.pos = "outer",
                            cat.dist = 0.03,
                            #cat.pos = 0,
                            print.mode = c("raw","percent"), # to add percentages of the shared taxa
                            margin = 0.08, 
                            lwd = 3,
                            lty ='blank',
                            auto_scale = FALSE, # it resize the circle
                            imagetype = "svg")
# Create a new plot
plot.new()

# Plot
venn.16s.st.p <- grid.arrange(venn.16s.st) # to view the plot properly; it worked before but not it's not working
# try this
library(grid); packageVersion("grid")
library(gridExtra); packageVersion("gridExtra")

# Combine the grobs into one
venn_grob <- grobTree(venn.16s.st)

# Now arrange it
venn.16s.st.p <- grid.arrange(venn_grob)

## D. LEfSe----
# Keep it at last, but during combining, put it first
library(microbiomeMarker); packageVersion("microbiomeMarker")
# https://yiluheihei.github.io/microbiomeMarker/articles/microbiomeMarker-vignette.html#introduction

### Genus rank, lda 3.5----
lefse.st3.5 <- run_lefse(
  pseq,
  group = "Sample_type",
  taxa_rank = "Genus", # to plot cladogram, we need to avoid this
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 3.5)
lefse.st3.5 # 25 microbiome markers with 5 variables (with p<0.01, 14 microbiome markers with 5 variables) 2240 taxa and  298 samples
#saveRDS(lefse.st3.5, "lefse3.5_microbiomeMarker_Sample_type_p0.05.rds")
lefse.st3.5 <- readRDS("lefse3.5_microbiomeMarker_Sample_type_p0.05.rds")
lefse.st3.5 # 25 microbiome markers with 5 variables 2240 taxa and  298 samples

# if want to adjust p value
#marker_table(lefse.st)$padj <- p.adjust(marker_table(lefse.st)$pvalue, method = "BH")

#### Bar plot for effect size----
da.st3.5 <- plot_ef_bar(lefse.st3.5) # this plot is not bad
da.st3.5 + scale_fill_manual(values = c("Gill_swab"= "#4363d8", "Skin_swab" = "#8F7700FF", "Pond_water"= "#5EC747"))

#### Make the barplot better with ggplot----
# Extract the LEfSe results and convert to a tibble
lefse.st3.5_results <- as_tibble(marker_table(lefse.st3.5))

# Check the structure of the LEfSe results
head(lefse.st3.5_results)

# Convert enrich_group to a factor and specify the order
lefse.st3.5_results$enrich_group <- factor(lefse.st3.5_results$enrich_group, levels = c("Gill_swab", "Skin_swab", "Pond_water"))

# Reorder features within each enriched group based on LDA score
lefse.st3.5_results <- lefse.st3.5_results %>%
  arrange(enrich_group, ef_lda) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

# Custom plot using ggplot2
lefse.st3.5.plot <- ggplot(lefse.st3.5_results, aes(x = feature, y = ef_lda, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Gill_swab"= "#4363d8", "Skin_swab" = "#8F7700FF", "Pond_water"= "#5EC747"),
                    labels = c("Gill", "Skin", "Water")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 15)), # Center title and increase margin
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(size = 14, face = "italic", margin = margin(r = 5)), # Italicize y-axis labels
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14)) +
  labs(fill = "Sample Type",
       #title = "Sample Types",
       x = "Enriched taxa",
       y = "LDA Score (log10)") +
  guides(color = guide_legend(nrow = 1))
lefse.st3.5.plot

# Combine plots----
fig.2 <- cowplot::plot_grid(lefse.st3.5.plot,
                            st.c2,
                            st.s2,
                            st.pd2,
                            st.beta3 + theme(legend.position = "none"),
                            venn.16s.st.p,
                            #venn.core.plot,
                            labels = "AUTO")
fig.2
# Save as 2000*1200


## Calculation of ASV across samples
#skin <- pseq %>% 
#  ps_filter(Sample_type == "Skin_swab")
#skin # 29594 taxa and 120 samples
#
#gill <- pseq %>% 
#  ps_filter(Sample_type == "Gill_swab")
#gill # 31264 taxa and 118 samples
#
#water <- pseq %>% 
#  ps_filter(Sample_type == "Pond_water")
#water #  18017 taxa and 60 samples
#g <- 31264
#s <- 29594
#w <- 18017
#c <- 16403
#
#gs <- (16403+12410)
#gs # 28813
#gw <- (16403+1300)
#gw #17703
#sw <- (16403+232)
#sw #16635
#
#wp <- (c/w*100)
#wp #  91.04179% of water asv is shared across three sample type
#gp <- (gw/w*100)
#gp # 98.2572% of water asv is shared with gill
#sp <- (sw/w*100)
#sp # 92.32947% of water asv is shared with skin
#gsp1 <- (gs/g*100)
#gsp1 # 92.16031
#gsp2 <- (gs/s*100)
#gsp2 # 97.36095



