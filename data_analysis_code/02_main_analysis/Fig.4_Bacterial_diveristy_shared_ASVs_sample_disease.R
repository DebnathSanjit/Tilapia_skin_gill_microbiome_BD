# date: 20241101 
# Author: Sanjit Debnath

# with these code, I plotted alpha-beta diversity, shared taxa (asv) 
# using venn diagram for prokaryotes across different diseased and non-diseased sample types
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

#to make the variable using my desired color set variables desired color
sample.colors <- c("Gill_swab"= "#4363d8", "Skin_swab" = "#8F7700FF", "Pond_water"= "#5EC747",
                   "Non-diseased" = "#1ff8ff", "Diseased" = "#b249d5",
                   "Gill_Diseased" = "#F35E5A", "Gill_Non-diseased" = "#9DCC00",
                   "Skin_Diseased" = "#B68A06", "Skin_Non-diseased" = "#7cb9cb",
                   "Water_Non-diseased" = "#619CFF", "Water_Diseased" = "#FDBF6F")

## Setup working dictionary 
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

# Create a named vector for renaming the facet titles
facet_labels <- c("Gill_swab" = "Gill", "Skin_swab" = "Skin", "Pond_water" = "Water")


## Chao1 richness----
sd.c <- ggplot(alpha_estimates, aes(x = Sample_Disease, y = Chao1)) + # consider rare species also
  #geom_boxplot() +
  #facet_wrap(~ Sample_type, labeller = labeller(Sample_type = facet_labels)) + # Use labeller to rename facets
  geom_boxplot(aes(color = Sample_Disease), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_Disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  stat_summary(fun = mean, geom = "point", shape = 20, size = 6, color = "black", fill = "blue") + # to add mean
  scale_color_manual(values = sample.colors, labels = c(
                                                  "Gill_Diseased" = "GD", "Gill_Non-diseased" = "GND",
                                                  "Skin_Diseased" = "SD", "Skin_Non-diseased" = "SND",
                                                  "Water_Non-diseased" = "WND", "Water_Diseased" = "WD")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", # Legend position
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
    axis.line = element_line(color = "black")
    #strip.text = element_text(size = 18, face = "bold")  # Increase facet title size
  ) +
  scale_x_discrete(labels = c(
                            "Gill_Diseased" = "GD",
                            "Gill_Non-diseased" = "GND",
                            "Skin_Diseased" = "SD",
                            "Skin_Non-diseased" = "SND",
                            "Water_Non-diseased" = "WND",
                            "Water_Diseased" = "WD")) +
  labs(color = "Disease state",
       x = "Reported disease",
       y = "Chao1") +
  guides(color = guide_legend(nrow = 1))
sd.c

#### Perform Kw test----
alpha_estimates %>%
  group_by(Sample_type) %>%
  kruskal_test(Chao1 ~ Reported_disease)

## Shannon diversity----
sd.s <- ggplot(alpha_estimates, aes(x = Sample_Disease, y = Shannon)) + # consider rare species also
  #geom_boxplot() +
  #facet_wrap(~ Sample_type, labeller = labeller(Sample_type = facet_labels)) + # Use labeller to rename facets
  geom_boxplot(aes(color = Sample_Disease), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_Disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  stat_summary(fun = mean, geom = "point", shape = 20, size = 6, color = "black", fill = "blue") + # to add mean
  scale_color_manual(values = sample.colors, labels = c(
                                                      "Gill_Diseased" = "GD", "Gill_Non-diseased" = "GND",
                                                      "Skin_Diseased" = "SD", "Skin_Non-diseased" = "SND",
                                                      "Water_Non-diseased" = "WND", "Water_Diseased" = "WD")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", # Legend position
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
    axis.line = element_line(color = "black")
    #strip.text = element_text(size = 18, face = "bold")  # Increase facet title size
  ) +
  scale_x_discrete(labels = c(
                              "Gill_Diseased" = "GD",
                              "Gill_Non-diseased" = "GND",
                              "Skin_Diseased" = "SD",
                              "Skin_Non-diseased" = "SND",
                              "Water_Non-diseased" = "WND",
                              "Water_Diseased" = "WD")) +
  labs(color = "Disease state",
       x = "Reported disease",
       y = "Shannon") +
  guides(color = guide_legend(nrow = 1))
sd.s

#### Perform Kw test----
alpha_estimates %>%
  group_by(Sample_type) %>%
  kruskal_test(Shannon ~ Reported_disease)

## Phylogenetic diversity----
# we need to separately calculte the phylogenetic diversity
# follow: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html#diversities
library(microbiome) # data analysis and visualisation
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(picante); packageVersion("picante") # for PD

### Calculate phylogenetic diversity----
rarefy.asvtab <- as.data.frame(ps_rarefy@otu_table)
rarefy.tree <- ps_rarefy@phy_tree

# We first need to check if the tree is rooted or not 
ps_rarefy@phy_tree
# Rooted; includes branch lengths

df.pd <- pd(t(rarefy.asvtab), rarefy.tree, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
datatable(df.pd)

# to plot the PD, we need to add the PD result with the previous alpha divesity dataframe
# get the metadata out 
alpha_estimate_withPD <- meta(alpha_estimates)

# Add the results of PD to this dataframe and then plot.
alpha_estimate_withPD$Phylogenetic_Diversity <- df.pd$PD #ndis.meta got all the metadata and alpha divesity results

# Reorder the factor levels
alpha_estimate_withPD$Sample_type <- factor(alpha_estimate_withPD$Sample_type, 
                                            levels = c("Gill_swab", "Skin_swab", "Pond_water"))


# Faith's PD (after calculating PD diversity)
sd.pd <- ggplot(data = alpha_estimate_withPD, aes(y = Phylogenetic_Diversity, x = Sample_Disease)) +
  #geom_boxplot() +
  #facet_wrap(~ Sample_type, labeller = labeller(Sample_type = facet_labels)) + # Use labeller to rename facets
  geom_boxplot(aes(color = Sample_Disease), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Sample_Disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  stat_summary(fun = mean, geom = "point", shape = 20, size = 6, color = "black", fill = "blue") + # to add mean
  scale_color_manual(values = sample.colors, labels = c(
                                                      "Gill_Diseased" = "GD", "Gill_Non-diseased" = "GND",
                                                      "Skin_Diseased" = "SD", "Skin_Non-diseased" = "SND",
                                                      "Water_Non-diseased" = "WND", "Water_Diseased" = "WD")) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(), # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", # Legend position
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
    axis.line = element_line(color = "black")
    #strip.text = element_text(size = 18, face = "bold")  # Increase facet title size
  ) +
  scale_x_discrete(labels = c(
                            "Gill_Diseased" = "GD",
                            "Gill_Non-diseased" = "GND",
                            "Skin_Diseased" = "SD",
                            "Skin_Non-diseased" = "SND",
                            "Water_Non-diseased" = "WND",
                            "Water_Diseased" = "WD")) +
  labs(#title = "Faith's PD",
    color = "Disease state",  # Customize legend title
    y = "PD",  # Add y-axis title to the plot
    x = "Reported diseae") +
  guides(color = guide_legend(nrow = 1)) + # Set number of columns in the legend
  ylim(0, 120)
sd.pd

#### Perform Kw test----
alpha_estimate_withPD %>%
  group_by(Sample_type) %>%
  kruskal_test(Phylogenetic_Diversity ~ Reported_disease)

## B. Beta diversity----
### Gill ----
pseq_gill <- pseq %>% ps_filter(Sample_type == "Gill_swab")

gill.beta <- pseq_gill %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Sample_Disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Sample_Disease)
  ) +
  scale_color_manual(values = sample.colors, 
                     labels = c("Gill_Diseased" = "GD", 
                                "Gill_Non-diseased" = "GND")) +
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
    axis.line = element_line(color = "black")) +
  labs(color = "Disease state",  # Customize legend title
       title = "Gill") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
gill.beta

#### PERMANOVA Test----
metadata_gill <- sample_data(pseq_gill) %>%
  data.frame() %>%
  tibble()
bray_dist.gill = phyloseq::distance(pseq_gill, method="bray")
PERM_bray.gill <- adonis2(bray_dist.gill ~ Reported_disease, data = metadata_gill)

# Add test result on the plot
gill.beta2 <- gill.beta + annotate(geom = "label",
                                   label = paste("PERMANOVA: R² = ", round(PERM_bray.gill["Reported_disease","R2"], 3), 
                                                 ", p = ", PERM_bray.gill["Reported_disease", "Pr(>F)"], sep = ""),
                                   x=Inf, y=Inf, hjust = 1, vjust = 1)
gill.beta2

### Skin----
pseq_skin <- pseq %>% ps_filter(Sample_type == "Skin_swab")

skin.beta <- pseq_skin %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Sample_Disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Sample_Disease)
  ) +
  scale_color_manual(values = sample.colors, 
                     labels = c("Skin_Diseased" = "SD", 
                                "Skin_Non-diseased" = "SND")) +  # Ensure points and box outlines use the same colors
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
    axis.line = element_line(color = "black")) +
  labs(color = "Disease state",  # Customize legend title
       title = "Skin") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
skin.beta

#### PERMANOVA Test----
metadata_skin <- sample_data(pseq_skin) %>%
  data.frame() %>%
  tibble()
bray_dist.skin = phyloseq::distance(pseq_skin, method="bray")
PERM_bray.skin <- adonis2(bray_dist.skin ~ Reported_disease, data = metadata_skin)

# Add test result on the plot
skin.beta2 <- skin.beta + annotate(geom = "label",
                                   label = paste("PERMANOVA: R² = ", round(PERM_bray.skin["Reported_disease","R2"], 3), 
                                                 ", p = ", PERM_bray.skin["Reported_disease", "Pr(>F)"], sep = ""),
                                   x=Inf, y=Inf, hjust = 1, vjust = 1)
skin.beta2

### Water---- 
pseq_water <- pseq %>% ps_filter(Sample_type == "Pond_water")

water.beta <- pseq_water %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Sample_Disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Sample_Disease)
  ) +
  scale_color_manual(values = sample.colors, 
                     labels = c("Water_Non-diseased" = "WND", "Water_Diseased" = "WD")) +  # Ensure points and box outlines use the same colors
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
    axis.line = element_line(color = "black")) +
  labs(color = "Disease state",  # Customize legend title
       title = "Water") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
water.beta

#### PERMANOVA Test----
metadata_water <- sample_data(pseq_water) %>%
  data.frame() %>%
  tibble()
bray_dist.water = phyloseq::distance(pseq_water, method="bray")
PERM_bray.water <- adonis2(bray_dist.water ~ Reported_disease, data = metadata_water)

# Add test result on the plot
water.beta2 <- water.beta + annotate(geom = "label",
                                     label = paste("PERMANOVA: R² = ", round(PERM_bray.water["Reported_disease","R2"], 3), 
                                                   ", p = ", PERM_bray.water["Reported_disease", "Pr(>F)"], sep = ""),
                                     x=Inf, y=Inf, hjust = 1, vjust = 1)
water.beta2

## C.Shared ASVs----
# Load packages
library(UpSetR); packageVersion("UpSetR")
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess")
library(dplyr); packageVersion("dplyr")

### Diseased and non-diseased across sample types----
# use the short name of the sample types
name_mapping <- c(
  "Gill_Diseased" = "GD",
  "Gill_Non-diseased" = "GND",
  "Skin_Diseased" = "SD",
  "Skin_Non-diseased" = "SND",
  "Water_Non-diseased" = "WND",
  "Water_Diseased" = "WD") # This will change the name in the metadata, be careful if the above codes need to rund again,

# Change the Levels of the Factor
sample_data(pseq)$Sample_Disease <- recode_factor(
  sample_data(pseq)$Sample_Disease,
  !!!name_mapping
)

# Set the color maping
sets.bar.color <- c(
  "GD" = "#F35E5A",
  "GND" = "#9DCC00",
  "SD" = "#B68A06",
  "SND" = "#7cb9cb",
  "WND" = "#619CFF",
  "WD" = "#FDBF6F")

# Prepare the data for upset plot
st3 <- get_upset(pseq, factorNames = "Sample_Disease")
# Check the sum of ASVs present in each group
colSums(st3 == 1)
#GD   GND    SD   SND   WND    WD 
#26885 25911 24129 22446 13963 13894 

### Plot shared ASVs
upset.st3 <- upset(st3, 
                   mainbar.y.label = "Shared ASVs between diseased and\nnon-diseased ponds across sample types", 
                   #mainbar.y.max = 8000,
                   nsets = 6, # for more than 5 set, need to specify this
                   nintersects = NA, # for more than 5 set, need to specify this
                   matrix.color = "#F740CE", 
                   main.bar.color = "#9370FF",
                   sets.bar.color = sets.bar.color,
                   sets.x.label = "ASVs per group",
                   att.pos = "bottom", 
                   att.color = "#9370FF",
                   order.by = c("freq", "degree"), 
                   decreasing = c(TRUE,FALSE),
                   show.numbers = "yes", number.angles = 45,
                   point.size = 5, 
                   line.size = 1, 
                   mb.ratio = c(0.7, 0.3),
                   shade.alpha = 0.5, 
                   matrix.dot.alpha = 1,
                   scale.sets = "identity",
                   text.scale = c(2.5,2,2,2.5,2,2.3),
                   set_size.angles = 90, # angle to rotate the set size plot x-axis text
                   set_size.show = TRUE, # display the set sizes on the set size bar chart
                   set_size.numbers_size = 3.5, #adjust the size of the numbers
                   set_size.scale_max = NULL)
upset.st3
# Save as 2000*1114 and use this in the paper
# Note: This Upset plot can't be directly combined with ggplots. So save the 
## combined alpha and beta diversity plots, save the upset plot with same width,
### Then in Inkscape, combine both plots

# Combine plots----
fig.4 <- cowplot::plot_grid(sd.c,
                   sd.s,
                   sd.pd,
                   gill.beta2 + theme(legend.position = "none"),
                   skin.beta2 + theme(legend.position = "none"),
                   water.beta2 + theme(legend.position = "none"),
                   labels = "AUTO",
                   align = "v",
                   ncol = 3)
fig.4
# Save as 2000*850 and then in a canvas 210*212 mm, load the upset plot and this one, 
## Transform 40% and then adjust accordingly. Don't forget to add ** on gill for PD


