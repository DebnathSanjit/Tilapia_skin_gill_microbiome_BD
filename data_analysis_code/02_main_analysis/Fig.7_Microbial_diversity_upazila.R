# date: 20241101 (updated from 20240903 (updated from 20240703))
# Author: Sanjit Debnath

# with these code, I plotted alpha and beta diversity for prokaryotic and 
## microeukaryotic for different upazila 
# These all works, just need to run these codes, save the final plot,
# just a little modification in inkscape. 
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

upazila.colors <- c("Daudkandi" = "#75C333", "Chandina" = "#EC823C", "Barura" = "#100AFF", 
                    "Lalmai" = "#1B9E77", "Laksam" = "#5c47b8", "Nangalkot" = "#02B7EB", "Faridganj" = "#FC81DF")

## Setup working dictionary first
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)
#theme_set(theme(panel.background = element_blank(), axis.line = element_line(color = "black")))
# this gives a clear background, but when I want to combine several plots, it shows problem
# but if i use patchwork, then it works. so i will use pathwork

# Effect of geographical locations----
# Prokaryotes----
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

## Alpha diversity----
ps_rarefy_pseq <- rarefy_even_depth(pseq, rngseed = 1234) #subsample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#8187oTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_pseq <- estimate_richness(ps_rarefy_pseq, measures = c("Chao1", "Shannon"))
alpha_estimates_pseq <- cbind(alpha_estimates_pseq, sample_data(pseq))

# Reorder the factor levels
alpha_estimates_pseq$Sample_type <- factor(alpha_estimates_pseq$Sample_type, 
                                           levels = c("Gill_swab", "Skin_swab", "Pond_water"))

# Define the correct order of Pond_name
pond_order <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20")

# Modify Pond_name factor levels in the sample data
alpha_estimates_pseq$Pond_name <- factor(alpha_estimates_pseq$Pond_name, levels = pond_order)

# Plot alpha diversity
# Reorder the factor levels
alpha_estimates_pseq$Upazila <- factor(alpha_estimates_pseq$Upazila, 
                                       levels = c("Daudkandi", "Chandina", "Barura", "Lalmai", "Laksam", "Nangalkot", "Faridganj"))
### Chao1----
up.c <- ggplot(alpha_estimates_pseq, aes(x = Upazila, y = Chao1)) +
  geom_boxplot() +
  #facet_wrap(~ Sample_type) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 1, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_y_continuous(breaks = seq(0, 8000, by = 500)) +  # Set y-axis breaks with a desired interval
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0.5, size = 18,  face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
                               ) +
  labs(color = "Upazila",
       title = "Chao1", 
       x = "Upazila",
       y = "Chao1") 
up.c

#### Kruskal-Wallis test----
kruskal.test(Chao1 ~ Upazila, data = alpha_estimates_pseq)


##### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
stat.test.up.c <- dunn_test(Chao1 ~ Upazila, data = alpha_estimates_pseq, p.adjust.method = "BH")

##### Add the stat----
stat.test.up.c <- stat.test.up.c %>% add_xy_position(x = "Upazila")
# Adjust y.position
stat.test.up.c2 <- stat.test.up.c %>%
  mutate(y.position = y.position * 0.7)  # Adjust this value as needed
# plot
up.c2 <- up.c + stat_pvalue_manual(stat.test.up.c2, label = "p.adj.signif", 
                                   step.increase = 0.05, tip.length = 0.005, 
                                   size = 8,  # Adjust the size of the asterisk here
                                   hide.ns = TRUE)
up.c2

### Shannon----
up.s <- ggplot(alpha_estimates_pseq, aes(x = Upazila, y = Shannon)) +
  geom_boxplot() +
  #facet_wrap(~ Sample_type) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  scale_y_continuous(breaks = seq(1, 12, by = 1)) +  # Set y-axis breaks from 1 to 12 at intervals of 1
  theme(
    plot.title = element_text(hjust = 0.5, size = 18,  face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
                ) +
  labs(color = "Upazila",
       title = "Shannon", 
       x = "Upazila",
       y = "Shannon") +
  guides(color = guide_legend(ncol = 1)) 
up.s

#### Kruskal-Wallis test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_pseq)

##### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
stat.test.up.s <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_pseq, p.adjust.method = "BH")

##### Add stat value----
stat.test.up.s <- stat.test.up.s %>% add_xy_position(x = "Upazila")
up.s2 <- up.s + stat_pvalue_manual(stat.test.up.s, label = "p.adj.signif", 
                                   step.increase = 0.05, tip.length = 0.005, 
                                   size = 8,  # Adjust the size of the asterisk here
                                   hide.ns = TRUE)
up.s2

### Phylogenetic diversity----
# we need to separately calculte the phylogenetic diversity
# follow: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html#diversities
library(microbiome); packageVersion("microbiome") # data analysis and visualisation
library(DT); packageVersion("DT") # interactive tables in html and markdown
library(data.table); packageVersion("data.table") # alternative to data.frame
library(picante); packageVersion("picante") # for PD

#### Calculate phylogenetic diversity----
rarefy.asvtab <- as.data.frame(ps_rarefy_pseq@otu_table)
rarefy.tree <- ps_rarefy_pseq@phy_tree

# We first need to check if the tree is rooted or not 
ps_rarefy_pseq@phy_tree
# Rooted; includes branch lengths

# t(ou_table) transposes the table for use in picante 
df.pd <- pd(t(rarefy.asvtab), rarefy.tree, include.root=T) # takes areound 2 minutes
#datatable(df.pd) # view pd table

## to plot the PD, we need to add the PD result with the previous alpha divesity dataframe
## get the metadata out 
alpha_estimate_withPD <- meta(alpha_estimates_pseq)

# Add the results of PD to this dataframe and then plot.
alpha_estimate_withPD$Phylogenetic_Diversity <- df.pd$PD #alpha_estimate_withPD got all the metadata and alpha divesity results

# Reorder the factor levels
alpha_estimate_withPD$Upazila <- factor(alpha_estimate_withPD$Upazila, 
                                        levels = c("Daudkandi", "Chandina", "Barura", "Lalmai", "Laksam", "Nangalkot", "Faridganj"))
# Modify Pond_name factor levels in the sample data
#alpha_estimate_withPD$Pond_name <- factor(alpha_estimate_withPD$Pond_name, levels = pond_order)

#### Plot Phylogenetic diversity----
up.pd <- ggplot(data = alpha_estimate_withPD, aes(y = Phylogenetic_Diversity, x = Upazila)) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  scale_y_continuous(breaks = seq(0, 250, by = 50)) +  # Set y-axis breaks with a desired interval
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
  ) +
  #scale_x_discrete(labels = c("Gill_swab" = "Gill", 
  #                            "Skin_swab" = "Skin", 
  #                            "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Upazila",  # Customize legend title
       y = "PD",  # Add y-axis title to the plot
       x = "Upazila",
       title = "Faith's PD") +
  guides(color = guide_legend(ncol = 1)) # Set number of columns in the legend
up.pd

##### Perform Kw test----
kruskal.test(Phylogenetic_Diversity ~ Upazila, data = alpha_estimate_withPD)

# perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.up.pd <- dunn_test(Phylogenetic_Diversity ~ Upazila, data = alpha_estimate_withPD, p.adjust.method = "BH")

##### Add stat value----
stat.test.up.pd <- stat.test.up.pd %>% add_xy_position(x = "Upazila")
up.pd2 <- up.pd + stat_pvalue_manual(stat.test.up.pd, label = "p.adj.signif", 
                                     step.increase = 0.05, tip.length = 0.005, 
                                     size = 8,  # Adjust the size of the asterisk here
                                     hide.ns = TRUE)
up.pd2

## B.Beta diversity (Bray-Curtis)----
# beta diversity using bray-curtis
up.beta.cal <- pseq %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") 

up.beta <- up.beta.cal %>%
  ord_calc(
    method = "PCoA") %>% 
  ord_plot(
    axes = c(1, 2),
    color = "Upazila", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(color = Upazila)
  ) +
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",#c(vjust = 0.9, hjust = 0.2),# need to be outside to get the legend
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
    ) +
  #scale_x_discrete(labels = c("Gill_swab" = "Gill", 
  #                            "Skin_swab" = "Skin", 
  #                            "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Upazila",  # Customize legend title
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(ncol = 1))  # Set number of columns in the legend
up.beta

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_pseq <- sample_data(pseq) %>%
  data.frame() %>%
  tibble()
bray_dist = phyloseq::distance(pseq, method="bray")
PERM_bray.up <- adonis2(bray_dist ~ Upazila, data = metadata_pseq)

# to add test result on the plot
up.beta2 <- up.beta + annotate(geom = "label",
                               label = paste("PERMANOVA: R² = ", round(PERM_bray.up["Upazila","R2"], 3), 
                                             ", p = ", PERM_bray.up["Upazila", "Pr(>F)"], sep = ""),
                               x=Inf, y=Inf, hjust = 1, vjust = 1)
up.beta2

#### Pairwise PERMANOVA (pairwiseAdonis):----
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

# Perform pairwise comparisons
pairwise_result.up <- pairwise.adonis(bray_dist, metadata_pseq$Upazila, p.adjust.m = "BH")
pairwise_result.up
# Save as csv file
#write.csv(pairwise_result.up, file = "pairwise_permanova_results_for_upazila_16s_20241101.csv", row.names = FALSE)


#load phyloseq object
ps.18s <- readRDS("phyloseq_FSS_BD_metadata_v7_18S_20240703.rds")
ps.18s # 2961 taxa and 57 samples

## Estimate alpha diversity----
ps_rarefy_ps.18s <- rarefy_even_depth(ps.18s, rngseed = 1234) #subsample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#295OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates_ps.18s <- estimate_richness(ps_rarefy_ps.18s, measures = c("Chao1", "Shannon"))
alpha_estimates_ps.18s <- cbind(alpha_estimates_ps.18s, sample_data(ps.18s))

# Microeukaryotes----
# Define the correct order of Pond_name
pond_order <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20")

# Modify Pond_name factor levels in the sample data
alpha_estimates_ps.18s$Pond_name <- factor(alpha_estimates_ps.18s$Pond_name, levels = pond_order)

# Plot alpha diversity
# Reorder the factor levels
alpha_estimates_ps.18s$Upazila <- factor(alpha_estimates_ps.18s$Upazila, 
                                         levels = c("Daudkandi", "Chandina", "Barura", "Lalmai", "Laksam", "Nangalkot", "Faridganj"))
### Chao1----
up.c.18s <- ggplot(alpha_estimates_ps.18s, aes(x = Upazila, y = Chao1)) +
  geom_boxplot() +
  #facet_wrap(~ Sample_type) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  scale_y_continuous(breaks = seq(0, 2000, by = 250)) +  # Set y-axis breaks with a desired interval
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
  ) +
  labs(color = "Upazila",
       title = "Chao1", 
       x = "Upazila",
       y = "Chao1") +
  ylim(0,800) # after setting ylim, p value doesn't add automatically
up.c.18s

## Kruskal-Wallis test----
kruskal.test(Chao1 ~ Upazila, data = alpha_estimates_ps.18s)

##### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
stat.test.up.c.18s <- dunn_test(Chao1 ~ Upazila, data = alpha_estimates_ps.18s, p.adjust.method = "BH")

##### Add stat value----
stat.test.up.c.18s <- stat.test.up.c.18s %>% add_xy_position(x = "Upazila")
# Adjust y.position
stat.test.up.c2.18s <- stat.test.up.c.18s %>%
  mutate(y.position = y.position * 0.7)  # Adjust this value as needed

up.c2.18s <- up.c.18s + stat_pvalue_manual(stat.test.up.c2.18s, label = "p.adj.signif", 
                                           step.increase = 0.05, tip.length = 0.005, 
                                           size = 8,  # Adjust the size of the asterisk here
                                           hide.ns = TRUE)
up.c2.18s # it doesn't add the p value. add it manually

### Shannon----
up.s.18s <- ggplot(alpha_estimates_ps.18s, aes(x = Upazila, y = Shannon)) +
  geom_boxplot() +
  #facet_wrap(~ Sample_type) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  #scale_y_continuous(breaks = seq(0, , by = 500)) +  # Set y-axis breaks with a desired interval
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
  ) +
  labs(color = "Upazila",
       title = "Shannon", 
       x = "Upazila",
       y = "Shannon")
up.s.18s

## Kruskal-Wallis test----
kruskal.test(Shannon ~ Upazila, data = alpha_estimates_ps.18s)

###### Perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)----
#stat.test.up.s.18s <- dunn_test(Shannon ~ Upazila, data = alpha_estimates_ps.18s, p.adjust.method = "BH")
#
#stat.test.up.s.18s <- stat.test.up.s.18s %>% add_xy_position(x = "Upazila")
#up.s2.18s <- up.s.18s + stat_pvalue_manual(stat.test.up.s.18s, label = "p.adj.signif", 
#                                           step.increase = 0.05, tip.length = 0.005, 
#                                           size = 8,  # Adjust the size of the asterisk here
#                                           hide.ns = TRUE)
#up.s2.18s

### Phylogenetic diversity----
rarefy.asvtab.18s <- as.data.frame(ps_rarefy_ps.18s@otu_table)
rarefy.tree.18s <- ps_rarefy_ps.18s@phy_tree

# We first need to check if the tree is rooted or not 
ps_rarefy_ps.18s@phy_tree
# Rooted; includes branch lengths

df.pd.18s <- pd(t(rarefy.asvtab.18s), rarefy.tree.18s, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
datatable(df.pd.18s)
## to plot the PD, we need to add the PD result with the previous alpha divesity dataframe
## get the metadata out 
alpha_estimate_withPD.18s <- meta(alpha_estimates_ps.18s)

# Add the results of PD to this dataframe and then plot.
alpha_estimate_withPD.18s$Phylogenetic_Diversity <- df.pd.18s$PD #alpha_estimate_withPD.18s got all the metadata and alpha divesity results

# Reorder the factor levels
alpha_estimate_withPD.18s$Upazila <- factor(alpha_estimate_withPD.18s$Upazila, 
                                        levels = c("Daudkandi", "Chandina", "Barura", "Lalmai", "Laksam", "Nangalkot", "Faridganj"))

#### Plot Phylogenetic diversity----
up.pd.18s <- ggplot(data = alpha_estimate_withPD.18s, aes(y = Phylogenetic_Diversity, x = Upazila)) +
  geom_boxplot(aes(color = Upazila), fill = NA) +  # Outline color, no fill for box
  geom_jitter(aes(color = Upazila), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  scale_y_continuous(breaks = seq(0, 400, by = 50)) +  # Set y-axis breaks with a desired interval
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
  ) +
  #scale_x_discrete(labels = c("Gill_swab" = "Gill", 
  #                            "Skin_swab" = "Skin", 
  #                            "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Upazila",  # Customize legend title
       y = "PD",  # Add y-axis title to the plot
       x = "Upazila",
       title = "Faith's PD") +
  guides(color = guide_legend(ncol = 1)) + # Set number of columns in the legend
  ylim(0, 250)
up.pd.18s

##### Perform Kw test----
kruskal.test(Phylogenetic_Diversity ~ Upazila, data = alpha_estimate_withPD.18s)

##### perform Dunn test (https://www.youtube.com/watch?v=pyLQmUfrel8)
stat.test.up.pd.18s <- dunn_test(Phylogenetic_Diversity ~ Upazila, 
                                 data = alpha_estimate_withPD.18s, 
                                 p.adjust.method = "BH")

###### Add the value on the plotBox plot----
stat.test.up.pd.18s2 <- stat.test.up.pd.18s %>% add_xy_position(x = "Upazila")
up.pd.18s2 <- up.pd.18s + stat_pvalue_manual(stat.test.up.pd.18s2, 
                                             label = "p.adj.signif", 
                                     step.increase = 0.05, tip.length = 0.005, 
                                     size = 8,  # Adjust the size of the asterisk here
                                     hide.ns = TRUE)
up.pd.18s2 # p value doesn't add, so add manually

## B.Beta diversity (Bray-Curtis)----
# beta diversity using bray-curtis
up.beta.18s <- ps.18s %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    color = "Upazila", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(color = Upazila)
  ) +
  scale_color_manual(values = upazila.colors) +  # Ensure points and box outlines use the same colors
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "right",#c(vjust = 0.9, hjust = 0.2),# need to be outside to get the legend
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove border
    panel.background = element_blank(),#rect(fill = "white"),  # Set background color to white
    axis.line.x = element_line(color = "black"),  # Add black line to x-axis
    axis.line.y = element_line(color = "black")   # Optional: Add border to panel
  ) +
  #scale_x_discrete(labels = c("Gill_swab" = "Gill", 
  #                            "Skin_swab" = "Skin", 
  #                            "Pond_water" = "Water"))  + # Insert line breaks in category names
  labs(color = "Upazila",  # Customize legend title
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(ncol = 1))  # Set number of columns in the legend
up.beta.18s

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_ps.18s <- sample_data(ps.18s) %>%
  data.frame() %>%
  tibble()
bray_dist.18s = phyloseq::distance(ps.18s, method="bray")
PERM_bray.18s <- adonis2(bray_dist.18s ~ Upazila, data = metadata_ps.18s)

# to add test result on the plot
up.beta2.18s <- up.beta.18s + annotate(geom = "label",
                                       label = paste("PERMANOVA: R² = ", round(PERM_bray.18s["Upazila","R2"], 3), 
                                                     ", p = ", PERM_bray.18s["Upazila", "Pr(>F)"], sep = ""),
                                       x=Inf, y=Inf, hjust = 1, vjust = 1)
up.beta2.18s

#### Pairwise PERMANOVA (pairwiseAdonis):----
pairwise_result.18s <- pairwise.adonis(bray_dist.18s, metadata_ps.18s$Upazila, p.adjust.m = "BH")
pairwise_result.18s

#write.csv(pairwise_result.18s, file = "pairwise_permanova_results_18s_20240805.csv", row.names = FALSE)

# Combine plots----
comb.up <- cowplot::plot_grid(up.c2 + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.s2 + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.pd2 + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.beta2 + theme(legend.position = "none"),
                         up.c2.18s  + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.s.18s  + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.pd.18s2  + theme(axis.title.x = element_blank(), legend.position = "none"),
                         up.beta2.18s +theme(legend.position = "none"),
                         align = "v",# this was was keeping lots of gaps between each plot
                         labels = c("A Prokaryotes","", "", "", "B Microeukaryotes"),
                         ncol = 4)
comb.up
# get legend
legend.up <- get_legend(up.s)
# final plot
comb.up.ab <- cowplot::plot_grid(comb.up,
                                 legend.up,
                                 rel_widths = c(2, 0.15),
                                 ncol = 2)
comb.up.ab
# saved as 2000*1200

# Significance heatmap instead of *----
## Make Significance heatmap for Chao1----
# Extract relevant columns
heatmap_data.up.c <- stat.test.up.c[, c("group1", "group2", "p.adj.signif")]

# Pivot data to wide format: rows as group1, columns as group2
heatmap_wide.up.c <- pivot_wider(
  heatmap_data.up.c,
  names_from = group2,
  values_from = p.adj.signif,
  values_fill = list(p.adj.signif = NA))  # Fill missing values with NA

# Set row names from 'group1' and remove the 'group1' column
rownames(heatmap_wide.up.c) <- heatmap_wide.up.c$group1

# Convert to a matrix
heatmap_matrix.up.c <- as.matrix(heatmap_wide.up.c)
heatmap_matrix.up.c <- heatmap_matrix.up.c[, -1]  # Remove 'group1' column as it's now the row names

library(reshape2): packageVersion("reshape2")

# Melt the matrix into long format for ggplot2
heatmap_melted.up.c <- melt(heatmap_matrix.up.c, na.rm = FALSE)
colnames(heatmap_melted.up.c) <- c("group1", "group2", "Significance")

# Load the dplyr package
library(dplyr); packageVersion("dplyr")

# Remove rows with NA in the Significance column
heatmap_melted.up.c <- heatmap_melted.up.c %>%
  filter(!is.na(Significance))

# Define a color scale based on significance levels
color_scale <- scale_fill_manual(
  values = c("ns" = "grey", "*" = "#ccf9ff", "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"),
  na.value = "white"
)

# plot
up.c.sig <- ggplot(data = heatmap_melted.up.c, aes(x = group2, y = group1, fill = Significance)) +
  geom_tile(color = "white") +  # Creates the tiles of the heatmap
  scale_fill_manual(values = c("ns" = "grey", "*" = "#ccf9ff", 
                               "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"))+
  labs(x = "Group 2", y = "Group 1", fill = "Significance") +  # Axis labels
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 25),
        legend.position = c(vjust = 0.3, hjust = 0.9),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove border
        axis.line = element_line(color = "black")) +
  labs(fill = "FDR")+
  scale_x_discrete(labels = c("Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan",
                              "Faridganj" = "Far")) +
  scale_y_discrete(labels = c("Daudkandi" = "Dau",
                              "Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan")) +
  guides(fill = guide_legend(nrow = 1))
up.c.sig
# save 350*200

## Make heatmap for shannon----
# Extract relevant columns
heatmap_data.up.s <- stat.test.up.s[, c("group1", "group2", "p.adj.signif")]

# Pivot data to wide format: rows as group1, columns as group2
heatmap_wide.up.s <- pivot_wider(
  heatmap_data.up.s,
  names_from = group2,
  values_from = p.adj.signif,
  values_fill = list(p.adj.signif = NA))  # Fill missing values with NA

# Set row names from 'group1' and remove the 'group1' column
rownames(heatmap_wide.up.s) <- heatmap_wide.up.s$group1

# Convert to a matrix
heatmap_matrix.up.s <- as.matrix(heatmap_wide.up.s)
heatmap_matrix.up.s <- heatmap_matrix.up.s[, -1]  # Remove 'group1' column as it's now the row names

# Melt the matrix into long format for ggplot2
heatmap_melted.up.s <- melt(heatmap_matrix.up.s, na.rm = FALSE)
colnames(heatmap_melted.up.s) <- c("group1", "group2", "Significance")

# Remove rows with NA in the Significance column
heatmap_melted.up.s <- heatmap_melted.up.s %>%
  filter(!is.na(Significance))

# Define a color scale based on significance levels
color_scale <- scale_fill_manual(
  values = c("ns" = "grey", "*" = "#ccf9ff", "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"),
  na.value = "white"
)

# plot
up.s.sig <- ggplot(data = heatmap_melted.up.s, aes(x = group2, y = group1, fill = Significance)) +
  geom_tile(color = "white") +  # Creates the tiles of the heatmap
  scale_fill_manual(values = c("ns" = "grey", "*" = "#ccf9ff", 
                               "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"))+
  labs(x = "Group 2", y = "Group 1", fill = "Significance") +  # Axis labels
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 25),
        legend.position = c(vjust = 0.3, hjust = 0.9),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove border
        axis.line = element_line(color = "black")) +
  labs(fill = "FDR")+
  scale_x_discrete(labels = c("Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan",
                              "Faridganj" = "Far")) +
  scale_y_discrete(labels = c("Daudkandi" = "Dau",
                              "Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan")) +
  guides(fill = guide_legend(nrow = 1))
up.s.sig
# save 350*200

## Make heatmap for PD----
# Extract relevant columns
heatmap_data.up.pd <- stat.test.up.pd[, c("group1", "group2", "p.adj.signif")]

# Pivot data to wide format: rows as group1, columns as group2
heatmap_wide.up.pd <- pivot_wider(
  heatmap_data.up.pd,
  names_from = group2,
  values_from = p.adj.signif,
  values_fill = list(p.adj.signif = NA))  # Fill missing values with NA

# Set row names from 'group1' and remove the 'group1' column
rownames(heatmap_wide.up.pd) <- heatmap_wide.up.pd$group1

# Convert to a matrix
heatmap_matrix.up.pd <- as.matrix(heatmap_wide.up.pd)
heatmap_matrix.up.pd <- heatmap_matrix.up.pd[, -1]  # Remove 'group1' column as it's now the row names

# Melt the matrix into long format for ggplot2
heatmap_melted.up.pd <- melt(heatmap_matrix.up.pd, na.rm = FALSE)
colnames(heatmap_melted.up.pd) <- c("group1", "group2", "Significance")

# Remove rows with NA in the Significance column
heatmap_melted.up.pd <- heatmap_melted.up.pd %>%
  filter(!is.na(Significance))

# Define a color scale based on significance levels
color_scale <- scale_fill_manual(
  values = c("ns" = "grey", "*" = "#ccf9ff", "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"),
  na.value = "white"
)

# plot
up.pd.sig <- ggplot(data = heatmap_melted.up.pd, aes(x = group2, y = group1, fill = Significance)) +
  geom_tile(color = "white") +  # Creates the tiles of the heatmap
  scale_fill_manual(values = c("ns" = "grey", "*" = "#ccf9ff", 
                               "**" = "#7ce8ff", "***" = "#55d0ff", "****" = "#00acdf"))+
  labs(x = "Group 2", y = "Group 1", fill = "Significance") +  # Axis labels
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 25),
        legend.position = c(vjust = 0.3, hjust = 0.9),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove border
        axis.line = element_line(color = "black")) +
  labs(fill = "FDR")+
  scale_x_discrete(labels = c("Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan",
                              "Faridganj" = "Far")) +
  scale_y_discrete(labels = c("Daudkandi" = "Dau",
                              "Chandina" = "Cha", 
                              "Barura" = "Bar", 
                              "Lalmai" = "Lal", 
                              "Laksam" = "Lak", 
                              "Nangalkot" = "Nan")) +
  guides(fill = guide_legend(nrow = 1))
up.pd.sig
# save 350*200 and add this to the comb.up.ab plot in inkscape. I tried using R but doesn't work as i want
## As the legend keys are not square, transform the height by 70%
