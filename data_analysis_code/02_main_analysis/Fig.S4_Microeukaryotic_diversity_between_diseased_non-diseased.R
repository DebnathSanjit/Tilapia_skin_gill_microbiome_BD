# date: 20240825 (fine tuned from 20240610)
# Author: Sanjit Debnath

# with these code, I plotted alpha-beta diversity, shared taxa (asv) 
# using venn diagram between diseased and non-diseased microeukaryotes
# just a little modification in inkscape. not using the chao1 result but plooting in case needed

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(microbiome); packageVersion("microbiome") # data analysis and visualisation
library(DT); packageVersion("DT") # interactive tables in html and markdown
library(data.table); packageVersion("data.table") # alternative to data.frame
library(picante); packageVersion("picante") # for PD
library(gridExtra); packageVersion("gridExtra")
# Libraries for tests
#library(ggpubr); packageVersion("ggpubr")
#library(rstatix); packageVersion("rstatix")
#library(dunn.test); packageVersion("dunn.test")

## Setup working dictionary 
setwd("C:/Users/scd226/OneDrive - University of Exeter/Experiment & Research/Data_analysis/2.Field_study_samples_BD/R_scripts/statistical_analysis")

#load theme
theme_set(theme_bw())
set.seed(1234)

# Microeukaryotic----
#load phyloseq object
ps.18s <- readRDS("phyloseq_FSS_BD_metadata_v7_18S_20240703.rds")
ps.18s # 2961 taxa and 57 samples

## A. Alpha diversity----
set.seed(1234)#The set. seed() function in R is used to create reproducible results when writing code that involves creating variables that take on random values. By using the set. seed() function, you guarantee that the same random values are produced each time you run the code
ps_rarefy.18s <- rarefy_even_depth(ps.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#295OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates.18s <- estimate_richness(ps_rarefy.18s, measures = c("Chao1", "Shannon"))
alpha_estimates.18s <- cbind(alpha_estimates.18s, sample_data(ps.18s))

### Calculate phylogenetic diversity----
# follow: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html#diversities
rarefy.asvtab.18s <- as.data.frame(ps_rarefy.18s@otu_table)
rarefy.tree.18s <- ps_rarefy.18s@phy_tree

# We first need to check if the tree is rooted or not 
ps_rarefy.18s@phy_tree
# Rooted; includes branch lengths

df.pd.18s <- pd(t(rarefy.asvtab.18s), rarefy.tree.18s, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
#datatable(df.pd.18s)

# to plot the PD, we need to add the PD result with the previous alpha divesity dataframe
# get the metadata out 
alpha_estimate_withPD.18s <- meta(alpha_estimates.18s)

# Add the results of PD to this dataframe and then plot.
alpha_estimate_withPD.18s$Phylogenetic_Diversity <- df.pd.18s$PD 

# Reshape data to long format
alpha_long <- alpha_estimate_withPD.18s %>%
  pivot_longer(cols = c(Chao1, Shannon, Phylogenetic_Diversity),
               names_to = "Index",
               values_to = "Value")

# Reorder the factor levels
alpha_long$Index <- factor(alpha_long$Index, 
                           levels = c("Chao1", "Shannon", "Phylogenetic_Diversity"))

# Create a named vector for renaming the facet titles
facet_labels.18s <- c("Chao1" = "A Chao1", "Shannon" = "B Shannon", "Phylogenetic_Diversity" = "C PD")

# use customised colour for water
water.color <- c("Non-diseased" = "#619CFF", "Diseased" = "#FDBF6F")

# Create the plot
# Instead of plotting seperate, then can be plotted using the below funciton
library(patchwork); packageVersion("patchwork")

# Create a list of unique Index levels
indices <- unique(alpha_long$Index)

# Create a list to hold the individual plots
plots <- list()

# Define custom y-axis labels for each Index
y_axis_labels <- c(
  "Chao1"="Chao1", "Shannon"="Shannon", "Phylogenetic_Diversity" = "PD"
)

# Loop through each Index level and create a separate plot
for (i in indices) {
  p <- ggplot(subset(alpha_long, Index == i), aes(x = Reported_disease, y = Value)) +
    geom_boxplot(aes(color = Reported_disease), fill = NA) +
    geom_jitter(aes(color = Reported_disease), width = 0.2, alpha = 0.8, size = 2, shape = 16) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 6, color = "black", fill = "blue") +
    labs(y = y_axis_labels[i], title = i) +  # Set custom y-axis label and title for each plot
    scale_color_manual(values = water.color,
                       labels = c("Diseased" = "WD", "Non-diseased" = "WND")) +  # Ensure points and box outlines use the same colors
    theme(
      plot.title = element_text(hjust = 0, size = 16, margin = margin(b = 15)),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
      axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
      axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_blank(),
      legend.spacing.y = unit(1, "cm"),
      legend.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      strip.background = element_blank(),  # Remove background shade for facet titles
      strip.text = element_text(size = 18, face = "bold"))
  
  # Add the plot to the list
  plots[[i]] <- p
}

# Combine all the plots using patchwork
combined_plot <- wrap_plots(plots, ncol = 3)  # Adjust `ncol` as needed

# Display the combined plot
combined_plot

##### Perform Kw test----
## Chao1
kruskal.test(Chao1 ~ Reported_disease, data = alpha_estimates.18s)
#data:  Chao1 by Sample_type
#Kruskal-Wallis chi-squared = 2.189, df = 1, p-value = 0.139
#
## Shannon
kruskal.test(Shannon ~ Reported_disease, data = alpha_estimates.18s)
##data:  Shannon by Sample_type
##Kruskal-Wallis chi-squared = 0.31681, df = 1, p-value = 0.5735
#
## PD
kruskal.test(Phylogenetic_Diversity ~ Reported_disease, data = alpha_estimate_withPD.18s)
##data:  Phylogenetic_Diversity by Sample_type
##Kruskal-Wallis chi-squared = 0.45621, df = 1, p-value = 0.4994


## B. Beta diversity----
# Water
water.beta.18s <- ps.18s %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "PCoA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Reported_disease", 
    #fill = "Reported_disease",
    shape = 16, 
    alpha = 1,
    size = 3
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse( #comment out to remove too many ellipses
    ggplot2::aes(colour = Reported_disease)
  ) +
  scale_color_manual(values = water.color, labels = c("Non-diseased" = "WND", 
                                                      "Diseased" = "WD")) +  # Ensure points and box outlines use the same colors
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
       title = "PCoA") +
  guides(color = guide_legend(nrow = 1))  # Set number of columns in the legend
water.beta.18s

### PERMANOVA Test----
metadata_water.18s <- sample_data(ps.18s) %>%
  data.frame() %>%
  tibble()
bray_dist.water.18s = phyloseq::distance(ps.18s, method="bray")
PERM_bray.water.18s <- adonis2(bray_dist.water.18s ~ Reported_disease, data = metadata_water.18s)

# to add test result on the plot
water.beta.18s2 <- water.beta.18s + annotate(geom = "label",
                                             label = paste("PERMANOVA: R² = ", round(PERM_bray.water.18s["Model","R2"], 3), 
                                                           ", p = ", PERM_bray.water.18s["Model", "Pr(>F)"], sep = ""),
                                             x=Inf, y=Inf, hjust = 1, vjust = 1)
water.beta.18s2

## Venn diagram----
# load libraries
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess") # needed to make the vennlist
library(VennDiagram); packageVersion("VennDiagram") # for venn diagram
library(gridExtra); packageVersion("gridExtra")# need to visualise venn diagram properly

### Diseased and non-diseased----
vennlist.18s.dis <- get_vennlist(obj = ps.18s, factorNames = "Reported_disease")

# to rename the name of the lists
names(vennlist.18s.dis) <- c("WD", "WND")

# Plot
venn.18s.dis <- venn.diagram(vennlist.18s.dis,
                             filename = NULL, # can't make it work
                             disable.logging = TRUE, # stop log file output in dictionary
                             main = "Shared microeukaryotic ASVs", 
                             main.pos = c(0.5, 1.0), 
                             main.fontface = "bold",
                             main.fontfamily = "serif", main.col = "black",
                             main.cex = 1.5, main.just = c(0.5, 1),
                             fill = c("#FDBF6F","#619CFF"),
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
venn.18s.dis.p <- grid.arrange(venn.18s.dis) # to view the plot properly (it worked before but not it's not working)

# try this
library(grid); packageVersion("grid")
library(gridExtra); packageVersion("gridExtra")

# Combine the grobs into one
venn_grob.18s <- grobTree(venn.18s.dis)

# Now arrange it
venn.18s.dis.p <- grid.arrange(venn_grob.18s)

# Combine 18s----
comb.bv.18s <- cowplot::plot_grid(water.beta.18s2 + theme(legend.position = "none"),
                                  venn.18s.dis.p,
                                  labels = c("D", "E"))

combine.18s <- cowplot::plot_grid(combined_plot + theme(legend.position = "none"),
                                  comb.bv.18s,# + theme(legend.position = "none"),
                                  ncol = 1,
                                  rel_heights = c(4.5,5.5),
                                  align = "V")
combine.18s
# save as 1500*900

