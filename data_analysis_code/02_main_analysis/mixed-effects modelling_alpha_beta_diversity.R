# date: 20250616
# Author: Sanjit Debnath

# This script is for modelling alpha and beta diversity, not for individual factor but for mixed factors
## When I used Kruskal-Wallis, it doesn't cover the repeated measure from each pond, also it does not consider the effect of several factors
## such as disease condition or location, which could be a confounding factor. so it is necessary to run mixed-effct model to cover all these
#### So I made these figure and model, just a little modification in inkscape. 

# Load Libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(microViz); packageVersion("microViz")
library(vegan); packageVersion("vegan") # needed for PERMANOVA test
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(stringr); packageVersion("stringr") # to wrap text
library(gridExtra); packageVersion("gridExtra")
library(microbiome); packageVersion("microbiome") # data analysis and visualisation
library(DT); packageVersion("DT")# interactive tables in html and markdown
library(data.table); packageVersion("data.table") # alternative to data.frame
library(picante); packageVersion("picante") # for PD
# Libraries for tests
library(ggpubr); packageVersion("ggpubr")
library(rstatix); packageVersion("rstatix")
library(dunn.test); packageVersion("dunn.test")
library(writexl); packageVersion("writexl")
library(lmerTest); packageVersion("lmerTest")
library(emmeans); packageVersion("emmeans")
library(tibble); packageVersion("tibble")
library(performance); packageVersion("performance")    # For model diagnostics (optional)
library(AICcmodavg); packageVersion("AICcmodavg")     # Alternative AIC model comparison
library(glmmTMB); packageVersion("glmmTMB")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

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

## To compare gill and skin sample collected from same fish, add a column named Fish_ID where same fish will have same name
# Step 1: Copy metadata
sample_df <- data.frame(sample_data(pseq))

# Step 2: Extract Pond ID (e.g., "P5")
sample_df$Pond_name <- sub("_.*", "", sample_df$Sample)

# Step 3: Extract Fish number (e.g., "2" from SSw2 or GSw2)
sample_df$Fish_num <- ifelse(grepl("Sw[0-9]+", sample_df$Sample),
                             sub(".*Sw([0-9]+)", "\\1", sample_df$Sample),
                             NA)

# Step 4: Combine to make Fish_ID (e.g., "P5_F2")
sample_df$Fish_ID <- ifelse(!is.na(sample_df$Fish_num),
                            paste0(sample_df$Pond_name, "_F", sample_df$Fish_num),
                            NA)

# Step 5: Assign Fish_ID back to phyloseq object
sample_data(pseq)$Fish_ID <- sample_df$Fish_ID

# Prokaryotes----
# Define the correct order of Pond_name
pond_order <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20")

# Modify Pond_name factor levels in the sample data
sample_data(pseq)$Pond_name <- factor(sample_data(pseq)$Pond_name, levels = pond_order)

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

# Reshape data into long format
alpha_long <- alpha_estimate_withPD %>%
  pivot_longer(cols = c(Chao1, Shannon, Phylogenetic_Diversity), 
               names_to = "Diversity_metric", 
               values_to = "Diversity_value")

#### Plot alpha diversity####
plot_alpha <- ggplot(alpha_long, aes(x = Pond_name, y = Diversity_value)) +
  geom_boxplot(aes(color = Sample_type,), fill = NA) +
  #geom_jitter(aes(color = Sample_type), width = 0.2, alpha = 0.8, size = 2, shape = 16) +  # Increase size of jitter points, set shape to round
  facet_wrap(~Diversity_metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = sample.colors, labels = c("Gill", "Skin", "Water")) +  # Ensure points and box outlines use the same colors
  #theme_bw() +
  theme(
    #plot.title = element_blank(),#text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)), # Increase margin for plot title
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, margin = margin(t = 5)),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14, margin = margin(r = 5)),
    axis.title.x = element_blank(),#text(size = 18, face = "bold", margin = margin(t = 10)),  # Increase margin for x-axis title
    axis.title.y = element_text(angle = 90, size = 18, face = "bold", margin = margin(r = 10)),  # Increase margin for y-axis title
    legend.position = "bottom", #c(vjust = 0.9, hjust = 0.75),
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank(),
    legend.spacing.y = unit(1, "cm"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold")  # For facet strip
  ) + 
  labs(x = "Pond", y = "Alpha Diversity", 
       color = "Sample Type", fill = "Sample Type") 
plot_alpha
# This if alternate plot of figure 2, showing variation across sample types for each pond. But I am not using this
## Statistical modelling for this and also diseased vs non-diseased are done below

### Mixed-effect model####
## Does alpha diversity differ across sample types, and within each sample type, does disease affect alpha diversity?

#### Chao1####
model_c_lmm <- lmerTest::lmer(Chao1 ~ Sample_type * Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimate_withPD)
summary(model_c_lmm) # fits better

##### Check residuals#####
plot(model_c_lmm)
qqnorm(resid(model_c_lmm)); qqline(resid(model_c_lmm)) # If points mostly fall along the line, residuals are approximately normal.
hist(resid(model_c_lmm), breaks = 30, main = "Residual Histogram (Chao1)", xlab = "Residuals") #Look bell shaped
plot(density(alpha_estimate_withPD$Chao1)) # Lmm acceptable
shapiro.test(resid(model_c_lmm)) # W = 0.98941, p-value = 0.02898 # mild deviation, LMM maybe acceptable
check_model(model_c_lmm)

##### Save fixed effects summary----
fixed_df_c <- summary(model_c_lmm)$coefficients %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Coefficient")
# Save to Excel
#write_xlsx(fixed_df_c, "fixed_effects_summary_chao1_interaction_lmm_20250616.xlsx")

##### Save overall model summary----
### Use Anova(model_c_lmm, type = "III") as the main result, and the whole model summary as supplementary
chao1_summary <- Anova(model_c_lmm, type = "III") # main result

# Extract fixed effects summary
fixed_df_chao2 <- as.data.frame(chao1_summary)
# Create a new column 'Coefficient' with row names
fixed_df_chao2$Coefficient <- rownames(fixed_df_chao2)  
# Reorder columns to have 'Coefficient' first:
fixed_df_chao2 <- fixed_df_chao2[, c("Coefficient", setdiff(names(fixed_df_chao2), "Coefficient"))]

# Add significance stars manually
fixed_df_chao2$Signif <- cut(fixed_df_chao2$`Pr(>Chisq)`,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                             labels = c("***", "**", "*", ".", ""))
# Save to Excel
#write_xlsx(fixed_df_chao2, "fixed_effects_summary_chao2_interaction_lmm_20250616.xlsx")

##### Pairwise comparison----
# If interaction is significant, it tells whether the effect of disease varies by sample type
Anova(model_c_lmm, type = "III") # not significant, only interpret main effects directly.

# If Sample_type is significant (yes) → Do pairwise comparisons between sample types (ignoring disease)
pairwise_c_st <- emmeans(model_c_lmm, pairwise ~ Sample_type, adjust = "BH", lmer.df = "satterthwaite")
pairwise_c_st_df <- as.data.frame(pairwise_c_st$contrasts) # Extract the pairwise comparisons dataframe
#write_xlsx(pairwise_c_st_df, "pairwise_comparisons_chao1_st_20250617.xlsx")

# All pair are significant
## no need for disease, as only one pair
# Since Upazila (location) is significant and has >2 levels, follow this 
pairwise_c_up <- emmeans(model_c_lmm, pairwise ~ Upazila, adjust = "BH", lmer.df = "satterthwaite")
pairwise_c_up_df <- as.data.frame(pairwise_c_up$contrasts) # Extract the pairwise comparisons dataframe
#write_xlsx(pairwise_c_up_df, "pairwise_comparisons_chao1_up_20250617.xlsx")

#### Shannon####
#model_s_lmm <- lmerTest::lmer(Shannon ~ Sample_type * Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimate_withPD)
#summary(model_s_lmm)
##### Check residuals####
#plot(model_s_lmm)
#qqnorm(resid(model_s_lmm)); qqline(resid(model_s_lmm)) # If points mostly fall along the line, residuals are approximately normal.
#hist(resid(model_s_lmm), breaks = 30, main = "Residual Histogram (Shannon)", xlab = "Residuals") #Look left skewed
#plot(density(alpha_estimate_withPD$Shannon))
#shapiro.test(resid(model_s_lmm)) # W = 0.9295, p-value = 1.121e-10 # Strong deviation, LMM not acceptable
#check_model(model_s_lmm)
## This needs glmm

#### PD####
model_pd_lmm <- lmerTest::lmer(Phylogenetic_Diversity ~ Sample_type * Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimate_withPD)
summary(model_pd_lmm) # fits better

##### Check residuals#####
plot(model_pd_lmm) # seems normally distributed
qqnorm(resid(model_pd_lmm)); qqline(resid(model_pd_lmm)) # If points mostly fall along the line, residuals are approximately normal.
hist(resid(model_pd_lmm), breaks = 30, main = "Residual Histogram (PD)", xlab = "Residuals") #Look bell shapred
plot(density(alpha_estimate_withPD$Phylogenetic_Diversity)) # Lmm acceptable
shapiro.test(resid(model_pd_lmm)) # W = 0.97991, p-value = 0.0003397 # mild deviation, LMM not acceptable
check_model(model_pd_lmm)
# Shapiro is very sensitive to sample size. without that, Chao1 and PD seems normally distributed but not shannon

##### Save fixed effects summary----
fixed_df_pd <- summary(model_pd_lmm)$coefficients %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Coefficient")
# Save to Excel
#write_xlsx(fixed_df_pd, "fixed_effects_summary_pd_interaction_lmm_20250616.xlsx")

##### Save overall model summary----
pd_summary <- Anova(model_pd_lmm, type = "III") # main result

# Extract fixed effects summary
fixed_df_pd2 <- as.data.frame(pd_summary)
# Create a new column 'Coefficient' with row names
fixed_df_pd2$Coefficient <- rownames(fixed_df_pd2)  
# Reorder columns to have 'Coefficient' first:
fixed_df_pd2 <- fixed_df_pd2[, c("Coefficient", setdiff(names(fixed_df_pd2), "Coefficient"))]

# Add significance stars manually
fixed_df_pd2$Signif <- cut(fixed_df_pd2$`Pr(>Chisq)`,
                           breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                           labels = c("***", "**", "*", ".", ""))
# Save to Excel
#write_xlsx(fixed_df_pd2, "fixed_effects_summary_pd2_interaction_lmm_20250616.xlsx")

##### Pairwise comparison####
# If interaction is significant, it tells whether the effect of disease varies by sample type
pd_summary <- Anova(model_pd_lmm, type = "III") # significant, pairwise comparison for interaction and fixed effect are necessary
pd_summary # Sample_type:Reported_disease is significant (p < 0.05): Do pairwise comparisons within levels of disease or sample type

###### Pairwise comparison for interaction----
emm_pd <- emmeans(model_pd_lmm, ~ Sample_type * Reported_disease, # for pairwise compariosn, no need to add upazila as usually pairwise comparisons focus on the interaction factors (Sample_type and Disease)
                lmer.df = "satterthwaite") # # between diseased and non-diseased within each sample type
pairwise_pd_int <- pairs(emm_pd, by = "Sample_type", adjust = "BH")
pairwise_pd_int

# Convert to data frame
pairwise_pd_int_df <- as.data.frame(pairwise_pd_int)
#write_xlsx(pairwise_pd_int_df, "pairwise_comparisons_pd_interaction_20250617.xlsx")

###### Pairwise comparison for fixed effect----
# Sample_type is significant (yes) → Do pairwise comparisons between sample types (ignoring disease)
pairwise_pd_st <- emmeans(model_pd_lmm, pairwise ~ Sample_type, adjust = "BH", lmer.df = "satterthwaite")
pairwise_pd_st # All pair are significant
pairwise_pd_st_df <- as.data.frame(pairwise_pd_st$contrasts) # Extract the pairwise comparisons dataframe
#write_xlsx(pairwise_pd_st_df, "pairwise_comparisons_pd_st_20250617.xlsx")

## no need for disease, as only one pair
# Since Upazila (location) is significant and has >2 levels, follow this 
pairwise_pd_up <- emmeans(model_pd_lmm, pairwise ~ Upazila, adjust = "BH", lmer.df = "satterthwaite")
pairwise_pd_up_df <- as.data.frame(pairwise_pd_up$contrasts) # Extract the pairwise comparisons dataframe
#write_xlsx(pairwise_pd_up_df, "pairwise_comparisons_pd_up_20250617.xlsx")

### GLMM####
# Use Gamma distribution if Shannon is left-skewed
#### Shannon----
model_s_glmm <- glmmTMB(Shannon ~ Sample_type * Reported_disease + Upazila + (1 | Pond_name),
                        family = Gamma(link = "log"), data = alpha_estimate_withPD)
summary(model_s_glmm)

##### Check residuals----
hist(residuals(model_s_glmm), breaks = 30, main = "Residuals Histogram", xlab = "Residuals")
#hist(residuals(model_shannon_glmm))

# Fitted values and residuals
fitted_vals2 <- fitted(model_s_glmm)
residuals_vals2 <- resid(model_s_glmm, type = "pearson")  # or "deviance", depending on your diagnostic goal
check_model(model_s_glmm)

# Basic plot
plot(fitted_vals2, residuals_vals2,
     xlab = "Fitted values",
     ylab = "Pearson Residuals",
     main = "Fitted vs Residuals")
abline(h = 0, col = "red", lty = 2)

# Check residuals
qqnorm(residuals(model_s_glmm)); qqline(residuals(model_s_glmm))

##### Save fixed effect results----
# Extract the structure of the coefficient matrix
coef_shannon2 <- summary(model_s_glmm)$coefficients
# extract conditional fixed effects matrix
fixed_cond_shannon2 <- coef_shannon2$cond  
# save as dataframe
fixed_df_shannon2 <- data.frame(
  Term = rownames(fixed_cond_shannon2),
  fixed_cond_shannon2,
  row.names = NULL)
# Save to Excel
#write_xlsx(fixed_df_shannon2, "fixed_effects_summary_shannon_interaction_glmm_20250616.xlsx")

##### Save overall model summary----
shannon_summary <- Anova(model_s_glmm, type = "III") # main result
shannon_summary # upazila significant

# Extract fixed effects summary
fixed_df_shannon2 <- as.data.frame(shannon_summary)
# Create a new column 'Coefficient' with row names
fixed_df_shannon2$Coefficient <- rownames(fixed_df_shannon2)  
# Reorder columns to have 'Coefficient' first:
fixed_df_shannon2 <- fixed_df_shannon2[, c("Coefficient", setdiff(names(fixed_df_shannon2), "Coefficient"))]

# Add significance stars manually
fixed_df_shannon2$Signif <- cut(fixed_df_shannon2$`Pr(>Chisq)`,
                                breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                                labels = c("***", "**", "*", ".", ""))
# Save to Excel
#write_xlsx(fixed_df_shannon2, "fixed_effects_summary_shannon2_interaction_lmm_20250616.xlsx")

##### Pairwise comparison----
# If interaction is significant, it tells whether the effect of disease varies by sample type
shannon_summary <- Anova(model_s_glmm, type = "III") # main result
shannon_summary # only upazila significant, pairwise comparison only for upazila is necessary

###### Pairwise comparison for fixed effect----
# Since Upazila (location) is significant and has >2 levels, follow this 
pairwise_s_up <- emmeans(model_s_glmm, pairwise ~ Upazila, adjust = "BH", lmer.df = "satterthwaite")
pairwise_s_up_df <- as.data.frame(pairwise_s_up$contrasts) # Extract the pairwise comparisons dataframe
#write_xlsx(pairwise_s_up_df, "pairwise_comparisons_shannon_up_20250617.xlsx")


## B. Beta diversity (Bray-Curtis)----
st_bray <- pseq %>%
  tax_transform(rank = "unique", trans = "compositional") %>%
  dist_calc(dist = "bray") 

# Sampling_month
st.beta <- st_bray %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    axes = c(1, 2),
    colour = "Sample_type",
    fill = "Sample_type",
    shape = "Reported_disease",
    alpha = 1,
    size = 3
  ) +
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(group = interaction(Sample_type, Reported_disease), colour = Sample_type, linetype = Reported_disease),
    type = "t"
  ) +
  facet_wrap(~Upazila)+
  scale_color_manual(values = sample.colors, labels = c("Gill", "Water", "Skin")) +  # Ensure points and box outlines use the same colors
  scale_fill_manual(values = sample.colors, labels = c("Gill", "Water", "Skin")) +  # Ensure points and box outlines use the same colors
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
       fill = "Sample type",
       #y = "Chao1",  # Add y-axis title to the plot
       #x = "Reported disease",
       title = "PCoA") +
  guides(color = guide_legend(nrow = 1))
st.beta

### PERMANOVA Test----
#Plot PERMANOVA with phyloseq
metadata_st <- sample_data(pseq) %>%
  data.frame() %>%
  tibble()
bray_dist.st = phyloseq::distance(pseq, method="bray")
#### PERMANOVA with all factors----
PERM_all_bray2 <- adonis2(bray_dist.st ~ Sample_type + Reported_disease + Upazila,
                          data = metadata_st, permutations = 999, 
                          strata = metadata_st$Pond_name, 
                          by = "margin") # To get significance of individual terms
PERM_all_bray2

# Extract fixed effects summary
PERM_all_bray2 <- as.data.frame(PERM_all_bray2)
# Create a new column 'Coefficient' with row names
PERM_all_bray2$Coefficient <- rownames(PERM_all_bray2)  
# Reorder columns to have 'Coefficient' first:
PERM_all_bray2 <- PERM_all_bray2[, c("Coefficient", setdiff(names(PERM_all_bray2), "Coefficient"))]

# Save to Excel
#write_xlsx(PERM_all_bray2, "permanova_results_for_all_facotrs_16s_20250616.xlsx")
#write.csv(PERM_all_bray2, file = "permanova_results_for_all_facotrs_16s_20250616.csv", row.names = FALSE)

##### Pairwise adonis----
library(pairwiseAdonis); packageVersion("pairwiseAdonis")
pair.st <- pairwise.adonis(bray_dist.st, metadata_st$Sample_type,  p.adjust.m = "BH")
pair.st
pair.dis <- pairwise.adonis(bray_dist.st, metadata_st$Reported_disease,  p.adjust.m = "BH")
pair.dis
pair.up <- pairwise.adonis(bray_dist.st, metadata_st$Upazila,  p.adjust.m = "BH")
pair.up
pair.sd <- pairwise.adonis(bray_dist.st, metadata_st$Sample_Disease,  p.adjust.m = "BH")
pair.sd

## Save all together
library(openxlsx); packageVersion("openxlsx")

wb <- createWorkbook()
addWorksheet(wb, "Sample_type")
writeData(wb, "Sample_type", as.data.frame(pair.st))

addWorksheet(wb, "Reported_disease")
writeData(wb, "Reported_disease", as.data.frame(pair.dis))

addWorksheet(wb, "Upazila")
writeData(wb, "Upazila", as.data.frame(pair.up))

addWorksheet(wb, "Sample_disease")
writeData(wb, "Sample_disease", as.data.frame(pair.sd))

#saveWorkbook(wb, "pairwise_adonis_all_factors_20250617.xlsx", overwrite = TRUE)

#### Different way to plot beta divestiy ----
ordination <- cmdscale(bray_dist.st, k = 2, eig = TRUE)
ordination_df <- as.data.frame(ordination$points)
ordination_df$Sample_type <- metadata_st$Sample_type
ordination_df$Reported_disease <- metadata_st$Reported_disease
ordination_df$Upazila <- metadata_st$Upazila
# plot
ggplot(ordination_df, aes(x = V1, y = V2, color = Sample_type, shape = Reported_disease)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", level = 0.95) +
  facet_wrap(~Upazila) +  # or ~Reported_disease
  labs(x = "PCoA1", y = "PCoA2") +
  theme_minimal()

# Microeukaryotes----
#load phyloseq object
ps.18s <- readRDS("phyloseq_FSS_BD_metadata_v7_18S_20240703.rds")
ps.18s # 2961 taxa and 57 samples

## A. Alpha diversity----
ps_rarefy.18s <- rarefy_even_depth(ps.18s, rngseed = 1234) #sub-sample to an even sequencing depth (important for alpha diversity, but there is a debate to its importance)
#295OTUs were removed because they are no longer 
#present in any sample after random subsampling
alpha_estimates.18s <- estimate_richness(ps_rarefy.18s, measures = c("Chao1", "Shannon"))
alpha_estimates.18s <- cbind(alpha_estimates.18s, sample_data(ps.18s))

### Calculate phylogenetic diversity----
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

### Statistical modelling for alpha diversity----
# Consider residuals from a simple LMM

## Chao1
model_c.18s <- lmerTest::lmer(Chao1 ~ Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimates.18s)
summary(model_c.18s)
Anova(model_c.18s, type = "III")

## save results as excel file
# Extract fixed effects summary
fixed_df_chao1.18s <- summary(model_c.18s)$coefficients %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Coefficient")
# Save to Excel
#write_xlsx(fixed_df_chao1.18s, "fixed_effects_summary_chao1_18s_lmm_20250612.xlsx")

##### Check residuals#####
plot(model_c.18s) # no specific shape
qqnorm(resid(model_c.18s)); qqline(resid(model_c.18s)) # points mostly fall along the line, residuals are approximately normal.
hist(resid(model_c.18s), breaks = 30, main = "Residual Histogram (Chao1)", xlab = "Residuals") #Look bell shaped
plot(density(alpha_estimates.18s$Chao1)) # Lmm acceptable
shapiro.test(resid(model_c.18s)) # W = 0.97699, p-value = 0.3471 # Normal, LMM maybe acceptable

## Shannon
model_s.18s <- lmerTest::lmer(Shannon ~ Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimates.18s)
summary(model_s.18s)

##### Check residuals#####
plot(model_s.18s)
qqnorm(resid(model_s.18s)); qqline(resid(model_s.18s)) # points mostly fall along the line, residuals are approximately normal.
hist(resid(model_s.18s), breaks = 30, main = "Residual Histogram (Shannon)", xlab = "Residuals") # left skewed
plot(density(alpha_estimates.18s$Shannon)) # Lmm not acceptable
shapiro.test(resid(model_s.18s)) # W = 0.88163, p-value = 4.563e-05 # strongly derivated, LMM not acceptable, try glmm

## PD
model_pd.18s <- lmerTest::lmer(Phylogenetic_Diversity ~ Reported_disease + Upazila + (1 | Pond_name), data = alpha_estimate_withPD.18s)
summary(model_pd.18s)

Anova(model_pd.18s, type = "III")

## save results as excel file
# Extract fixed effects summary
fixed_df_pd.18s <- summary(model_pd.18s)$coefficients %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Coefficient")

# Save to Excel
#write_xlsx(fixed_df_pd.18s, "fixed_effects_summary_pd_18s_lmm_20250612.xlsx")

##### Check residuals#####
plot(model_pd.18s)
qqnorm(resid(model_pd.18s)); qqline(resid(model_pd.18s)) # If points mostly fall along the line, residuals are approximately normal.
hist(resid(model_pd.18s), breaks = 30, main = "Residual Histogram (Chao1)", xlab = "Residuals") #Look bell shaped
plot(density(alpha_estimate_withPD.18s$Phylogenetic_Diversity)) # Lmm acceptable
shapiro.test(resid(model_pd.18s)) # W = 0.95476, p-value = 0.03249 # mild derivated, LMM maybe acceptable

### Try GLMMM for shannon####
#library(glmmTMB); packageVersion("glmmTMB")
#library(DHARMa); packageVersion("DHARMa")

###### Shannon----
model_shannon_glmm.18s <- glmmTMB(
  Shannon ~ Reported_disease + Upazila + (1 | Pond_name),
  family = Gamma(link = "log"), data = alpha_estimate_withPD.18s)
summary(model_shannon_glmm.18s)

Anova(model_shannon_glmm.18s, type = "III")

## save results as excel file
# Extract the structure of the coefficient matrix
coef_shannon.18s <- summary(model_shannon_glmm.18s)$coefficients

# extract conditional fixed effects matrix
fixed_cond_shannon.18s <- coef_shannon.18s$cond  
# save as dataframe
fixed_df_shannon.18s <- data.frame(
  Term = rownames(fixed_cond_shannon.18s),
  fixed_cond_shannon.18s,
  row.names = NULL)

# Save to Excel
#write_xlsx(fixed_df_shannon.18s, "fixed_effects_summary_shannon_18s_glmm_20250612.xlsx")

# Histogram of residuals
hist(residuals(model_shannon_glmm.18s), breaks = 30, main = "Residuals Histogram", xlab = "Residuals")
#hist(residuals(model_shannon_glmm))

# Fitted values and residuals
fitted_vals2.18s <- fitted(model_shannon_glmm.18s)
residuals_vals2.18s <- resid(model_shannon_glmm.18s, type = "pearson")  # or "deviance", depending on your diagnostic goal

# Basic plot
plot(fitted_vals2.18s, residuals_vals2.18s,
     xlab = "Fitted values",
     ylab = "Pearson Residuals",
     main = "Fitted vs Residuals")
abline(h = 0, col = "red", lty = 2)

# Check residuals
qqnorm(residuals(model_shannon_glmm.18s)); qqline(residuals(model_shannon_glmm.18s))


#### PERMANOVA with all factors----
metadata_water.18s <- sample_data(ps.18s) %>%
  data.frame() %>%
  tibble()
bray_dist.water.18s = phyloseq::distance(ps.18s, method="bray")
PERM_all_bray.18s <- adonis2(bray_dist.water.18s ~ Reported_disease + Upazila,
                          data = metadata_water.18s, permutations = 999, 
                          strata = metadata_water.18s$Pond_name, 
                          by = "margin") # To get significance of individual terms
PERM_all_bray.18s # p > 0.05, no pairwise comparison needed

# Extract fixed effects summary
PERM_all_bray.18s <- as.data.frame(PERM_all_bray.18s)
# Create a new column 'Coefficient' with row names
PERM_all_bray.18s$Coefficient <- rownames(PERM_all_bray.18s)  
# Reorder columns to have 'Coefficient' first:
PERM_all_bray.18s <- PERM_all_bray.18s[, c("Coefficient", setdiff(names(PERM_all_bray.18s), "Coefficient"))]

# Save to Excel
#write_xlsx(PERM_all_bray.18s, "permanova_results_for_all_factors_18s_20250617.xlsx")







