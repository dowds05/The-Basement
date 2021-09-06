##################################################################################################################
# Load all applicable libraries for PERMANOVA, diversity, compostition analysis-----------------------------------
# libraries() include: tidyverse, vegan, ade4, cluster, FD, 
library(tidyverse)
library(vegan)
library(ade4)
library(cluster)
library(FD)
library(lattice)
library(permute)
library(readxl)
library(adespatial)
library(gclus)
library(plyr)
library(data.table)
library(phyloseq)
library(DECIPHER)
library(janitor)
library(devtools)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(microbiome)

BiocManager::install("microbiome")

##################################################################################################################
#-------------- Importing Datasets -------------------------------------------------------------------------------
##################################################################################################################
#------------import a csv file------------------------------------------------------------------------------------
GS2 <- read_csv("~/Desktop/Bob_Aim3/R/DATA/tot.csv", col_names = FALSE)
# Import environmental file for PERMANOVA. It's analogous to a metadata or map file in QIIME
env <- read_csv("~/Desktop/Bob_Aim3/R/DATA/env.csv", col_names = TRUE)
# use unite to make treatment column of effects
env <- env %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
glimpse(GS2)
# prior analysis ID 630 unique taxa
unique1 <- read_csv("~/Desktop/Bob_Aim3/R/DATA/unique.csv", col_names = TRUE)
colnames(unique1) <- c('Strain')
# use filter to subset data for network analysis
unique <- l %>% group_by(Strain) %>% filter(Strain %in% unique1$Strain) %>%
   distinct(.keep_all = TRUE)
unique <- unique %>% unite(Treatment, Activity, Diet, Sex, sep = "_", remove = FALSE)
##################################################################################################################
#-------------- PIVOT TABLES TO COUNT READS BY SAMPLE ------------------------------------------------------------
##################################################################################################################
pivot <-  unique %>% group_by(Phylum, Sample) %>% tally() %>% 
  pivot_wider(names_from = "Sample", values_from = n, names_sort=TRUE) %>% mutate_all(~replace(., is.na(.), 0))
dim(pivot)
# Total number of reads per sample
colSums(pivot[,2:227])
# Mean, Median, Range of 57 samples
summary(colSums(pivot[,2:227]))
#-------------------- From long-to-wide format -------------------------------------------------------------------
# wide format 
wide <- pivot %>% t() %>% .[-1,]
class(wide) <- "numeric"
# Make a data frame
wide <- as.data.frame(wide)
# check row sums, confirm rowSums=colSums
rowSums(wide)
# Change column names to Description (or 'species/OTU/BBH')
colnames(wide) <- pivot$Phylum
#----------- Transpose this file into long format, if needed ----------------------------------------------------
long <- as.data.frame(t(wide))
colSums(long)
##################################################################################################################
#-------------- Phyloseq Tutorial --------------------------------------------------------------------------------
##################################################################################################################
# https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
# we need three files:
# 1) long format (OTUs = rows, samples = columns) = otu_mat
# 2) OTU taxonomy (OTU, Kingdom, phylyum, etc) = tax_mat
# 3) Samples = samples_df
otu_mat <- long %>% rownames_to_column("Phylum")
colSums(otu_mat[,-1])
# make tax_otu file (BBHs as rows, taxnomoy as columns)
tax_mat1 <- l %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Phylum) %>% group_by(Phylum)
# This filters tax_mat1 by Strains in otu_mat. ?filter()
tax_mat <- tax_mat1 %>% group_by(Phylum) %>% filter(Phylum %in% otu_mat$Phylum) %>% distinct(.keep_all = FALSE)
# UPDATE WEEK! Make target for subsetting samples_df, 
#target <- c("D_0", "Week_3")
# Subset of the GS2.env file, subset Metadata file for conformable arrays
samples_df <- as.data.frame(env) %>% select(Sample, Sample, Activity, Diet, Sex, Time, Treatment)
#samples_df <- filter(samples_df, Time %in% target)
# use unite to make treatment column of effects
samples_df <- samples_df %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
row.names(samples_df) <- samples_df$Sample
#samples_df <- samples_df %>% select (-Sample)
# paste0 'sample' to eliminate numbers in sample name. ?paste0
rownames(samples_df) <- paste0("Sample",1:nrow(samples_df))
samples_df
# --------------------- Pre-processing for Phyloseq ---------------------------------------------------------------
# Assign Species as. row.names()
row.names(otu_mat) <- otu_mat$Phylum
otu_mat <- otu_mat %>% select(-Phylum)
# paste0 'sample' to eliminate numbers in sample name
colnames(otu_mat) <- paste0("Sample",1:ncol(otu_mat))
colSums(otu_mat)
# use dplyr distinct to remove rows with duplicate values, if needed. ?filter ?distinct
tax_mat <- tax_mat %>%  filter(Phylum %in% pivot$Phylum) %>% distinct(Phylum, .keep_all = TRUE)
# Assign Species as row.names()
row.names(tax_mat) <- tax_mat$Phylum
# Matrixes for otu_mat & tax_mat
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
# files for phyloseq(). ?sample_data
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
# Phyloseq uses the phyloseq() command to rapidly combine the otu_table, tax_table and samples
ps <- phyloseq(OTU, TAX, samples)
ps
# Visualize data
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)
##################################################################################################################
#-------------- SHINY PHYLOSEQ -----------------------------------------------------------------------------------
##################################################################################################################
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#references
# To save R objects (e.g. ps, pivot, etc.) as .RData file for Shiny-phyloseq
save(ps, otu_mat, tax_mat, samples, OTU, TAX, file = "unique.RData")
# This will open Shiny Phyloseq, a GUI that allows for rapid figure generation
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")
# REMOVE Shiny Phyloseq data
rm(datalist, distlist, env_psdata, pkg1, RasterRstudio, shiny_phyloseq_ggtheme_list, theme_blank_custom,
   animation_steps, d3DefaultDistance, d3DefaultLinkScaleFactor, d3NetworkColorVar, d3NodeLabelVar,
   default_netLabel, deppkgs, graphicFormats, i, interval, kovera_A, kovera_k, LinkDistThreshold,
   loop, netdist, netThreshColorVariableDefault, netThreshDistanceMethod, netThreshShapeVariableDefault,
   OTUSumDefault, R_min_version, R_version, rasterGraphicFormats, RstudioPNGsp, SampleSumDefault, vectorGraphicFormats,
   av, component_options, dist_to_edge_table, fail_gen, get_facet_grid, ggfilegen, ggsave2, 
   install_missing_packages, numericInputRow, output_phyloseq_print_html, scaled_distance, shiny_phyloseq_print,
   simpletime, tablify_phyloseq_component, textInputRow)

# --------------------- Pre-processing for Phyloseq --------------------------------------------------------------
?transform
ps.t <- microbiome::transform(ps, transform = "log10", target = "samples")
?ordinate

#------------import a csv file------------------------------------------------------------------------------------
# NETWORK analysis using 630 unique taxa. https://joey711.github.io/phyloseq/plot_network-examples.html
?plot_net
plot_net(ps, distance = "bray", color="Treatment", shape = "Treatment", maxdist = 0.25, type = 'samples', 
         point_label = 'Sample', point_size = 6)

plot_net(ps, distance = "bray", maxdist = 0.2, type = 'taxa')


plot_net(ps, maxdist = 0.2, type = 'taxa', color = 'Diet')

?make_network
ig <- make_network(ps, type = "taxa", dist.fun="bray", max.dist=0.25)

?plot_network
plot_network(ig, ps, type = 'taxa')

plot_network(ig, ps, color="Diet", shape="Diet", line_weight=0.4, label=NULL)
plot_net(ig, ps, )


ig <- make_network(ps, max.dist=0.2)

plot_network(ig, ps, color="Activity", shape="Activity", line_weight=0.4, label=NULL)



plot_heatmap(ps, sample.label="Treatment")



# Relative Abundance Histograms  ---------------------------------------------------------------------------------#

data(peerj32)
p <- boxplot_abundance(peerj32$phyloseq, x='time', y='Akkermansia',
                       line='subject')
p

boxplot_abundance(ps, x='Activity', y='Muribaculum_intestinale_strain_YL27', line='subject')
ps


# REMOVE Phyloseq data
rm(pivot, wide, long, otu_mat, pie, ps, samples, samples_df, tax_mat, tax_mat1, OTU, TAX)

# REMOVE Phyloseq data
rm(ig, q1, unique1, unique)
