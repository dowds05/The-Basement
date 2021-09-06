# Update to DESeq2 script on 2/17/21 
# UPDATE: contrast, species, factors
##################################################################################################################
#-------------- DESeq2 UPDATE ------------------------------------------------------------------------------------
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
library(gplots)
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)

BiocManager::install("apeglm")

##################################################################################################################
#-------------- Importing Datasets -------------------------------------------------------------------------------
##################################################################################################################
# UPDATE FACTOR!
merge <- bind_rows(baseline, w8)
# verify levels of merge: Activity, Sample, etc.
levels(as.factor(merge$Time))
##################################################################################################################
#-------------- PIVOT TABLES TO COUNT READS BY SAMPLE ------------------------------------------------------------
##################################################################################################################
pivot <-  merge %>% group_by(Strain, Sample) %>% tally() %>% 
  pivot_wider(names_from = "Sample", values_from = n, names_sort=TRUE) %>% mutate_all(~replace(., is.na(.), 0))
dim(pivot)
# Total number of reads per sample
colSums(pivot[,2:227])
# Mean, Median, Range of samples
summary(colSums(pivot[,2:227]))
#-------------------- From long-to-wide format -------------------------------------------------------------------
# wide format 
wide <- pivot %>% t() %>% .[-1,]
class(wide) <- "numeric"
# Make a data frame
wide <- as.data.frame(wide)
# check row sums, confirm rowSums=colSums
rowSums(wide)
# Change column names to Description (or 'Strain/OTU/BBH')
colnames(wide) <- pivot$Strain
#----------- Transpose this file into long format, if needed -----------------------------------------------------
long <- as.data.frame(t(wide))
colSums(long)
##################################################################################################################
#-------------- Phyloseq -----------------------------------------------------------------------------------------
##################################################################################################################
otu_mat <- long %>% rownames_to_column("Strain")
colSums(otu_mat[,-1])
# make tax_otu file (BBHs as rows, taxnomoy as columns)
tax_mat1 <- l %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Strain, Strain) %>% group_by(Strain)
# This filters tax_mat1 by Strains in otu_mat
tax_mat <- tax_mat1 %>% group_by(Strain) %>% filter(Strain %in% otu_mat$Strain) %>% distinct(.keep_all = FALSE)
# UPDATE FACTOR! Make target for subsetting samples_df.
# verify levels of merge for target
levels(as.factor(merge$Time))
target <- c("D_0")
# Subset of the GS2.env file, subset Metadata file for conformable arrays
samples_df <- as.data.frame(env) %>% select(OG_Sample, Sample, Activity, Diet, Sex, Time, Treatment)
samples_df <- filter(samples_df, Time %in% target)
# make row.names Sample #
row.names(samples_df) <- samples_df$Sample
samples_df <- samples_df %>% select (-Sample)
# paste0 'sample' to eliminate numbers in sample name
rownames(samples_df) <- paste0("Sample",1:nrow(samples_df))
samples_df
# --------------------- Pre-processing for Phyloseq --------------------------------------------------------------
# Assign Strain as. row.names()
row.names(otu_mat) <- otu_mat$Strain
otu_mat <- otu_mat %>% select(-Strain)
# paste0 'sample' to eliminate numbers in sample name
colnames(otu_mat) <- paste0("Sample",1:ncol(otu_mat))
colSums(otu_mat)
# use dplyr distinct to remove rows with duplicate values
tax_mat <- tax_mat %>%  filter(Strain %in% pivot$Strain) %>% distinct(Strain, .keep_all = TRUE)
# Assign Strain as row.names()
row.names(tax_mat) <- tax_mat$Strain
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
#############################################################################################
#-------------- DIFFERNTIAL ABUNDANCE (DA) TESTING ------------------------------------------
#############################################################################################
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html
# UPDATE FACTOR! Differential abundance testing. ?phyloseq_to_deseq2
diagdds = phyloseq_to_deseq2(ps, ~ Sex)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# new DESeq object for factor labeling. ?DESeq
dds <- DESeq(diagdds)
# UPDATE FACTOR! verify factor levels for setting up contrasts.
levels(as.factor(samples_df$Sex))
# UPDATE FACTOR! Create results df. contrast <- c("condition", "level_to_compare", "base_level")
contrast <- c("Sex", "Male", "Female")
# Extract results. ?results
res <- results(dds, contrast = contrast)
res
#res %>% data.frame() %>% view()

# Alpha to sort res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
dim(sigtab)
# UPDATE FACTOR! write sigtab to .csv to make figures in Excel. Rename with appropriate factor
write.csv(sigtab, file = "Female_vs_Male.baseline.da.csv")

# UPDATE FACTOR & CONTRAST! Note type='apeglm' & 'ashr' outperform OG 'normal' shrinkage estimator. ?lfcShrink
shrunk <- lfcShrink(dds, contrast = contrast, type = 'normal')
# MA plot shows the mean of the normalized counts versus the log2FoldChange for all genes tested
# MA plot unsrunken fold changes. ?plotMA
plotMA(res, ylim=c(-5,5), main = "res")
# MA plot  shrunken fold changes
plotMA(shrunk, ylim=c(-5,5), main = "shrunk")

## Volcano plot
#ggplot(res) + geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
#  ggtitle("Mov10 overexpression") + xlab("log2 fold change") + ylab("-log10 adjusted p-value")

# UPDATE FACTOR & CONTRAST! title = 'base_level vs. level_to_compare'
EnhancedVolcano(shrunk, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', 
                title = 'Female vs Male', pCutoff = 10e-3, FCcutoff = 1, 
                pointSize = 3.0, labSize = 4.0)

# Check factor labels = data
# Search for different Strain
check <- merge %>% filter(str_detect(Strain, "^Sphingobium_yanoikuyae_strain_CD09_2_CD09_2"))
dim(check)
# PIVOT TABLE on BBH
check.pivot <- check %>% group_by(Treatment, Sex) %>%
  tally() %>% pivot_wider(names_from = "Treatment", values_from = n, names_sort=TRUE) %>% 
  mutate_all(~replace(., is.na(.), 0))
dim(check.pivot)
view(check.pivot)

# REMOVE files to look at other factors
rm(diagdds, dds, contrast, res, alpha, sigtab, shrunk, check, check.pivot)


# REMOVE Phyloseq data for next analysis
rm(check, check.pivot, target, merge, ps, otu_mat, samples, samples_df, tax_mat, tax_mat1, OTU, TAX, long, pivot, wide)

