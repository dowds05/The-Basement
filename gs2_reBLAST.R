#############################################################################################
# Load all applicable libraries for PERMANOVA, diversity, compostition analysis--------------
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

# Use install.packages('package_name') to install packages with base R commands
install.packages('RColorBrewer')
# if install.packages() doesnt work, try this:
if (!requireNamespace("Stack", quietly = TRUE))
  install.packages("Stack")
# to bypass install.packages() install with BiocManager::install("package_name")
BiocManager::install("rgr")
#############################################################################################
#-------------- Importing Datasets ----------------------------------------------------------
#############################################################################################
#------------import a csv file---------------------------------------------------------------
GS2 <- read_csv("~/Desktop/Bob_Aim3/R/DATA/tot.csv", col_names = FALSE)
# Import environmental file for PERMANOVA. It's analogous to a metadata or map file in QIIME
env <- read_csv("~/Desktop/Bob_Aim3/R/DATA/env.csv", col_names = TRUE)
# use unite to make treatment column of effects
env <- env %>% unite(Treatment, Activity, Diet, Sex, sep = "_", remove = FALSE)
glimpse(GS2)
# delete columns. re-BLAST made two identical columns for subject & phylogeny. delete copy
GS2 = subset(GS2, select = -c(X2,X6,X7))
# Change column names from re-BLAST (10/2020)
colnames(GS2) <- c('Query','Subject', 'Percent','Length','Sample','Activity','Diet','Sex','Time', 'OG_Sample')
glimpse(GS2)
# Split phylogeny by ";" using separate
?separate()
GS2 <- GS2 %>% 
  separate(col=Subject, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
glimpse(GS2)
# Filter by sequences >= 1000 bp
l <- GS2 %>% filter(Length >= 1000) 
# Make Treatment column. Note without Sex
l <- l %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
# Double check how many samples and levels of each factor with levels
samples <- levels(as.factor(l$Sample)) 
levels(as.factor(l$OG_Sample))
levels(as.factor(l$Activity))
levels(as.factor(l$Diet))
levels(as.factor(l$Sex))
levels(as.factor(l$Time))
# remove Limnobacter from the analysis
search <- l %>% filter(str_detect(Strain, "^Limnobacter_"))
# strains of Limnobacter
levels(as.factor(search$Strain))
# Delete any reads that have Limnobacter from analysis
?filter
l <- l %>% filter(Strain != c('Limnobacter_sp._130_strain_Limnobacter_sp._130'))
l <- l %>% filter(Strain != c('Limnobacter_sp._MED105__________'))
l <- l %>% filter(Strain != c('Limnobacter_sp._SAORIC-580'))
l <- l %>% filter(Strain != c('Limnobacter_sp._SAORIC-690'))

#-------------------- df of each week ----------------------------------------------------------------------------
# Make Treatment column. Note with Sex
GS2 <- GS2 %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
l <- l %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
env <- env %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)

baseline <- l %>% filter(Time == "D_0")
w3 <- l %>% filter(Time == "Week_3")
w4 <- l %>% filter(Time == "Week_4")
w8 <- l %>% filter(Time == "Week_8")
w12 <- l %>% filter(Time == "Week_12")

#-------------------- write.csv for l ----------------------------------------------------------------------------
write.csv(l, file = "l.csv")
write.csv(baseline, file = "baseline.csv")
write.csv(w3, file = "w3.csv")
write.csv(w4, file = "w4.csv")
write.csv(w8, file = "w8.csv")
write.csv(w12, file = "w12.csv")
#############################################################################################
#-------------- PIVOT TABLES TO COUNT READS BY SAMPLE ---------------------------------------
#############################################################################################
# PIVOT TABLE on QA/QC reads
# Pivot on w12
pivot <- l %>% group_by(Strain, Sample) %>% tally() %>% 
  pivot_wider(names_from = "Sample", values_from = n, names_sort=TRUE) %>% mutate_all(~replace(., is.na(.), 0))
dim(pivot)
#-------------- EXPORT AS CSV --------------------------------------------------------------------------------------
write.csv(env, file = "env.prism.csv") # Use name from search

# Total number of reads per sample
abun <- colSums(pivot[,2:227])
hist(abun)
# dataset:
data=data.frame(value=rnorm(100))
# ggplot2 for total abundance histograms
p <- ggplot(abun, aes(x=value)) + 
  geom_histogram()

# fastQ
fastq <- colSums(pivot[,2:227])
# fastA
fasta <- colSums(pivot[,2:227])
# BBH
bbh <- colSums(pivot2[,2:227])
# QA/QC 
quality_alignments <- colSums(pivot[,2:227])
# Mean, Median, Range of 57 samples
summary(colSums(pivot[,2:225]))
#--------------    Write data files         ?write.csv()    --------------------------------------------
write.csv(fastq, file = "fastq_counts.csv")
write.csv(fasta, file = "fasta_counts.csv")
write.csv(quality_alignments, file = "quality_alignments_counts")
write.csv(bbh, file = "bbh_counts")
write.csv(pivot, file = "unique_counts.csv")
#-------------- DELETE IF USING ANOTHER R SCRIPT FILE --------------------------------------------------------------
rm(fastq, fasta, bbh, quality_alignments, pivot2)
#-------------------- From long-to-wide format ----------------------------------------------
# wide format 
wide <- pivot %>% t() %>% .[-1,]
class(wide) <- "numeric"
# Make a data frame
wide <- as.data.frame(wide)
# check row sums, confirm rowSums=colSums
rowSums(wide)
# Change column names to Description (or 'species/OTU/BBH')
colnames(wide) <- pivot$Strain
#----------- Transpose this file into long format, if needed -------------------------------
long <- as.data.frame(t(wide))
colSums(long)
#############################################################################################
#-------------- Phyloseq Tutorial -----------------------------------------------------------
#############################################################################################
# https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
# we need three files:
# 1) long format (OTUs = rows, samples = columns) = otu_mat
# 2) OTU taxonomy (OTU, Kingdom, phylyum, etc) = tax_mat
# 3) Samples = samples_df
otu_mat <- long %>% rownames_to_column("Strain")
colSums(otu_mat[,-1])
# make tax_otu file (BBHs as rows, taxnomoy as columns)
tax_mat1 <- l %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain) %>% 
  group_by(Strain)
# This filters tax_mat1 by Strains in otu_mat
?filter()
tax_mat <- tax_mat1 %>% group_by(Strain) %>% filter(Strain %in% otu_mat$Strain) %>%
  distinct(.keep_all = FALSE)
# Subset of the GS2.env file, subset Metadata file for conformable arrays
samples_df <- as.data.frame(env) %>% 
  select(OG_Sample, Sample, Activity, Diet, Sex, Time, Treatment)
# use unite to make treatment column of effects
samples_df <- samples_df %>% unite(Treatment, Activity, Diet, sep = "_", remove = FALSE)
row.names(samples_df) <- samples_df$Sample
samples_df <- samples_df %>% select (-Sample)
# paste0 'sample' to eliminate numbers in sample name
?paste0
rownames(samples_df) <- paste0("Sample",1:nrow(samples_df))
samples_df
# --------------------- Pre-processing for Phyloseq ------------------------------------------
# Assign Species as. row.names()
row.names(otu_mat) <- otu_mat$Strain
otu_mat <- otu_mat %>% select(-Strain)
# paste0 'sample' to eliminate numbers in sample name
colnames(otu_mat) <- paste0("Sample",1:ncol(otu_mat))
colSums(otu_mat)
# use dplyr distinct to remove rows with duplicate values, if needed
?filter()
?distinct()
tax_mat <- tax_mat %>%  filter(Strain %in% pivot$Strain) %>% distinct(Strain, .keep_all = TRUE)
# Assign Species as row.names()
row.names(tax_mat) <- tax_mat$Strain
# Matrixes for otu_mat & tax_mat
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
# files for phyloseq()
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
?sample_data
samples = sample_data(samples_df)
# Phyloseq uses the phyloseq() command to rapidly combine the otu_table, 
# tax_table and samples
ps <- phyloseq(OTU, TAX, samples)
ps
# Visualize data
sample_names(ps)
rank_names(ps)
sample_variables(ps)
otu_table(ps)
#############################################################################################
#-------------- SHINY PHYLOSEQ --------------------------------------------------------------
#############################################################################################
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#references
# To save R objects (e.g. ps, pivot, etc.) as .RData file for Shiny-phyloseq
save(ps, otu_mat, tax_mat, samples, OTU, TAX, file = "gs2.RData")
# This will open Shiny Phyloseq, a GUI that allows for rapid figure generation
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")
# remove Shiny Phyloseq data
rm(datalist, distlist, env_psdata, pkg1, RasterRstudio, shiny_phyloseq_ggtheme_list, theme_blank_custom,
   animation_steps, d3DefaultDistance, d3DefaultLinkScaleFactor, d3NetworkColorVar, d3NodeLabelVar,
   default_netLabel, deppkgs, graphicFormats, i, interval, kovera_A, kovera_k, LinkDistThreshold,
   loop, netdist, netThreshColorVariableDefault, netThreshDistanceMethod, netThreshShapeVariableDefault,
   OTUSumDefault, R_min_version, R_version, rasterGraphicFormats, RstudioPNGsp, SampleSumDefault, vectorGraphicFormats,
   av, component_options, dist_to_edge_table, fail_gen, get_facet_grid, ggfilegen, ggsave2, 
   install_missing_packages, numericInputRow, output_phyloseq_print_html, scaled_distance, shiny_phyloseq_print,
   simpletime, tablify_phyloseq_component, textInputRow)
#############################################################################################
#-------------- RELATIVE ABUNDANCE HISTOGRAMS -----------------------------------------------
#############################################################################################
# Absolute abundance plot
plot_bar(ps, levels = level_order)
# Remove taxonomic levels filled with NAs.  ?subset_taxa
subset <- subset_taxa(ps, Phylum != "NA")

# transform counts to relative abundance
ra <- transform_sample_counts(subset, function(x) x*100/sum(x))
# agglomerate counts. UPDATE TAXONOMY! ?tax_glom
tg <- tax_glom(ra, "Class")
# next get the otu_table. UPDATE TAXONOMY!
propData <- as.data.frame(t(otu_table(tg)))

# get melted dataframe. ?melt
plotDF <- propData %>%
  rownames_to_column(var="Sample") %>%
  melt() %>%
  magrittr::set_names(c("Sample", "taxa", "relab"))

## Using sample as id variables. add taxonomy to plotting DF. UPDATE TAXONOMY!
taxonomy.2 <- tax_mat %>%
  as.data.frame() %>%
  rownames_to_column(var="taxa")
plotDF <- left_join(plotDF, taxonomy.2, by="taxa") %>%
  select(Sample, taxa, relab, Class)

# Warning: Column `taxa` joining factor and character vector, coercing into character vector
# aggregate all unclassified taxa levels. UPDATE TAXONOMY! ?mutate
tot_unc <- plotDF %>%
  group_by(Sample) %>%
  filter(grepl("unclassified", Class)) %>%
  summarise(relab=sum(relab)) %>%
  mutate(taxa = "unclassified", Class = "unclassified") 

# add the unclassified aggregated rel ab to plotting data frame. UPDATE TAXONOMY!
plotDF <- plotDF %>%
  filter(!grepl("unclassified", Class)) %>%
  bind_rows(tot_unc)
#-------------- EXPORT AS CSV --------------------------------------------------------------------------------------
write.csv(plotDF, file = "fig1_data.csv") 

# add metadata to plotting DF
df <- samples_df %>%
  rownames_to_column(var="Sample")
plotDF <- left_join(plotDF, df, by="Sample")
# check sample levels
levels(as.factor(plotDF$Sample))
# lets plot stacked bars using ggplot2. UPDATE TAXONOMY!
mycolors <- rev(colorRampPalette(brewer.pal(12, "Paired"))(length(unique(plotDF$Class))))
#this vector might be useful for other plots/analyses
level_order <- c('Sample225', 'Sample226', 'Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6', 'Sample7', 'Sample8',
                 'Sample9', 'Sample10', 'Sample11', 'Sample12', 'Sample13', 'Sample14', 'Sample15', 'Sample16',
                 'Sample17', 'Sample18', 'Sample19', 'Sample20', 'Sample21', 'Sample22', 'Sample23', 'Sample24',
                 'Sample25', 'Sample26', 'Sample27', 'Sample28', 'Sample29', 'Sample30', 'Sample31', 'Sample32',
                 'Sample33', 'Sample34', 'Sample35', 'Sample36', 'Sample37', 'Sample38', 'Sample39', 'Sample40',
                 'Sample41', 'Sample42', 'Sample43', 'Sample44', 'Sample45', 'Sample46', 'Sample47', 'Sample48',
                 'Sample49', 'Sample50', 'Sample51', 'Sample52', 'Sample53', 'Sample54', 'Sample55', 'Sample56',
                 'Sample57', 'Sample58', 'Sample59', 'Sample60', 'Sample61', 'Sample62', 'Sample63', 'Sample64',
                 'Sample65', 'Sample66', 'Sample67', 'Sample68', 'Sample69', 'Sample70', 'Sample71', 'Sample72',
                 'Sample73', 'Sample74', 'Sample75', 'Sample76', 'Sample77', 'Sample78', 'Sample79', 'Sample80',
                 'Sample81', 'Sample82', 'Sample83', 'Sample84', 'Sample85', 'Sample86', 'Sample87', 'Sample88',
                 'Sample89', 'Sample90', 'Sample91', 'Sample92', 'Sample93', 'Sample94', 'Sample95', 'Sample96',
                 'Sample97', 'Sample98', 'Sample99', 'Sample100', 'Sample101', 'Sample102', 'Sample103', 'Sample104',
                 'Sample105', 'Sample106', 'Sample107', 'Sample108', 'Sample109', 'Sample110', 'Sample111', 'Sample112',
                 'Sample113', 'Sample114', 'Sample115', 'Sample116', 'Sample117', 'Sample118', 'Sample119', 'Sample120',
                 'Sample121', 'Sample122', 'Sample123', 'Sample124', 'Sample125', 'Sample126', 'Sample127', 'Sample128',
                 'Sample129', 'Sample130', 'Sample131', 'Sample132', 'Sample133', 'Sample134', 'Sample135', 'Sample136',
                 'Sample137', 'Sample138', 'Sample139', 'Sample140', 'Sample141', 'Sample142', 'Sample143', 'Sample144',
                 'Sample145', 'Sample146', 'Sample147', 'Sample148', 'Sample149', 'Sample150', 'Sample151', 'Sample152',
                 'Sample153', 'Sample154', 'Sample155', 'Sample156', 'Sample157', 'Sample158', 'Sample159', 'Sample160',
                 'Sample161', 'Sample162', 'Sample163', 'Sample164', 'Sample165', 'Sample166', 'Sample167', 'Sample168',
                 'Sample169', 'Sample170', 'Sample171', 'Sample172', 'Sample173', 'Sample174', 'Sample175', 'Sample176',
                 'Sample177', 'Sample178', 'Sample179', 'Sample180', 'Sample181', 'Sample182', 'Sample183', 'Sample184',
                 'Sample185', 'Sample186', 'Sample187', 'Sample188', 'Sample189', 'Sample190', 'Sample191', 'Sample192',
                 'Sample193', 'Sample194', 'Sample195', 'Sample196', 'Sample197', 'Sample198', 'Sample199', 'Sample200',
                 'Sample201', 'Sample202', 'Sample203', 'Sample204', 'Sample205', 'Sample206', 'Sample207', 'Sample208',
                 'Sample209', 'Sample210', 'Sample211', 'Sample212', 'Sample213', 'Sample214', 'Sample215', 'Sample216',
                 'Sample217', 'Sample218', 'Sample219', 'Sample220', 'Sample221', 'Sample222', 'Sample223', 'Sample224')
# Relative Abundance Histogram. UPDATE TAXONOMY
ggplot(plotDF, aes(x = factor(Sample, levels = level_order), y = relab,  fill=Phylum)) + 
  geom_bar(stat="identity", position = "stack") + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "bottom") +
  ylab("Relative Abundance") +
  xlab("Samples") +
  ggtitle("Relative Abundance") +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))

# Baseline Samples
female <- plotDF %>% filter(Sample == "Sample1")
male <- plotDF %>% filter(Sample == "Sample2")
baseline_plot <- bind_rows(female, male)
# Baseline Relative Abudnance Histogram
bp <- ggplot(baseline_plot, aes(x = factor(Sample), y = relab,  fill=Class)) + 
  geom_bar(stat="identity", position = "stack") + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "bottom") +
  ylab("Relative Abundance") +
  xlab("Samples") +
  ggtitle("Baseline Relative Abundance") +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))
bp
# Phyla level prevlance figure
# Create table, number of features for each phyla
percent <- table(tax_table(ps)[, "Phylum"], exclude = NULL)
# write.csv for percentages
write.csv(percent, file = "percent.csv")
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(subset),
               MARGIN = ifelse(taxa_are_rows(subset), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(subset),
                    tax_table(subset))
# Low prevalance phyla?
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Total Abundance vs. Prevalence[Frac. Samples]
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(subset),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")



# REMOVE relative abundance data:
rm(subset, ra, tg, plotDF, propData, taxonomy.2, tot_unc, level_order, df,
   mycolors, prevdf, male, female, baseline, bp)
# -------------- HEATMAPS -----------------------------------------------------------------

# plot_heatmap to make heatmaps
?plot_heatmap
plot_heatmap(ps, method = "NMDS", distance = "bray")

plot_heatmap(ps, method = "NMDS", distance = "bray", 
             taxa.label = "Strain", taxa.order = "Strain", 
             low="beige", high="red", na.value="beige")

# Plot Chao1 richness estimator and Shannon diversity estimator
plot_richness(subset, measures=c("Chao1", "Shannon"), x="sample", color = "Treatment") + 
  geom_point(size=3) + theme_gray()

#-------------- ORDINATION -----------------------------------------------------------------
nmds <- ordinate(ps, "NMDS", "bray")
p1 = plot_ordination(ps, w12.ord, type="taxa", color="Phylum", title="Week 3 taxa")
p1

# noisy plots, use facet_wrap() to split data to make more intuitive
p1 + facet_wrap(~Phylum, 3)


# There are many built-in distances that can be used
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# PCoA NOTE: TRY different distances (e.g. "bray", "jaccard", "kulczynski", "gower", etc.)
pcoa <- ordinate(ps, method = "PCoA", distance = "bray")
pcoa
# Plot BBHs of above distance
?plot_ordination
?geom_point()
#


pcoa.plot <- plot_ordination(ps, pcoa, type="samples", color="Diet", shape="Treatment",
                             title="Week 3 PCoA, Bray Distance")
pcoa.plot + geom_point(size=3) + facet_wrap(~Sex, 3) 


# new pcoa plot
?scale_shape_manual()
?stroke
pcoa.plot2 <- plot_ordination(ps, pcoa, type = "samples", color = "Treatment", shape = "Treatment",
                              title = "Week 3 PCoA")
pcoa.plot2 + geom_point(size=4) +
  scale_shape_manual(values=c(0, 16, 17, 18)) +
  scale_color_manual(values=c('#33FF00','#CCFF99', '#CC9900', '#FF3300')) +
  theme(legend.position="right") +
  facet_wrap(~Sex, 3) +
  theme_gray()

pcoa.plot2 + geom_point(size=4, fill) +
  scale_shape_manual(values=c(20, 21, 22, 23)) +
  scale_color_manual(values=c('#33FF00','#CCFF99', '#CC9900', '#FF3300')) +
  facet_wrap(~Sex, 3) +
  ggtitle("20, 21, 22, 23")
# PERMANOVA for testing significance
# With respect to "Behavior"
# Compute distance matrix using vegdist
dist.pcoa <- vegdist(wide, method = 'bray')
head(dist.pcoa)
dist.pcoa
perm <- adonis(dist.pcoa ~ Activity + Diet + Sex + Treatment, data = samples_df, permutations = 999,
               method = 'bray')
perm
#-----------Check assumption of homogeneity of multivariate dispersion-----------
bd <- betadisper(dist.pcoa, samples_df$Activity)
bd <- betadisper(dist.pcoa, samples_df$Diet)
bd <- betadisper(dist.pcoa, samples_df$Sex)
# View differences with boxplot
boxplot(bd)
# boxplot of Average distance to median shows that there might be lower dispersion 
# F-Test
# Testing Homogeneity of dispersion btw each condition. Thus, we do NOT want p<0.05. 
anova(bd)
# permutation test
permutest(bd)
# Remove Ordination & PERMANOVA testing
rm(bd, perm, dist.pcoa)
#############################################################################################
#-------------- DIFFERNTIAL ABUNDANCE (DA) TESTING ------------------------------------------
#############################################################################################
# This is an expanding part of microbiome analysis. If you wish to find the 'heavy hitters'
# or determine through statistical testing specific bacteria significantly different abundant
# between two ecosystems, this is called differential abundance testing

# See Weiss_S_2017_Normalization-microbial-differential-abundances, 
# McMurdie_PJ_Holmes_S_2014_rarefying-microbiome-data-is-inadmissible
# Love_MI_2014_Moderated-estimation-of-fold-change-and-dispersion-for-RNA-seq-data-with-DESeq2
# Phylyoseq recommends DESeq2, but also edgeR, Voom, analysis of composition of microbiomes
# (ANCOM) & metagenomeSeq
# The following two lines actually do all the complicated DESeq2 work. The function 
# phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with 
# dispersions estimated, using the experimental design formula, also shown (the ~Activity term).
# We need to remove baseline samples to DA testing
ps = subset_samples(ps, Activity != "Baseline")
diagdds = phyloseq_to_deseq2(ps, ~ Treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# new DESeq object for factor labeling
dds <- DESeq(diagdds)
# Use factor() command to define levels to see in results() print summary, below
dds$Activity <- factor(dds$Activity, levels = c("Baseline", "Exercise", "Sedentary"))
contrast <- c("Activity", "Exercise", "Sedentary")
# The results function call reates a table of results of the tests, very fast. 
res = results(dds, contrast=contrast, cooksCutoff = FALSE) # Note the "Wald test p-value: Activity Sedentary
res
dim(res)
# vs Exercise. Sedentary is our reference point
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
dim(sigtab)
View(sigtab)
sigtab
# write sigtab to .csv to make figures in Excel. Rename with appropriate factor
write.csv(sigtab, file = "full.data.da.activity.csv")
# Rm files to look at other factors
rm(diagdds, dds, contrast, res, alpha, sigtab)

# DA Testing looking at Diet
diagdds = phyloseq_to_deseq2(ps, ~ Diet)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# new DESeq object for factor labeling
dds <- DESeq(diagdds)
# set up contrasts for hypothesis testing
contrast <- c("Diet", "VeryHighFat", "Control")
# The results function call reates a table of results of the tests, very fast. 
res = results(dds, contrast=contrast, cooksCutoff = FALSE) # Note the "Wald test p-value: Activity Sedentary
res
# vs Exercise. Sedentary is our reference point
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
dim(sigtab)
View(sigtab)
sigtab
# write sigtab to .csv to make figures in Excel. Rename with appropriate factor
write.csv(sigtab, file = "full.data.da.diet.csv")
# Rm files to look at other factors
rm(diagdds, dds, contrast, res, alpha, sigtab)

# DA Testing looking at Sex
diagdds = phyloseq_to_deseq2(ps, ~ Sex)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# new DESeq object for factor labeling
dds <- DESeq(diagdds)
# set up contrasts for hypothesis testing
contrast <- c("Sex", "Female", "Male")
# The results function call reates a table of results of the tests, very fast. 
res = results(dds, contrast=contrast, cooksCutoff = FALSE) # Note the "Wald test p-value: Activity Sedentary
res
# vs Exercise. Sedentary is our reference point
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
dim(sigtab)
View(sigtab)
sigtab
# write sigtab to .csv to make figures in Excel. Rename with appropriate factor
write.csv(sigtab, file = "full.data.da.sex.csv")
# Rm files to look at other factors
rm(diagdds, dds, contrast, res, alpha, sigtab)

# Let's look at the OTUs that were significantly different between the two tissues. 
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Free up space for Exploratory Data Analysis
rm(ps, otu_mat, res, samples, samples_df, tax_mat, tax_mat1, OTU, TAX, long, pivot, wide)

#############################################################################################
#-------------- LEfSe -----------------------------------------------------------------------
#############################################################################################
# These lines of code create relative abundance counts for each OTU by each sample. Create a
# pivot table first and then use 

# Make relative abundance values for otu table (i.e. phyloseq otu_table(ps) = w12.p = week 3 pivot table)
w12r  = transform_sample_counts(ps, function(x) x / sum(x) )
w12fr = filter_taxa(w12r, function(x) mean(x) > 1e-5, TRUE)

dim(otu_table(ps))
dim(otu_table(w12r))
# save 
w12.ra <- otu_table(w12r)
#another way to export .csv filesas csv to local directory
write.csv(w12.ra, file = "w12.ra.csv")

#############################################################################################
#------------- if Phyloseq doesn't work, use vegan for beta diversity -----------------------
#############################################################################################
#-------------- make GS2.wide data as a matrix with numeric values for lm() ----------------
# this of this as, 'just the data'
dim(wide)
spe <- as.matrix(wide[,1:4620])
str(spe)
rowSums(spe)

#############################################################################################
#-------------Exploratory Data Analysis------------------------------------------------------
#############################################################################################

#--------------checking for normality qqnorm-------------------------------------------------
qqnorm(spe)

wide[1:5, 1:10] # Display only 5 lines and 10 columns
head(spe) # Display only the first 6 lines
tail(spe) # Display only the last 6 rows
nrow(spe) # Number of rows (samples)
ncol(spe) # Number of columns (species/OTU/BBH)
colnames(spe) # Column labels (descriptors = species)
rownames(wide) # Row labels (objects = samples)

# Search for different BBHs (i.e. ^ means beginning with) -----------------------------------

test <- pivot %>% filter(str_detect(Species, "^Muribacu"))
View(test)
colSums(test[,2:57])

#--------------------------------------------------------------------------------------------
# Count the cases for each abundance class
(ab <- table(unlist(wide)))
# Barplot of the distribution, all species confounded 
?barplot()
barplot(ab,
        las = 1,
        xlab = "Abundance class", ylab = "Frequency",
        col = gray(5 : 0 / 5))
# How many zeros?
sum(spe == 0)
# Proportion of zeros in the spe data set. Staticians would call this 'zero-inflated'
sum(spe == 0) / (nrow(spe) * ncol(spe))

#################################################################################################
## Compare species: number of occurrences. Compute the number of sites where each species is 
# present. To sum by columns, the second argument of apply(), MARGIN, is set to 2.
pres <- apply(wide > 0, 2, sum)
# Sort the results in increasing order
sort(pres)
# Compute percentage frequencies
relf <- 100 * pres/nrow(wide)
# Round the sorted output to 1 digit
round(sort(relf), 1)
# Plot the histograms
par(mfrow = c(1,2))
hist(pres, main = "Species Occurrences", right = FALSE, las = 1, 
     xlab = "Number of occurrences", ylab = "Number of species",
     col = "bisque")
hist(relf, main = "Species Relative Frequencies", right = FALSE, las = 1, 
     xlab = "Frequency of occurrences (%)", ylab = "Number of species", 
     breaks = seq(0, 100, by = 10), col = "bisque")
# Remove to save space
rm(ab, pres, relf, spe, ab)

########################################################################################
#------------- Alpha Diversity Measures ------------------------------------------------
########################################################################################
# Get help on the diversity() function
?diversity
# Compute alpha diversity indices of the fish communities
N0 <- rowSums(wide > 0) # species richness
N0 <- specnumber(wide) #species richness (alternate)
H <- diversity(wide) # Shannon entropy (base e)
Hb2 <- diversity(wide, base = 2) # Shannon entropy (base 2)
N1 <- exp(H) # Shannon diversity (base e)
# (number of abundant species)

N1b2 <- 2^Hb2 # Shannon diversity (base2)
N2 <- diversity(wide, "inv") # Simpson diversity
# (number of dominant species)

J <- H / log(N0)  # Pielou evenness
E10 <- N1 / N0 # Shannon evenness (Hill's ratio)
E20 <- N2 / N0 # Simpson evenness (Hill's ratio)
alpha <- (div <- data.frame(N0, H, Hb2, N1, N1b2, N2, E10, E20, J))

#############################################################################################
#-------------- NMDS FOR ORDINATION ---------------------------------------------------------
#############################################################################################

#-------------STRESS VALUE FOR PLOT------------------------------------------------
meta.nmds <- metaMDS(wide) #no transformation of species data is made here prior to bray curtis dissimilarities being calculated. (Bray Curtis is the default in R).
str(meta.nmds) # gives stress value for plot
stressplot(meta.nmds) # To gain the stress plot for stress values for your MDS

#-------------dim() for dimensions-------------------------------------------------
dim(wide)
dim(env)

#-----------Distance matirx--------------------------------------------------------
# Calculate the range 0.25 route range()
range(wide)
range(wide^0.25)

# Compute distance matrix using Bray-Curtis on 0.25 route
?vegdist()
bray <- vegdist(wide^0.25, method = 'bray')
kul <- vegdist(wide, method = 'kulczynski')

head(bray)
bray

#-----------NMDS------------------------------------------------------------------
# Run NMDS
?metaMDS
nmds <- metaMDS(bray)
nmds <- metaMDS(kulczynski)
plot(nmds, type = 'text')
nmds
# You can view the coordinates of each sample with
points <- as.data.frame(nmds[["points"]])
nmds[["nobj"]]
#-----------NMDS PLOT WITH GGPLOT2-----------------------------------------------
qplot(MDS1, MDS2, data = points)

# make a dataframe from the env file
# Subset the env file, and use bind_cols to make df for ggplot
samples_df <- as.data.frame(env) %>% 
  select(OG_Sample, Sample, Activity, Diet, Sex, Time, Treatment)
row.names(samples_df) <- samples_df$Sample
samples_df <- samples_df %>% select (-Sample)
?bind_cols()
data <- bind_cols(samples_df, points)
# use ggplot() to create more complex graphics
nmds.plot <- ggplot(data, aes(x=MDS1, y=MDS2, col=Diet)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "NMDS, Kulczynski, Raw Data")
nmds.plot
# PERMANOVA for testing significance
?adonis
perm <- adonis(bray ~ Activity + Diet + Sex + Treatment + Time, data = data, permutations = 999,
               method = 'kulczynski')
perm

#-----------Check assumption of homogeneity of multivariate dispersion-----------
bd <- betadisper(bray, data$Activity)
bd <- betadisper(bray, data$Diet)
bd <- betadisper(bray, data$Sex)
bd <- betadisper(bray, data$Treatment)
bd <- betadisper(bray, data$Time)
# View differences with boxplot. Look at distance from centroid. Any outiers? 
boxplot(bd)
# Use F-Test Testing Homogeneity of dispersion btw each condition. 
# Thus, we do NOT want p<0.05. 
anova(bd)
# permutation test, # Thus, we do NOT want p<0.05
permutest(bd)
#we cannot find a statistically different dispersion, therefore we pass

# REMOVE to save RAM
rm(bd, data, nmds, nmds.plot, perm, points, samples_df, bray, kulczynski, meta.nmds)

#############################################################################################
#-------------- PCOA FOR ORDINATION ---------------------------------------------------------
#############################################################################################
#-------------dim() for dimensions-----------------------------------------------------------
dim(wide)
dim(env)

#-----------Distance matrix------------------------------------------------------------------
# Calculate the range 0.25 route range()
range(wide)
range(wide^0.25)

# Compute distance matrix using Bray-Curtis on 0.25 route
?vegdist()
bray <- vegdist(wide, method = 'bray')
kul <- vegdist(wide, method = 'kulczynski')
#-----------PCoA----------------------------------------------------------------------------
?pcoa
pcoa <- pcoa(bray, correction="none", rn=TRUE)

# S3 method for pcoa
biplot(pcoa, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1, rn=NULL, main="PCoA Bray Distance")
# You can view the coordinates of each sample with
vectors <- as.data.frame(pcoa[["vectors"]])
points <- vectors[,1:2]
#-----------PCoA PLOT WITH GGPLOT2-----------------------------------------------
qplot(Axis.1, Axis.2, data = points)

# make a dataframe from the env file
# Subset the env file, and use bind_cols to make df for ggplot
samples_df <- as.data.frame(env) %>% 
  select(OG_Sample, Sample, Activity, Diet, Sex, Time, Treatment)
row.names(samples_df) <- samples_df$Sample
samples_df <- samples_df %>% select (-Sample)
?bind_cols()
data <- bind_cols(samples_df, points)
# use ggplot() to create more complex graphics
pcoa.plot <- ggplot(data, aes(x=Axis.1, y=Axis.2, col=Diet)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "PCoA, Bray, Raw Data")
pcoa.plot
# PERMANOVA for testing significance
?adonis
perm <- adonis(bray ~ Activity + Diet + Sex + Treatment + Time, data = data, permutations = 999,
               method = 'bray')
perm

#-----------Check assumption of homogeneity of multivariate dispersion-----------
bd <- betadisper(kulczynski, data$Activity)
bd <- betadisper(kulczynski, data$Diet)
bd <- betadisper(kulczynski, data$Sex)
bd <- betadisper(kulczynski, data$Treatment)
bd <- betadisper(kulczynski, data$Time)
# View differences with boxplot. Look at distance from centroid. Any outiers? 
boxplot(bd)
# Use F-Test Testing Homogeneity of dispersion btw each condition. 
# Thus, we do NOT want p<0.05. 
anova(bd)
# permutation test, # Thus, we do NOT want p<0.05
permutest(bd)
#we cannot find a statistically different dispersion, therefore we pass

# REMOVE PCoA data frames and matrixes to save space
rm(pcoa, pcoa.plot, vectors, bray)

#----------- Remove files for next week ------------------------------------------
rm(long, pivot, wide)


