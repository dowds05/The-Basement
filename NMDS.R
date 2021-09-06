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

# to bypass install.packages() install with BiocManager::install("package_name")
BiocManager::install("DESeq2")
##################################################################################################################
#-------------- Importing Datasets -------------------------------------------------------------------------------
##################################################################################################################
# Make Treatment column. Note with Sex
l <- l %>% unite(Treatment, Activity, Diet, Time, sep = "_", remove = FALSE)
env <- env %>% unite(Treatment, Activity, Diet, Time, sep = "_", remove = FALSE)

baseline <- l %>% filter(Time == "D_0")
w3 <- l %>% filter(Time == "Week_3")
w4 <- l %>% filter(Time == "Week_4")
w8 <- l %>% filter(Time == "Week_8")
w12 <- l %>% filter(Time == "Week_12")
# UPDATE FACTOR!
merge <- bind_rows(baseline, w3)
# verify levels of merge: Time, Sample, etc.
levels(as.factor(merge$Sample))
##################################################################################################################
#-------------- PIVOT TABLES TO COUNT READS BY SAMPLE ------------------------------------------------------------
##################################################################################################################
pivot <-  l %>% group_by(Strain, Sample) %>% tally() %>% 
  pivot_wider(names_from = "Sample", values_from = n, names_sort=TRUE) %>% mutate_all(~replace(., is.na(.), 0))
dim(pivot)
# Total number of reads per sample
colSums(pivot[,2:227])
# Mean, Median, Range of 57 samples
summary(colSums(pivot[,2:59]))
#-------------------- From long-to-wide format -------------------------------------------------------------------
# wide format 
wide <- pivot %>% t() %>% .[-1,]
class(wide) <- "numeric"
# Make a data frame
wide <- as.data.frame(wide)
# check row sums, confirm rowSums=colSums
rowSums(wide)
# Change column names to Description (or 'species/OTU/BBH')
colnames(wide) <- pivot$Strain
##################################################################################################################
#-------------- NMDS FOR ORDINATION ------------------------------------------------------------------------------
##################################################################################################################
#-------------STRESS VALUE FOR PLOT-------------------------------------------------------------------------------
# No transformation of species data is made here prior to bray curtis dissimilarities being 
# calculated. (Bray Curtis is the default in R).
meta.nmds <- metaMDS(wide)
str(meta.nmds) # gives stress value for plot
stressplot(meta.nmds) # To gain the stress plot for stress values for your MDS
#-------------dim() for dimensions--------------------------------------------------------------------------------
dim(wide)
dim(env)
#-----------Distance matirx---------------------------------------------------------------------------------------
# Calculate the range 0.25 route range()
range(wide)
range(wide^0.25)
# ?decostand
t <- decostand(wide, "log")
rowSums(t)
rowSums(wide)
range(t)
# Compute distance matrix using Bray-Curtis on 0.25 route. ?vegdist() --------------------------------------------
bray <- vegdist(t, method = 'bray')
#-----------NMDS--------------------------------------------------------------------------------------------------
# Run NMDS. ?metaMDS
nmds <- metaMDS(bray)
plot(nmds, type = 'text')
nmds
# You can view the coordinates of each sample with
points <- as.data.frame(nmds[["points"]])
nmds[["nobj"]]

# UPDATE WEEK! Make target for subsetting samples_df, 
target <- c("D_0", "Week_3")
# Subset of the GS2.env file, subset Metadata file for conformable arrays
df <- filter(env, Time %in% target)
# make a dataframe from the env file with bind_cols() samples_df & points. ?bind_cols()
data <- bind_cols(env, points)
# Make Treatment column for activity, diet, sex, time ?unite()
data <- data %>% unite(Treatment, Activity, Diet, Sex, Time, sep = "_", remove = FALSE)

#-------------- EXPORT AS CSV --------------------------------------------------------------------------------------
write.csv(data, file = "nmds_tot.csv") # Use name from search

# UPDATE WEEK! use ggplot() to create more complex graphics. ?ggplot
treatment <- ggplot(data, aes(x=MDS1, y=MDS2, col=Treatment, shape=Treatment)) +
  geom_point(size = 5) +
  scale_shape_manual(values=c(17, 1, 0, 16, 15)) +
  scale_color_manual(values=c('red1','black', 'deeppink1', 'purple1', 'darkorange1')) +
  #scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(cols = vars(Sex)) +
  theme(legend.position="bottom") +
  labs(title = "Week_8, Log Transformation, Bray")
treatment # screen capture this plot/export for publication or presentation
# UPDATE WEEK!
activity <- ggplot(data, aes(x=MDS1, y=MDS2, col=Activity, shape=Activity)) +
  geom_point(size = 5) +
  scale_shape_manual(values=c(17, 15, 16)) +
  scale_color_manual(values=c('red1', 'purple2', 'darkgoldenrod2')) +
  #scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(cols = vars(Sex)) +
  theme(legend.position="bottom") +
  labs(title = "Week_8, Log Transformation, Bray")
activity # screen capture this plot/export for publication or presentation
# UPDATE WEEK!
diet <- ggplot(data, aes(x=MDS1, y=MDS2, col=Diet, shape=Diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values=c(17, 15, 16)) +
  scale_color_manual(values=c('red1', 'royalblue1', 'darkseagreen3')) +
  #scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(cols = vars(Sex)) +
  theme(legend.position="bottom") +
  labs(title = "Week_8, Log Transformation, Bray")
diet # screen capture this plot/export for publication or presentation

# PERMANOVA for testing significance. ?adonis --------------------------------------------------------------------
colnames(env)
head(bray)
bray
# PERMANOVA using ?adonis()
perm <- adonis2(bray ~ Activity + Diet + Sex + Treatment + Time,
               data = data, permutations = 999, method = 'bray')
perm
#-----------Check assumption of homogeneity of multivariate dispersion. ?betadisper ------------------------------
bd.activity <- betadisper(bray, data$Activity)
bd.diet <- betadisper(bray, data$Diet)
bd.sex <- betadisper(bray, data$Sex)
bd.time <- betadisper(bray, data$Time)
bd.treat <- betadisper(bray, data$Treatment)
# View differences with boxplot
boxplot(bd.activity)
boxplot(bd.diet)
boxplot(bd.sex)
boxplot(bd.time)
boxplot(bd.treat)
# boxplot of Average distance to median shows that there might be lower dispersion 
# F-Test Testing Homogeneity of dispersion btw each condition. Thus, we do NOT want p<0.05. 
aov.activity <- anova(bd.activity)
aov.diet <- anova(bd.diet)
aov.sex <- anova(bd.sex)
aov.time <- anova(bd.time)
aov.treat <- anova(bd.treat)
# permutation test
permutest(bd.activity)
permutest(bd.diet)
permutest(bd.sex)
permutest(bd.time)
permutest(bd.treat)
# UPDATE WEEK! bind_rows() perm & assumptions
perm_w8_assumptions <- bind_rows(perm$aov.tab, aov.activity, aov.diet, aov.sex, aov.time, aov.treat)

# UPDATE WEEK! -- EXPORT AS CSV ----------------------------------------------------------------------------------
write.csv(perm_w8_assumptions, file = "perm_w8_assumptions.csv")

#  UPDATE WEEK! REMOVE Ordination & PERMANOVA testing
rm(merge, pivot, wide, t, meta.nmds, bray, nmds, treatment, activity, diet, df, points, target, data, perm, 
   bd.activity, bd.diet, bd.sex, bd.treat, bd.time, aov.activity, aov.diet, aov.sex, aov.time, aov.treat, 
   perm_w8_assumptions)



