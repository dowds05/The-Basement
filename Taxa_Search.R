##################################################################################################################
#---------------------------- Searching for BBHs for MUSCLE Alignment --------------------------------------------
##################################################################################################################
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

# Search for different Species. ?filter
species <- l %>% filter(str_detect(Species, "^Muri"))
dim(species)
# Look how many unique Species/Strains of BBH
levels(as.factor(species$Strain))
levels(as.factor(species$Activity))
levels(as.factor(species$Diet))
levels(as.factor(species$Sex))
levels(as.factor(species$Sample))
View(species)
# Search for different Strain
strain <- l %>% filter(str_detect(Strain, "^Muribaculum_intestinale_strain_YL27"))
dim(strain)
# Look how many unique Species/Strains of BBH
levels(as.factor(strain$Strain))
levels(as.factor(strain$Activity))
levels(as.factor(strain$Diet))
levels(as.factor(strain$Sex))
levels(as.factor(strain$Sample))

# Search for different Family
family <- l %>% filter(str_detect(Family, "^De"))
dim(family)
# Look how many unique Species/Strains of BBH
levels(as.factor(family$Species))
levels(as.factor(family$Activity))
levels(as.factor(family$Diet))
levels(as.factor(family$Sex))
levels(as.factor(family$Sample))

# Search for different Genus
genus <- l %>% filter(str_detect(Genus, "^Muri"))
dim(genus)
# Look how many unique Species/Strains of BBH
levels(as.factor(genus$Strain))
levels(as.factor(genus$Activity))
levels(as.factor(genus$Diet))
levels(as.factor(genus$Sex))
levels(as.factor(genus$Sample))

# Search for different Phylum
phylum <- l %>% filter(str_detect(Phylum, "^Bacteroidetes"))
dim(phylum)
# Look how many unique Species/Strains of BBH
levels(as.factor(phylum$Strain))
levels(as.factor(phylum$Activity))
levels(as.factor(phylum$Diet))
levels(as.factor(phylum$Sex))
levels(as.factor(phylum$Sample))
#-------------- SCATTER PLOT -------------------------------------------------------------------------------------
# basic scatterplot
sp <- ggplot(species, aes(x=Length, y=Percent)) + geom_point() + 
  labs(title = "Angelakisella_massiliensis_strain_Marseille-P3217 | 340404 reads")
sp

# limit percent?
p <- strain %>% filter(Percent >= 85.00) %>% filter(Length >= 3500)
dim(p)

sp2 <- ggplot(p, aes(x=Length, y=Percent)) + geom_point() + 
  labs(title = "Duncaniella_sp._B8 | 409 filtered reads")
sp2

##################################################################################################################
#------------------------------- PIVOT TABLES TO COUNT READS BY SAMPLE -------------------------------------------
##################################################################################################################
# PIVOT TABLE on BBH
search.pivot <- p %>% group_by(Treatment, Sex, Time) %>% tally() %>% 
  pivot_wider(names_from = "Treatment", values_from = n, names_sort=TRUE) %>% mutate_all(~replace(., is.na(.), 0))
dim(search.pivot)
view(search.pivot)
# Total number of reads per sample
colSums(search.pivot[,4:48])
hist(colSums(search.pivot[,4:8]))
# Mean, Median, Range of 57 samples. Look to see sample read distribution
summary(colSums(search.pivot[,4:8]))
summary(q1$Percent)
summary(q1$Length)

#-------------- EXPORT AS CSV ------------------------------------------------------------------------------------
write.csv(strain, file = "D_sp_B8.csv") # Use name from search
write.csv(species, file = "muri.species.total.csv")
#-------------- EXPORT Query by INDIVIDUAL TREATMENT to search.txt for filtering on Amarel -----------------------
# UPDATE TREATMENT, TIME, SEX
levels(as.factor(species$Time))
levels(as.factor(species$Treatment))
q1 <- species %>% filter(str_detect(Treatment, "^Exercise_Control")) %>% 
  filter(str_detect(Time, "^Week_8")) %>% filter(str_detect(Sex, "^Male"))
dim(q1)
query <- q1 %>% select(Query)
write.table(query, file = "/Users/lab/Desktop/Bob_Aim3/Linux/DATA/query.fa", sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)

#-------------- EXPORT Query by TAXA on Amarel -------------------------------------------------------------------
# UPDATE TREATMENT, TIME, SEX
levels(as.factor(p$Strain))
levels(as.factor(p$Treatment))
query <- p %>% select(Query)
dim(query)
write.table(query, file = "/Users/lab/Desktop/Bob_Aim3/Linux/DATA/FASTA/taxa.txt", sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)

#-------------- DELETE IF USING ANOTHER R SCRIPT FILE ------------------------------------------------------------
rm(species, strain, family, genus, phylum, search.pivot, sp, sp2, q1, query, p)
