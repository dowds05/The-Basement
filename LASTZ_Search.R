#############################################################################################
#-------------- Searching for BBHs for MUSCLE Alignment -------------------------------------
#############################################################################################
# Following Differential Abundance testing (DESeq2, edgeR, LEfSe, etc.) use filter() to 
# Search for different BBHs (i.e. ^ means beginning with) -----------------------------------
search <- w3 %>% filter(str_detect(Species, "^Clostridium disporicum"))
View(search)
# Look how many unique Species/Strains of BBH
levels(as.factor(search$Strain))
#############################################################################################
#-------------- PIVOT TABLES TO COUNT READS BY SAMPLE ---------------------------------------
#############################################################################################
# PIVOT TABLE on BBH
search.pivot <- search %>% 
  group_by(Species, Strain, Sample, Query, Percent, Length) %>% 
  tally() %>% 
  pivot_wider(names_from = "Sample", values_from = n, names_sort=TRUE) %>% 
  mutate_all(~replace(., is.na(.), 0))
dim(search.pivot)
# Total number of reads per sample
colSums(search.pivot[,5:40])
# Mean, Median, Range of 57 samples
summary(colSums(search.pivot[,5:40]))
summary(search.pivot$Percent)
summary(search.pivot$Length)
# Look how many unique Strains of BBH
levels(as.factor(search.pivot$Strain))

#-------------- EXPORT AS CSV -----------------------------------------------------------------
write.csv(search.pivot, file = "C_disporicum.csv")


# Select query from each sample to 'search by index in Geneious
# Filter week by above search
bug.list <- search%>% 
  filter(Species %in% search3$Species) %>%
  distinct(.keep_all = FALSE)
# Save query as list and save as file. Then, use grep on command line to search "l" or "GS2" data file
query <- bug.list$Query
query

?pull

?heatmap
c.dis <- read_csv("~/Desktop/c_dis.csv",col_names = FALSE)
c.dis <- as.numeric(c.dis)
heatmap(c.dis)
search %>% group_by(Sample, Query) %>% pull(Sample) %>% pull(var = "1")
