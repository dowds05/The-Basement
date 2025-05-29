#####################################################################################################################################################
#-------------- LOAD & INSTALL PACKAGES -------------------------------------------------------------------------------------------------------------
#####################################################################################################################################################
sessionInfo() # To check R & MacOS versions
# Old install code  ---------------------------------------------------------------------------------------------------------------------------------
install.packages('ggrepel') # If need base R install code
BiocManager::install('ggpubr') # Install additional packages
BiocManager::valid() # For problems with old packages after R re-install, use this and follow the instructions
# Install   -----------------------------------------------------------------------------------------------------------------------------------------
packages <- c("NBZIMM", "remotes", "R.rsp", "plyr", "tidyverse", "vegan", "cluster", "FD", "lattice", "readxl", "adespatial", "DECIPHER", "janitor", "phyloseq", "ggpubr", "DESeq2", "gplots", "EnhancedVolcano", "apeglm", "microeco", "devtools", "ggrepel")
for(x in packages){
  if(!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
  }
}
# REMOVE Load & Install values ----------------------------------------------------------------------------------------------------------------------
rm(packages, x)
