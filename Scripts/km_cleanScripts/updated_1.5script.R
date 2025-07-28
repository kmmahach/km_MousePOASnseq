# R-4.3.1
# Cell counts output from CellRanger v7.0.1, details found in >txt file 
# Souporcell output from https://github.com/wheaton5/souporcell, 
# details found in >txt file 

root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))


# load functions
source("./functions/QC_filtering_fun.R")

# libraries
lapply(c("tidyverse","Seurat", "SeuratExtend", "collapse", "clustree",
         "scCustomize", "HGNChelper","openxlsx"), 
       library, character.only = T)


#### Seurat Object Integration ####
load("./data/filtered_counts_withScType.rda")

# feature selection
features <- SelectIntegrationFeatures(object.list = l.dfs,
                                      nfeatures = 3000)
# identify anchors
prep.ldfs <- PrepSCTIntegration(object.list = l.dfs,
                                anchor.features = features)

# find anchors - may be time-consuming
anchors <- FindIntegrationAnchors(object.list = prep.ldfs,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction = "rpca", # runs faster integration
                                  k.anchor = 20) # increase strength of alignment



# integrate - retain all genes
int.ldfs <- IntegrateData(anchorset = anchors,
                          normalization.method = "SCT",
                          dims = 1:30)