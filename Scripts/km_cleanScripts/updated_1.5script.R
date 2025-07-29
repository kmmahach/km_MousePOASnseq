# R-4.3.1
# Cell counts output from CellRanger v7.0.1, details found in >txt file 
# Souporcell output from https://github.com/wheaton5/souporcell, 
# details found in >txt file 

root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))

# path to save plots
qc_plots_path <- paste0(root.dir, "/QC_filtering")

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
  # rm(l.dfs)

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
  # save(int.ldfs, file = "./data/integrated_seurat.rda",
  #      compress = "xz") # definitely too big to store for free!

int.ldfs@project.name = "integrated.data"

# dimensionality reduction
int.ldfs %>%
  RunPCA() %>% 
  RunUMAP(reduction = "pca", 
          dims = 1:30) %>% 
  PrepSCTFindMarkers(assay = "SCT") %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) -> int.ldfs

# check umap
DimPlot(int.ldfs,
        reduction = "umap",
        group.by = "orig.ident")

# find best cluster resolution?
res.1 = seq(0,1,0.2)
res.2 = seq(0,0.8,0.2) # reduced range

  find.cluster.range(int.ldfs,
                     qc_plots_path,
                     int_range1 = res.1,
                     int_range2 = res.2)
# best res = 0.4 
int.ldfs <- FindClusters(int.ldfs, resolution = 0.4)
