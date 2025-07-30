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

int.ldfs@project.name = "integrated.data"

  # save(int.ldfs, file = "./data/integrated_seurat.rda",
  #      compress = "xz") # definitely too big to store for free!

# dimensionality reduction
int.ldfs %>%
  RunPCA() %>% 
  RunUMAP(reduction = "pca", 
          dims = 1:30) %>% 
  PrepSCTFindMarkers(assay = "SCT") %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) -> int.ldfs

#### Check Results of Integration ####
# check umap
DimPlot(int.ldfs,
        reduction = "umap",
        group.by = "orig.ident")

ggsave(paste0(qc_plots_path, 
              "/integrated/dimplot.ident.integrated.png"),
       height = 10,
       width = 10)

DimPlot(int.ldfs,
        reduction = "umap",
        label = TRUE,
        repel = TRUE)

ggsave(paste0(qc_plots_path, 
              "/integrated/clusters.dimplot.integrated.png"),
       height = 10,
       width = 10)

# find best cluster resolution?
res.1 = seq(0,1,0.2)
res.2 = seq(0,0.8,0.2) # reduced range

  find.cluster.range(int.ldfs,
                     qc_plots_path,
                     int_range1 = res.1,
                     int_range2 = res.2)
# best res = 0.4 
int.ldfs <- FindClusters(int.ldfs, resolution = 0.4)

# check cells per group per cluster
table(int.ldfs@active.ident, 
      int.ldfs@meta.data$orig.ident) %>%
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>% 
  pivot_longer(cols = unique(int.ldfs$orig.ident),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total) -> cells.per.cluster.table

# graph
cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()

ggsave(paste0(qc_plots_path, 
              "/integrated/cells.per.cluster.integrated.png"),
       width = 5,
       height = 4)

# make figures - presentation
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
# empty cluster
DimPlot(int.ldfs,
        reduction = "umap",
        cols = c(rep('grey',
                     43)),
        pt.size = 1) +
  theme(legend.position = 'none')

ggsave(paste0(qc_plots_path, 
              "/integrated/empty.dimplot.comparison.all.png"),
       width = 10,
       height = 10)

# presentation simple cluster
DimPlot(int.ldfs,
        reduction = "umap",
        label = T,
        pt.size = 1,
        label.size = 10,
        repel = T) +
  theme(legend.position = 'none')

ggsave(paste0(qc_plots_path, 
              "/integrated/clusters.dimplot.comparison.all.png"),
       width = 10,
       height = 10)

# presentation across samples cluster
DimPlot(int.ldfs,
        reduction = "umap",
        split.by = 'orig.ident',
        pt.size = 1) +
  theme(legend.position = 'none')

ggsave(paste0(qc_plots_path, 
              "/integrated/domVsub.clusters.dimplot.comparison.all.png"),
       width = 20,
       height = 10)

# presentation across samples genotype
DimPlot(int.ldfs,
        reduction = "umap",
        split.by = 'orig.ident',
        group.by = 'indiv_genotype',
        pt.size = 1) 

ggsave(paste0(qc_plots_path, 
              "/integrated/domVsub.genotype.dimplot.comparison.all.png"),
       width = 20,
       height = 10)

# presentation across samples nfeature
FeaturePlot(int.ldfs,
            reduction = "umap",
            split.by = 'orig.ident',
            features = 'nFeature_SCT',
            pt.size = 1)

ggsave(paste0(qc_plots_path, 
              "/integrated/domVsub.nfeature.dimplot.comparison.all.png"),
       width = 20,
       height = 10)

#### ScType on Integrated SeuratObj ####

