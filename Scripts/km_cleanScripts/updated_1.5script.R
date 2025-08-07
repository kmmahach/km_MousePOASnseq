# R-4.3.1
# Seurat-4.4.0
# Cell counts output from CellRanger-7.0.1, details found in Scripts/cellranger_poa_IMC.txt
# Souporcell output from https://github.com/wheaton5/souporcell, details found in Scripts/souporcell_poa_IMC.txt 
# ScType wrapper function adapted from https://github.com/IanevskiAleksandr/sc-type 

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# path to save plots
qc_plots_path <- paste0(root.dir, "/QC_filtering")

# load functions
source("./functions/QC_filtering_fun.R")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "Seurat", "SeuratExtend", "collapse", "clustree",
                "scCustomize", "HGNChelper","openxlsx", "ggalluvial"),
              out_prefix = "1.5")

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

  save(int.ldfs, file = "./data/integrated_seurat.rda",
       compress = "xz") # definitely too big to store for free!

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

 find.cluster.range2(int.ldfs,
                     qc_plots_path,
                     int_range1 = res.1,
                     int_range2 = res.2)
# best res = 0.4 
 l.dfs <- lapply(l.dfs, \(x) {
   x = FindClusters(x, resolution = 0.4)
   return(x)
 }
)

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

# load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# gene list
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# set tissue type
tissue = "Brain" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

l.dfs <- annotate.with.sctype(int.ldfs, 
                              "sctype.integrated", 
                              qc_plots_path) %>% 
  set_names("int.ldfs")

### compare individual to integrated
## graph alluvial plot
l.dfs$int.ldfs@meta.data %>%
  select(c(sctype.integrated,
           sctype.ind)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(sctype.integrated, -Freq.scale),
             axis2 = reorder(sctype.ind,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.integrated)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.integrated",
                              "sctype.ind"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()

ggsave(paste0(qc_plots_path, 
              "/integrated/alluvial.cells.by.sctype.indVintegrated.png"),
       width = 15,
       height = 10)

## graph by individual 
for (i in unique(l.dfs$int.ldfs$orig.ident)) {
  l.dfs$int.ldfs@meta.data %>%
    filter(orig.ident == i) %>% 
    select(c(sctype.integrated,
             sctype.ind)) %>%
    table() %>%
    as.data.frame() %>%
    mutate(Count = sum(Freq)) %>%
    ungroup() %>%
    mutate(Freq.scale = 100*Freq/Count) %>%
    ggplot(aes(axis1 = reorder(sctype.integrated, -Freq.scale),
               axis2 = reorder(sctype.ind, -Freq.scale),
               y = Freq.scale)) +
    geom_alluvium(aes(fill = sctype.integrated)) +
    geom_stratum() +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("sctype.all",
                                "sctype.ind"),
                     expand = c(0.15, 0.05)) +
    scale_fill_viridis_d() +
    theme_classic() +
    ggtitle(paste(i))
  
  
  ggsave(paste0(qc_plots_path, 
                "/integrated/alluvial.cells.by.sctype.indVintegrated.",
               i, ".png"),
         width = 15,
         height = 10)
}

## graph alluvial plot with clusters
l.dfs$int.ldfs@meta.data %>% 
  select(c(integrated_snn_res.0.4,
           sctype.integrated,
           sctype.ind)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(sctype.integrated, -Freq.scale),
             axis2 = reorder(integrated_snn_res.0.4, -Freq.scale),
             axis3 = reorder(sctype.ind,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.integrated)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all",
                              "sctype.ind"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()

ggsave(paste0(qc_plots_path, 
              "/integrated/alluvial.cells.by.sctype.indVintegrated.clusters.png"),
       width = 22,
       height = 12)

## counts per sample
l.dfs$int.ldfs@meta.data %>% 
  select(c(orig.ident,
           sctype.integrated)) %>%
  table() %>%
  as.data.frame() %>% 
  ggplot(aes(x = reorder(sctype.integrated, -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  xlab('')+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 4))

ggsave(paste0(qc_plots_path, 
              "/integrated/cell.count.per.sample.sctype.integrated.png"),
       height = 5,
       width = 5)

#### Integrated ScType vs. HypoMap [stopped here] ####
# line 1952 in seurat_snseq_mouse_IMC.R
# save integrated seurat object with sctype annotation

int.ldfs <- l.dfs$int.ldfs

  save(int.ldfs, file = "./data/integrated_seurat_withScType.rda",
       compress = "xz")
  
# run HypoMap and continue in updated_1.75script.R 
  # mostly exploratory comparisons between celltype annotation methods