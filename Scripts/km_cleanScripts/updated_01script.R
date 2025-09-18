# R-4.3.1
# Cell counts output from CellRanger-7.0.1, details found in Scripts/cellranger_poa_IMC.txt
# Souporcell output from https://github.com/wheaton5/souporcell, details found in Scripts/souporcell_poa_IMC.txt 

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# Seurat-4.4.0; results may vary slightly from other versions of v4/v5
# remove.packages("Seurat")
# remotes::install_version("Seurat", version = "4.4.0")
library(Seurat)
# check
packageVersion("Seurat")
# [1] ‘4.4.0’

# requires SeuratObject install (is compatible with all Seurat versions)
packageVersion("SeuratObject") 
# [1] ‘5.1.0’

# load functions
source("./functions/QC_filtering_fun.R")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "patchwork", "SeuratExtend", "scCustomize", 
                "collapse", "clustree", "ggrepel", "HGNChelper","openxlsx"),
              out_prefix = "01")


#### Read in data ####
## count data
temp <- list.files(root.dir, "*filtered_feature_bc_matrix.h5", 
                   recursive = TRUE, full.names = TRUE)

# automated naming:
lapply(temp, \(x) { Read10X_h5(x) }) %>% 
  set_names(
    str_extract(temp, "(?<=/)[^/]+(?=/outs/)") %>% 
    str_replace_all("([a-z])([A-Z])", "\\1.\\2") %>%
    str_to_lower() %>%
    paste0(".data")) -> l.df

# l.df <- l.df[order(names(l.df))] # already in abc order

  save(l.df, file = "./data/rawdata.rda",
       compress = compression) # too big for github

# souporcell output
path <- list.files(root.dir, "*clusters.tsv",
                   recursive = TRUE, full.names = TRUE)

# automated naming: 
lapply(path, \(x) { as.data.frame(read_tsv(x)) }) %>%
  set_names(
    str_extract(path, "(?<=/)[^/]+(?=/souporcell/)") %>% 
    str_replace_all("([a-z])([A-Z])", "\\1.\\2") %>%
    str_to_lower() %>%
    paste0(".clusters")) -> l.clust

# organize
l.clust[order(names(l.clust))] %>% 
  lapply(\(x) {
    x = x %>%
      rename(cell.id = barcode,
             doublet = status,
             indiv_genotype = assignment) %>% 
      select(cell.id, doublet, indiv_genotype) 
  }
) -> l.clust

  save(l.clust, file = "./data/souporcell_output.rda")


#### Inital QC ####
## Check range of genes by min cell count that expr
steps = c(3,
          seq(from = 10, 
              to = 100, 
              by = 10),
          seq(from = 200, 
              to = 1000, 
              by = 100))

l.df2 <- gene.min.cells.range(l.df, steps) # may take a while

  save(l.df2, file = "./data/rawdata_withGeneByCell_counts.rda",
       compress = compression) # still too big for github

## Generate graphs
qc_plots_path <- paste0(root.dir, "/QC_filtering")

  gene.min.cells.plot(l.df2, qc_plots_path) 

## Seurat QC steps
l.dfs <- make.seurat.obj(l.df2) # should work on both l.df and l.df2

# add metadata
map2(l.dfs, l.clust, ~ {
  # df <- column_to_rownames(.y, "cell.id")
  # df = .y
  AddMetaData(.x, metadata = .y)
  }
) -> l.dfs

  rm(l.df, l.df2, l.clust) # save space in env

# label mitochondrial genes
lapply(l.dfs, \(x) {
  PercentageFeatureSet(x, pattern = "^mt-",
                       col.name = "percent.mt") 
  } 
) -> l.dfs


## Graph QC metrics
# violin plots 
  plot.qc.metrics(l.dfs, qc_plots_path)
  plot.doublets.qc(l.dfs, qc_plots_path)

# choose filter cutoffs 
filter.categories <- c("nFeature_min", "nFeature_max", "mito.max")
  fdf <- c(700, 1400, 5)
  fsf <- c(700, 1500, 5)
  mdf <- c(700, 3000, 5)
  msf <- c(700, 2500, 5)

filters <- matrix(rbind(fdf,fsf,mdf,msf),
                  nrow = length(l.dfs),
                  ncol = length(filter.categories),
                  dimnames = list(names(l.dfs),
                                  filter.categories))


# visualize
  plot.feature.scatter(l.dfs, filters, qc_plots_path)
  check.filters(l.dfs, filters, qc_plots_path)
  plot.filters(l.dfs, filters, qc_plots_path)
  
fltr.ldfs <- make.filtered.seurat(l.dfs, filters)

# dimensionality reduction
l.dfs <- lapply(fltr.ldfs, \(x, dims = 1:15) {
  x = x %>% 
    SCTransform() %>% 
    RunPCA() %>%
    FindNeighbors(dims = dims) %>%
    RunUMAP(dims = dims)
  
  return(x)
  }
)

# plot PCA, UMAP, var features, and cluster
  plot.dim.clust(l.dfs, qc_plots_path)

# find best cluster resolution?
  find.cluster.range(l.dfs, qc_plots_path) # use auto resolution ranges

# may not replicate exact number of clusters depending on which Seurat version 
l.dfs <- lapply(l.dfs, \(x) {
  x = FindClusters(x, resolution = 1)
  return(x)
  }
)

  save(l.dfs, file = "./data/filtered_counts.rda",
       compress = compression) # still too big for github

## cell type annotation with ScType

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


# add cell types to Seurat object clusters
l.dfs <- annotate.with.sctype(l.dfs, "sctype.ind", qc_plots_path)

  save(l.dfs, file = "./data/filtered_counts_withScType.rda",
       compress = compression) # still too big for github
