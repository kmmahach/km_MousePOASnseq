# R-4.3.1
# Cell counts output from CellRanger v(?), details found in >txt file 
# Souporcell output from https://github.com/wheaton5/souporcell, 
# details found in >txt file 

root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))


# load functions
source("./functions/QC_filtering_fun.R")

# libraries
lapply(c("tidyverse","Seurat", "SeuratExtend", "collapse", "clustree",
         "scCustomize"), 
       library, character.only = T)

#### Read in data ####
## count data
# path <- list.dirs(path = paste0(root.dir, "count"), 
#                   recursive = TRUE, full.names = TRUE)
# l.df <- lapply(ls(pattern="*.data"), \(x) get(x)) # could use for indiv dfs

# .h5 files on github (easiest), or
# if accessing netdir, then 
# filtered_feature_bc_matrix directories contain the same data

temp <- list.files(root.dir, "*filtered_feature_bc_matrix.h5", 
                   recursive = TRUE, full.names = TRUE)

l.df <- lapply(temp, \(x) { Read10X_h5(x) })
names(l.df) = c("female.dom.data",
                "female.sub.data",
                "male.dom.data",
                "male.sub.data")

# l.df <- l.df[order(names(l.df))] # already in abc order

# save(l.df, file = "./data/rawdata.rda",
#      compress = "xz") # too big for github


# souporcell output
path <- list.files(paste0(net.dir, "/souporcell"), "*clusters.tsv",
                   recursive = TRUE, full.names = TRUE)

# automated naming:
  lapply(path, \(x) { read_tsv(x) }) %>% 
    set_names(
      str_extract(
        path, "[^/]+(?=_soupercell)" ) %>% 
        paste(str_sub(.,4),
              str_extract(., "[a-z]{1,3}"), 
              "clusters",
              sep = ".") %>%
        sub("[a-z]*.", "", .)) -> l.clust # whew that's not very pretty 
                                           # but it works

# base r version is even worse:
  # path |>
  #   lapply(function(x) read.delim(x, header = TRUE, sep = "\t")) |>
  #   (\(lst) {
  #     names(lst) <- paste(
  #       substring(regmatches(path,
  #                            regexpr("[^/]+(?=_soupercell)", 
  #                                    path, perl = TRUE)), 4),
  #       substr(regmatches(path, regexpr("[^/]+(?=_soupercell)", 
  #                                       path, perl = TRUE)), 1, 3),
  #       "clusters", sep = ".")
  #     lst 
  #     })() -> l.clust

# manual naming:
  # (easy to read but also to make mistakes...)

  # names(l.clust) = c("female.dom.clusters",
  #                    "male.dom.clusters",
  #                    "female.sub.clusters",
  #                    "male.sub.clusters")

l.clust[order(names(l.clust))] %>% lapply(\(x) {
    x = x %>%
      rename(cell.id = barcode,
             doublet = status,
             indiv_genotype = assignment,
             indiv1 = cluster0,
             indiv2 = cluster1,
             indiv3 = cluster2) 
    }) -> l.clust

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

l.df2 <- gene.min.cells.range(l.df, steps)
  # save(l.df2, file = "./data/rawdata_withGeneByCell_counts.rda",
  #      compress = "xz") # still too big for github

## Generate graphs
qc_plots_path <- paste0(root.dir, "/QC_filtering")

  gene.min.cells.plot(l.df2, qc_plots_path) 

## Seurat QC steps
l.dfs <- make.seurat.obj(l.df) # should work on both l.df and l.df2

# add metadata
l.dfs <- mapply(AddMetaData, l.dfs, l.clust)
   # rm(l.df, l.df2, l.clust) # save space in env

# label mitochondrial genes
l.dfs <- lapply(l.dfs, \(x) {
  PercentageFeatureSet(x, pattern = "^mt-",
                       col.name = "percent.mt") 
  } 
)

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
fltr.ldfs <- make.filtered.seurat(l.dfs, filters)

# dimensionality reduction
l.dfs <- lapply(fltr.ldfs, \(x, dims = c(1:15)) {
  x = x %>% 
        SCTransform() %>% 
        RunPCA() %>% 
        FindNeighbors(dims = dims) %>% 
        RunUMAP(dims = dims)
  }
)
  
# find best cluster resolution?
find.cluster.range(l.dfs, qc_plots_path)

# plot PCA, UMAP, var features
plot.dim.clust(l.dfs, 1, qc_plots_path)

                 