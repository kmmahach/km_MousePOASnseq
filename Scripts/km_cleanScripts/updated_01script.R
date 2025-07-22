# R-4.3.1
# Cell counts output from CellRanger v(?), details found in >txt file 
# Souporcell output from https://github.com/wheaton5/souporcell, 
# details found in >txt file 

root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq/"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))


# load functions
source("./functions/QC_filtering_fun.R")

# libraries
lapply(c("tidyverse","Seurat", "SeuratExtend", "collapse", "clustree"), 
       library, character.only = T)

#### Read in data ####
## count data
# path <- list.dirs(path = paste0(root.dir, "count"), 
#                   recursive = TRUE, full.names = TRUE)
# l.df <- lapply(ls(pattern="*.data"), \(x) get(x)) # could use for indiv dfs

# .h5 files on github, or
# if accessing netdir, then 
# filtered_feature_bc_matrix directories contain the same data

temp <- list.files(root.dir, pattern="*filtered_feature_bc_matrix.h5", 
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
path <- list.files(path = paste0(net.dir, "souporcell"), pattern = "*clusters.tsv",
                  recursive = TRUE, full.names = TRUE)

l.clust <- lapply(path, \(x) { read_tsv(x) })
names(l.clust) = c("female.dom.clusters",
                   "male.dom.clusters",
                   "female.sub.clusters",
                   "male.sub.clusters")

l.clust <- l.clust[order(names(l.clust))] %>% 
  lapply(\(x) {
    x = x %>%
      rename(cell.id = barcode,
             doublet = status,
             indiv_genotype = assignment,
             indiv1 = cluster0,
             indiv2 = cluster1,
             indiv3 = cluster2)
    return(x)
  })


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
gene.min.cells.plot(l.df2, "") # autosaves .png in ./QC_filtering if no outdir

## Seurat QC steps
l.dfs <- make.seurat.obj(l.df) # should work on both l.df and l.df2

# add metadata
l.dfs <- mapply(AddMetaData, l.dfs, l.clust)
  # rm(l.df, l.df2, l.clust) # save space in env

# label mitochondrial genes
l.dfs <- lapply(l.dfs, \(x) {
    PercentageFeatureSet(x, pattern = "^mt-",
                         col.name = "percent.mt") } )

# graph QC metrics
# violin plots 
plot.qc.metrics(l.dfs, "")
plot.doublets.qc(l.dfs, "")
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
plot.feature.scatter(l.dfs, filters, "")
  make.filtered.seurat(l.dfs, filters) -> fltr.ldfs

# find best cluster resolution?
find.cluster.range(fltr.ldfs, "")


                 