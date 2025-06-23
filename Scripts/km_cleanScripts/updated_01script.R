# R-4.3.1
# Cell counts output from CellRanger v(?), details found in >txt file 
# Souporcell output from https://github.com/wheaton5/souporcell, 
# details found in >txt file 

root.dir <- "/stor/home/kmm7552/KMmPOA_snRNAseq/"
setwd(paste0(root.dir, "/Scripts/KM_cleanScripts"))


# load functions
source("./functions/QC_preprocessing_fun.R")

# libraries
lapply(c("tidyverse","Seurat", "collapse"), 
       library, character.only = T)

#### Read in data ####
# count data
path <- list.dirs(path = paste0(root.dir, "count"), 
                  recursive = TRUE, full.names = TRUE)

temp <- list.files(path, pattern="*_matrix", full.names = TRUE)
l.df <- lapply(temp, \(x) { Read10X(x) })
names(l.df) = c("female.dom.data",
                "female.sub.data",
                "male.dom.data",
                "male.sub.data")

# l.df <- l.df[order(names(l.df))] # already in abc order
save(l.df, file = "./data/rawdata.rda",
     compress = "xz") # too big for github

# l.df <- lapply(ls(pattern="*.data"), \(x) get(x)) # could use for indiv dfs

# souporcell output
path <- list.files(path = paste0(root.dir, "souporcell"), pattern = "*.tsv",
                  recursive = TRUE, full.names = TRUE)

l.clust <- lapply(path, \(x) { read_tsv(x) })
names(l.clust) = c("female.dom.clusters",
                   "male.dom.clusters",
                   "female.sub.clusters",
                   "male.sub.clusters")

l.clust <- l.clust[order(names(l.clust))] %>% 
  lapply(\(x) {
    x = x %>%
      rename(Cell.id = barcode,
             Doublet = status,
             Genotype = assignment)
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
  save(l.df2, file = "./data/rawdata_withGeneByCell_counts.rda",
       compress = "xz") # still too big for github

## Generate graphs
  # could probably do this in a prettier way but it works
mapply("c",l.df2,names(l.df2),SIMPLIFY=FALSE) %>% 
  lapply(\(x) {
    names(x)[[4]] = "group"
    return(x)
  } 
    ) -> plotlist
gene.min.cells.plot(plotlist)

## Seurat QC steps
l.dfs <- make.seurat.obj(plotlist) # should also work on l.df & l.df2
  rm(l.df, l.df2, plotlist)

