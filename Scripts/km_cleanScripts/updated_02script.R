# R-4.3.1, Seurat v.5.3.0
# HypoMap reference data from https://github.com/lsteuernagel/mapscvi

root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))

# libraries and functions
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github("lsteuernagel/mapscvi")

lapply(c("tidyverse","Seurat", "ggalluvial", "mapscvi", "SeuratDisk"), 
       library, character.only = T)


#### load data ####
# load integrated Seurat object
# load('./data/mouse.snseq.combined.sct.RData') # IMC final copy
load('./data/integrated_seurat_withScType.rda')

# load HypoMap reference data
reference_hypoMap_downsample <- mapscvi::reference_hypoMap_downsample

# check intersection of data with reference
# 25982 overlap out of 26007 
reference_hypoMap_downsample@assays[["RNA"]]@data@Dimnames[[1]] %>% 
  intersect(dimnames(int.ldfs[["RNA"]])[[1]]) %>% 
  length()

## trying to add back in genes/cell names (works but not necessary?)
  int.ldfs[["RNA"]]@layers$counts@Dimnames = dimnames(int.ldfs[["RNA"]])


# someone else had similar prob w/ RNA assay: 
  # https://github.com/mojaveazure/seurat-disk/issues/27
# here was their solution
DefaultAssay(int.ldfs) <-  'RNA' # temporarily making 'RNA' active assay
int.ldfs <- FindVariableFeatures(int.ldfs)
DefaultAssay(int.ldfs) <-  'SCT' # returning 'SCT' as the default assay

  save(int.ldfs, file = "./data/integrated_seurat_withScType.rda",
       compress = "xz") # save to be safe, although following script may run fine regardless
  
#### Map data to HypoMap ####
  
## Run bash script ./data/km_rscvi_docker/run_scvi_docker.sh
## Requires docker install! If running on lambcomp02 or other rental POD, comes pre-installed
## Do not try to run docker desktop on Windows if using VPN to access data files - can't bind mount to remote server
## If running on local machine (not recommended w/o large RAM), check all filepaths or clone this repo

  # This pulls the docker image from lsteuernagel/mapscvi and modifies it to pre-load all dependencies
  # (see ./data/km_rscvi_docker/Dockerfile),
  # then uses a modified map_new_seurat_hypoMap wrapper function 
  # (./data/km_rscvi_docker/mapscvi_wrapper.R) to perform appropriate data prep and mapping,
  # and finally saves the new Seurat object: ./data/integrated_seurat_withHypoMap.rda
  
load('./data/integrated_seurat_withHypoMap.rda')
  


