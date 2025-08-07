# load only necessary libraries
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(mapscvi)

# try test data first!
# query_romanov = mapscvi::map_new_seurat_hypoMap(mapscvi::query_romanov,
#                                                 suffix="query_romanov",
#                                                 label_col = "C66_named",
#                                                 max_epochs=20)
# 
# query_romanov

# test data works fine, now try real data
load('./data/integrated_seurat_withScType.rda')

# problems with mapscvi and Seurat v5 compatibility...
# also, assay argument not passed from prepare_query to predict_query? 
# edited wrapper function to prepare+predict+project onto reference: 
source("./data/km_rscvi_docker/mapscvi_wrapper.R")

# requires "Batch_ID" metadata 
int.ldfs$Batch_ID <- paste(int.ldfs$orig.ident)

int.ldfs = UPDATEDmap_new_seurat_hypoMap(int.ldfs,
                                         suffix = "integrated.data",
                                         label_col = "C66_named",
                                         max_epochs = 20,
                                         assay = "SCT")
print(int.ldfs)


save(int.ldfs, file = "./data/integrated_seurat_withHypoMap.rda",
     compress = "xz")



