#### Mouse snseq seurat analysis
### Assigning cells with HypoMap
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
## https://github.com/lsteuernagel/mapscvi

### set working directory
setwd("/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/")

#### load libraries ####
##install libraries
# devtools::install_github("lsteuernagel/mapscvi")
# remotes::install_github("mojaveazure/seurat-disk")
#load libraries
library(tidyverse)
library(mapscvi)
library(Seurat)
library(ggalluvial)

#### load data ####
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

##load reference
load('reference_hypoMap_downsample.RData')

## check intersection of data with reference
# 25982 overlap out of 26007 
reference_hypoMap_downsample@assays[["RNA"]]@data@Dimnames[[1]] %>% intersect(mouse.snseq.combined.sct@assays[["RNA"]]@data@Dimnames[[1]]) %>% 
  length()

#### functions ####
### plot_query_labels_test

#' Plot labels on query and reference to compare visually.
#'
#' Is a wrapper around Seurat::DimPlot functions that allows to plot qeury and reference side-by-side or on top of each other
#' The overlay parameters allow to change the appearance of query points and the alpha of the reference points.
#'
#' @param query_seura_object query seurat object
#' @param reference_seurat reference seurat object
#' @param label_col the column name in reference_map metadata with labels to propagate
#' @param label_col_query he column name in query_seura_object metadata with labels to propagate
#' @param overlay overlay query onto label
#' @param bg_col if non null, the reference will take this color and the query ill be overlayed colored by label_col_query. This will ignore overlay_color and overlay_alpha
#' @param overlay_color color of overlay points
#' @param overlay_alpha = alpha for overlay plot
#' @param query_pt_size numeric for pt size of query cells. defaults to NULL which will inherit the pt.size from the reference DimPlot
#' @param query_umap name of query umap. defaults to "umap_scvi"
#' @param reference_umap name of reference umap. defaults to "umap_scvi"
#' @param labelonplot put labels on plot. defaults to TRUE
#' @param cols_plot passed to Seurat::DimPlot cols argument. Defaults to NULL
#' @param noaxes don't plot axes on UMAP. defaults to TRUE
#' @param nolegend don't plot legend. defaults to TRUE
#' @param ... additional arguments to Seurat::DimPlot
#'
#' @return plot
#'
#' @export
#'
#' @import SeuratObject Seurat cowplot ggplot2


plot_query_labels_test = function(query_seura_object,reference_seurat,label_col,label_col_query = "predicted", overlay = FALSE, bg_col = "grey80", overlay_color = "red", overlay_alpha = 0.5,query_pt_size=NULL, query_umap = "umap_scvi",reference_umap="umap_scvi",labelonplot=TRUE,cols_plot=NULL,noaxes=TRUE,nolegend=TRUE,...){
  
  # check
  if(is.null(reference_seurat)){stop("Please provide reference seurat with latent space, umap and metadata")}
  if(! (reference_umap %in% names(reference_seurat@reductions))){stop("Cannot find '",reference_umap,"' in provided reference_seurat.") }
  if(! (query_umap %in% names(query_seura_object@reductions))){ stop("Cannot find '",reference_umap,"' in provided reference_seurat.")}
  if(! label_col %in% colnames(reference_seurat@meta.data)){stop("Cannot find '",label_col,"' in reference_seurat to label data.") }
  
  # overlay mode
  if(overlay){
    # extract data for overlay from query
    plot_data = cbind(query_seura_object@reductions[[query_umap]]@cell.embeddings,query_seura_object@meta.data)
    # plot reference UMAP
    p_full=Seurat::DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap,label = labelonplot, raster = FALSE,...)
    # save and remove geom_text layer
    if(labelonplot){
      save_geom_text = p_full$layers[[2]]
      p_full$layers[[2]] =NULL
    }
    # recreate plot if all points are bg col
    if(!is.null(bg_col)){
      reference_seurat$dummy = NA
      p_full=Seurat::DimPlot(reference_seurat,group.by = "dummy", raster = FALSE)+scale_color_manual(values = bg_col,na.value=bg_col)
    }
    if(noaxes){p_full = p_full+Seurat::NoAxes()}
    # adjust alpha
    p_full[[1]]$layers[[1]]$aes_params$alpha = min(overlay_alpha,1)
    # get pt size
    if(is.null(query_pt_size)){
      pt_size = p_full[[1]]$layers[[1]]$aes_params$size
    }else{
      pt_size = query_pt_size
    }
    # plot query points on top
    if(!is.null(bg_col)){
      # if bg color is set we are plotting the color with the query points
      # do the label_col and label_col query overlap ? then use color scale from reference
      if(length(intersect(unique(query_seura_object@meta.data[,label_col_query]),unique(reference_seurat@meta.data[,label_col])))>0){
        test_plot = DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap, raster = FALSE,cols=cols_plot) # testplot from full data
        testplot_build=ggplot_build(test_plot)$data[1][[1]] # dataframe with colors
        color_mapping_df=as.data.frame(cbind(testplot_build[,"colour"],reference_seurat@meta.data[,c(label_col)])) %>% dplyr::distinct(V1,V2) # make a df with label_col and colours
        color_mapping <- as.character(color_mapping_df$V1) # convert to a named vector for scale_color_manual
        names(color_mapping) <- color_mapping_df$V2
        # add points to plot
        p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2],color=label_col_query),size=pt_size)+
          ggplot2::scale_color_manual(values = color_mapping,na.value= bg_col)
      }else{# if not use default mapping
        p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2],color=label_col_query),size=pt_size)
      }
      p_full = p_full+ggtitle(label_col_query)
    }else{
      # if no bg color use overlay color instead
      p_full=p_full+ggplot2::geom_point(data=plot_data,ggplot2::aes_string(x=colnames(plot_data)[1],y=colnames(plot_data)[2]),size=pt_size,color=overlay_color)
    }
    # add geom_labels back
    if(labelonplot){p_full$layers[[3]] = save_geom_text}
  }else{
    # need labels in query
    if(! label_col_query %in% colnames(query_seura_object@meta.data)){
      stop("Cannot find '",label_col_query,"' in query_seura_object to label data.")
    }
    # browser()
    # side-by-side
    xlims = c(min(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,1])-0.5,max(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,1])+0.5)
    ylims = c(min(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,2])-0.5,max(reference_seurat@reductions[[reference_umap]]@cell.embeddings[,2])+0.5)
    p1 = Seurat::DimPlot(reference_seurat,group.by = label_col,reduction = reference_umap,label = labelonplot, raster = FALSE,...)+xlim(xlims)+ylim(ylims)
    # get pt size
    if(is.null(query_pt_size)){
      pt_size = p1[[1]]$layers[[1]]$aes_params$size
    }else{
      pt_size = query_pt_size
    }
    p2 = Seurat::DimPlot(query_seura_object,group.by = label_col_query,reduction = query_umap,label = labelonplot,pt.size = pt_size, raster = FALSE,...)+xlim(xlims)+ylim(ylims)
    if(noaxes){
      p1 = p1+Seurat::NoAxes()
      p2 = p2+Seurat::NoAxes()
    }
    
    p_full = cowplot::plot_grid(p1,p2)
    
  }
  
  p_full
}

#### Map data to HypoMap (not working) ####
### prepare mouse data
test = prepare_query(mouse.snseq.combined.sct,
                     suffix="mouse.snseq",
                     normalize=FALSE,
                     batch_var = 'orig.ident',
                     assay = 'SCT',
                     global_seed = 12345)


names(test@reductions)

### mapping mouse data to hypomap reference
test = mapscvi::map_new_seurat_hypoMap(mouse.snseq.combined.sct,
                                                suffix="mouse.snseq",
                                                max_epochs=20,
                                       assay = 'SCT')

names(test@reductions)
#> [1] "scvi"      "umap_scvi"
#> 
# Plotting results:
#   We can take a look at the top clusters that were found in the query:
  
  head(sort(table(query_romanov@meta.data$predicted),decreasing = TRUE),n = 10)
#> 
#>                  C66-1: GLU-1          C66-28: Mixed.GABA-2 
#>                           138                            79 
#>              C66-4: Trh.GLU-2         C66-22: Caprin2.GLU-6 
#>                            69                            63 
#>            C66-43: Nts.GABA-1 C66-55: Ermn.Oligodendrocytes 
#>                            60                            38 
#>          C66-29: Vipr2.GABA-2          C66-12: Gpr149.GLU-3 
#>                            33                            29 
#>            C66-18: Rfx4.GLU-4           C66-5: Tent5a.GLU-2 
#>                            29                            26
# The package provides plotting functions to visualize the query cells on the reference:
  
  # We can plot query and reference side by side. Here we set labelonplot to False to prevent cluster labels from being plotted and we use the object ‘reference_hypoMap_full’ that is automatically loaded when using the wrapper above.

plot_query_labels(query_seura_object=query_romanov,reference_seurat=mapscvi::reference_hypoMap_downsample,label_col="C66_named",overlay = FALSE,labelonplot = FALSE)


# Overlay them query over the reference. The overlay parameters allow to change the behavior of query points. We can use the Seurat::DimPlot parameters to further adjust the plots. E.g. by decreasing the size of the labels.

plot_query_labels(query_seura_object=query_romanov,reference_seurat=mapscvi::reference_hypoMap_downsample,label_col="C66_named",overlay = TRUE,query_pt_size = 0.4,labelonplot = FALSE,label.size=1)
#> Scale for 'colour' is already present. Adding another scale for 'colour',
#> which will replace the existing scale.


# The mapping also returns a prediction probability based on the similarity to the the neighbors in there reference which an indicate how well different cells mapped to their assigned celltypes.

Seurat::FeaturePlot(query_romanov,features = "prediction_probability")+Seurat::NoAxes()

#### Map data to HypoMap detailed (not working) (fail?) ####
# The prepare_query function is able to load Seurat, SingleCellExperiment or matrix objects for mapping.

query_seurat_object = prepare_query(mouse.snseq.combined.sct,
                     suffix="mouse.snseq",
                     batch_var = 'orig.ident',
                     assay = 'SCT',
                     global_seed = 12345)
# This new seurat object is compatible with the downstream functions for mapping the data.

# Next, predict_query can be used to embed the query data into the latent space of scvi. We have to specify a model path and the number of epochs for training during mapping (10-20 should be sufficient.).

# The mapscvi package comes with two models for the neuron and full hypothalamus map in the extdata which can be found using system.file as below.

model_path = paste0(system.file('extdata/models/hypoMap_harmonized_scVI_model/"', package = 'mapscvi'),"/")
max_epochs = 20
query_seurat_object = predict_query(query_seurat_object,
                     model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", 
                                              package = 'mapscvi'),
                     max_epochs = max_epochs,
                     assay = 'SCT',
                     global_seed = 12345)
names(test@reductions)
#> [1] "scvi"
# The scvi reduction is a pca-like low dimensional space that can be used to embed the data into the same UMAP as the reference object. This requires an existing UMAP model in the reference Seurat object that was calculated based on the same scvi latent space.
# 
# We use the mapscvi::reference_hypoMap_downsample object that comes with the package and contains for each cell the scvi and umap reduction values as well as metadata.

mapscvi::reference_hypoMap_downsample
#> An object of class Seurat 
#> 57362 features across 104925 samples within 1 assay 
#> Active assay: RNA (57362 features, 0 variable features)
#>  2 dimensional reductions calculated: scvi, umap_scvi
# Then we can calculate nearest neighbors and UMAP based on the reference. Additionally we can project labels (any categorical metadata column from the reference) using the same nearest neighbor graph. This helps with consistency between label propagation and UMAP. However there might be more accurate ways to propagate labels using other classifiers such as a random forest or scANVI.
# 
# To propagate labels with the project_query function we can provide a vector of the same length as the reference cells (and same order!). Preferably this is a column from the metadata of the reference seurat object.

cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C66_named
reference_reduction = "scvi"
test = project_query(query_seurat_object = test,
                                      reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                      reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                      query_reduction = "scvi",
                                      label_vec =cluster_labels)
# This can then be used to plot the results side-by side:
  
  plot_query_labels(query_seura_object=lamanno_seurat_object,reference_seurat=mapscvi::reference_hypoMap_downsample,label_col="C66_named",overlay = FALSE,labelonplot = FALSE)


plot_query_labels(query_seura_object=lamanno_seurat_object,reference_seurat=mapscvi::reference_hypoMap_downsample,label_col="C66_named",overlay = TRUE,query_pt_size = 0.4,labelonplot = FALSE,label.size=1)

head(sort(table(lamanno_seurat_object@meta.data$predicted),decreasing = TRUE),n = 10)

Seurat::FeaturePlot(lamanno_seurat_object,features = "prediction_probability")+Seurat::NoAxes()




#### function breakdown predict query (not working) ####
query_seurat_object = mouse.snseq.combined.sct
model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", 
                         package = 'mapscvi')
max_epochs = 20
assay = 'SCT'
global_seed = 12345
query_reduction = "scvi"
var_names = NULL
use_reticulate = FALSE


# predict_query = function (query_seurat_object, model_path, query_reduction = "scvi", 
#                           var_names = NULL, max_epochs = 30, assay = "RNA", use_reticulate = FALSE, 
#                           global_seed = 12345) 
# {
  # if (!file.exists(paste0(model_path, "model.pt"))) {
  #   stop("Error: Please provide a valid model_path to an scvi model at ", 
  #        model_path)
  # }
  # if (!is.null(var_names)) {
  #   var_df = data.frame(var_names = rownames(SeuratObject::GetAssayData(query_seurat_object, 
  #                                                                       slot = "counts", assay = assay)))
  #   rownames(var_df) = var_df$var_names
  #   included_var_features = intersect(var_features, var_df$var_names)
  #   matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object, 
  #                                                             slot = "counts", assay = assay)[included_var_features, 
  #                                                             ])
  #   var_df = data.frame(var_names = rownames(matrix_for_anndata))
  #   rownames(var_df) = var_df$var_names
  #   message("Matrix for anndata dim ", dim(matrix_for_anndata)[1], 
  #           " ", dim(matrix_for_anndata)[2])
  # }
  # else {
    matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object, 
                                                              slot = "counts", assay = assay))
    message("Matrix for anndata dim ", dim(matrix_for_anndata)[1], 
            " ", dim(matrix_for_anndata)[2])
  # }
  # if (use_reticulate) {
  #   if (!requireNamespace("reticulate", quietly = TRUE)) {
  #     warning("The reticulate package must be installed to use this function when use_reticulate is set to TRUE.")
  #     return(NULL)
    # }
    # pd <- reticulate::import("pandas", convert = FALSE)
    # sc <- reticulate::import("scanpy", convert = FALSE)
    # scvi <- reticulate::import("scvi", convert = FALSE)
    # adata_query <- sc$AnnData(X = t(matrix_for_anndata), 
    #                           obs = query_seurat_object@meta.data, var = var_df)
    # scvi$model$SCVI$prepare_query_anndata(adata_query, model_path)
    # vae_q = scvi$model$SCVI$load_query_data(adata = adata_query, 
    #                                         reference_model = model_path)
    # message(max_epochs)
    # vae_q$train(max_epochs = as.integer(max_epochs), plan_kwargs = list(weight_decay = 0), 
    #             progress_bar_refresh_rate = 0)
    # scvi_prediction = vae_q$get_latent_representation()
    # scvi_prediction = as.matrix(scvi_prediction)
    # colnames(scvi_prediction) = paste0("scVI_", 1:ncol(scvi_prediction))
    # rownames(scvi_prediction) = colnames(matrix_for_anndata)
  # }
  # else {
    temp_dir = paste0(tempdir(), "/")
    temp_seurat = SeuratObject::CreateSeuratObject(counts = matrix_for_anndata, 
                                                   meta.data = query_seurat_object@meta.data, 
                                                   project = query_seurat_object@project.name,
                                                   assay = assay) ### add assay name to match 
    
    h5Seurat_filename = paste0(temp_dir, "temp_", temp_seurat@project.name, 
                               ".h5Seurat")
    SeuratDisk::SaveH5Seurat(object = temp_seurat, filename = h5Seurat_filename, 
                             overwrite = TRUE)
    updated_name = gsub(".h5Seurat", paste0("_", assay, 
                                            ".h5ad"), h5Seurat_filename)
    SeuratDisk::Convert(h5Seurat_filename, dest = updated_name, 
                        assay = assay, verbose = FALSE, overwrite = TRUE)
    output_file = paste0(temp_dir, "predicted_", temp_seurat@project.name, 
                         ".txt")

    
#### error with old CPU    
        system(paste0("python3 -u ", system.file("python/map_scvi2.py", 
                                             package = "mapscvi"), " ", updated_name, " ", model_path, 
                  " ", output_file, " ", max_epochs))
#system("python3 -u /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/python/map_scvi2.py /tmp/Rtmp41iNxH/temp_SeuratProject_SCT.h5ad /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/extdata/models/hypoMap_harmonized_scVI_model/ /tmp/Rtmp41iNxH/predicted_SeuratProject.txt 20")   
        
## Move files from tmp to new location?
#Need to change 'updated_name' and 'output_file'
        
#source("python3 -u /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/python/map_scvi2.py /tmp/Rtmp41iNxH/temp_SeuratProject_SCT.h5ad /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/extdata/models/hypoMap_harmonized_scVI_model/ /tmp/Rtmp41iNxH/predicted_SeuratProject.txt 20")   
      
# temp_dir.2 = paste0(tempfile(tmpdir='/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp'), "/")

temp_dir.2 = '/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/'

temp_seurat = SeuratObject::CreateSeuratObject(counts = matrix_for_anndata, 
                                               meta.data = query_seurat_object@meta.data %>%  
                                                 rownames_to_column('Cell_tmp') %>% 
                                                 mutate(Cell_ID = Cell_tmp,
                                                        Batch_ID = orig.ident) %>% 
                                                 column_to_rownames('Cell_tmp'), ### need to add 'Batch_ID' and 'Cell_ID' column to metadata 
                                               project = query_seurat_object@project.name,
                                               assay = assay) ### add assay name to match 



h5Seurat_filename = paste0(temp_dir.2, "temp_", temp_seurat@project.name, 
                           ".h5Seurat")
SeuratDisk::SaveH5Seurat(object = temp_seurat, filename = h5Seurat_filename, 
                         overwrite = TRUE)
updated_name = gsub(".h5Seurat", paste0("_", assay, 
                                        ".h5ad"), h5Seurat_filename)
SeuratDisk::Convert(h5Seurat_filename, dest = updated_name, 
                    assay = assay, verbose = FALSE, overwrite = TRUE)
output_file = paste0(temp_dir.2, "predicted_", temp_seurat@project.name, 
                     ".txt")


#### error with old CPU 
## need to run command on ccbbcomp02 command line
# system(paste0("python3 -u ", system.file("python/map_scvi2.py", 
#                                          package = "mapscvi"), " ", updated_name, " ", model_path, 
#               " ", output_file, " ", max_epochs))
        

system("python3 -u /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/python/map_scvi2.py /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp/temp_SeuratProject_SCT.h5ad /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/extdata/models/hypoMap_harmonized_scVI_model/ /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/predicted_SeuratProject.txt 20")

system("python3 -u /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/python/map_scvi2.py /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/temp_SeuratProject_SCT.h5ad /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.0/mapscvi/extdata/models/hypoMap_harmonized_scVI_model/ /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/predicted_SeuratProject.txt 20")

### run in R
    scvi_prediction = data.table::fread(output_file, header = TRUE, 
                                        data.table = F)
    rownames_x = as.character(scvi_prediction[, 1])
    scvi_prediction = scvi_prediction[, 2:ncol(scvi_prediction)]
    scvi_prediction = as.matrix(apply(scvi_prediction, 2, 
                                      as.numeric))
    rownames(scvi_prediction) = rownames_x
  # }
  query_dimred <- Seurat::CreateDimReducObject(embeddings = as.matrix(scvi_prediction), 
                                               stdev = as.numeric(apply(scvi_prediction, 2, stats::sd)), 
                                               assay = assay, key = query_reduction)
  query_seurat_object@reductions[[query_reduction]] = query_dimred
  return(query_seurat_object)
# }

#### function breakdown predict query (working) ####
  query_seurat_object = mouse.snseq.combined.sct
  query_seurat_object = prepare_query(query_seurat_object,
                       suffix="mouse.snseq",
                       normalize=FALSE,
                       batch_var = 'orig.ident',
                       assay = 'SCT',
                       global_seed = 12345)
  
  model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", 
                           package = 'mapscvi')
  max_epochs = 20
  assay = 'SCT'
  global_seed = 12345
  query_reduction = "scvi"
  var_names = NULL
  use_reticulate = FALSE
  
  
  # predict_query = function (query_seurat_object, model_path, query_reduction = "scvi", 
  #                           var_names = NULL, max_epochs = 30, assay = "RNA", use_reticulate = FALSE, 
  #                           global_seed = 12345) 
  # {
  # if (!file.exists(paste0(model_path, "model.pt"))) {
  #   stop("Error: Please provide a valid model_path to an scvi model at ", 
  #        model_path)
  # }
  # if (!is.null(var_names)) {
  #   var_df = data.frame(var_names = rownames(SeuratObject::GetAssayData(query_seurat_object, 
  #                                                                       slot = "counts", assay = assay)))
  #   rownames(var_df) = var_df$var_names
  #   included_var_features = intersect(var_features, var_df$var_names)
  #   matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object, 
  #                                                             slot = "counts", assay = assay)[included_var_features, 
  #                                                             ])
  #   var_df = data.frame(var_names = rownames(matrix_for_anndata))
  #   rownames(var_df) = var_df$var_names
  #   message("Matrix for anndata dim ", dim(matrix_for_anndata)[1], 
  #           " ", dim(matrix_for_anndata)[2])
  # }
  # else {
  matrix_for_anndata = as.matrix(SeuratObject::GetAssayData(query_seurat_object, 
                                                            slot = "counts", assay = assay))
  message("Matrix for anndata dim ", dim(matrix_for_anndata)[1], 
          " ", dim(matrix_for_anndata)[2])
  # }
  # if (use_reticulate) {
  #   if (!requireNamespace("reticulate", quietly = TRUE)) {
  #     warning("The reticulate package must be installed to use this function when use_reticulate is set to TRUE.")
  #     return(NULL)
  # }
  # pd <- reticulate::import("pandas", convert = FALSE)
  # sc <- reticulate::import("scanpy", convert = FALSE)
  # scvi <- reticulate::import("scvi", convert = FALSE)
  # adata_query <- sc$AnnData(X = t(matrix_for_anndata), 
  #                           obs = query_seurat_object@meta.data, var = var_df)
  # scvi$model$SCVI$prepare_query_anndata(adata_query, model_path)
  # vae_q = scvi$model$SCVI$load_query_data(adata = adata_query, 
  #                                         reference_model = model_path)
  # message(max_epochs)
  # vae_q$train(max_epochs = as.integer(max_epochs), plan_kwargs = list(weight_decay = 0), 
  #             progress_bar_refresh_rate = 0)
  # scvi_prediction = vae_q$get_latent_representation()
  # scvi_prediction = as.matrix(scvi_prediction)
  # colnames(scvi_prediction) = paste0("scVI_", 1:ncol(scvi_prediction))
  # rownames(scvi_prediction) = colnames(matrix_for_anndata)
  # }
  # else {
  # needed to create folder
  temp_dir.2 = '/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/'
  
  temp_seurat = SeuratObject::CreateSeuratObject(counts = matrix_for_anndata, 
                                                 meta.data = query_seurat_object@meta.data %>%  
                                                   rownames_to_column('Cell_tmp') %>% 
                                                   mutate(Cell_ID = Cell_tmp,
                                                          Batch_ID = orig.ident) %>% 
                                                   column_to_rownames('Cell_tmp'), ### need to add 'Batch_ID' and 'Cell_ID' column to metadata 
                                                 project = query_seurat_object@project.name,
                                                 assay = assay) ### add assay name to match 
  
  
  
  h5Seurat_filename = paste0(temp_dir.2, "temp_", temp_seurat@project.name, 
                             ".h5Seurat")

## save file
  SeuratDisk::SaveH5Seurat(object = temp_seurat, filename = h5Seurat_filename, 
                           overwrite = TRUE)
  updated_name = gsub(".h5Seurat", paste0("_", assay, 
                                          ".h5ad"), h5Seurat_filename)
  SeuratDisk::Convert(h5Seurat_filename, dest = updated_name, 
                      assay = assay, verbose = FALSE, overwrite = TRUE)
  output_file = paste0(temp_dir.2, "predicted_", temp_seurat@project.name, 
                       ".txt")
  
  
  #### error with old CPU! 
  ## need to run command on ccbbcomp02 command line
  # system(paste0("python3 -u ", system.file("python/map_scvi2.py",
  #                                          package = "mapscvi"), " ", updated_name, " ", model_path,
  #               " ", output_file, " ", max_epochs))
  
  # paste0("python3 -u ", system.file("python/map_scvi2.py",
  #                                                                            package = "mapscvi"), " ", updated_name, " ", model_path,
  #                                                 " ", output_file, " ", max_epochs)
  
  # ##same thing written out
  # python3 -u /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.3/mapscvi/python/map_scvi2.py /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/temp_mouse.snseq_SCT.h5ad /stor/home/imc/R/x86_64-pc-linux-gnu-library/4.3/mapscvi/extdata/models/hypoMap_harmonized_scVI_model/ /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/tmp2/predicted_mouse.snseq.txt 20
  # 
  ### run in R
  scvi_prediction = data.table::fread(output_file, header = TRUE, 
                                      data.table = F)
  rownames_x = as.character(scvi_prediction[, 1])
  scvi_prediction = scvi_prediction[, 2:ncol(scvi_prediction)]
  scvi_prediction = as.matrix(apply(scvi_prediction, 2, 
                                    as.numeric))
  rownames(scvi_prediction) = rownames_x
  # }
  query_dimred <- Seurat::CreateDimReducObject(embeddings = as.matrix(scvi_prediction), 
                                               stdev = as.numeric(apply(scvi_prediction, 2, stats::sd)), 
                                               assay = assay, key = query_reduction)
  query_seurat_object@reductions[[query_reduction]] = query_dimred
  return(query_seurat_object)
  # }
  
#### project_query ####
### specify labels and reduction  
cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C66_named
reference_reduction = "scvi"

### need to add 'Batch_ID' and 'Cell_ID' column to metadata 
query_seurat_object@meta.data = query_seurat_object@meta.data %>% 
rownames_to_column('Cell_tmp')%>% 
  mutate(Cell_ID = Cell_tmp,
         Batch_ID = orig.ident) %>% 
  column_to_rownames('Cell_tmp')


### project data into UMAP space from reference
## predict cell types
query_seurat_object = project_query(query_seurat_object = query_seurat_object,
                                        reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                        reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                        query_reduction = "scvi",
                                        label_vec =cluster_labels)  
  
### save to original object
mouse.snseq.combined.sct = query_seurat_object

#### graph ####
### projection UMAP and Predicted UMAP
plot_query_labels(query_seura_object=mouse.snseq.combined.sct,
                    reference_seurat=mapscvi::reference_hypoMap_downsample,
                    label_col="C66_named",
                    overlay = FALSE,
                    labelonplot = FALSE)
ggsave('HypoMap/Projection UMAP and predicted.png',
       height = 10,
       width = 10)


### Just projection onto total
plot_query_labels(query_seura_object=mouse.snseq.combined.sct,
                  # reference_seurat=mapscvi::reference_hypoMap_downsample,
                  reference_seurat=reference_hypoMap_downsample,
                  label_col="C66_named",
                  overlay = TRUE,
                  query_pt_size = 0.4,
                  labelonplot = TRUE,
                  label.size=1)
ggsave('HypoMap/Predicted UMAP projection.png',
       height = 10,
       width = 10)




### Prediction probability
Seurat::FeaturePlot(mouse.snseq.combined.sct,
                    features = "prediction_probability")+
  Seurat::NoAxes()
ggsave('HypoMap/Predicted probability UMAP.png',
       height = 10,
       width = 10)
  
### Prediction probability
Seurat::DimPlot(mouse.snseq.combined.sct,
                    group.by = "predicted")+
  Seurat::NoAxes() 
ggsave('HypoMap/Predicted UMAP.png',
       height = 10,
       width = 10)  
  
  
  
  

#### compare cell type results ####
### create broad cell type category for C66_labeled
cell.type.hypomap.data = readxl::read_excel('42255_2022_657_MOESM3_ESM_HypoMap_data.xlsx',
                                            sheet = 6)

## subset into the cell types used for hypomap
cell.type.hypomap.data = cell.type.hypomap.data %>% 
  dplyr::select(c(cluster_name,
                  parent_id)) %>% 
  droplevels() %>% 
  distinct()

## pull out cell types of interest
# create list of cells
mouse.snseq.cell.hypomap.list = mouse.snseq.combined.sct@meta.data %>% 
  pull(predicted) %>% 
  unique() %>% 
  c()

# filter down to those cells
cell.type.hypomap.data = cell.type.hypomap.data %>% 
  filter(cluster_name %in% mouse.snseq.cell.hypomap.list)

### load parent id information
cell.type.hypomap.data.parent.id = read.csv('HypoMap_cluster_parent_ids.csv')

## add parent id
cell.type.hypomap.data = cell.type.hypomap.data %>% 
  left_join(cell.type.hypomap.data.parent.id) %>% 
  dplyr::select(-c(parent_id))
  
## create metadata of broad cell type
metadata.mouse.snseq = mouse.snseq.combined.sct@meta.data
# add broad cell data
metadata.mouse.snseq = metadata.mouse.snseq %>% 
  left_join(cell.type.hypomap.data %>% 
              dplyr::rename(predicted = cluster_name))

## select relevant variables
metadata.mouse.snseq = metadata.mouse.snseq %>% 
  column_to_rownames('Cell_ID') %>% 
  dplyr::select(c(predicted,
           parent_id.exp,
           parent_id.broad))

### add metadata to seurat object
mouse.snseq.combined.sct = AddMetaData(object = mouse.snseq.combined.sct,
                                       metadata = metadata.mouse.snseq)





#### add unknown cells ####
#### project_query prob cell type score 
### specify labels and reduction  
cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C66_named
reference_reduction = "scvi"

### need to add 'Batch_ID' and 'Cell_ID' column to metadata 
query_seurat_object = mouse.snseq.combined.sct
query_seurat_object@meta.data = query_seurat_object@meta.data %>% 
  rownames_to_column('Cell_tmp')%>% 
  mutate(Cell_ID = Cell_tmp,
         Batch_ID = orig.ident) %>% 
  column_to_rownames('Cell_tmp')


### project data into UMAP space from reference
## predict cell types
query_seurat_object = project_query(query_seurat_object = query_seurat_object,
                                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                    query_reduction = "scvi",
                                    label_vec =cluster_labels,
                                    result_type="all")  

## add unknown 
### use prediction probability to assign cell type
predicted.prob.df = query_seurat_object@meta.data

# predicted tmp
tmp = project_query(query_seurat_object = query_seurat_object,
                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                    query_reduction = "scvi",
                    label_vec =cluster_labels)


tmp.df = tmp@meta.data


tmp.df = tmp.df %>% 
  mutate(predicted.prob.tmp = ifelse(prediction_probability >= 0.75,
                                     predicted,
                                     'Unknown'))

## combine dataframes
predicted.prob.df = predicted.prob.df %>% 
  rownames_to_column('cell.id') %>% 
  full_join(tmp.df %>% 
              dplyr::select(c(predicted.prob.tmp,
                              prediction_probability,
                              predicted)) %>% 
              rownames_to_column('cell.id'))



## create cell type score
test.df = predicted.prob.df %>%
  mutate(neuron = rowSums(dplyr::select(., contains(c('GABA', 'GLU'))))) %>%
  mutate(astrocytes = rowSums(dplyr::select(., contains(c('Astrocytes', 'Tanycytes','Ependymal'))))) %>%
  mutate(oligodendrocytes = rowSums(dplyr::select(., contains(c('OPC', 'Oligodendrocytes'))))) %>%
  mutate(immune = rowSums(dplyr::select(., contains(c('Immune'))))) %>%
  mutate(vascular = rowSums(dplyr::select(., contains(c('ParsTuber'))))) %>%
  mutate(parstuber = rowSums(dplyr::select(., contains(c('Fibroblasts', 'Endothelial', 'Mural')))))

## create column for comparison
test.df = test.df %>% 
  mutate(neuron.keep = ifelse(neuron >= 0.75, 1, 0)) %>% 
  mutate(astrocytes.keep = ifelse(astrocytes >= 0.75, 1, 0)) %>% 
  mutate(oligodendrocytes.keep = ifelse(oligodendrocytes >= 0.75, 1, 0)) %>% 
  mutate(immune.keep = ifelse(immune >= 0.75, 1, 0)) %>% 
  mutate(vascular.keep = ifelse(vascular >= 0.75, 1, 0)) %>% 
  mutate(parstuber.keep = ifelse(parstuber >= 0.75, 1, 0)) 

## check counts
test.df = test.df %>% 
  mutate(unknown.count = rowSums(dplyr::select(., contains(c('.keep'))))) %>% 
  mutate(unknown.count = as.numeric(unknown.count))

##graph
test.df %>% 
  dplyr::select(c(predicted,
                  prediction_probability)) %>% 
  droplevels() %>% 
  ggplot(aes(x = prediction_probability,
             fill = prediction_probability > 0.75)) +
  geom_vline(xintercept = 0.75) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(predicted~.) +
  theme_classic() +
  ggtitle('Hypomap predicted probability')
ggsave("HypoMap/Cell type predicted prob.png")

test.df %>% 
  dplyr::select(c(neuron,
                  predicted.prob.tmp,
                  prediction_probability,
                  predicted)) %>% 
  droplevels() %>% 
  ggplot(aes(x = neuron,
             color = neuron > 0.75,
             y = prediction_probability)) +
  geom_vline(xintercept = 0.75) +
  geom_hline(yintercept = 0.75) +
  geom_point() +
  facet_wrap(predicted~.) +
  theme_classic() +
  ylim(0,1)+
  ggtitle('Hypomap neuron score')
ggsave("HypoMap/Neuron score predicted prob.png")


test.df %>% 
  dplyr::select(c(neuron,
                  predicted.prob.tmp,
                  prediction_probability,
                  predicted)) %>% 
  filter(prediction_probability < 0.75) %>% 
  droplevels() %>% 
  ggplot(aes(x = neuron,
             fill = neuron > 0.75)) +
  geom_vline(xintercept = 0.75) +
  geom_hline(yintercept = 0.75) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(predicted~.) +
  theme_classic()+
  ggtitle('Hypomap Unknown Neuron score')
ggsave("HypoMap/Neuron score predicted prob histogram.png")

### use prediction probability to assign cell type
# celltype call
test.df.meta = test.df %>% 
  mutate(predicted.prob = ifelse(unknown.count > 0,
                                 predicted,
                                 'Unknown'),
         parent_id.exp.prob = ifelse(unknown.count > 0,
                                     parent_id.exp,
                                 'Unknown'),
         parent_id.broad.prob = ifelse(unknown.count > 0,
                                       parent_id.broad,
                                 'Unknown'))
## add metadata to object
mouse.snseq.combined.sct = AddMetaData(mouse.snseq.combined.sct,
                                       metadata = test.df.meta$predicted.prob ,
                                       col.name = 'predicted.prob')
mouse.snseq.combined.sct = AddMetaData(mouse.snseq.combined.sct,
                                       metadata = test.df.meta$parent_id.exp.prob ,
                                       col.name = 'parent_id.exp.prob')
mouse.snseq.combined.sct = AddMetaData(mouse.snseq.combined.sct,
                                       metadata = test.df.meta$parent_id.broad.prob ,
                                       col.name = 'parent_id.broad.prob')





#### save data point ####
##just single cell object
# save(mouse.snseq.combined.sct,
#      file = "mouse.snseq.combined.sct.RData")
# load('mouse.snseq.combined.sct.RData')

#### graph comparison ####
### graph parent umaps
## Prediction broad
Seurat::DimPlot(mouse.snseq.combined.sct,
                group.by = "parent_id.broad.prob",
                cols = c('#1b9e77',
                  '#d95f02',
                  '#7570b3',
                  '#e7298a',
                  '#66a61e',
                  '#e6ab02',
                  '#a6761d',
                  'grey'))+
  Seurat::NoAxes() 
ggsave('HypoMap/Predicted broad UMAP.png',
       height = 10,
       width = 10)  


##
Seurat::DimPlot(mouse.snseq.combined.sct,
                group.by = "parent_id.broad.prob.label",
                cols = c('#1b9e77',
                         '#d95f02',
                         '#7570b3',
                         '#e7298a',
                         '#66a61e',
                         '#e6ab02',
                         '#a6761d',
                         'grey'))+
  Seurat::NoAxes() +
  ggtitle('') +
  theme(text = element_text(size = 15,
                            face = "bold"))
ggsave('HypoMap/Predicted broad UMAP poster legend.png',
       height = 5.25,
       width = 5.25)  

## predicted exp
Seurat::DimPlot(mouse.snseq.combined.sct,
                group.by = "parent_id.exp.prob")+
  Seurat::NoAxes() 
ggsave('HypoMap/Predicted expanded UMAP.png',
       height = 10,
       width = 10)  

### graph cell counts across samples
## broad
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(orig.ident,
                  parent_id.broad.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = parent_id.broad.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts broad.png',
       height = 10,
       width = 10)  

## expanded
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(orig.ident,
                  parent_id.exp.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = parent_id.exp.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts expanded.png',
       height = 10,
       width = 10)  

#reduce 
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(orig.ident,
                  parent_id.exp.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  ggplot(aes(x = reorder(parent_id.exp.prob,
                         -Freq),
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts expanded reduce.png',
       height = 10,
       width = 10)  

## predicted
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = predicted.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts predicted.png',
       height = 10,
       width = 10)

#reduce 
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  ggplot(aes(x = reorder(predicted.prob,
                         -Freq),
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts predicted reduce.png',
       height = 10,
       width = 10)

#reduce 
#neurons
#genotype
##get list of cell types with at least 50 in one sample
neuron.cell.types.predicted.50 = mouse.snseq.combined.sct@meta.data %>% 
  filter(parent_id.broad.prob %in% c("C7-2: GABA",
                                     "C7-1: GLU")) %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob,
                  Genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  pull(predicted.prob) %>% 
  unique()
#graph
mouse.snseq.combined.sct@meta.data %>% 
  filter(parent_id.broad.prob %in% c("C7-2: GABA",
                                   "C7-1: GLU")) %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob,
                  Genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(predicted.prob %in% neuron.cell.types.predicted.50) %>%
  ggplot(aes(x = reorder(predicted.prob,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_hline(yintercept = 50) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = Genotype,
                 group=orig.ident),
             size = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts predicted reduce neurons.png',
       height = 10,
       width = 10)



#scale
#reduce 
#neurons
#genotype
##get list of cell types with at least 50 in one sample
neuron.cell.types.predicted.50 = mouse.snseq.combined.sct@meta.data %>% 
  filter(parent_id.broad.prob %in% c("C7-2: GABA",
                                     "C7-1: GLU")) %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob,
                  Genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  pull(predicted.prob) %>% 
  unique()
#dataframe of scale
neuron.genotype.scale = mouse.snseq.combined.sct@meta.data %>% 
  filter(parent_id.broad.prob %in% c("C7-2: GABA",
                                     "C7-1: GLU")) %>%
  dplyr::select(c(orig.ident,
                  Genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.neuron.count = Freq)

#dataframe of cluster counts
neuron.cluster.genotype.count = mouse.snseq.combined.sct@meta.data %>% 
  filter(parent_id.broad.prob %in% c("C7-2: GABA",
                                     "C7-1: GLU")) %>% 
  dplyr::select(c(orig.ident,
                  predicted.prob,
                  Genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(predicted.prob %in% neuron.cell.types.predicted.50) %>% 
  as.data.frame()

#combine counts and scale
neuron.cluster.genotype.count = neuron.cluster.genotype.count %>% 
  full_join(neuron.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.neuron.count)

#graph
neuron.cluster.genotype.count %>% 
  ggplot(aes(x = reorder(predicted.prob,
                         -Percent),
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = Genotype,
                 group=orig.ident),
             size = 3) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('HypoMap/Predicted cell counts predicted reduce neurons scale.png',
       height = 10,
       width = 10)

# ### graph alluvial
# ### hypomap predicted vs parent ids 
# mouse.snseq.combined.sct@meta.data %>%
#   select(c(parent_id.broad.prob,
#            parent_id.exp.prob,
#            predicted.prob)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis1 = reorder(predicted.prob, -Freq.scale),
#              axis2 = reorder(parent_id.exp.prob,-Freq.scale),
#              axis3 = reorder(parent_id.broad.prob,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.broad.prob)) +
#   geom_stratum() +
#   scale_x_discrete(limits = c("predicted.prob",
#                               "parent_id.exp.prob",
#                               "parent_id.broad.prob"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial hypomap vs parent_id broad and expanded.png',
#        height = 10,
#        width = 10)


# ### hypomap sctype expanded vs hypomap broad
# mouse.snseq.combined.sct@meta.data %>%
#   select(c(sctype.all.hypomap.exp,
#            parent_id.broad.prob)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.broad.prob,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.broad.prob)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#                   aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.broad.prob"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.broad.png',
#        height = 10,
#        width = 10)
# 
# ## predicted
# mouse.snseq.combined.sct@meta.data %>%
#   mutate(parent_id.broad.pred = ifelse(prediction_probability > 0.50, 
#                                        parent_id.broad,
#                                        "Unknown")) %>% 
#   select(c(sctype.all.hypomap.exp,
#            parent_id.broad.pred)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.broad.pred,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.broad.pred)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.broad.pred"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.broad predicted.png',
#        height = 10,
#        width = 10)
# 
# ### hypomap sctype expanded vs hypomap expanded
# mouse.snseq.combined.sct@meta.data %>%
#   select(c(sctype.all.hypomap.exp,
#            parent_id.exp)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.exp,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.exp)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.exp"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.exp.png',
#        height = 10,
#        width = 10)
# 
# ## predicted
# mouse.snseq.combined.sct@meta.data %>%
#   mutate(parent_id.exp.pred = ifelse(prediction_probability > 0.50, 
#                                      parent_id.exp,
#                                        "Unknown")) %>% 
#   select(c(sctype.all.hypomap.exp,
#            parent_id.exp)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.exp,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.exp)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.exp"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.exp predicted.png',
#        height = 10,
#        width = 10)
# 
# ### hypomap sctype expanded vs hypomap 
# mouse.snseq.combined.sct@meta.data %>%
#   select(c(sctype.all.hypomap.exp,
#            parent_id.exp)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.exp,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.exp)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.exp"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.exp.png',
#        height = 10,
#        width = 10)
# 
# ## predicted
# mouse.snseq.combined.sct@meta.data %>%
#   mutate(parent_id.exp.pred = ifelse(prediction_probability > 0.50, 
#                                      parent_id.exp,
#                                      "Unknown")) %>% 
#   select(c(sctype.all.hypomap.exp,
#            parent_id.exp)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
#              axis3 = reorder(parent_id.exp,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = parent_id.exp)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctype.all.hypomap.exp",
#                               "parent_id.exp"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('HypoMap/Alluvial sctype.hypomap.exp vs parent_id.exp predicted.png',
#        height = 10,
#        width = 10)

###  sctype expanded vs hypomap broad
mouse.snseq.combined.sct@meta.data %>%
  select(c(sctype.all,
           parent_id.broad.prob)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis2 = reorder(sctype.all, -Freq.scale),
             axis3 = reorder(parent_id.broad.prob,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = parent_id.broad.prob)) +
  geom_stratum() +
  geom_text(stat = "stratum",
                  aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all",
                              "parent_id.broad.prob"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('HypoMap/Alluvial sctype.all vs parent_id.broad.png',
       height = 10,
       width = 10)

### for presentation
# projection with broad cell types
plot_query_labels(query_seura_object=mouse.snseq.combined.sct,
                  reference_seurat=reference_hypoMap_downsample,
                  label_col_query = "parent_id.broad",
                  label_col="C7_named",
                  overlay = TRUE,
                  query_pt_size = 0.4,
                  labelonplot = TRUE,
                  label.size=5
)
ggsave('HypoMap/Predicted UMAP projection broad label.png',
       height = 10,
       width = 10)
#without label
plot_query_labels_test(query_seura_object=mouse.snseq.combined.sct,
                       reference_seurat=reference_hypoMap_downsample,
                       label_col_query = "parent_id.broad.prob",
                       label_col="C7_named",
                       overlay = TRUE,
                       query_pt_size = 0.4,
                       labelonplot = FALSE,
                       label.size=5
)
ggsave('HypoMap/Predicted UMAP projection broad.png',
       height = 10,
       width = 10)

#without label
# add unknown 
plot_query_labels_test(query_seura_object=mouse.snseq.combined.sct,
                       reference_seurat=reference_hypoMap_downsample,
                       label_col_query = "parent_id.broad.prob",
                       label_col="C7_named",
                       overlay = TRUE,
                       query_pt_size = 0.4,
                       labelonplot = FALSE,
                       label.size=5,
                       cols = c('#1b9e77',
                                '#d95f02',
                                '#7570b3',
                                '#e7298a',
                                '#66a61e',
                                '#e6ab02',
                                '#a6761d',
                                'grey')
)
ggsave('HypoMap/Predicted UMAP projection broad unknown.png',
       height = 10,
       width = 10)


#### check unknown cell types ####
### compare known to unknown cell counts
## broad cell types
mouse.snseq.combined.sct@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  dplyr::select(c(parent_id.broad,
                  orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(mouse.snseq.combined.sct@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              dplyr::select(c(parent_id.broad,
                              orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(parent_id.broad,
                         -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  theme_classic()  +theme(axis.text.x = element_text(angle = 45, 
                                                               vjust = 0.5,
                                                               hjust=0.5)) +
  xlab('')
ggsave('HypoMap/Known vs unknown cell counts broad.png',
       height = 10,
       width = 10)  

## expanded cell types
mouse.snseq.combined.sct@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  dplyr::select(c(parent_id.exp,
                  orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(mouse.snseq.combined.sct@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              dplyr::select(c(parent_id.exp,
                              orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(parent_id.exp,
                         -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  theme_classic()  +theme(axis.text.x = element_text(angle = 45, 
                                                     vjust = 0.5,
                                                     hjust=0.5)) +
  xlab('')
ggsave('HypoMap/Known vs unknown cell counts exp.png',
       height = 10,
       width = 10)  

## expanded cell types reduce
mouse.snseq.combined.sct@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  dplyr::select(c(parent_id.exp,
                  orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(mouse.snseq.combined.sct@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              dplyr::select(c(parent_id.exp,
                              orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  filter(Freq > 50) %>% 
  ggplot(aes(x = reorder(parent_id.exp,
                         -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  theme_classic()  +theme(axis.text.x = element_text(angle = 45, 
                                                     vjust = 0.5,
                                                     hjust=0.5)) +
  xlab('')
ggsave('HypoMap/Known vs unknown cell counts exp reduce.png',
       height = 10,
       width = 10)  

## predicted cell types
mouse.snseq.combined.sct@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  dplyr::select(c(predicted,
                  orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(mouse.snseq.combined.sct@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              dplyr::select(c(predicted,
                              orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(predicted,
                         -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  theme_classic()  +theme(axis.text.x = element_text(angle = 90, 
                                                     vjust = 1,
                                                     hjust=0)) +
  xlab('')
ggsave('HypoMap/Known vs unknown cell counts.png',
       height = 10,
       width = 10)  

## predicted cell types reduce
mouse.snseq.combined.sct@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  dplyr::select(c(predicted,
                  orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(mouse.snseq.combined.sct@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              dplyr::select(c(predicted,
                              orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  filter(Freq > 50) %>% 
  ggplot(aes(x = reorder(predicted,
                         -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75)) +
  theme_classic()  +theme(axis.text.x = element_text(angle = 90, 
                                                     vjust = 1,
                                                     hjust=0)) +
  xlab('')
ggsave('HypoMap/Known vs unknown cell counts reduce.png',
       height = 10,
       width = 10)  


### check distribution of probability for cell types
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::select(c(predicted,
                  orig.ident,
                  prediction_probability,
                  parent_id.broad)) %>% 
  droplevels() %>% 
  ggplot(aes(x = prediction_probability,
             fill = prediction_probability > 0.25)) +
  geom_vline(xintercept = 0.25) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(parent_id.broad~.) +
  theme_classic()
ggsave('HypoMap/histogram predicted probability.png',
       height = 10,
       width = 10)  



#### project_query spatial ####
### specify labels and reduction  
cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C286_named
reference_reduction = "scvi"

### need to add 'Batch_ID' and 'Cell_ID' column to metadata 
query_seurat_object = mouse.snseq.combined.sct
query_seurat_object@meta.data = query_seurat_object@meta.data %>% 
  rownames_to_column('Cell_tmp')%>% 
  mutate(Cell_ID = Cell_tmp,
         Batch_ID = orig.ident) %>% 
  column_to_rownames('Cell_tmp')


### project data into UMAP space from reference
## predict cell types
query_seurat_object = project_query(query_seurat_object = query_seurat_object,
                                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                    query_reduction = "scvi",
                                    label_vec =cluster_labels)  



# # graph
# plot_query_labels(query_seura_object=query_seurat_object,
#                   reference_seurat=mapscvi::reference_hypoMap_downsample,
#                   label_col="C286_named",
#                   overlay = FALSE,
#                   labelonplot = FALSE)
# 
# Seurat::DimPlot(query_seurat_object,
#                 group.by = "predicted") +
#   Seurat::NoLegend()

## add unknown 
### use prediction probability to assign cell type
predicted.prob.df = query_seurat_object@meta.data
# celltype call
predicted.prob.df = predicted.prob.df %>% 
  mutate(predicted.prob = ifelse(prediction_probability >= 0.25,
                                 predicted,
                                 'Unknown'))

# combine with spatial data
data.spatial = read.csv('Clustser C286 brain annotation hypomap.csv')
predicted.prob.df.spatial = predicted.prob.df %>% 
  dplyr::select(c(predicted)) %>% 
  left_join(data.spatial %>% 
              dplyr::select(c(cluster_name,
                              Region_summarized)) %>% 
              dplyr::rename(predicted = cluster_name))

## add metadata to object
query_seurat_object = AddMetaData(query_seurat_object,
                                       metadata = predicted.prob.df.spatial$Region_summarized,
                                       col.name = 'Region_summarized')

# Seurat::DimPlot(query_seurat_object,
#                 group.by = "predicted.prob") +
#   Seurat::NoLegend()
# 
# 
# Seurat::DimPlot(query_seurat_object,
#                 group.by = "Region_summarized")

## graph counts
query_seurat_object@meta.data %>% 
  dplyr::select(c(orig.ident,
                  Region_summarized)) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(y = reorder(Region_summarized,
                         -Freq),
             x = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')
ggsave('HypoMap/Brain region count all.png',
       height = 10,
       width = 10)  

## graph counts
# known
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  dplyr::select(c(orig.ident,
                  Region_summarized)) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(y = reorder(Region_summarized,
                         -Freq),
             x = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')
ggsave('HypoMap/Brain region count known.png',
       height = 10,
       width = 10)  

# known
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  dplyr::select(c(orig.ident,
                  Region_summarized)) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(y = reorder(Region_summarized,
                         -Freq),
             x = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')
ggsave('HypoMap/Brain region count known.png',
       height = 10,
       width = 10) 

# known
# reduce
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  dplyr::select(c(orig.ident,
                  Region_summarized,
                  Genotype)) %>% 
  table() %>% 
  data.frame() %>% 
  filter(Region_summarized %in% c('Medial preoptic area',
                                  'Lateral preoptic area',
                                  'Lateral hypotalamic area',
                                  'Zona incerta',
                                  'Paraventricular hypothalamic nucleus')) %>% 
  ggplot(aes(x = reorder(Region_summarized,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = Genotype,
                 group=orig.ident),
             size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('') 
ggsave('HypoMap/Brain region count known reduce genotype.png',
       height = 10,
       width = 10) 

# known
# reduce
#scale
nuclei.genotype.count = query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  dplyr::select(c(orig.ident,
                  Region_summarized,
                  Genotype)) %>% 
  table() %>% 
  data.frame() %>% 
  filter(Region_summarized %in% c('Medial preoptic area',
                                  'Lateral preoptic area',
                                  'Lateral hypotalamic area',
                                  'Zona incerta',
                                  'Paraventricular hypothalamic nucleus')) 

#dataframe of scale
nuclei.genotype.scale = query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  dplyr::select(c(orig.ident,
                  Genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.count = Freq)

#combine counts and scale
nuclei.genotype.count = nuclei.genotype.count %>% 
  full_join(nuclei.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.count)

#graph
nuclei.genotype.count %>% 
  ggplot(aes(x = reorder(Region_summarized,
                         -Percent),
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = Genotype,
                 group=orig.ident),
             size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('') 
ggsave('HypoMap/Brain region count known reduce genotype scale.png',
       height = 10,
       width = 10) 





