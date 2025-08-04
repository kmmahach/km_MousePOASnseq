UPDATEDmap_new_seurat_hypoMap = function(query_seurat_object,
                                  suffix = "query",
                                  assay = "RNA",
                                  subset_col = "",
                                  label_col = "",
                                  subset_values = NULL,
                                  max_epochs = 20,
                                  model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", 
                                                         package = 'mapscvi'),
                                  reference_seurat = mapscvi::reference_hypoMap_downsample,
                                  reference_reduction = "scvi",
                                  inferred_sex_varname = "inferred_sex" ,
                                  pythonEnv = NULL, 
                                  use_reticulate = FALSE,
                                  global_seed = 12345){
  
  # prepare
  query_seurat_object = prepare_query(query_seurat_object,
                                      suffix = suffix,
                                      assay = assay,
                                      subset_col = subset_col,
                                      subset_values = subset_values,
                                      normalize = FALSE,
                                      batch_var = "Batch_ID",
                                      global_seed = global_seed)
  
  # check if model path is valid:
  if(!file.exists(paste0(model_path,"/","model.pt"))){
    stop("Error: Please provide a valid model_path to an scvi model at ",model_path)
  }
  
  # predict with scvi
  query_seurat_object = predict_query(query_seurat_object,
                                      model_path,
                                      max_epochs = max_epochs,
                                      assay=assay,
                                      use_reticulate = use_reticulate,
                                      pythonEnv = pythonEnv,
                                      global_seed = global_seed)
  # check
  if(is.null(reference_seurat)){
    stop("Error: Please provide reference seurat with latent space, umap and metadata")
  }
  if(! (reference_reduction %in% names(reference_seurat@reductions) & paste0("umap_",reference_reduction) %in% names(reference_seurat@reductions))){
    stop("Error: Cannot find '",reference_reduction,"' or 'umap_",reference_reduction,"'  in provided reference_seurat.")
  }
  
  # make label_vec
  reference_map_metadata = reference_seurat@meta.data
  label_vec = NULL
  if(!is.null(reference_map_metadata)){
    if(label_col %in% colnames(reference_map_metadata)){
      label_vec = reference_map_metadata[,label_col]
    }else{
      message("label_col does not exist in provided metadata. Cannot propagate labels!")
    }
  }else{
    message("No metadata provided. Cannot propagate labels!")
  }
  # project onto reference
  query_seurat_object = project_query(query_seurat_object,
                                      reference_seurat@reductions[[reference_reduction]],
                                      reference_seurat@reductions[[paste0("umap_",reference_reduction)]],
                                      label_vec = label_vec,
                                      global_seed = global_seed)
  
  # return
  return(query_seurat_object)
}