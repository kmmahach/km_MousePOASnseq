#### functions ####

# get names of list elements for plotting
namelist <- \(...) {
  group <- as.list(do.call("names", list(...))) 
  return(group)
}

# plot min cell range of genes by number of cells w/ expression
gene.min.cells.range <- function(list_of_dgCMatrices, 
                                 steps_vector, 
                                 min.features = 700) { 
  
  stopifnot("`list_of_dgCMatrices` must be a list" = is.list(list_of_dgCMatrices))
  stopifnot("`steps_vector` must be numeric" = is.numeric(steps_vector))
  stopifnot("`steps_vector` must be a vector" = is.vector(steps_vector))
  
  lapply(list_of_dgCMatrices, \(x) {
      temp.m <- matrix(0, nrow = length(steps), ncol = 2)
  
        for(i in 1:length(steps)) {
    
        temp <- CreateSeuratObject(counts = x,
                                   project = "temp",
                                   min.cells = steps[i],
                                   min.features = 700)
    
      temp.m[i,] = c(steps[i], nrow(temp[["RNA"]]))
        print(paste0("finishing step ", steps[i]))
    
        rm(temp) 
       }

      temp.df <- data.frame(temp.m)
      colnames(temp.df) <- c("min.cells", "gene.count")
      append(x, temp.df)
    }
  )
}

gene.min.cells.plot <- function(list_with_minCells, 
                                outdir = getwd()) {
  
  stopifnot("`list_with_minCells` must be a list" = is.list(list_with_minCells))
  
  newlist = mapply("c",list_with_minCells,
                   namelist(list_with_minCells),SIMPLIFY=FALSE)
  
  lapply(newlist, \(x) {
    
    title <- gsub('.{5}$', '', x[[length(x)]])
    df <- data.frame(cbind(x[["min.cells"]], x[["gene.count"]]))
    colnames(df) <- c("min.cells", "gene.count")
  
    ggplot(df,aes(x = min.cells,
                  y = gene.count,
                  label = paste(min.cells,
                                gene.count,
                                sep = ", "))) +
      geom_line() +
      geom_label() +
      ylim(0,max(df$gene.count)) +
      theme_classic() +
      ggtitle(paste(title)) 
         
    
    if(length(grep(title, list.files(outdir))) > 0) {
      plotname <- paste0(outdir, "/", title, 
                         "/min.cells.", 
                         title, ".png")
    } else {
      
      if(!is.null(outdir)) {
        message("outdir not specified; saving plots to working directory")
        }
      
      plotname <- paste0(outdir, "/min.cells.",
                         title, ".png")
      
      }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 10, 
           bg = "white")
    }
  )
}

# create seurat objects from dcCMatrices
make.seurat.obj <- function(list_of_dgCMatrices) { 
  stopifnot("`list_of_dgCMatrices` must be a list" = is.list(list_of_dgCMatrices))
    
    if(!require(collapse)) {
      message("trying to install collapse")
      install.packages("collapse")
      
      if(!require(collapse)) {
        stop("install collapse package before running")
      }
      
    } else {
      
      newlist = mapply("c",list_of_dgCMatrices,
                       namelist(list_of_dgCMatrices),SIMPLIFY=FALSE)
  
      lapply(newlist, \(x) {
      
        proj.name = x[[length(x)]]
        x = x[-length(x)]
        
        if ({class(x)=="dgCMatrix"}) {
          
          CreateSeuratObject(counts = x,
                             project = paste(proj.name),
                             min.cells = 3,
                             min.features = 200) -> x 
          return(x) 
          
        } else {
      
          get_elem(x, \(x){class(x)=="dgCMatrix"},
                   recursive = TRUE, keep.class = TRUE) -> cts
    
          CreateSeuratObject(counts = cts, 
                             project = paste(proj.name),
                             min.cells = 3,
                             min.features = 200) -> x
          return(x)
        }
      }
    )
  }
}

# violin plots of features, counts, and percent mito genes
plot.qc.metrics <- function(list_of_SeuratObj, 
                            outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  lapply(list_of_SeuratObj, \(x) {
    
    title <- gsub('.{5}$', '', x@project.name)
    
    VlnPlot(x, assay = "RNA", 
            layer = "counts",
            features = c("nFeature_RNA",
                         "nCount_RNA",
                         "percent.mt"),
            ncol = 3, combine = TRUE) +
      plot_annotation(paste(title)) &
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 20,
                                      face = "bold",
                                      hjust = 0.5)) 
    
    
    if(length(grep(title, list.files(outdir))) > 0) {
      plotname <- paste0(outdir, "/", title, 
                         "/feature.count.mito.", 
                         title, ".png")
    } else {
      
      if(!is.null(outdir)) {
        message("outdir not specified; saving plots to working directory")
      } 
      
      plotname <- paste0(outdir, "/feature.count.mito.",
                         title, ".png")
    }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 10, 
           bg = "white")
    }
  )
}

plot.doublets.qc <- function(list_of_SeuratObj, 
                             outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  if(!require(SeuratExtend)) {
    message("trying to install SeuratExtend")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("huayc09/SeuratExtend")
    
    if(!require(SeuratExtend)) {
      stop("install SeuratExtend package before running")
    }
  } else {
    
    lapply(list_of_SeuratObj, \(x) {
      
      title <- gsub('.{5}$', '', x@project.name)
      
      VlnPlot2(x,
               features = c("nFeature_RNA",
                            "nCount_RNA",
                            "percent.mt"),
               group.by = "doublet", 
               scales = "free_y", 
               ncol = 3, nrow = 3) +
        theme(legend.position = "none") +
        ggtitle(paste(title))
      
      if(length(grep(title, list.files(outdir))) > 0) {
        
        plotname <- paste0(outdir, "/", title,
                           "/doublets.",
                           title, ".png")
        
        } else {
          
          if(!is.null(outdir)) {
            message("outdir not specified; saving plots to working directory")
          } 
          
          plotname <- paste0(outdir, "/doublets.",
                         title, ".png")
          }
      
      ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 5, 
           bg = "white")
      }
    )
  }
}

plot.feature.scatter <- function(list_of_SeuratObj, 
                                 filters, 
                                 outdir = getwd()) {
  
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  stopifnot("`filters` must be numeric" = is.numeric(filters))
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(collapse)) {
      stop("install gridExtra package before running")
      }
    
    } else {
      
      lapply(list_of_SeuratObj, \(x) {
        
        my_filters <- filters[(rownames(filters) %in% x@project.name),]
        title <- gsub('.{5}$', '', x@project.name)
        
        scatter1 <-  FeatureScatter(x,
                                    feature1 = "nCount_RNA",
                                    feature2 = "nFeature_RNA") +
          theme(legend.position = "none")
        
        scatter2 <- FeatureScatter(subset(x,
                                          subset = nFeature_RNA < 6000),
                                   feature1 = "nCount_RNA",
                                   feature2 = "nFeature_RNA") +
          geom_hline(yintercept = my_filters["nFeature_min"],
                     linetype = "dotted") +
          geom_hline(yintercept = my_filters["nFeature_max"],
                     linetype = "dotted") +
          annotate("text", x=Inf, y=my_filters["nFeature_min"],
                   label=paste(my_filters["nFeature_min"]),
                   vjust = -0.5, hjust = 1) +
          annotate("text", x=Inf,
                   y=my_filters["nFeature_max"],
                   label=paste(my_filters["nFeature_max"]),
                   vjust = -0.5, hjust = 1) +
          theme(legend.position = "none")
      
      
        hist <- x@meta.data %>%
          subset(nFeature_RNA < 6000) %>%
          mutate(Subset = ifelse(nFeature_RNA > my_filters["nFeature_min"] 
                                 & nFeature_RNA < my_filters["nFeature_max"],
                                 'Keep',
                                 'Remove')) %>%
          ggplot(aes(nFeature_RNA,
                     fill = Subset))+
          geom_histogram(binwidth = 10) +
          theme_bw() +
          ggtitle("UMI counts") +
          scale_fill_manual(values = c('#88e788',
                                       'black'))
        
        
        if(length(grep(title, list.files(outdir))) > 0) {
          
          plotname <- paste0(outdir, "/", title,
                             "/feature.scatter.",
                             title, ".png")
          } else {
            
            if(!is.null(outdir)) {
              message("outdir not specified; saving plots to working directory")
              }
            
            plotname <- paste0(outdir, "/feature.scatter.",
                           title, ".png")
            }
        
        plots <- arrangeGrob(scatter1, scatter2, hist,
                             ncol = 2, widths = c(1,1.5),
                             layout_matrix = rbind(c(1, 3),
                                                   c(2, 3)),
                             
                             top = grid::textGrob(paste(title),
                                                  gp=grid::gpar(fontsize=24)))
        
        ggsave(plotname,
               plots,
               units = "in", 
               width = 10, 
               height = 6.5, 
               bg = "white")
        }
      )
    }
  }

check.filters <- function(list_of_SeuratObj, 
                          filters, 
                          outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(gridExtra)) {
      stop("install gridExtra package before running")
    }
    
  } else {
  
    lapply(list_of_SeuratObj, \(x) {
      
      my_filters <- filters[(rownames(filters) %in% x@project.name),]
      title <- gsub('.{5}$', '', x@project.name)
      
      x1 <- subset(x, subset =
                     nFeature_RNA >= my_filters["nFeature_min"] &
                     nFeature_RNA <= my_filters["nFeature_max"] &
                     percent.mt < my_filters["mito.max"] & 
                     doublet == 'singlet')
      
      xsum <- x1@meta.data %>% 
        summarize(max = max(nFeature_RNA),
                  min = min(nFeature_RNA),
                  median = median(nFeature_RNA),
                  count = n())
      show(xsum)
      
      scatter <- FeatureScatter(x1, 
                                feature1 = "nCount_RNA",
                                feature2 = "nFeature_RNA") +
        theme(legend.position = "none")
      
      hist <- x1@meta.data %>%
        ggplot(aes(nFeature_RNA))+
        geom_histogram(binwidth = 10) +
        theme_bw() +
        ggtitle("filtered UMI counts")
      
      
      if(length(grep(title, list.files(outdir))) > 0) {
        
        plotname <- paste0(outdir, "/", title,
                           "/filtered.feature.scatter.",
                           title, ".png")
      } else {
        
        if(!is.null(outdir)) {
          message("outdir not specified; saving plots to working directory")
        }
        
        plotname <- paste0(outdir, "/filtered.feature.scatter.",
                           title, ".png")
      } 
      
      plots <- arrangeGrob(scatter, hist,
                           ncol = 2, widths = c(1,1.5),
                           
                           top = grid::textGrob(paste(title), 
                                                gp=grid::gpar(fontsize=24)))
      
      ggsave(plotname,
             plots,
             units = "in", 
             width = 10, 
             height = 5, 
             bg = "white")
      }
    )
  }
}

make.filtered.seurat <- function(list_of_SeuratObj, 
                                 filters) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  lapply(list_of_SeuratObj, \(x) {
    
    my_filters <- filters[(rownames(filters) %in% x@project.name),]
    
    subset(x, subset = 
             nFeature_RNA >= my_filters["nFeature_min"] &
             nFeature_RNA <= my_filters["nFeature_max"] &
             percent.mt < my_filters["mito.max"] & 
             doublet == 'singlet')
    }
  )
}

find.cluster.range <- function(list_of_SeuratObj, 
                               outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  if(!require(clustree)) {
    message("trying to install clustree")
    install.packages("clustree")
    
    if(!require(clustree)) {
      stop("install clustree package before running")
    }
  }
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(gridExtra)) {
      stop("install gridExtra package before running")
    }
    
  } else {
    
    lapply(list_of_SeuratObj, \(x) {
      
      title <- gsub('.{5}$', '', x@project.name)
      
      x.clustree <- FindClusters(object = x,
                                 resolution = seq(0,2,0.2))
      
      ctx1 <- clustree(x.clustree,
                       prefix = "SCT_snn_res.") +
        ggtitle('clusters across resolutions') +
        guides(edge_colour = "none", edge_alpha = "none") +
        theme(legend.position = "bottom") +
        scale_edge_color_continuous(low = "black",
                                    high = "black")
      
      x.clust.reduced <- FindClusters(object = x,
                                      resolution = seq(0,1.4,0.2))
      
      ctx2 <- clustree(x.clust.reduced,
                       prefix = "SCT_snn_res.") +
        ggtitle('reduced resolution range') +
        guides(edge_colour = "none", edge_alpha = "none") +
        theme(legend.position = "bottom") +
        scale_edge_color_continuous(low = "black",
                                    high = "black")
      
      if(length(grep(title, list.files(outdir))) > 0) {
        
        plotname <- paste0(outdir, "/", title,
                           "/clustree.resolution.",
                           title, ".png")
      } else {
        
        if(!is.null(outdir)) {
          message("outdir not specified; saving plots to working directory")
        }
        
        plotname <- paste0(outdir, "/clustree.resolution.",
                           title, ".png")
      } 
      
      plots <- arrangeGrob(ctx1, ctx2,
                           ncol = 2, widths = c(1,1),
                           
                           top = grid::textGrob(paste(title), 
                                                gp=grid::gpar(fontsize=24)))
      
      message("**saving plots to outdir**")
      
      ggsave(plotname,
             plots,
             units = "in", 
             width = 15, 
             height = 5, 
             bg = "white")
      } 
    )
  }
}

plot.dim.clust <- function(list_of_SeuratObj, 
                           res = 1, 
                           outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(gridExtra)) {
      stop("install gridExtra package before running")
    }
    
  } else {
  
    lapply(list_of_SeuratObj, \(x) {
      
      title <- gsub('.{5}$', '', x@project.name)
      
      x <- FindClusters(x, resolution = res)
      
      ebp <- ElbowPlot(x)
      vfp <- VariableFeaturePlot_scCustom(x) 
      dmp <- DimPlot(x, 
                     reduction = "umap",
                     label = TRUE,
                     repel = TRUE)
      ftp <- FeaturePlot(x, 
                         reduction = "umap",
                         features = "nFeature_RNA",
                         label = TRUE,
                         repel = TRUE) +
        labs(color = "nFeature_RNA")
      
      plots <- arrangeGrob(ebp,dmp,ftp,vfp,
                           ncol = 2, nrow = 2,
                           
                           top = grid::textGrob(paste(title), 
                                                gp=grid::gpar(fontsize=24)))
      
      if(length(grep(title, list.files(outdir))) > 0) {
        plotname <- paste0(outdir, "/", title, 
                           "/elbow.dim.varfeat.", 
                           title, ".png")
      } else {
      
        if(!is.null(outdir)) {
          message("outdir not specified; saving plots to working directory")
        } 
        
        plotname <- paste0(outdir, "/elbow.dim.varfeat.",
                           title, ".png")
      }
      
      ggsave(plotname,
             plots,
             units = "in", 
             width = 12, 
             height = 10, 
             bg = "white")
      }
    )
  }
}

annotate.with.sctype <- function(list_of_SeuratObj,
                              name = "sctype.ind",
                              outdir = getwd()) {
  stopifnot("`list_of_SeuratObj` must be a list" = is.list(list_of_SeuratObj))
  
  lapply(list_of_SeuratObj, \(x) {
    title <- gsub('.{5}$', '', x@project.name)
    
    # from sctype wrapper function
    # source(https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R)
    message(paste0("using ScType wrapper to annotate brain cell types for ", x@project.name))
    
    es.max <-  sctype_score(scRNAseqData = x[["SCT"]]$scale.data,
                            scaled = TRUE, 
                            gs = gs_list$gs_positive, 
                            gs2 = gs_list$gs_negative)
    
    # Extract top cell types for each cluster
    cL_resutls = do.call("rbind", lapply(unique(x@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(x@meta.data[x@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(x@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    
    x@meta.data[name] = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      x@meta.data[x@meta.data$seurat_clusters == j,name] = as.character(cl_type$type[1])
    }
    
    #graph umap
    DimPlot(x, reduction = "umap",
            group.by = name) +
      ggtitle(paste(title))
    
    if(length(grep(title, list.files(outdir))) > 0) {
      plotname <- paste0(outdir, "/", title, 
                         "/sctype.dimplot.", 
                         title, ".png")
    } else {
      
      if(!is.null(outdir)) {
        message("outdir not specified; saving plots to working directory")
      }
      
      plotname <- paste0(outdir, "/sctype.dimplot.",
                         title, ".png")
      
    }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 8, 
           bg = "white")
    
    return(x)
    }
  )
}


#### from seurat_snseq_mouse_IMC.R ####
### modify gene_sets_prepare function to integrate with user added data
## do not require xlsx for gene_sets_prepare
## remove checkGeneSymbols function 
# remove forcing genes to uppercase
gene_sets_prepare.df = function(db_file, cell_type){
  
  cell_markers = db_file #change to dataframe
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}
