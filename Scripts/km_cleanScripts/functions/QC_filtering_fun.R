#### accessory functions ####
# package version and dependency log
load_packages <- function(pkgs, 
                          out_prefix = "package_log", 
                          auto_install = TRUE) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- "package_logs"
  if (!dir.exists(out_dir)) dir.create(out_dir)
  report_file <- file.path(out_dir, paste0(out_prefix, "_package_report_", timestamp, ".txt"))
  
  # install missing from CRAN/Bioconductor
  install_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("🔍 Installing missing package: ", pkg)
      tryCatch({
        install.packages(pkg)
        if (!requireNamespace(pkg, quietly = TRUE)) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }
      }, error = function(e) {
        message("✖ Could not install ", pkg, ": ", e$message)
      })
    }
  }
  
  if (auto_install) {
    missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
    if (length(missing_pkgs) > 0) invisible(lapply(missing_pkgs, install_missing))
  }
  
  con <- file(report_file, open = "wt")
  on.exit(close(con), add = TRUE)  
  
  writeLines(c(
    "=== PACKAGE LOG REPORT ===",
    paste("Generated on:", Sys.time()),
    paste("Packages:", paste(pkgs, collapse = ", ")), "",
    "---- Package Load Status ----"
  ), con)
  
  for (pkg in pkgs) {
    cat("\n📦 ", pkg, "\n")          # console
    writeLines(paste0("\n📦 ", pkg), con)  # file
    
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      cat("✔ Loaded successfully\n")        # console
      writeLines("✔ Loaded successfully", con) # file
    }, error = function(e) {
      cat("✖ Failed to load: ", e$message, "\n")   # console
      writeLines(paste0("✖ Failed to load: ", e$message), con) # file
    })
  }
  
  cat("\n---- Package Versions ----\n")
  writeLines("\n---- Package Versions ----", con)
  versions <- sapply(pkgs, function(pkg) {
    tryCatch(as.character(packageVersion(pkg)), error = function(e) "NOT INSTALLED")
  })
  print(versions) 
  writeLines(capture.output(print(versions)), con) 
  
  writeLines("\n---- Package Dependencies (Recursive) ----", con)
  deps <- tools::package_dependencies(pkgs, recursive = TRUE)
  writeLines(capture.output(print(deps)), con)
  
  writeLines("\n---- Unique Dependencies ----", con)
  unique_deps <- sort(unique(unlist(deps)))
  writeLines(capture.output(print(unique_deps)), con)
  
  message("✅ Combined package log saved to: ", report_file)
}

# get names of list elements for plotting
namelist <- \(...) {
  group <- as.list(do.call("names", list(...))) 
  return(group)
}

# make list with automatic name propagation 
as_named_list <- function(...) {
  dots <- rlang::enquos(...)
  
  # name and evaluate
  named_list <- purrr::imap(dots, ~ rlang::eval_tidy(.x) |> 
                              `names<-`(rlang::as_name(.x)))
  # flatten
  setNames(purrr::map(dots, rlang::eval_tidy), purrr::map_chr(dots, rlang::as_name))
}

#### Seurat QC and plotting functions ####
# note: will work with Seurat v5 and Seurat v4;
# but may produce different outputs due to changes in Seurat functions!

# compatible version of scCustomize variable feature plot
VariableFeaturePlot_patched <- function(object,
                                        assay = "SCT",
                                        nfeatures = 3000,
                                        label = TRUE,
                                        repel = TRUE) {
  sct_model <- object[[assay]]@SCTModel.list$model1
  
  if (is.null(sct_model)) {
    stop("No SCT model found. Make sure SCTransform() has been run.")
  }
  
  df <- sct_model@feature.attributes
  df$feature <- rownames(df)
  
  # Order by residual variance and tag top variable features
  df <- df[order(df$residual_variance, decreasing = TRUE), ]
  df$variable <- "Non-variable"
  df$variable[1:nfeatures] <- "Variable"
  df$variable <- factor(df$variable, levels = c("Non-variable", "Variable"))
  
  # Split data
  df_var <- df[df$variable == "Variable", ]
  df_non <- df[df$variable == "Non-variable", ]
  
  # Calculate dynamic axis limits
  y_min <- min(df$residual_variance, na.rm = TRUE)
  y_max <- max(df$residual_variance, na.rm = TRUE)
  x_vals <- df$gmean[df$gmean > 0]  # exclude zero or negatives for log scale
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  
  # Build plot
  p <- ggplot(df, aes(x = gmean, y = residual_variance)) +
    geom_point(data = df_non, aes(color = variable), size = 0.5) +
    geom_point(data = df_var, aes(color = variable), size = 0.5) +
    scale_color_manual(values = c("Non-variable" = "black", "Variable" = "red"),
                       labels = c(paste0("Non-variable (", nrow(df_non), ")"),
                                  paste0("Variable (", nrow(df_var), ")"))) +
    scale_x_log10() +  
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    xlab("Geometric mean of expression") +
    ylab("Residual variance") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.title = element_blank(), legend.position = "top",
          legend.text = element_text(size = 15))
  
  # Add gene labels if desired
  if (label) {
    label_data <- df_var[1:min(10, nrow(df_var)), ]  # Cap at 10 labels
    if (repel) {
      p <- p + geom_text_repel(data = label_data,
                               aes(label = feature),
                               size = 3,
                               max.overlaps = Inf)
    } else {
      p <- p + geom_text(data = label_data,
                         aes(label = feature),
                         size = 3, vjust = 1)
    }
  }
  
  return(p)
}

# plot min cell range of genes by number of cells w/ expression
gene.min.cells.range <- function(list_of_dgCMatrices, 
                                 steps_vector = NULL, 
                                 min.features = 700) { 
  if(is.null(steps_vector)) {
    stop("input required for `steps_vector`")
  }
  stopifnot("`steps_vector` must be numeric" = is.numeric(steps_vector))
  
  if(!is.vector(steps_vector)) {
    steps_vector <- as.vector(steps_vector)
  }
  
  if (is.list(list_of_dgCMatrices) && !is.null(names(list_of_dgCMatrices))) {
    list_of_dgCMatrices <- list_of_dgCMatrices
  } else {
    list_of_dgCMatrices <- as_named_list(!!rlang::enquo(list_of_dgCMatrices))
  }
  
  steps.loop <- \(x) {
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
  
  lapply(list_of_dgCMatrices, steps.loop)
}

gene.min.cells.plot <- function(list_with_minCells, 
                                outdir = NULL) {
  
  if (is.list(list_with_minCells) && !is.null(names(list_with_minCells))) {
    list_with_minCells <- list_with_minCells
  } else {
    list_with_minCells <- as_named_list(!!rlang::enquo(list_with_minCells))
  }
  
  newlist = mapply("c",
                   list_with_minCells,
                   namelist(list_with_minCells),
                   SIMPLIFY=FALSE)
  
  generate.graphs <- \(x) {
    
    suffix = x[[length(x)]]
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
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
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/min.cells.", 
                         title, ".png")
    } else {
      
      plotname <- paste0(outdir, "/min.cells.",
                         title, ".png")
      
    }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 10, 
           bg = "white")
  }
  
  lapply(newlist, generate.graphs)
}

# create seurat objects from dcCMatrices
make.seurat.obj <- function(list_of_dgCMatrices) { 
  
  if (is.list(list_of_dgCMatrices) && !is.null(names(list_of_dgCMatrices))) {
    list_of_dgCMatrices <- list_of_dgCMatrices
  } else {
    list_of_dgCMatrices <- as_named_list(!!rlang::enquo(list_of_dgCMatrices))
  }
  
  if(!require(collapse)) {
    message("trying to install collapse")
    install.packages("collapse")
    
    if(!require(collapse)) {
      stop("install collapse package before running")
    }
  }
  
  newlist = mapply("c",
                   list_of_dgCMatrices,
                   namelist(list_of_dgCMatrices),
                   SIMPLIFY=FALSE)
  
  make.obj <- \(x) {
    
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
  
  lapply(newlist, make.obj)
}

# violin plots of features, counts, and percent mito genes (all)
plot.qc.metrics <- function(list_of_SeuratObj, 
                            outdir = NULL) {
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  generate.graphs <- \(x) {
    
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
    VlnPlot(x, assay = "RNA", 
            slot = "counts",
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
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/feature.count.mito.", 
                         title, ".png")
    } else {
      
      plotname <- paste0(outdir, "/feature.count.mito.",
                         title, ".png")
    }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 10, 
           bg = "white")
  }
  
  lapply(list_of_SeuratObj, generate.graphs)
}

# violin plots of features, counts, and percent mito genes by doublet status
plot.doublets.qc <- function(list_of_SeuratObj, 
                             outdir = NULL) {
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  if(!require(SeuratExtend)) {
    message("trying to install SeuratExtend")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("huayc09/SeuratExtend")
    
    if(!require(SeuratExtend)) {
      stop("install SeuratExtend package before running")
    }
  }
  
  generate.graphs <- \(x) {
    
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
    VlnPlot2(x,
             features = c("nFeature_RNA",
                          "nCount_RNA",
                          "percent.mt"),
             group.by = "doublet", 
             scales = "free_y", 
             ncol = 3, nrow = 3) +
      theme(legend.position = "none") +
      ggtitle(paste(title))
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/doublets.", 
                         title, ".png")
      
    } else {
      
      plotname <- paste0(outdir, "/doublets.",
                         title, ".png")
    }
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 5, 
           bg = "white")
  }
  
  lapply(list_of_SeuratObj, generate.graphs)
}

# nCount_RNA x nFeature_RNA with chosen filter thresholds
plot.feature.scatter <- function(list_of_SeuratObj, 
                                 filters, 
                                 outdir = NULL) {
  
  stopifnot("`filters` must be numeric" = is.numeric(filters))
  
  if(!is.matrix(filters)) {
    filters <- as.matrix(filters)
  }
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(collapse)) {
      stop("install gridExtra package before running")
    }
  }
  generate.graphs <- \(x) {
    
    my_filters <- filters[(rownames(filters) %in% Seurat::Project(x)),]
    
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
    scatter1 <- FeatureScatter(x,
                               feature1 = "nCount_RNA",
                               feature2 = "nFeature_RNA") +
      geom_hline(yintercept = my_filters["nFeature_min"],
                 linetype = "dotted") +
      geom_hline(yintercept = my_filters["nFeature_max"],
                 linetype = "dotted") +
      annotate("text", x=Inf, y=my_filters["nFeature_min"],
               label=paste(my_filters["nFeature_min"]),
               vjust = -0.25, hjust = 1) +
      annotate("text", x=Inf,
               y=my_filters["nFeature_max"],
               label=paste(my_filters["nFeature_max"]),
               vjust = -0.25, hjust = 1) +
      theme(legend.position = "none")
    
    scatter2 <- FeatureScatter(subset(x,
                                      subset = nFeature_RNA < 6000),
                               feature1 = "nCount_RNA",
                               feature2 = "nFeature_RNA") +
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
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/feature.scatter.", 
                         title, ".png")
      
    } else {
      
      plotname <- paste0(outdir, "/feature.scatter.",
                         title, ".png")
    }
    
    plots <- arrangeGrob(scatter1, scatter2, hist,
                         ncol = 2, widths = c(1,2),
                         layout_matrix = rbind(c(1, 3),
                                               c(2, 3)),
                         
                         top = grid::textGrob(paste(title),
                                              gp=grid::gpar(fontsize=24)))
    
    ggsave(plotname,
           plots,
           units = "in", 
           width = 12, 
           height = 7.5, 
           bg = "white")
  }
  
  lapply(list_of_SeuratObj, generate.graphs)
}

# nCount_RNA x nFeature_RNA and nFeature_RNA histo after applying filters
check.filters <- function(list_of_SeuratObj, 
                          filters, 
                          outdir = NULL) {
  
  stopifnot("`filters` must be numeric" = is.numeric(filters))
  
  if(!is.matrix(filters)) {
    filters <- as.matrix(filters)
  }
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(gridExtra)) {
      stop("install gridExtra package before running")
    }
  }
  
  generate.graphs <- \(x) {
    
    my_filters <- filters[(rownames(filters) %in% Seurat::Project(x)),]
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
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
    
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/filtered.feature.scatter.", 
                         title, ".png")
    } else {
      
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
  
  lapply(list_of_SeuratObj, generate.graphs)
}

# filter data in Seurat obj with specified thresholds
make.filtered.seurat <- function(list_of_SeuratObj, 
                                 filters) {
  
  stopifnot("`filters` must be numeric" = is.numeric(filters))
  
  if(!is.matrix(filters)) {
    filters <- as.matrix(filters)
  }
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  apply.filters <- \(x) {
    
    my_filters <- filters[(rownames(filters) %in% Seurat::Project(x)),]
    
    subset(x, subset =
             nFeature_RNA >= my_filters["nFeature_min"] &
             nFeature_RNA <= my_filters["nFeature_max"] &
             percent.mt < my_filters["mito.max"] & 
             doublet == 'singlet')
    
  }
  
  lapply(list_of_SeuratObj, apply.filters)
}

# test different clustering resolutions
find.cluster.range <- function(list_of_SeuratObj,
                               outdir = NULL,
                               int_range1 = seq(0,2,0.2),
                               int_range2 = seq(0,1.4,0.2)) {
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
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
  }
  
  generate.graphs <- \(x) {
    
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
    x.clustree <- Seurat::FindClusters(object = x,
                                       # graph.name = paste0(DefaultAssay(x), "_snn"),
                                       resolution = int_range1)
    
    ctx1 <- clustree(x.clustree,
                     prefix = paste0(DefaultAssay(x), "_snn_res.")) +
      labs(title = 'clusters across resolutions') +
      guides(edge_colour = "none", edge_alpha = "none") +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_edge_color_continuous(low = "black",
                                  high = "black")
    
    ctx1
    
    x.clust.reduced <- Seurat::FindClusters(object = x,
                                            resolution = int_range2)

    ctx2 <- clustree(x.clust.reduced,
                     prefix = paste0(DefaultAssay(x), "_snn_res.")) +
      labs(title = 'reduced resolution range') +
      guides(edge_colour = "none", edge_alpha = "none") +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_edge_color_continuous(low = "black",
                                  high = "black")
    
    ctx2

    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")

    } else {
      outdir = outdir
    }

    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {

      plotname <- paste0(outdir, "/", title,
                         "/clustree.resolution.",
                         title, ".png")

    } else {
      plotname <- paste0(outdir, "/clustree.resolution.",
                         title, ".png")
    }

    plots <- arrangeGrob(ctx1, ctx2,
                         ncol = 2, widths = c(1,1),

                         top = grid::textGrob(paste(title),
                                              gp=grid::gpar(fontsize=24)))

    message("**saving plots**")

    ggsave(plotname,
           plots,
           units = "in",
           width = 17,
           height = 6.5,
           bg = "white")

    return(x)
  }
  
  lapply(list_of_SeuratObj, generate.graphs)
}


# UMAPs with cluster assignment, nFeature_RNA count; elbow plot; variable features
plot.dim.clust <- function(list_of_SeuratObj,
                           outdir = NULL) {
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  if(!require(gridExtra)) {
    message("trying to install gridExtra")
    install.packages("gridExtra")
    
    if(!require(gridExtra)) {
      stop("install gridExtra package before running")
    }
  }
  
  generate.graphs <- \(x) {
    
    suffix <- Seurat::Project(x)
    
    if(substr(suffix, nchar(suffix) - 4, nchar(suffix)) == ".data") {
      title <- gsub('.{5}$', '', suffix)
      
    } else {
      title <- suffix
    }
    
    ebp <- ElbowPlot(x) + 
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    vfp <- VariableFeaturePlot_patched(x) + 
      theme(plot.margin = unit(c(0.25, 1, 1, 1), "cm"))
    dmp <- DimPlot(x, 
                   reduction = "umap",
                   label = TRUE,
                   repel = TRUE) + 
      labs(color = "cluster") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    ftp <- FeaturePlot(x, 
                       reduction = "umap",
                       features = "nFeature_RNA",
                       label = TRUE,
                       repel = TRUE) +
      labs(title = NULL, color = "nFeature_RNA") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    plots <- arrangeGrob(ebp,dmp,vfp,ftp,
                         widths = c(0.75, 1),
                         ncol = 2, nrow = 2,
                         
                         top = grid::textGrob(paste(title), 
                                              gp=grid::gpar(fontsize=24)))
    
    
    if(is.null(outdir)) {
      outdir = paste(getwd())
      message("outdir not specified; saving plots to working directory")
      
    } else {
      outdir = outdir
    }
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/elbow.dim.varfeat.", 
                         title, ".png")
      
    } else {
      
      plotname <- paste0(outdir, "/elbow.dim.varfeat.",
                         title, ".png")
    }
    
    ggsave(plotname,
           plots,
           units = "in", 
           width = 14, 
           height = 12, 
           bg = "white")
    
    }
  
  lapply(list_of_SeuratObj, generate.graphs)
}

# use ScType annotation to assign cell types (save to 'name' in metadata)
annotate.with.sctype <- function(list_of_SeuratObj,
                                 name = "sctype.ind",
                                 outdir = NULL) {
  
  if (is.list(list_of_SeuratObj) && !is.null(names(list_of_SeuratObj))) {
    list_of_SeuratObj <- list_of_SeuratObj
  } else {
    list_of_SeuratObj <- as_named_list(!!rlang::enquo(list_of_SeuratObj))
  }
  
  mod.sctype.wrapper <- \(x) {
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
    
    if(length(grep(paste0("/",title,"$"), list.dirs(outdir))) > 0) {
      
      plotname <- paste0(outdir, "/", title, 
                         "/sctype.dimplot.", 
                         title, ".png")
      
    } else {
      
      plotname <- paste0(outdir, "/sctype.dimplot.",
                         title, ".png")
    
    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 8, 
           bg = "white")
    
    return(x)
    }
  }
  
  lapply(list_of_SeuratObj, mod.sctype.wrapper)
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

### from HypoMap_snseq_mouse_IMC.R ####
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

