#### accessory functions ####
# package version and dependency log
load_packages <- function(pkgs, 
                          out_prefix = "package_log", 
                          auto_install = TRUE,
                          folder = getwd()) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- paste0(folder, "/package_logs")
  if (!dir.exists(out_dir)) dir.create(out_dir)
  report_file <- file.path(out_dir, paste0(out_prefix, "_package_report_", timestamp, ".txt"))
  
  # install missing from CRAN/Bioconductor
  install_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing missing package: ", pkg)
      tryCatch({
        install.packages(pkg)
        if (!requireNamespace(pkg, quietly = TRUE)) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }
      }, error = function(e) {
        message("Could not install ", pkg, ": ", e$message)
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
    cat("\n ", pkg, "\n")          # console
    writeLines(paste0("\n ", pkg), con)  # file
    
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      cat("Loaded successfully\n")        # console
      writeLines("✔ Loaded successfully", con) # file
    }, error = function(e) {
      cat("Failed to load: ", e$message, "\n")   # console
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
  
  message("Combined package log saved to: ", report_file)
}

# subset Seurat obj by gene list
subset_by_gene <- function(seurat_obj, 
                           subset_genes, 
                           slot = "counts",
                           min_count = 2) {
  
  gene_list <- lapply(subset_genes, function(gene) {
    
    counts <- GetAssayData(seurat_obj, slot = slot)[gene, ]
    keep_cells <- names(counts[counts >= min_count])
    
    cat(paste0(gene, ": ", length(keep_cells), " cells retained\n"))
    
    # Subset object
    obj <- subset(seurat_obj, cells = keep_cells)
    obj@project.name <- gene
    obj 
    
  } )
  
  names(gene_list) <- subset_genes
  return(gene_list)
}

#### network analysis and graphing functions ####

# (old) function for limma_trend
run_limmatrend_neuron_cluster <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatmale.dom.data + treatfemale.dom.data)/2 - (treatmale.sub.data + treatfemale.sub.data)/2, 
                               DvsSM = treatmale.dom.data - treatmale.sub.data, 
                               DvsSF = treatfemale.dom.data - treatfemale.sub.data,
                               MvsF = (treatmale.dom.data + treatmale.sub.data)/2 - (treatfemale.dom.data + treatfemale.sub.data)/2,
                               MvsFD = treatmale.dom.data - treatfemale.dom.data,
                               MvsFS = treatmale.sub.data - treatfemale.sub.data,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    
    fit <- lmFit(y, design = design)
    fit <- contrasts.fit(fit, contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    ttDvsS <- topTable(fit, 
                       n = Inf,
                       coef = "DvsS",
                       adjust.method = "BH")
    ttDvsSM <- topTable(fit, 
                        n = Inf,
                        coef = "DvsSM",
                        adjust.method = "BH")
    ttDvsSF <- topTable(fit, 
                        n = Inf,
                        coef = "DvsSF",
                        adjust.method = "BH")
    ttMvsF <- topTable(fit, 
                       n = Inf,
                       coef = "MvsF",
                       adjust.method = "BH")
    ttMvsFD <- topTable(fit, 
                        n = Inf,
                        coef = "MvsFD",
                        adjust.method = "BH")
    ttMvsFS <- topTable(fit, 
                        n = Inf,
                        coef = "MvsFS",
                        adjust.method = "BH")
  })
  
  # Open pdf file
  pdf(file= paste0('./neurons/stats/limma_trend/cluster_',
                   j, '/', j,
                   '_limmatrend.histograms.pdf'))
  
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(ttDvsS$P.Value, 50)
  hist(ttDvsS$adj.P.Val, 50)
  hist(ttDvsSM$P.Value, 50)
  hist(ttDvsSM$adj.P.Val, 50)
  hist(ttDvsSF$P.Value, 50)
  hist(ttDvsSF$adj.P.Val, 50)
  hist(ttMvsF$P.Value, 50)
  hist(ttMvsF$adj.P.Val, 50)
  hist(ttMvsFD$P.Value, 50)
  hist(ttMvsFD$adj.P.Val, 50)
  hist(ttMvsFS$P.Value, 50)
  hist(ttMvsFS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  # pdf(file= "./neuron_cluster/neuropeptides/limma_trend/limmatrend.MDS.pdf" )
  # # create a 2X1 grid
  # par( mfrow= c(1,1) )
  # limma::plotMDS(dge,
  #                top = 20,
  #                col = as.numeric(as.factor(L$condt)),
  #                pch = 19)
  # plotMD(fit)
  # dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF,
       ttMvsFD = ttMvsFD,
       ttMvsFS = ttMvsFS)
}

# function to get dotplot data
DotPlot.data = function (object, assay = NULL, 
                         features, 
                         cols = c("lightgrey", "blue"), 
                         col.min = -2.5, col.max = 2.5, 
                         dot.min = 0, dot.scale = 6, 
                         idents = NULL, 
                         group.by = NULL, 
                         split.by = NULL, 
                         cluster.idents = FALSE, 
                         scale = TRUE, scale.by = "radius", 
                         scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  return(data.plot)
}

# function for bigger network graphs
ModuleUMAPPlot.size = function (seurat_obj, sample_edges = TRUE, edge_prop = 0.2, 
                                label_hubs = 5, edge.alpha = 0.25, vertex.label.cex = 0.5, 
                                label_genes = NULL, return_graph = FALSE, keep_grey_edges = TRUE, dot.size = 3, edge.size = 0.5,
                                wgcna_name = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  TOM <- GetTOM(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$module)
  mods <- mods[mods != "grey"]
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 
                                                 "hub"]]
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(label_hubs) %>% 
      .$gene_name
  })
  names(hub_list) <- mods
  hub_labels <- as.character(unlist(hub_list))
  print("hub labels")
  print(hub_labels)
  print(label_genes)
  if (is.null(label_genes)) {
    label_genes <- hub_labels
  }
  else {
    if (!any(label_genes %in% umap_df$gene)) {
      stop("Some genes in label_genes not found in the UMAP.")
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  print(label_genes)
  selected_modules <- modules[umap_df$gene, ]
  selected_modules <- cbind(selected_modules, umap_df[, c("UMAP1", 
                                                          "UMAP2", "hub", "kME")])
  selected_modules$label <- ifelse(selected_modules$gene_name %in% 
                                     label_genes, selected_modules$gene_name, "")
  selected_modules$fontcolor <- ifelse(selected_modules$color == 
                                         "black", "gray50", "black")
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% 
                                          label_genes, "black", selected_modules$color)
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), 
                                               function(i) {
                                                 gene1 = as.character(edge_df[i, "Var1"])
                                                 gene2 = as.character(edge_df[i, "Var2"])
                                                 col1 <- selected_modules[selected_modules$gene_name == 
                                                                            gene1, "color"]
                                                 col2 <- selected_modules[selected_modules$gene_name == 
                                                                            gene2, "color"]
                                                 if (col1 == col2) {
                                                   col = col1
                                                 }
                                                 else {
                                                   col = "grey90"
                                                 }
                                                 col
                                               })
  if (!keep_grey_edges) {
    edge_df <- edge_df %>% subset(color != "grey90")
  }
  groups <- unique(edge_df$color)
  if (sample_edges) {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * 
                                              edge_prop))
      cur_df[cur_sample, ]
    }))
  }
  else {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), 
                              wt = value)
    }))
  }
  edge_df <- temp
  print(dim(edge_df))
  edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(subset(edge_df, color == "grey90"), subset(edge_df, 
                                                              color != "grey90"))
  edge_df$color_alpha <- ifelse(edge_df$color == "grey90", 
                                alpha(edge_df$color, alpha = edge_df$value/2), alpha(edge_df$color, 
                                                                                     alpha = edge_df$value))
  selected_modules <- rbind(subset(selected_modules, hub == 
                                     "other"), subset(selected_modules, hub != "other"))
  selected_modules <- rbind(subset(selected_modules, label == 
                                     ""), subset(selected_modules, label != ""))
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, 
                                     vertices = selected_modules)
  print("making net")
  print(head(edge_df))
  print(head(selected_modules))
  if (return_graph) {
    return(g)
  }
  plot(g, layout = as.matrix(selected_modules[, c("UMAP1", 
                                                  "UMAP2")]), edge.color = adjustcolor(E(g)$color_alpha, 
                                                                                       alpha.f = edge.alpha), vertex.size = V(g)$kME * dot.size, edge.curved = 0, 
       edge.width = edge.size, vertex.color = V(g)$color, vertex.label = V(g)$label, 
       vertex.label.dist = 1.1, vertex.label.degree = -pi/4, 
       vertex.label.family = "Helvetica", vertex.label.font = 3, 
       vertex.label.color = V(g)$fontcolor, vertex.label.cex = 0, 
       vertex.frame.color = V(g)$framecolor, margin = 0)
}


#### From Neuropeptide_network_snseq_mouse_IMC.R ####

### create remove_nonexp
## from: https://almeidasilvaf.github.io/BioNERO/reference/remove_nonexp.html 
#' Remove genes that are not expressed based on a user-defined threshold
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param method Criterion to filter non-expressed genes out.
#' One of "mean", "median", "percentage", or "allsamples". Default is "median".
#' @param min_exp If method is 'mean', 'median', or 'allsamples',
#' the minimum value for a gene to be considered expressed.
#' If method is 'percentage', the minimum value each gene must have in
#' at least n percent of samples to be considered expressed.
#' @param min_percentage_samples In case the user chooses 'percentage' as method,
#' expressed genes must have expression >= min_exp in at least this percentage.
#' Values must range from 0 to 1.
#'
#' @return Filtered gene expression data frame or `SummarizedExperiment` object.
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom matrixStats rowMedians
#' @seealso
#'  \code{\link[matrixStats]{rowMedians}}
#'  \code{\link[WGCNA]{goodSamplesGenes}}
#' @rdname remove_nonexp
#' @examples
#' data(zma.se)
#' filt_exp <- remove_nonexp(zma.se, min_exp = 5)
remove_nonexp <- function(exp, method="median", min_exp=1, min_percentage_samples=0.25) {
  # fexp <- handleSE(exp)
  fexp <- exp
  
  if(method == "median") {
    final_exp <- fexp[matrixStats::rowMedians(as.matrix(fexp)) >= min_exp,]
  } else if (method == "mean") {
    final_exp <- fexp[rowMeans(fexp) >= min_exp,]
  } else if (method == "percentage") {
    min_n <- ncol(fexp) * min_percentage_samples
    final_exp <- fexp[rowSums(fexp >= min_exp) >= min_n, ]
  } else if (method == "allsamples") {
    final_exp <- fexp[rowSums(fexp >= min_exp) == ncol(fexp), ]
  } else {
    stop("No method specified. Please, choose a filtering method - mean, median or percentage")
  }
  
  # if(is(exp, "SummarizedExperiment")) {
  #   final_exp <- exp2SE(final_exp, exp)
  # }
  
  return(final_exp)
}


#### add pvalue to mantel test figure
#' Mantel test for a set of correlation matrices
#' @description
#' `r badge('stable')`
#'
#' This function generate a pairwise matrix of plots to compare the similarity
#' of two or more correlation matrices. In the upper diagonal are presented the
#' plots and in the lower diagonal the result of Mantel test based on
#' permutations.
#'
#'
#' @param ... The input matrices. May be an output generated by the function
#'   `lpcor` or a coerced list generated by the function `as.lpcor`
#' @param type The type of correlation if an obect generated by the function
#'   `lpcor` is used. 1 = Linear correlation matrices, or 2 = partial
#'   correlation matrices.
#' @param nrepet The number of permutations. Default is 1000
#' @param names An optional vector of names of the same length of `...` .
#' @param prob The error probability for Mantel test.
#' @param diag Logical argument. If `TRUE`, the Kernel density is shown in
#'   the diagonal of plot.
#' @param export Logical argument. If `TRUE`, then the plot is exported to
#'   the current directory.
#' @param main The title of the plot, set to 'auto'.
#' @param file.type The format of the file if `export = TRUE`.  Set to
#'   `'pdf'`. Other possible values are `*.tiff` using `file.type
#'   = 'tiff'`.
#' @param file.name The name of the plot when exported. Set to `NULL`,
#'   i.e., automatically.
#' @param width The width of the plot, set to `8`.
#' @param height The height of the plot, set to `7`.
#' @param resolution The resolution of the plot if `file.type = 'tiff'` is
#'   used. Set to `300` (300 dpi).
#' @param size.point The size of the points in the plot. Set to `0.5`.
#' @param shape.point The shape of the point, set to ` 19`.
#' @param alpha.point The value for transparency of the points: 1 = full color.
#' @param fill.point The color to fill the points. Valid argument if points are
#'   between 21 and 25.
#' @param col.point The color for the edge of the point, set to `black`.
#' @param minsize The size of the letter that will represent the smallest
#'   correlation coefficient.
#' @param maxsize The size of the letter that will represent the largest
#'   correlation coefficient.
#' @param signcol The colour that indicate significant correlations (based on
#'   the `prob` value.), set to 'green'.
#' @param alpha The value for transparency of the color informed in
#'   `signcol`, when 1 = full color. Set to 0.15.
#' @param diagcol The color in the kernel distribution. Set to 'gray'.
#' @param col.up.panel,col.lw.panel,col.dia.panel The color for the opper, lower
#'   and diagonal pannels. Set to 'gray', 'gray', and 'gray', respectively.
#' @param pan.spacing The space between the pannels. Set to 0.15.
#' @seealso [mantel_test()]
#' @param digits The number of digits to show in the plot.
#' @return An object of class `gg, ggmatrix`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # iris dataset
#' lpc <- iris %>%
#'        group_by(Species) %>%
#'        lpcor() %>%
#'        pairs_mantel(names = c('setosa', 'versicolor', 'virginica'))
#'
#'
#' # mtcars dataset
#' mt_num <- select_numeric_cols(mtcars)
#' lpdata <- as.lpcor(cor(mt_num[1:5]),
#'                    cor(mt_num[1:5]),
#'                    cor(mt_num[2:6]),
#'                    cor(mt_num[4:8])) %>%
#'           pairs_mantel()
#'}
pairs_mantel.pvalue <- function (..., type = 1, nrepet = 1000, names = NULL, prob = 0.05, 
                                 diag = FALSE, export = FALSE, main = "auto", file.type = "pdf", 
                                 file.name = NULL, width = 8, height = 7, resolution = 300, 
                                 size.point = 0.5, shape.point = 19, alpha.point = 1, fill.point = NULL, 
                                 col.point = "black", minsize = 2, maxsize = 3, signcol = "green", 
                                 alpha = 0.15, diagcol = "gray", col.up.panel = "gray", col.lw.panel = "gray", 
                                 col.dia.panel = "gray", pan.spacing = 0.15, digits = 2) 
{
  class <- list(...)
  if (!type %in% c(1:2)) {
    stop("The argument type must be 1 (linear correlation) or 2 (partial correlation).")
  }
  if (sum(lapply(class, function(x) !any(class(x) %in% c("lpcor_group", 
                                                         "lpcor", "mahala_group", "covcor_design", "group_clustering", 
                                                         "clustering") == TRUE)) > 0)) {
    stop("The object must be of the class lpcor. Please use 'as.lpcorr' to convert correlation matrices into the correct format.")
  }
  if (any(class(...) == "lpcor_group")) {
    data <- lapply(...[[2]], function(x) {
      x[["linear.mat"]]
    })
  }
  if (any(class(...) == "group_clustering")) {
    data <- lapply(...[[2]], function(x) {
      x$distance
    })
  }
  if (!any(class(...) %in% c("lpcor_group", "group_clustering"))) {
    data <- lapply(..., function(x) {
      x
    })
  }
  w <- c(21:25)
  if (is.null(fill.point) == TRUE && any(w == shape.point)) {
    stop(call. = FALSE, "If 'shape.point' is a value between 21 and 25, you must provide a color for fill the shape using the argument 'fill.point.'")
  }
  for (i in 1:length(data)) {
    if (i == 1) {
      Dataset <- data.frame(var = as.vector(t(data[[1]])[lower.tri(data[[1]], 
                                                                   diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix 1")
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[1]
      }
    }
    if (i >= 2) {
      Dataset <- mutate(Dataset, var = as.vector(t(data[[i]])[lower.tri(data[[i]], 
                                                                        diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix ", 
                                                                    i)
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[i]
      }
    }
  }
  dim <- nrow(data[[1]])
  my_custom_cor <- function(data, mapping, color = I("black"), 
                            sizeRange = c(minsize, maxsize), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    sig <- symnum(ct[[3]], corr = FALSE, na = FALSE, cutpoints = c(0, 
                                                                   0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
                                                                                                      "*", ""))
    r <- ct[[1]]
    rt <- paste(format(r, digits = 2)[1],
                format(ct[[3]], digits = 2)[1],
                sep = ' , ')
    cex <- max(sizeRange)
    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(label = as.character(rt), mapping = aes(), 
                        xP = 0.5, yP = 0.5, size = I(percent_of_range(cex * 
                                                                        abs(r), sizeRange)), color = color, ...) + geom_text(aes_string(x = 0.8, 
                                                                                                                                        y = 0.8), label = sig, size = I(cex), color = color, 
                                                                                                                             ...) + theme_classic() + theme(panel.background = element_rect(color = col.lw.panel), 
                                                                                                                                                            axis.line = element_blank(), axis.ticks = element_blank(), 
                                                                                                                                                            axis.text.y = element_blank(), axis.text.x = element_blank())
  }
  my_custom_smooth <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    pval <- ct[[3]]
    p <- ggplot(data = data, mapping = mapping)
    if (is.null(fill.point) == FALSE) {
      p <- p + geom_point(color = I(col.point), fill = fill.point, 
                          shape = shape.point, size = size.point, alpha = alpha.point)
    }
    else {
      p <- p + geom_point(color = I(col.point), shape = shape.point, 
                          size = size.point, alpha = alpha.point)
    }
    p <- p + theme_classic() + theme(panel.background = element_rect(fill = "white", 
                                                                     color = col.up.panel), axis.line = element_blank(), 
                                     axis.ticks = element_blank(), axis.text.y = element_blank(), 
                                     axis.text.x = element_blank())
    if (pval < prob) {
      p <- p + theme(panel.background = element_rect(fill = ggplot2::alpha(signcol, 
                                                                           alpha)))
    }
    p
  }
  ggally_mysmooth <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) + geom_density(fill = alpha(diagcol, 
                                                                       1)) + theme_classic() + theme(panel.background = element_rect(fill = alpha("white", 
                                                                                                                                                  1), color = col.dia.panel))
  }
  if (main == "auto") {
    title <- paste0("Mantel's test with ", nrepet, " resamples")
  }
  else {
    title <- main
  }
  if (diag == TRUE) {
    diag <- list(continuous = ggally_mysmooth)
  }
  else {
    diag <- NULL
  }
  p1 <- GGally::ggpairs(Dataset, title = title, diag = diag, 
                        lower = list(continuous = my_custom_cor), upper = list(continuous = my_custom_smooth), 
                        axisLabels = "none") + theme(panel.spacing = grid::unit(pan.spacing, 
                                                                                "lines"))
  if (export == FALSE) {
    return(p1)
  }
  else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("Pairs of Mantel's test.pdf", width = width, 
          height = height)
    }
    else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    print(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "Pairs of Mantel's test.tiff", width = width, 
           height = height, units = "in", compression = "lzw", 
           res = resolution)
    }
    else tiff(filename = paste0(file.name, ".tiff"), width = width, 
              height = height, units = "in", compression = "lzw", 
              res = resolution)
    print(p1)
    dev.off()
  }
}

#### add confidence interval to mantel test figure
#' Mantel test for a set of correlation matrices
#' @description
#' `r badge('stable')`
#'
#' This function generate a pairwise matrix of plots to compare the similarity
#' of two or more correlation matrices. In the upper diagonal are presented the
#' plots and in the lower diagonal the result of Mantel test based on
#' permutations.
#'
#'
#' @param ... The input matrices. May be an output generated by the function
#'   `lpcor` or a coerced list generated by the function `as.lpcor`
#' @param type The type of correlation if an obect generated by the function
#'   `lpcor` is used. 1 = Linear correlation matrices, or 2 = partial
#'   correlation matrices.
#' @param nrepet The number of permutations. Default is 1000
#' @param names An optional vector of names of the same length of `...` .
#' @param prob The error probability for Mantel test.
#' @param diag Logical argument. If `TRUE`, the Kernel density is shown in
#'   the diagonal of plot.
#' @param export Logical argument. If `TRUE`, then the plot is exported to
#'   the current directory.
#' @param main The title of the plot, set to 'auto'.
#' @param file.type The format of the file if `export = TRUE`.  Set to
#'   `'pdf'`. Other possible values are `*.tiff` using `file.type
#'   = 'tiff'`.
#' @param file.name The name of the plot when exported. Set to `NULL`,
#'   i.e., automatically.
#' @param width The width of the plot, set to `8`.
#' @param height The height of the plot, set to `7`.
#' @param resolution The resolution of the plot if `file.type = 'tiff'` is
#'   used. Set to `300` (300 dpi).
#' @param size.point The size of the points in the plot. Set to `0.5`.
#' @param shape.point The shape of the point, set to ` 19`.
#' @param alpha.point The value for transparency of the points: 1 = full color.
#' @param fill.point The color to fill the points. Valid argument if points are
#'   between 21 and 25.
#' @param col.point The color for the edge of the point, set to `black`.
#' @param minsize The size of the letter that will represent the smallest
#'   correlation coefficient.
#' @param maxsize The size of the letter that will represent the largest
#'   correlation coefficient.
#' @param signcol The colour that indicate significant correlations (based on
#'   the `prob` value.), set to 'green'.
#' @param alpha The value for transparency of the color informed in
#'   `signcol`, when 1 = full color. Set to 0.15.
#' @param diagcol The color in the kernel distribution. Set to 'gray'.
#' @param col.up.panel,col.lw.panel,col.dia.panel The color for the opper, lower
#'   and diagonal pannels. Set to 'gray', 'gray', and 'gray', respectively.
#' @param pan.spacing The space between the pannels. Set to 0.15.
#' @seealso [mantel_test()]
#' @param digits The number of digits to show in the plot.
#' @return An object of class `gg, ggmatrix`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # iris dataset
#' lpc <- iris %>%
#'        group_by(Species) %>%
#'        lpcor() %>%
#'        pairs_mantel(names = c('setosa', 'versicolor', 'virginica'))
#'
#'
#' # mtcars dataset
#' mt_num <- select_numeric_cols(mtcars)
#' lpdata <- as.lpcor(cor(mt_num[1:5]),
#'                    cor(mt_num[1:5]),
#'                    cor(mt_num[2:6]),
#'                    cor(mt_num[4:8])) %>%
#'           pairs_mantel()
#'}
pairs_mantel.ci <- function (..., type = 1, nrepet = 1000, names = NULL, prob = 0.05, 
                             diag = FALSE, export = FALSE, main = "auto", file.type = "pdf", 
                             file.name = NULL, width = 8, height = 7, resolution = 300, 
                             size.point = 0.5, shape.point = 19, alpha.point = 1, fill.point = NULL, 
                             col.point = "black", minsize = 2, maxsize = 3, signcol = "green", 
                             alpha = 0.15, diagcol = "gray", col.up.panel = "gray", col.lw.panel = "gray", 
                             col.dia.panel = "gray", pan.spacing = 0.15, digits = 2) 
{
  class <- list(...)
  if (!type %in% c(1:2)) {
    stop("The argument type must be 1 (linear correlation) or 2 (partial correlation).")
  }
  if (sum(lapply(class, function(x) !any(class(x) %in% c("lpcor_group", 
                                                         "lpcor", "mahala_group", "covcor_design", "group_clustering", 
                                                         "clustering") == TRUE)) > 0)) {
    stop("The object must be of the class lpcor. Please use 'as.lpcorr' to convert correlation matrices into the correct format.")
  }
  if (any(class(...) == "lpcor_group")) {
    data <- lapply(...[[2]], function(x) {
      x[["linear.mat"]]
    })
  }
  if (any(class(...) == "group_clustering")) {
    data <- lapply(...[[2]], function(x) {
      x$distance
    })
  }
  if (!any(class(...) %in% c("lpcor_group", "group_clustering"))) {
    data <- lapply(..., function(x) {
      x
    })
  }
  w <- c(21:25)
  if (is.null(fill.point) == TRUE && any(w == shape.point)) {
    stop(call. = FALSE, "If 'shape.point' is a value between 21 and 25, you must provide a color for fill the shape using the argument 'fill.point.'")
  }
  for (i in 1:length(data)) {
    if (i == 1) {
      Dataset <- data.frame(var = as.vector(t(data[[1]])[lower.tri(data[[1]], 
                                                                   diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix 1")
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[1]
      }
    }
    if (i >= 2) {
      Dataset <- mutate(Dataset, var = as.vector(t(data[[i]])[lower.tri(data[[i]], 
                                                                        diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix ", 
                                                                    i)
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[i]
      }
    }
  }
  dim <- nrow(data[[1]])
  my_custom_cor <- function(data, mapping, color = I("black"), 
                            sizeRange = c(minsize, maxsize), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    ct.2 <- ecodist::mantel(as.dist(D) ~ as.dist(D2), nboot = nrepet) # use mantel test from ecodist to get confidence interval
    sig <- symnum(ct[[2]], #change to get p-value
                  corr = FALSE, na = FALSE, cutpoints = c(0, 
                                                          0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
                                                                                             "*", ""))
    r <- ct[[1]]
    rt <- paste(format(ct.2[[2]], digits = 2)[1],
                " , ",
                format(ct.2[[1]], digits = 2)[1],
                ' , ',
                format(ct.2[[5]], digits = 2)[1], # add lower confidence interval
                '-',
                format(ct.2[[6]], digits = 2)[1], # add upper confidence interval
                sep = '')
    cex <- max(sizeRange)
    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(label = as.character(rt), mapping = aes(), 
                        xP = 0.5, yP = 0.5, size = I(percent_of_range(cex * 
                                                                        abs(r), sizeRange)), color = color, ...) + geom_text(aes_string(x = 0.8, 
                                                                                                                                        y = 0.8), label = sig, size = I(cex), color = color, 
                                                                                                                             ...) + theme_classic() + theme(panel.background = element_rect(color = col.lw.panel), 
                                                                                                                                                            axis.line = element_blank(), axis.ticks = element_blank(), 
                                                                                                                                                            axis.text.y = element_blank(), axis.text.x = element_blank())
  }
  my_custom_smooth <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    pval <- ct[[3]]  
    p <- ggplot(data = data, mapping = mapping)
    if (is.null(fill.point) == FALSE) {
      p <- p + geom_point(color = I(col.point), fill = fill.point, 
                          shape = shape.point, size = size.point, alpha = alpha.point)
    }
    else {
      p <- p + geom_point(color = I(col.point), shape = shape.point, 
                          size = size.point, alpha = alpha.point)
    }
    p <- p + theme_classic() + theme(panel.background = element_rect(fill = "white", 
                                                                     color = col.up.panel), axis.line = element_blank(), 
                                     axis.ticks = element_blank(), axis.text.y = element_blank(), 
                                     axis.text.x = element_blank())
    if (pval < prob) {
      p <- p + theme(panel.background = element_rect(fill = ggplot2::alpha(signcol, 
                                                                           alpha)))
    }
    p
  }
  ggally_mysmooth <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) + geom_density(fill = alpha(diagcol, 
                                                                       1)) + theme_classic() + theme(panel.background = element_rect(fill = alpha("white", 
                                                                                                                                                  1), color = col.dia.panel))
  }
  if (main == "auto") {
    title <- paste0("Mantel's test with ", nrepet, " resamples")
  }
  else {
    title <- main
  }
  if (diag == TRUE) {
    diag <- list(continuous = ggally_mysmooth)
  }
  else {
    diag <- NULL
  }
  p1 <- GGally::ggpairs(Dataset, title = title, diag = diag, 
                        lower = list(continuous = my_custom_cor), upper = list(continuous = my_custom_smooth), 
                        axisLabels = "none") + theme(panel.spacing = grid::unit(pan.spacing, 
                                                                                "lines"))
  if (export == FALSE) {
    return(p1)
  }
  else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("Pairs of Mantel's test.pdf", width = width, 
          height = height)
    }
    else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    print(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "Pairs of Mantel's test.tiff", width = width, 
           height = height, units = "in", compression = "lzw", 
           res = resolution)
    }
    else tiff(filename = paste0(file.name, ".tiff"), width = width, 
              height = height, units = "in", compression = "lzw", 
              res = resolution)
    print(p1)
    dev.off()
  }
}

#' @name testCov
#' @aliases testCov
#' @title Testing the equality of two sample covariance matrices in high dimension.
#' @description Testing the equality of two sample covariance matrices in high dimension using different methods.
#'
#' @param X the n x p training data, could be a \code{matrix} or a \code{data.frame} object.
#' @param Y the n x p training data matrix, could be a \code{matrix} or a \code{data.frame} object.
#' @param method a string incidating the method for the test. The current available
#'        methods are \code{ALL}, \code{HD}, \code{LC}, \code{CLX}, \code{Scott}.
#' @param J the number of repetition in the test
#' @param alpha the significant level of the test.
#' @param n.core the number of cores to be used in parallel when \code{HD}
#'        is called.
#'
#'
#' @return
#'
#' For any single method, the function returns an \code{htest} object.
#'
#' For method \code{ALL}: A list of four \code{htest} objects.
#'
#' HD refers to "Chang, J., Zhou, W., Zhou, W.-X., and Wang, L. (2016). Comparing large covariance
#' matrices under weak conditions on the dependence structure and its application to gene
#'  clustering. Biometrics. To appear"#'
#'
#' CLX refers to "Cai, T. T., Liu, W., and Xia, Y. (2013).
#' Two-sample covariance matrix testing and support recovery in high-dimensional
#' and sparse settings. Journal of the American Statistical Association 108, 265-277."
#'
#' Sc refers to "Schott, J. R. (2007). A test for the equality of covariance
#' matrices when the dimension is large relative to the sample size.
#' Computational Statistics and Data Analysis 51, 6535-6542."
#'
#' @examples
#' data(GO54)
#' testCov(GO54$X, GO54$Y, method = "ALL", J = 100)
#' data(GO26)
#' testCov(GO26$X, GO26$Y, method = "ALL", J = 100)
#'
#' @author Tong He
#' @export
#'
#'
testCov = function(X, Y, method = "ALL", J = 2500, alpha = 0.05, n.core = 1) {
  DNAME = deparse(substitute(X))
  DNAME = paste(DNAME, "and", deparse(substitute(Y)))
  X = data.matrix(X)
  Y = data.matrix(Y)
  if (method == "HD") {
    res = testCovHD(X = X, Y = Y, J = J, alpha = alpha, DNAME = DNAME,
                    n.core = n.core)
    res = res[[1]]
  } else if (method == "LC") {
    res = equalCovs(X = X, Y = Y, DNAME = DNAME)
  } else if (method == "CLX") {
    res = testCovHD(X = X, Y = Y, J = 1, alpha = alpha, DNAME = DNAME,
                    n.core = 1)
    res = res[[2]]
  } else if (method == "Scott") {
    message("The method \'Scott\' does not support skewed data.")
    res = testCovHD(X = X, Y = Y, J = 1, alpha = alpha, DNAME = DNAME,
                    n.core = 1)
    res = res[[3]]
  } else if (method == "ALL") {
    res = testCovHD(X = X, Y = Y, J = J, alpha = alpha, DNAME = DNAME,
                    n.core = n.core)
    # nms = colnames(res)
    res2 = equalCovs(X = X, Y = Y, alpha = alpha, DNAME = DNAME)
    res = c(res, list(res2))
    # colnames(res) = c(nms, 'LC')
    message("The method \'Scott\' does not support skewed data.")
  } else {
    msg = paste0("method \'", method, "\' not available.")
    stop(msg)
  }
  return(res)
}


testCovHD = function(X, Y, J = 2500, alpha = 0.05, DNAME, n.core = 1) {
  checkmate::checkMatrix(X)
  checkmate::checkMatrix(Y)
  
  p = ncol(X)
  n1 = nrow(X)
  n2 = nrow(Y)
  if (ncol(Y) != p) {
    stop('Different dimensions of X and Y.')
  }
  
  tmp = rep(c(rep(1, n1)/n1, rep(1, n2)/n2), J)
  scalev = matrix(tmp, ncol = 1)
  
  qalpha = -log(8*pi) - 2*log(log(1/(1-alpha)))
  cri = 4*log (p)-log (log (p)) + qalpha
  
  
  Sx = cov(X)*(n1-1)/n1
  Sy = cov(Y)*(n2-1)/n2
  
  xa = t(t(X) - colMeans(X))
  ya = t(t(Y) - colMeans(Y))
  
  vx = t(xa^2)%*%(xa^2)/n1 - 2/n1 * (t(xa)%*% xa) * Sx + Sx^2
  vy = t(ya^2)%*%(ya^2)/n2 - 2/n2 * (t(ya)%*% ya) * Sy + Sy^2
  
  numo = abs(Sx-Sy)
  deno = sqrt(vx/n1 + vy/n2)
  Tnm = max(numo/deno)
  
  xat = t(xa)/n1
  yat = t(ya)/n2
  # g = rnorm((n1+n2)*J)*scalev
  
  # ts = matrix(0,J,1)
  scalev = c(rep(1, n1)/n1, rep(1, n2)/n2)
  if (n.core > 1) {
    # for (j in 1:J) {
    cl = makeCluster(n.core)
    registerDoParallel(cl)
    ts = foreach(j = 1:J, .combine = rbind) %dopar% {
      # ind1 = ((j-1)*(n1+n2)+1):((j-1)*(n1+n2)+n1)
      # ind2 = ((j-1)*(n1+n2)+n1+1):((j-1)*(n1+n2)+n2+n1)
      
      g = rnorm(n1+n2)*scalev
      atmp = sum(g[1:n1])
      btmp = sum(g[(n1+1):(n1+n2)])
      
      ts1 = (t(xa*g[1:n1]) - xat*atmp)%*% xa
      ts2 = (t(ya*g[(n1+1):(n1+n2)]) - yat*btmp)%*% ya
      
      # atmp = sum(g[ind1])
      # btmp = sum(g[ind2])
      
      # ts1 = (t(xa*g[ind1]) - xat*atmp)%*% xa
      # ts2 = (t(ya*g[ind2]) - yat*btmp)%*% ya
      #ts[j] = max(abs(ts1-ts2)/deno)
      max(abs(ts1-ts2)/deno)
    }
    stopCluster(cl)
  } else {
    ts = matrix(0,J,1)
    for (j in 1:J) {
      g = rnorm(n1+n2)*scalev
      atmp = sum(g[1:n1])
      btmp = sum(g[(n1+1):(n1+n2)])
      
      ts1 = (t(xa*g[1:n1]) - xat*atmp)%*% xa
      ts2 = (t(ya*g[(n1+1):(n1+n2)]) - yat*btmp)%*% ya
      
      # ind1 = ((j-1)*(n1+n2)+1):((j-1)*(n1+n2)+n1)
      # ind2 = ((j-1)*(n1+n2)+n1+1):((j-1)*(n1+n2)+n2+n1)
      #
      # atmp = sum(g[ind1])
      # btmp = sum(g[ind2])
      #
      # ts1 = (t(xa*g[ind1]) - xat*atmp)%*% xa
      # ts2 = (t(ya*g[ind2]) - yat*btmp)%*% ya
      ts[j] = max(abs(ts1-ts2)/deno)
    }
  }
  
  ZCZt = Tnm > quantile(ts, 1-alpha)
  CLX = max((Sx-Sy)^2/(vx/n1+vy/n2))
  # CLX = CLX-(4*log(p)-log(log(p)))
  
  Sxx = Sx*n1/(n1-1)
  Syy = Sy*n2/(n2-1)
  
  SsS = (Sxx*n1 + Syy*n2)/(n1+n2)
  eta1 = ((n1-1)+2)*((n1-1)-1)
  eta2 = ((n2-1)+2)*((n2-1)-1)
  d1 = (1-(n1-1-2)/eta1)*sum(diag(Sxx%*%Sxx))
  d2 = (1-(n2-1-2)/eta2)*sum(diag(Syy%*%Syy))
  d3 = 2*sum(diag(Sxx %*% Syy))
  d4 = (n1-1)/eta1*sum(diag(Sxx))^2
  d5 = (n2-1)/eta2*sum(diag(Syy))^2
  th = 4*(((n1+n2-2)/((n1-1)*(n2-1)))^2)*
    ((n1+n2-2)^2/((n1+n2)*(n1+n2-2-1))*
       (sum(diag(SsS %*% SsS))-(sum(diag(SsS)))^2/(n1+n2-2)))^2
  Sc = (d1+d2-d3-d4-d5)/sqrt(th)
  
  # Res = matrix(0,3,3)
  # Res[1, ] = c(ZCZt, CLX>cri, abs(Sc)>qnorm(1-alpha/2))
  # Res[3, ] = c(Tnm, CLX, abs(Sc))
  # CLX = CLX-(4*log(p)-log(log(p)))
  # Res[2, ] = c(sum(ts>=Tnm)/J,
  #              1 - exp(-exp(-CLX/2)/(sqrt(8*pi))),
  #              (1-pnorm(abs(Sc)))*2)
  #
  # rownames(Res) = c('decision', 'p-value', 'test statistic')
  # colnames(Res) = c('HDtest', 'CLX', 'Scott')
  # return(Res)
  
  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))
  
  names(Tnm) = "Statistic"
  hd.res = list(statistics = Tnm, p.value = sum(ts>=Tnm)/J, alternative = "two.sided",
                method = "Two-Sample HD test", data.name = DNAME)
  class(hd.res) = "htest"
  
  names(CLX) = "Statistic"
  CLX.rev = CLX-(4*log(p)-log(log(p)))
  clx.p = 1 - exp(-exp(-CLX.rev/2)/(sqrt(8*pi)))
  clx.res = list(statistics = CLX, p.value = clx.p, alternative = "two.sided",
                 method = "Two-Sample CLX test", data.name = DNAME)
  class(clx.res) = "htest"
  
  names(Sc) = "Statistic"
  sc.p = (1-pnorm(abs(Sc)))*2
  sc.res = list(statistics = abs(Sc), p.value = sc.p, alternative = "two.sided",
                method = "Two-Sample Scott test", data.name = DNAME)
  class(sc.res) = "htest"
  
  res = list(HD = hd.res, CLX = clx.res, Scott = sc.res)
  return(res)
}

#' ' @name equalCovs
#' #' @aliases equalCovs
#' @title LC-test for equality of high dimensional covariances
#' @description Testing the equality of two high dimensional covariance matrices using the testing procedure by Li and Chen (2012).
#'
#' @param X The n x p data matrix from the sample 1
#' @param Y The n x p data matrix from the sample 2.
#' @param alpha The prescribed level of significance
#' @param DNAME Default input.
#'
#' @details Implementing testing procedure proposed by Li and Chen (2012) to test the equality of two sample high dimensional covariance matrices.
#'
#' @return Value of testing statistic, p-value, alternative hypothesis, and the name of testing procedure.
#'
#' @references J. Li and S. Chen (2012). Two sample tests for high-dimensional covariance matrices. Ann. Statist. 40, 908--940
#' @author Tong He
#' @export
#'


equalCovs <-function(X, Y, alpha, DNAME){
  # equalCovs <-function(sam1,sam2,size1,size2){
  ########################################################
  # obtain test statistic given in eqn (2.1) of the paper
  size1 = nrow(X)
  size2 = nrow(Y)
  A_mat <- X %*% t(X)
  out <- 0
  storage.mode(A_mat) <- "double"
  storage.mode(out) <- "double"
  nr <- as.integer(size1)
  
  # find A1
  #dyn.load("one1.dll")
  result1 <- .Fortran("code1", nr, A_mat, out=out, PACKAGE = "HDtest")
  A1 <- 2/(size1*(size1-1))*result1[[3]]
  
  # find A2
  #dyn.load("two2.dll")
  result2 <- .Fortran("code2", nr, A_mat, out=out, PACKAGE = "HDtest")
  A2 <- 4/(size1*(size1-1)*(size1-2))*result2[[3]]
  
  # find A3
  #dyn.load("three3.dll")
  result3 <- .Fortran("code3", nr, A_mat, out=out, PACKAGE = "HDtest")
  A3 <- 8/(size1*(size1-1)*(size1-2)*(size1-3))*result3[[3]]
  
  # obtain the statistic given by eqn (2.1) in our paper
  A_n1 <- A1-A2+A3
  
  # consider the sample 2
  
  B_mat <- Y %*% t(Y)
  out <- 0
  storage.mode(B_mat) <- "double"
  storage.mode(out) <- "double"
  nrr <- as.integer(size2)
  # find B1
  #dyn.load("one1.dll")
  result4 <- .Fortran("code1", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B1 <- 2/(size2*(size2-1))*result4[[3]]
  
  # find B2
  #dyn.load("two2.dll")
  result5 <- .Fortran("code2", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B2 <- 4/(size2*(size2-1)*(size2-2))*result5[[3]]
  
  # find B3
  #dyn.load("three3.dll")
  result6 <- .Fortran("code3", nrr, B_mat, out=out, PACKAGE = "HDtest")
  B3 <- 8/(size2*(size2-1)*(size2-2)*(size2-3))*result6[[3]]
  
  B_n2 <- B1-B2+B3
  
  ############################################################
  # obtain the test statistic given in eqn (2.2) of the paper
  ############################################################
  C_mat1 <- X %*% t(Y)
  C_mat2 <- Y %*% t(X)
  out <- 0
  storage.mode(C_mat1) <- "double"
  storage.mode(C_mat2) <- "double"
  storage.mode(out) <- "double"
  nrrr <- as.integer(size1)
  nl <- as.integer(size2)
  # find C1
  #dyn.load("four4.dll")
  result7 <- .Fortran("code4", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C1 <- -2/(size1*size2)*result7[[4]]
  
  # find C2
  #dyn.load("five5.dll")
  result8 <- .Fortran("code5", nl, nrrr, C_mat2, out=out, PACKAGE = "HDtest")
  C2 <- 4/(size1*size2*(size1-1))*result8[[4]]
  
  # find C3
  #dyn.load("five5.dll")
  result9 <- .Fortran("code5", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C3 <- 4/(size1*size2*(size2-1))*result9[[4]]
  
  # find C4
  #dyn.load("six6.dll")
  result10 <- .Fortran("code6", nrrr, nl, C_mat1, out=out, PACKAGE = "HDtest")
  C4 <- -4/(size1*size2*(size1-1)*(size2-1))*result10[[4]]
  
  # test statistic given by eqn (2.2) of the paper
  C_n <- C1+C2+C3+C4
  
  # the estimator
  T_n <- A_n1+B_n2+C_n
  
  # the standard deviation
  Sd_prime <- 2*(1/size1+1/size2)*((size1/(size1+size2))*A_n1+(size2/(size1+size2))*B_n2)
  
  test_stat <- T_n/Sd_prime
  pvalue <- 1-pnorm(test_stat)
  decision <- as.numeric(pvalue < alpha)
  # test <- matrix(0, 3, 1)
  # test[,1] <- c(decision, pvalue, test_stat)
  # colnames(test) = "LC"
  # rownames(test) = c('decision', 'p-value', 'test statistic')
  # return(test)
  # DNAME = deparse(substitute(X))
  # DNAME = paste(DNAME, "and", deparse(substitute(Y)))
  
  names(test_stat) = "Statistic"
  res = list(statistics = test_stat, p.value = pvalue, alternative = "two.sided",
             method = "Two-Sample LC test", data.name = DNAME)
  class(res) = "htest"
  
  return(res)
}
