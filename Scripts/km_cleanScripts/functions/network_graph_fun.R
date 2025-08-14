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

