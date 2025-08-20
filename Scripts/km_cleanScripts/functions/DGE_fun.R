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

# get names of list elements for plotting
namelist <- function(...) {
  group <- as.list(do.call("names", list(...))) 
  return(group)
}

# make list with automatic name propagation 
as_named_list <- function(...) {
  dots <- rlang::enquos(...)
  
  # name and evaluate
  values <- purrr::map(dots, rlang::eval_tidy)
  
  # names from quosures or fallback to default "unnamed"
  names <- purrr::map_chr(dots, function(dot) {
    if (rlang::quo_is_symbol(dot) || rlang::quo_is_call(dot)) {
      rlang::as_name(dot)
    } else {
      "unnamed"
    }
  })
  
  if (!is.null(named <- names(dots))) {
    names <- ifelse(named == "" | is.na(named), names, named)
  }
  setNames(values, names)
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

# subset Seurat obj by idents (set Idents(seurat_obj) prior)
subset_by_ident <- function(seurat_obj, 
                            subset_idents,
                            cluster = TRUE) {
  
  ident_list <- lapply(subset_idents, function(ident) {
    
    # subset to whatever metadata ident chosen 
    keep_cells <- colnames(seurat_obj)[Idents(seurat_obj) == ident]
    
    if(cluster) {
      cat(paste0("cluster ", ident, ": ", length(keep_cells), " cells retained\n"))
    } else {
      cat(paste0(ident, ": ", length(keep_cells), " cells retained\n")) 
    }
    
    # Subset object
    obj <- subset(seurat_obj, cells = keep_cells)
    
    if(cluster) {
      obj@project.name <- paste0("cluster_", ident)
    } else {
      obj@project.name <- as.character(ident)
    }
    
    obj
  })
  
  
  if(cluster) {
    names(ident_list) <- paste0("cluster_", subset_idents)
  } else {
    names(ident_list) <- subset_idents
  }

  return(ident_list)
}

# separate function to make melt matrix for RR
melt.matrix <- function (data,
                         na.rm = FALSE, 
                         value.name = "value",
                         ...) {
  Var1 <- Var2 <- NULL
  
  dt <- as.data.table(data)
  colnames(dt) <- as.character(1:ncol(dt))
  dt[, rownames := 1:nrow(dt)]
  
  melted_dt <- data.table::melt(dt, id.vars = "rownames", na.rm = na.rm, value.name = value.name)
  colnames(melted_dt)  <- c("Var1", "Var2", "value")
  melted_dt[, Var1 := as.double(Var1)]
  melted_dt[, Var2 := as.double(Var2)]
  
  
  
  return(melted_dt)
}

# Function to plot MERs faceted by predictor
plot_mer_facet <- function(mer_df, group) {
  
  p_est <- ggplot(mer_df, aes(x = interaction(Sex, Status), y = estimate, color = Status, shape = Sex)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.2,
                  position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ term, scales = "free_x") +
    labs(
      x = "Sex × Status",
      y = "Marginal Effect (MER)",
      title = paste0("MER: ", group, " neurons")
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_pred <- ggplot(mer, aes(x = interaction(Sex, Status), y = predicted, color = Status, shape = Sex)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Predicted Neuron Proportion", x = "Sex × Status",
         title = paste0("Predicted values: ", group, " neurons")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(grid.arrange(p_est, p_pred,
                      widths = c(1, 0.6),
                      nrow = 1))
}

# pull overlapping genes from RRHO2/RedRibbon results
get.overlaps <- function(x) {
  
  # overlaps per quadrant (RRHO2)
  RRHO_genes <- x[grep("^genelist", names(x))]
  
  overlap_RRHO2 <- map_dfr(names(RRHO_genes), function(nm) {
    overlap_name <- paste0("gene_list_overlap_", str_remove(nm, "genelist_"))
    tibble(quadrant = str_remove(nm, "genelist_"),
           genes = x[[nm]][[overlap_name]])
  })
  
  # overlaps per quadrant (RedRibbon)
  RedRibbon_genes <-x$RedRibbon.quads[names(x$RedRibbon.quads)]
  
  overlap_RedRibbon <- map_dfr(names(RedRibbon_genes), function(nm) {
    overlap_name <- paste(nm)
    tibble(quadrant = overlap_name,
           genes = x$df[x$RedRibbon.quads[[overlap_name]]$positions,1])
  })
  
  quad_genes <- as_named_list(overlap_RRHO2, overlap_RedRibbon)
  return(quad_genes)
  
}

get.percent.overlap <- function(dat) {
  
  sapply(dat, \(x) {
    
    if(is.null(quadrants_to_check)) {
      quadrants_to_check = unique(x$quadrant)
    }
    
    overlap_pct = data.frame()
    
    for(quad in quadrants_to_check) {
      x.quad = as.data.frame(x) %>% 
        filter(quadrant == quad) 
      
      pct.overlap = 100 * length(intersect(gene_list_input, x.quad$genes)) / length(x.quad$genes)
      
      total = data.frame(quad, 
                         pct.overlap, 
                         total.genes = length(x.quad$genes))
      overlap_pct = rbind(overlap_pct, total)
    }
    
    return(as_named_list(overlap_pct))
    
  })
}

#### processing for DGE/RRHO ####
prep.for.DGE <- function(list_of_SeuratObj,
                         assay = 'integrated',
                         SCTransformed = TRUE,
                         selection.method = 'vst',
                         pseudo_bulk = FALSE,
                         group.by = NULL) {
  
  if(!is.list(list_of_SeuratObj)) {
    list_of_SeuratObj <- as_named_list(list_of_SeuratObj)
  }
  
  # find union of variable features across all objects
  var_feats_list <- lapply(list_of_SeuratObj, \(x) {
    suppressWarnings(FindVariableFeatures(x, assay = assay,
                                          selection.method = selection.method)) -> x
    VariableFeatures(x, assay = assay)
  })
  
  union_feats <- unique(unlist(var_feats_list))
  message("Union of variable features: ", length(union_feats), " genes")
  
  get.data <- \(x) {
    
    gene <- Seurat::Project(x)
    cat(paste0("Subset is set to: ", gene, 
               "; if this is not what you intended, set project.name to subset\n"))
    
    if (is.null(x@meta.data$Cell_ID)) {
      x@meta.data$Cell_ID <- rownames(x@meta.data)
    }
    
    full_join(rownames_to_column(data.frame(x@reductions$umap@cell.embeddings), "Cell_ID"), 
              x@meta.data,
              join_by("Cell_ID")) -> x.umap
    
    get_assay_mat <- function(seurat_obj, assay_name, slot_name = "scale.data") {
      mat <- try(GetAssayData(seurat_obj, slot = slot_name, assay = assay_name), silent = TRUE)
      if(inherits(mat, "try-error") || nrow(mat) == 0 || ncol(mat) == 0) {
        message(paste0(slot_name, " slot empty or missing for assay '", assay_name, "', using 'data' slot instead"))
        mat <- GetAssayData(seurat_obj, slot = "data", assay = assay_name)
      }
      return(mat)
    }
    
    # Normalized & scaled counts for differential expression
    if (!SCTransformed) {
      
      message("Normalizing and scaling RNA assay")
      DefaultAssay(x) <- "RNA"
      
      x %>%
        NormalizeData(assay = "RNA", 
                      normalization.method = "LogNormalize", 
                      scale.factor = 1e4, 
                      verbose = FALSE) %>%
        FindVariableFeatures(assay = "RNA", 
                             verbose = FALSE) %>%
        ScaleData(assay = "RNA", 
                  verbose = FALSE) -> x
      
      norm_mat <- get_assay_mat(x, "RNA", "scale.data")
      norm_counts <- norm_mat[rownames(norm_mat) %in% union_feats, , drop = FALSE]
      
    } else {
      
      norm_mat <- get_assay_mat(x, "SCT", "scale.data")
      norm_counts <- norm_mat[rownames(norm_mat) %in% union_feats, , drop = FALSE]
    }
    
    # Raw counts from RNA assay data slot (un-normalized)
    raw_mat <- get_assay_mat(x, "RNA", "data")
    raw_counts <- raw_mat[rownames(raw_mat) %in% union_feats, , drop = FALSE]
    
    # Keep cells consistent
    common_cells <- intersect(colnames(norm_counts), x.umap$Cell_ID)
    norm_counts <- norm_counts[, common_cells, drop = FALSE]
    raw_counts <- raw_counts[, common_cells, drop = FALSE]
    x.umap <- x.umap %>% filter(Cell_ID %in% common_cells)
    
    # Reorder to match metadata cells
    norm_counts <- norm_counts[, x.umap$Cell_ID, drop = FALSE]
    raw_counts <- raw_counts[, x.umap$Cell_ID, drop = FALSE]
    
    # pseudo-bulk averaging per subset if requested
    if (pseudo_bulk && !is.null(group.by)) {
      if (all(group.by %in% colnames(x.umap))) {
        message("Performing pseudo-bulk averaging by groups: ", paste(group.by, collapse = ", "))
        grouping <- apply(x.umap[, group.by, drop = FALSE], 1, paste, collapse = "_")
        
        norm_counts <- sapply(split(seq_len(ncol(norm_counts)), grouping), function(idx) {
          Matrix::rowMeans(norm_counts[, idx, drop = FALSE])
        }) %>% as.matrix()
        
        raw_counts <- sapply(split(seq_len(ncol(raw_counts)), grouping), function(idx) {
          Matrix::rowMeans(raw_counts[, idx, drop = FALSE])
        }) %>% as.matrix()
        
        x.umap <- data.frame(group = unique(grouping))
        rownames(x.umap) <- x.umap$group
      } else {
        stop("One or more group.by columns not found in metadata")
      }
    }
    
    list(norm.counts = norm_counts, 
         raw.counts = raw_counts,
         meta.data = x.umap)
  }
  
  results <- lapply(list_of_SeuratObj, get.data)
  
  return(list(results = results,
              variable_features = union_feats))
}


run_limmatrend <- function(prep_results_list, 
                           outdir = "./limma_trend") {
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  subset_names <- names(prep_results_list)
  if (is.null(subset_names)) {
    subset_names <- paste0("subset", seq_along(prep_results_list))
  }
  
  results_list <- list()
  
  for (i in seq_along(prep_results_list)) {
    L <- prep_results_list[[i]]
    subset_name <- subset_names[i]
    
    message("Running limma-trend on subset '", subset_name, "' with ", ncol(L$raw.counts), " samples")
    
    if (!"orig.ident" %in% colnames(L$meta.data)) {
      stop("Metadata for subset '", subset_name, "' does not have 'orig.ident' column required for design")
    }
    
    orig_idents <- as.character(L$meta.data$orig.ident)
    clean_idents <- sub("\\.data$", "", orig_idents)  # remove trailing ".data"
    group <- factor(clean_idents)
    
    if (length(levels(group)) < 2) {
      warning(paste0("Subset '", subset_name, "' has fewer than 2 groups in orig.ident; skipping."))
      next
    }
    
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
    contrasts <- makeContrasts(
      Female_vs_Male = (female.dom + female.sub)/2 - (male.dom + male.sub)/2,
      Dom_vs_Sub = (female.dom + male.dom)/2 - (female.sub + male.sub)/2,
      Female_Dom_vs_Female_Sub = female.dom - female.sub,
      Male_Dom_vs_Male_Sub = male.dom - male.sub,
      Female_Dom_vs_Male_Dom = female.dom - male.dom,
      Female_Sub_vs_Male_Sub = female.sub - male.sub,
      levels = design
    )
    
    dge <- edgeR::DGEList(L$raw.counts)
    dge <- edgeR::calcNormFactors(dge)
    
    y <- methods::new("EList")
    y$E <- L$norm.counts
    
    fit <- limma::lmFit(y, design)
    fit <- limma::contrasts.fit(fit, contrasts)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    
    tt_to_df <- function(tt) {
      if (is.null(tt) || nrow(tt) == 0) {
        message("Empty or NULL topTable result, returning empty dataframe")
        return(data.frame(Gene = character(0)))
      }
      if (is.null(rownames(tt))) {
        message("No rownames in topTable result, creating a Gene column from row numbers")
        tt$Gene <- paste0("Gene_", seq_len(nrow(tt)))
      } else {
        tt <- tt %>% tibble::rownames_to_column("Gene")
      }
      return(tt)
    }
    
    
    ttFemale_vs_Male <- limma::topTable(fit, n = Inf, coef = "Female_vs_Male", adjust.method = "BH")
    ttDom_vs_Sub <- limma::topTable(fit, n = Inf, coef = "Dom_vs_Sub", adjust.method = "BH")
    ttFemale_Dom_vs_Female_Sub <- limma::topTable(fit, n = Inf, coef = "Female_Dom_vs_Female_Sub", adjust.method = "BH")
    ttMale_Dom_vs_Male_Sub <- limma::topTable(fit, n = Inf, coef = "Male_Dom_vs_Male_Sub", adjust.method = "BH")
    ttFemale_Dom_vs_Male_Dom <- limma::topTable(fit, n = Inf, coef = "Female_Dom_vs_Male_Dom", adjust.method = "BH")
    ttFemale_Sub_vs_Male_Sub <- limma::topTable(fit, n = Inf, coef = "Female_Sub_vs_Male_Sub", adjust.method = "BH")
    
    combined_df <- tt_to_df(ttFemale_vs_Male) %>%
      rename_with(~paste0(., "_Female_vs_Male"), -Gene) %>%
      dplyr::full_join(tt_to_df(ttDom_vs_Sub) %>% rename_with(~paste0(., "_Dom_vs_Sub"), -Gene), by = "Gene") %>%
      dplyr::full_join(tt_to_df(ttFemale_Dom_vs_Female_Sub) %>% rename_with(~paste0(., "_Female_Dom_vs_Female_Sub"), -Gene), by = "Gene") %>%
      dplyr::full_join(tt_to_df(ttMale_Dom_vs_Male_Sub) %>% rename_with(~paste0(., "_Male_Dom_vs_Male_Sub"), -Gene), by = "Gene") %>%
      dplyr::full_join(tt_to_df(ttFemale_Dom_vs_Male_Dom) %>% rename_with(~paste0(., "_Female_Dom_vs_Male_Dom"), -Gene), by = "Gene") %>%
      dplyr::full_join(tt_to_df(ttFemale_Sub_vs_Male_Sub) %>% rename_with(~paste0(., "_Female_Sub_vs_Male_Sub"), -Gene), by = "Gene")
    
    
    pdf(file = file.path(outdir, paste0("limmatrend_histograms_", subset_name, ".pdf")))
    par( mfrow= c(2,2) )
    hist(ttFemale_vs_Male$P.Value, 50, main = "F_vs_M_P.Values")
    hist(ttFemale_vs_Male$adj.P.Val, 50, main = "F_vs_M_adj.P.Values")
    hist(ttDom_vs_Sub$P.Value, 50, main = "D_vs_S_P.Values")
    hist(ttDom_vs_Sub$adj.P.Val, 50, main = "D_vs_S.adj P.Values")
    hist(ttFemale_Dom_vs_Female_Sub$P.Value, 50, main = "FD_vs_FS.P.Values")
    hist(ttFemale_Dom_vs_Female_Sub$adj.P.Val, 50, main = "FD_vs_FS.adj.P.Values")
    hist(ttMale_Dom_vs_Male_Sub$P.Value, 50, main = "MD_vs_MS.P.Values")
    hist(ttMale_Dom_vs_Male_Sub$adj.P.Val, 50, main = "MD_vs_MS.adj.P.Values")
    hist(ttFemale_Dom_vs_Male_Dom$P.Value, 50, main = "FD_vs_MD.P.Values")
    hist(ttFemale_Dom_vs_Male_Dom$adj.P.Val, 50, main = "FD_vs_MD.adj.P.Values")
    hist(ttFemale_Sub_vs_Male_Sub$P.Value, 50, main = "FS_vs_MS.P.Values")
    hist(ttFemale_Sub_vs_Male_Sub$adj.P.Val, 50, main = "FS_vs_MS.adj.P.Values")
    dev.off()
    
    pdf(file = file.path(outdir, paste0("MDS_", subset_name, ".pdf")))
    
    plotMDS(dge, 
            col = as.numeric(group), 
            pch = 19, 
            main = paste("MDS plot", subset_name))
    
    par(mfrow = c(2, 1))
    for(coef_name in colnames(contrasts)) {
      plotMD(fit, column = coef_name, main = paste("MD Plot -", coef_name))
    }
    
    dev.off()
    
    results_list[[subset_name]] <- list(
      ttFemale_vs_Male = ttFemale_vs_Male,
      ttDom_vs_Sub = ttDom_vs_Sub,
      ttFemale_Dom_vs_Female_Sub = ttFemale_Dom_vs_Female_Sub,
      ttMale_Dom_vs_Male_Sub = ttMale_Dom_vs_Male_Sub,
      ttFemale_Dom_vs_Male_Dom = ttFemale_Dom_vs_Male_Dom,
      ttFemale_Sub_vs_Male_Sub = ttFemale_Sub_vs_Male_Sub,
      combined_results = combined_df,
      meta.data = L$meta.data
    )
  }
  
  return(results_list)
}

# set colors to same scale for all plots
ggRedRibbon.rrho.scale <- function (self, 
                                    n = NULL, 
                                    labels = c("a", "b"), 
                                    show.quadrants = TRUE, 
                                    quadrants = NULL, 
                                    show.pval = TRUE,
                                    repel.force = 150, 
                                    base_size = 20, 
                                    .log10 = FALSE,
                                    new.max.log, # add value
                                    ...) {
  len <- length(self$data$a)
  
  if ( is.null(n) )
    n <- min(max(sqrt(len), 500), len)
  
  n.i <- n
  n.j <- n
  
  rrho <- rrho_rectangle(1, 1, len, len, n.i, n.j, self$data$a, self$data$b,  mode=self$enrichment_mode, LOG=TRUE)
  if (.log10)
  {
    rrho <- rrho / log(10)
  }
  log.label <- if (.log10) "log10" else "log"
  
  # set top of log scale
  max.log <- new.max.log
  
  
  if (0 == max.log)
  {
    min.log  <- -0.001
    max.log <- 0.001
  } else
    min.log <- - max.log
  
  # remove negative p-value scale
  ticks <- c(min.log,0, max.log)
  
  len.colors <- length(self$ggplot_colours)
  half.len.colors <- len.colors %/% 2
  colors.values <- seq(0, len.colors) /  len.colors
  
  
  ## Suppress warning RRHO: no visible binding for global variable ‘gg’
  Var1 <- Var2 <- value <- i <- j <- pvalue <- NULL
  
  gg <-  ggplot2::ggplot(melt.matrix(rrho), ggplot2::aes(Var1,Var2, fill=value)) +
    ggplot2::geom_raster() +
    ## ggplot2::scale_fill_gradientn(colours=self$ggplot_colours, name="-log p.val") +
    ggplot2::scale_fill_gradientn(colors = self$ggplot_colours,
                                  breaks = ticks,
                                  labels = format(ticks),
                                  limits=ticks[c(1,3)],
                                  ##limits=b[c(1,length(colors))],
                                  name=paste0("-", log.label, " p.val"),
                                  values=colors.values) +
    ggplot2::xlab(labels[1]) + ggplot2::ylab(labels[2]) +
    ## scale_x_continuous(labels = label_percent(accuracy = 1, scale = 100/n.i)) +
    ## scale_y_continuous(labels = label_percent(accuracy = 1, scale = 100/n.j) ) +
    ggplot2::scale_x_continuous(breaks = c(0 + n * 0.1, n - n * 0.1), labels = c("down", "up"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(0 + n * 0.1, n - n * 0.1), labels = c("down", "up"), expand = c(0, 0)) +
    ## ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(size=base_size),
                   legend.title = ggplot2::element_text(size = base_size * 7 / 10),
                   legend.text = ggplot2::element_text(size = base_size * 1 / 2),
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=base_size* 7 / 10,),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=base_size * 7 / 10, angle=90),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.ticks.length = ggplot2::unit(0, "pt"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.spacing = ggplot2::unit(0, "cm"),
                   plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
  
  
  ## find the middle of the plots
  if (show.quadrants || show.pval)
  {
    a_ltzero <- sum(self$data$a < 0)
    x.ind <- a_ltzero
    if ( x.ind == 0 )
      x.ind = len / 2
    
    b_ltzero <- sum(self$data$b < 0)
    y.ind <- b_ltzero
    if ( y.ind == 0 )
      y.ind = len / 2
    
    ## plot dotted quadrant lines
    if (show.quadrants)
    {
      gg  <- gg +
        ggplot2::geom_vline(ggplot2::aes(xintercept = x.ind * n.i / len), 
                            linetype = "dotted", colour = "gray10",size = 1) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = y.ind * n.j / len), 
                            linetype = "dotted", colour = "gray10",size = 1)
    }
    
    ## plot pvalue
    if (show.pval)
    {
      if (! is.null(quadrants) )
      {
        pval_size  <- as.integer(base_size * 1/5)
        quadrants_df <- as.data.frame(
          do.call(rbind, lapply(quadrants,
                                function (quadrant)
                                {
                                  if ( quadrant$pvalue > 0.05 || (! is.null(quadrant$padj) && quadrant$padj > 0.05) )
                                    return(NULL)
                                  
                                  pvalue <- if (.log10) quadrant$log_pvalue / log(10) else quadrant$log_pvalue
                                  pvalue.formatted <-  formatC(pvalue, format = "f", digits = 1)
                                  
                                  if (! is.null(quadrant$padj) )
                                  {
                                    padj <- if (.log10) quadrant$log_padj / log(10) else quadrant$log_padj
                                    pvalue.formatted <-  paste(pvalue.formatted,
                                                               "(padj =", formatC(padj, format = "f", digits = 1), ")")
                                  }
                                  data.frame(i=quadrant$i, j=quadrant$j,
                                             pvalue=pvalue.formatted, value=pvalue)
                                })))
        
        if ( nrow(quadrants_df) > 0 )
          gg <- gg +
          ggrepel::geom_text_repel(data=quadrants_df,
                                   ggplot2::aes(x=i * n.i / len, y=j * n.j / len,
                                                label=pvalue,
                                                colour = "gray"),
                                   hjust=1, vjust=1, colour = "black",
                                   force = repel.force, show.legend = FALSE, size = pval_size)
        
      }
    }
  }
  
  return(gg)
}

# invert quadrants to match pattern: upup = top right, downdown = bottom left
RRHO2_heatmap_axis.flip <- function (RRHO_obj, 
                                     maximum = NULL, 
                                     minimum = NULL, 
                                     colorGradient = NULL, 
                                     labels = NULL, 
                                     ...) {
  hypermat <- RRHO_obj$hypermat
  method <- RRHO_obj$method
  if (is.null(labels)) {
    labels <- RRHO_obj$labels
  }
  if (!is.null(maximum)) {
    hypermat[hypermat > maximum] <- maximum
  }
  else {
    maximum <- max(hypermat, na.rm = TRUE)
  }
  if (!is.null(minimum)) {
    hypermat[hypermat < minimum] <- minimum
  }
  else {
    minimum <- min(hypermat, na.rm = TRUE)
  }
  if (minimum > maximum) {
    stop("minimum > maximum, please check these function arguments!")
  }
  color.bar <- function(lut, min, max = -min, nticks = 11, 
                        ticks = seq(min, max, len = nticks), title = "") {
    scale <- (length(lut) - 1)/(max - min)
    plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", 
         xlab = "", yaxt = "n", ylab = "")
    mtext(title, 2, 2.3, cex = 0.8)
    axis(2, round(ticks, 0), las = 1, cex.lab = 0.8)
    for (i in 1:(length(lut) - 1)) {
      y <- (i - 1)/scale + min
      rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
    }
  }
  if (is.null(colorGradient)) {
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                     "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorGradient <- jet.colors(101)
  }
  layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  breaks <- seq(minimum, maximum, length.out = length(colorGradient) + 
                  1)
  image(hypermat[nrow(hypermat):1, ncol(hypermat):1], # flip the axes to match RedRibbon plots
        col = colorGradient, breaks = breaks, axes = FALSE, 
        ...)
  if (!is.null(labels)) {
    mtext(labels[2], 2, 0.5)
    mtext(labels[1], 1, 0.5)
  }
  if (method == "hyper") {
    atitle <- ifelse(RRHO_obj$log10.ind, "-log10(P-value)", 
                     "-log(P-value)")
    color.bar(colorGradient, min = minimum, max = maximum, 
              nticks = 6, title = atitle)
  }
  else if (method == "fisher") {
    atitle <- "log Odds"
    color.bar(colorGradient, min = minimum, max = maximum, 
              nticks = 6, title = atitle)
  }
  else {
    stop("internal error (1), please report this error to https://github.com/RRHO2/RRHO2/issues")
  }
  invisible(hypermat)
}


get.RRHO <- function(results_list,
                     compare.across,
                     group.by,
                     new.max.log = NULL,
                     outdir = "./RRHO") {
  
  rrho_results <- list()
  
  filter_comparisons <- \(x) {
    clean_names <- sub("^tt", "", names(x)[grepl("^tt", names(x))])
    
    compare.across.values <- compare.across
    group.by.values <- group.by
    
    # split into before vs after
    parts <- strsplit(clean_names, "_vs_")
    keep <- vapply(parts, function(x) {
      left  <- unlist(strsplit(x[1], "_"))
      right <- unlist(strsplit(x[2], "_"))
      
      # extract the compare.across and group.by values from each side
      left_compare  <- left[left %in% compare.across.values]
      right_compare <- right[right %in% compare.across.values]
      left_group  <- left[left %in% group.by.values]
      right_group <- right[right %in% group.by.values]
      
      # comparisons where compare.across values are different, AND group.by values are the same
      length(left_compare)  == 1 &&
        length(right_compare) == 1 &&
        left_compare != right_compare &&
        length(left_group)  == 1 &&
        length(right_group) == 1 &&
        left_group == right_group
    }, logical(1))
    clean_names[keep]
  }
  
  comparisons <- lapply(results_list, filter_comparisons)
  
  for (dataset_name in names(results_list)) {
    dataset <- results_list[[dataset_name]]
    comps <- comparisons[[dataset_name]]
    
    # skip datasets with no comparisons
    if (length(comps) < 2) next
    
    # pull comparison dfs
    left_df_raw  <- dataset[[paste0("tt", comps[1])]]
    right_df_raw <- dataset[[paste0("tt", comps[2])]]
    
    # ensure same genes
    common_genes <- intersect(rownames(left_df_raw), rownames(right_df_raw))
    left_df_raw  <- left_df_raw[common_genes, ]
    right_df_raw <- right_df_raw[common_genes, ]
    
    # flip signs if compare.across order differs
    left_compare_valUP <- strsplit(comps[1], "_vs_")[[1]][1]  # first part
    right_compare_valUP <- strsplit(comps[2], "_vs_")[[1]][1]
    
    if (grepl(compare.across[1], left_compare_valUP)) {
      left_sign <- 1
    } else {
      left_sign <- -1
    }
    if (grepl(compare.across[1], right_compare_valUP)) {
      right_sign <- 1
    } else {
      right_sign <- -1
    }
    
    left_df  <- data.frame(gene = rownames(left_df_raw),
                           value = -log10(left_df_raw$P.Value) * sign(left_df_raw$t) * left_sign,
                           stringsAsFactors = FALSE)
    
    right_df <- data.frame(gene = rownames(right_df_raw),
                           value = -log10(right_df_raw$P.Value) * sign(right_df_raw$t) * right_sign,
                           stringsAsFactors = FALSE)
    
    # create RRHO2 object
    rrho_obj <- RRHO2_initialize(right_df, left_df, log10.ind = TRUE,
                                 labels = c(comps[1], comps[2]))
    
    # Now create df for RedRibbon
    # needs to have an id (gene) col and one called 'a' and one called 'b'
    df <- right_df %>%
      dplyr::select(gene, value) %>%
      dplyr::rename('a' = value) %>%
      full_join(left_df %>%
                  dplyr::select(gene, value) %>%
                  dplyr::rename('b' = value),
                by = 'gene')
    
    # Create RedRibbon object
    rr <- RedRibbon(df, enrichment_mode="hyper-two-tailed")
    
    # Run the overlap using evolutionary algorithm,
    # computing permutation adjusted p-value for the four quadrants
    RedRibbon.quads <- quadrants(rr,
                                 algorithm = "ea",
                                 permutation = TRUE,
                                 whole = FALSE)
    
    # store rrho object and RedRibbon.quads in output list
    rrho_results[[dataset_name]] <- append(rrho_obj,
                                           as_named_list(RedRibbon.quads, df))
    
    # RedRibbon plot
    if(!is.null(new.max.log)) {
      p1 <- ggRedRibbon.rrho.scale(rr,
                                   quadrants = RedRibbon.quads,
                                   new.max.log = new.max.log) +
        coord_fixed(ratio = 1,
                    clip = "off") +
        xlab(gsub("_", " ", comps[1])) +
        ylab(gsub("_", " ", comps[2])) +
        ggtitle(dataset_name) +
        theme(plot.title = element_text(size = 25,
                                        face = "bold"),
              plot.margin = unit(c(0,0.5,0,0.5), units = "cm")) # top, right, bottom, left
      
    } else {
      
      p1 <- ggRedRibbon(rr, 
                        quadrants = RedRibbon.quads) +
        coord_fixed(ratio = 1,
                    clip = "off") +
          xlab(gsub("_", " ", comps[1])) +
          ylab(gsub("_", " ", comps[2])) +
          ggtitle(dataset_name) +
          theme(plot.title = element_text(size = 25,
                                          face = "bold"),
                plot.margin = unit(c(0,0.5,0,0.5), units = "cm")) # top, right, bottom, left
    }
    
    # Save plots
    out_dir <- file.path(outdir, dataset_name)
    
    if(!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    ggsave(file.path(out_dir, paste0("RedRibbon_concordance_by_", deparse(substitute(group.by)), ".png")),
           plot = p1, height = 10, width = 10, dpi = 300)
    
    png(file.path(out_dir, paste0("RRHO2_concordance_by_", deparse(substitute(group.by)), ".png")),
        height = 10, width = 10, units = "in", res = 300)

    plot.new()
    RRHO2_heatmap_axis.flip(rrho_obj)
    dev.off()
    
    cat("saving plots for subset:", dataset_name, "\n")
  }
  return(rrho_results)
}

# check for overlapping genes in RRHO2/RedRibbon quadrants 
check.quadrants <- function(rrho_results_list,
                            gene_list,
                            quadrants_to_check = NULL,
                            outdir = NULL) {
  
  if(!is.list(rrho_results_list)) {
    stop("input must be a results list output by get.RRHO")
  }
  
  if(!is.vector(gene_list) & !is.character(gene_list)) {
    gene_list_input = as.vector(gene_list$gene, mode = "character")
  } else {
    gene_list_input = gene_list
  }
  
  overlaps_list <- lapply(rrho_results_list, get.overlaps)
  
  overlaps_pct <- lapply(overlaps_list, get.percent.overlap)
  
  if(!is.null(outdir)) {
    
    for(subset_name in names(rrho_results_list)) {
      
      pctRRHO <- overlaps_pct[[subset_name]]$overlap_RRHO2.overlap_pct
      pctRedRibbon <- overlaps_pct[[subset_name]]$overlap_RedRibbon.overlap_pct
      
      overlaps_list[[subset_name]]$overlap_RRHO2 %>% 
        append(., as_named_list(pctRRHO)) -> overlaps_list[[subset_name]]$overlap_RRHO2
      
      overlaps_list[[subset_name]]$overlap_RedRibbon %>% 
        append(., as_named_list(pctRedRibbon)) -> overlaps_list[[subset_name]]$overlap_RedRibbon
      
      overlaps_list = overlaps_list %>% 
        append(., as_named_list(gene_list))
      
      if(length(grep(paste0("/",subset_name,"$"), list.dirs(outdir))) > 0) {

        dir <- paste(outdir, subset_name, subset_name, sep = "/")
        file <- paste(deparse(substitute(gene_list)), "overlaps.rda", sep = "_")
        
        filename <- paste(dir, file, sep = "_")
        
      } else {
        
        dir <- paste(outdir, subset_name, sep = "/")
        file <- paste(deparse(substitute(gene_list)), "overlaps.rda", sep = "_")
        
        filename <- paste(dir, file, sep = "_")
      }
      save(overlaps_list, file = filename)
    }
  }
  
  return(overlaps_pct)
}


plot.RRHO.counts <- function(rrho_results_list,
                             outdir) {
  
  for (set in names(rrho_results)) {
    
    gene_lists <- rrho_results[[set]]$RedRibbon.quads[ names(rrho_results[[set]]$RedRibbon.quads) ]
    
    map_dfr(names(gene_lists), function(q) {
      data.frame(quadrant = q,
                 overlap_length = length(rrho_results[[set]]$df[ gene_lists[[q]]$positions, 1 ] )) 
      }
    ) %>% 
      ggplot(aes(x = reorder(quadrant, -overlap_length), 
                 y = overlap_length,
                 label = overlap_length)) +
      geom_label() +
      theme_classic() +
      ylab('Number of genes per quadrant') +
      xlab('') +
      ggtitle(paste0(set,' RedRibbon quadrants'))
    
    cat("Saving RedRibbon counts plot for", set, "\n")
    
    ggsave(paste(outdir, set,'RedRibbon_quadGenes_count.png', sep = "/"),
           width = 5, height = 4)
    
    gene_lists <- rrho_results[[set]][ grep("^genelist", names(rrho_results[[set]])) ]
    
    map_dfr(names(gene_lists), function(q) {
      nm <- paste0("gene_list_overlap_", str_remove(q, "genelist_"))
      data.frame(quadrant = str_remove(q, "genelist_"),
                 overlap_length = length(gene_lists[[q]][[nm]]))
      }
    ) %>% 
      ggplot(aes(x = reorder(quadrant, -overlap_length), 
                 y = overlap_length,
                 label = overlap_length)) +
      geom_label() +
      theme_classic() +
      ylab('Number of genes per quadrant') +
      xlab('') +
      ggtitle(paste0(set,' RRHO2 quadrants'))
    
    ggsave(paste(outdir, set,'RRHO2_quadGenes_count.png', sep = "/"),
           width = 5, height = 4)
    
    cat("Saving RRHO2 counts plot for", set, "\n")
  }
}

plot.overlaps <- function(rrho_results_list,
                          gene_list,
                          quadrants_to_check = NULL,
                          subtitle = NULL,
                          group.by,
                          outdir) {
  
  if(!is.list(rrho_results_list)) {
    stop("input must be a results list output by get.RRHO")
  }
  
  if(!is.null(subtitle)) {
    subtitle = subtitle
  } else {
    subtitle = NA
  }
  
  gene_list_name <- deparse(substitute(gene_list))
    
  plot.dat <- function(dat) {
    dat <- mapply("c", dat, namelist(dat),
                SIMPLIFY = FALSE)
    
    lapply(dat, \(x) {
      
      analysis <- sub(".*?_", "", x[[3]])
      
      if(is.null(quadrants_to_check)) {
        quadrants_to_check = unique(x$quadrant)
      }
      
      for(quad in quadrants_to_check) {
        x.quad = as.data.frame(x) %>% 
          filter(quadrant == quad) 
        
        gene_list_q <- gene_list %>% 
          filter(gene %in% x.quad$genes) 
        
        col1 <- paste(group.by)
        col2 <- "count"
        
        gene_list_q %>% 
          group_by(get(group.by)) %>% 
          summarize(count = n()) %>% 
          rename_with(.cols = 1, ~paste(group.by)) %>% 
          ggplot(aes(reorder(.data[[col1]], -.data[[col2]]),
                     .data[[col2]])) +
          geom_point() +
          theme_classic() +
          xlab(col1) +
          ylab('Number of genes') +
          ggtitle(paste(analysis, quad, "genes in", gene_list_name, sep = " "),
                         subtitle = subtitle)
        
        ggsave(paste0(paste(outdir, set, analysis, sep = "/"), "_",
                      paste(quad, "Genes_in", gene_list_name, sep = "_"),
                      ".png"), height = 7, width = 7)
        
      }
      
    })
  }
  
  # plot number of cluster markers in RRHO2 each quadrant
  for (set in names(rrho_results_list)) {
    
    overlaps_list <- lapply(rrho_results, get.overlaps)
    
    cat("Saving RRHO2 plot(s) for", set, "\n")
    
    lapply(overlaps_list, plot.dat)
    
  }
}



get.GOtop15.RRHO <- function(rrho_results_list,
                             quadrants_to_check = NULL,
                             outdir) {
  
  if(!is.list(rrho_results_list)) {
    stop("input must be a results list output by get.RRHO")
  }
  
  plot.go <- function(dat) {
    
    # extract names
    raw_name <- dat[[3]]  
    set_name <- sub("\\..*$", "", raw_name)                
    analysis <- sub(".*_", "", sub(".*\\.", "", raw_name))  
    
    if(startsWith(outdir, "./")) {
      outdir <- sub("^\\./", "", outdir)
    } 
    
    cat("starting GO enrichment analysis on", set_name, "( genes from:", analysis, ")\n")
    
    x <- as.data.frame(dat[-length(dat)])  
    
    if(is.null(quadrants_to_check)) {
      quadrants_to_check <- unique(x$quadrant)
    }
    
    for(quad in quadrants_to_check) {
      x.quad <- x %>% 
        filter(quadrant == quad) 
      
      cat("generating GO results for", set_name, "quadrant", quad, "\n")
      
      GO.obj <- enrichGO(gene = unique(x.quad$genes), 
                         OrgDb = "org.Mm.eg.db",
                         keyType = "SYMBOL", 
                         ont = "BP")
      
      out_dir <- file.path(outdir, set_name)
      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      
      if(sum(GO.obj@result$p.adjust < 0.05) > 5) {
        
        # barplot
        png(file.path(out_dir,
                      paste(analysis,
                            quad,
                            "GO_results.png", 
                            sep = "_")),
            res = 300, width = 10, height = 10, units = 'in')
        
        print(barplot(GO.obj, showCategory = 15))
        
        dev.off()
        cat("saved barplot for quadrant", quad, "\n")
        
        # treeplot
        png(file.path(out_dir, 
                      paste(analysis, 
                            quad, 
                            "GO_resultsTree.png", 
                            sep = "_")),
            res = 300, width = 13.5, height = 10, units = 'in')
        
        suppressMessages({  print(GO.obj %>% 
                                    pairwise_termsim() %>% 
                                    treeplot())  })
        
        dev.off()
        cat("saved treeplot for quadrant", quad, "\n")
        
      } else {
        message("not enough significant GO terms; skipping quadrant ", quad)
      }
    }
  }
  
  # build overlaps list
  overlaps_list <- lapply(rrho_results_list, get.overlaps) %>% 
    unlist(., recursive = FALSE)
  
  overlaps_list <- mapply("c", overlaps_list, 
                          namelist(overlaps_list),
                          SIMPLIFY = FALSE)
  
  for(i in seq_along(overlaps_list)) {
    dat <- overlaps_list[[i]]
    plot.go(dat)
  }
}


# 
# gene_lists <- rrho_results[[set]]$RedRibbon.quads[ names(rrho_results[[set]]$RedRibbon.quads) ]
# cat("Saving RedRibbon plot(s) for", set, "\n")
# 
# for (q in quadrants_to_check) {
#   
#   map_dfr(names(gene_lists), function(q) {
#     data.frame(quadrant = q,
#                overlaps = rrho_results[[set]]$df[ gene_lists[[q]]$positions, 1 ])
#     }
#   ) -> quad_genes
#   
#   gene_list_q = subset(gene_list_input, intersect(quad_genes))
#   
#   gene_list_q %>% 
#     group_by(group.by) %>% 
#     summarize(count = n()) %>% 
#     ggplot(aes(x = reorder(as.character(cluster), -count),
#                y = count)) +
#     geom_point() +
#     theme_classic() +
#     xlab(paste(group.by)) +
#     ylab('Number of genes') +
#     ggtitle(paste0('Overlap of RedRibbon', q, "genes and ", deparse(substitute(gene_list)),
#                    subtitle = subtitle))
#   
#   ggsave(paste0(paste(outdir, set, "RedRibbon", sep = "/"), 
#                 paste(q, "Genes_in", deparse(substitute(gene_list)), sep = "_"),
#                 ".png"), height = 7, width = 7)
# }

# # create list of RRHO outcomes and deal with empty quadrants for graphing
# graphRRHO <- function(rr, df) {
#   
#   safe_quadrant_df <- function(quad_data, label, df) {
#     if (is.null(quad_data) || is.null(quad_data$positions) || length(quad_data$positions) == 0) {
#       message(dataset_name, ": quadrant '", label, "' is empty — skipping.")
#       return(data.frame(Gene = NA,
#                         RRquadrant = label,
#                         stringsAsFactors = FALSE))
#     } else {
#       return(data.frame(Gene = df[quad_data$positions, "gene"],
#                         RRquadrant = label,
#                         stringsAsFactors = FALSE))
#     }
#   }
#   
#   # Try to get quadrants, return empty df if it fails
#   quad <- tryCatch(
#     quadrants <- rr,
#     error = function(e) {
#       message("Error getting quadrants: ", e$message)
#       return(NULL)
#     }
#   )
#   
#   if (is.null(quad)) {
#     return(data.frame(Gene = character(0),
#                       RRquadrant = character(0),
#                       stringsAsFactors = FALSE))
#   }
#   
#   RedRibbon.quads <- rbind(
#     safe_quadrant_df(quad$upup, "upup", df),
#     safe_quadrant_df(quad$downdown, "downdown", df),
#     safe_quadrant_df(quad$updown, "updown", df),
#     safe_quadrant_df(quad$downup, "downup", df)
#   )
#   
#   return(RedRibbon.quads)
# }
