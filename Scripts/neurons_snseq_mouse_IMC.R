#### Mouse snseq seurat analysis
### hdWGCNA neurons
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###hdWGCNA
##https://smorabit.github.io/hdWGCNA/


### set working directory
setwd("/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/")

#### installation ####

### hdwgcna
# # create new conda environment for R
# conda create -n hdWGCNA -c conda-forge r-base r-essentials
# 
# # activate conda environment
# conda activate hdWGCNA

# install BiocManager
# install.packages("BiocManager")

# install Bioconductor core packages
# BiocManager::install()

# # install additional packages:
# install.packages(c("Seurat", "WGCNA", "igraph", "devtools"))

#Now you can install the hdWGCNA package using devtools.
# devtools::install_github('smorabit/hdWGCNA', ref='dev')

#load gene overlap
# BiocManager::install("GeneOverlap")

#### load libraries ####

#load libraries
library(Seurat)
library(GeneOverlap)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(cowplot)
library(patchwork)
library(corrplot)
library(ggrepel)
library(ggpattern)
library(multcomp)
library(emmeans)
library(limma)
library(edgeR)
library(RedRibbon)
library(tidyverse)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

#### load data ####
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

#### functions ####
### create function to get dotplot data
DotPlot.data = function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                         scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
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

### create function for bigger network graphs
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

#### neurons ####
mouse.snseq.combined.sct.neurons = mouse.snseq.combined.sct

# #set idents
# Idents(object = mouse.snseq.combined.sct.neurons) <- "sctype.all.hypomap.exp"
# 
# #subset to neurons
# mouse.snseq.combined.sct.neurons = subset(mouse.snseq.combined.sct.neurons,
#                                             idents = c("C25-3: GLU-3",
#                                                        "C25-6: GLU-6",
#                                                        "C25-11: GABA-1",  
#                                                        "C25-14: GABA-5",
#                                                        "C25-8: GLU-8",
#                                                        "C25-5: GLU-5",
#                                                        "C25-10: GABA-2"))

#set idents
Idents(object = mouse.snseq.combined.sct.neurons) <- "parent_id.broad.prob"

#subset to neurons
mouse.snseq.combined.sct.neurons = subset(mouse.snseq.combined.sct.neurons,
                                          idents = c("C7-2: GABA",
                                                     "C7-1: GLU"))

# need to set to integrated for clustering
DefaultAssay(mouse.snseq.combined.sct.neurons) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
mouse.snseq.combined.sct.neurons.recluster = mouse.snseq.combined.sct.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.4)

## graph 
# idents to new clusters
Idents(object = mouse.snseq.combined.sct.neurons.recluster) <- "seurat_clusters"

#graph umap
#clusters
mouse.snseq.combined.sct.neurons.recluster %>%
  DimPlot(label = TRUE) +
  theme_classic()
ggsave('./neurons/UMAP/UMAP neurons recluster clusters.png',
       height = 10,
       width = 10)
# #cell type
# mouse.snseq.combined.sct.neurons.recluster %>%
#   DimPlot(group.by = 'parent_id.broad.prob')+
#   theme_classic()
# ggsave('./neurons/UMAP/UMAP neurons recluster parent_id.broad.png',
#        height = 10,
#        width = 10)
# #cell type
# mouse.snseq.combined.sct.neurons.recluster %>%
#   DimPlot(group.by = 'orig.ident',
#           size = 3)+
#   theme_classic() +
#   scale_color_manual(values = c("red",
#                                 'cyan',
#                                 'black',
#                                 'green'))
# ggsave('./neurons/UMAP/UMAP neurons recluster orig.ident.png',
#        height = 10,
#        width = 10)
# 
# #cell type
# mouse.snseq.combined.sct.neurons.recluster %>%
#   DimPlot(group.by = 'orig.ident',
#           size = 3)+
#   theme_classic() +
#   scale_color_manual(values = c("black",
#                                 'red',
#                                 'black',
#                                 'red'))
# ggsave('./neurons/UMAP/UMAP neurons recluster status.png',
#        height = 10,
#        width = 10)
# 
# ## dom candidate genes
# DotPlot(mouse.snseq.combined.sct.neurons.recluster,
#         assay = 'SCT',
#             features = c("Lhx8",
#                          "Crym",
#                          "Ache",
#                          "Myb",
#                          "Mog")) +
#   theme_classic()
# ggsave('./neurons/Dom genes dotplot.png',
#        height = 10,
#        width = 10)


## scale counts
# celltype
# genotype
# dataframe of scale
neuron.genotype.scale = mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::select(c(orig.ident,
                  Genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.neuron.count = Freq)

#dataframe of cluster counts
neuron.cluster.genotype.count = mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::select(c(orig.ident,
                  seurat_clusters,
                  Genotype)) %>% 
  table() %>%
  as.data.frame()

#combine counts and scale
neuron.cluster.genotype.count = neuron.cluster.genotype.count %>% 
  full_join(neuron.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.neuron.count)

# #graph
# neuron.cluster.genotype.count %>%
#   ggplot(aes(x = seurat_clusters,
#              y = Percent,
#              color = orig.ident)) +
#   geom_boxplot(width = 0,
#                position = position_dodge(0.75)) +
#   geom_point(position = position_dodge(0.75),
#              aes(shape = Genotype,
#                  group = orig.ident))+
#   facet_grid(~seurat_clusters,
#              scales = "free_x") +
#   theme_classic()
# ggsave('./neurons/UMAP/Clusters vs orig.ident neurons recluster genotype scaled.png',
#        height = 10,
#        width = 10)

# graph
# poster
neuron.cluster.genotype.count %>%
  ggplot(aes(x = seurat_clusters,
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = Genotype,
                 group = orig.ident),
             size = 3)+
  facet_grid(~seurat_clusters,
             scales = "free_x") +
  theme_classic()+
  scale_color_manual(values = c("#CC5500",
    "#FF6A00",
    "#0077CC",
    "#0095FF"))
ggsave('./neurons/UMAP/Clusters vs orig.ident neurons recluster genotype scaled poster.png',
       height = 10,
       width = 10)

#### Neuron cluster stats ####
#### cluster names
neuron.cluster.names = mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  select(parent_id.exp.prob,
         parent_id.broad.prob,
         seurat_clusters)  %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/total)

## save neuron cluster names
write.csv(neuron.cluster.names,
          'neurons/neuron.cluster.names.csv')

### add variable for status and sex
neuron.cluster.genotype.count = neuron.cluster.genotype.count %>% 
  mutate(Sex = ifelse(grepl("Female",
                            orig.ident),
                      "Female",
                      "Male")) %>% 
  mutate(Status = ifelse(grepl("Dom",
                            orig.ident),
                      "Dom",
                      "Sub"))

# #### Anova
# ### run anova on every cluster
# # create empty data frame
# neuron.cluster.aov = data.frame(matrix(ncol = 6, nrow = 0))
# 
# #provide column names
# colnames(neuron.cluster.aov) <- c("Df",
#                   "Sum Sq",
#                   "Mean Sq",
#                   "F value",
#                   "Pr(>F)",
#                   "seurat_clusters")
# 
# ## loop through each cluster
# for (i in unique(neuron.cluster.genotype.count$seurat_clusters)) {
#   ## run two way anova on sex and status
#   tmp <- aov(Percent ~ Sex*Status, 
#                   data = neuron.cluster.genotype.count %>% 
#                     filter(seurat_clusters == i))
#   
#   # get results
#   tmp2 = summary(tmp)[[1]] %>% 
#     as.data.frame() %>% 
#     mutate(seurat_clusters = i)
#   
#   # combine in data frame
#   neuron.cluster.aov = neuron.cluster.aov %>% 
#     rbind(tmp2)
# }
# 
# # drop residuals
# neuron.cluster.aov = neuron.cluster.aov %>% 
#   na.omit()
# 
# ## adjust pvalue
# neuron.cluster.aov = neuron.cluster.aov %>% 
#   mutate(p_adj = p.adjust(`Pr(>F)`,
#                           method = 'fdr',
#                           n = nrow(neuron.cluster.aov)))
# 
# ## create column for comparison
# neuron.cluster.aov = neuron.cluster.aov %>% 
#   rownames_to_column('Type')
# 
# # clean up 
# neuron.cluster.aov = neuron.cluster.aov %>% 
#   mutate(Comp.type = case_when(grepl("Sex", Type) ~ "Sex",
#                                grepl("Status", Type) ~ "Status")) %>% 
#   mutate(Comp.type = ifelse(grepl("Sex:Status", Type),
#                             "Sex:Status",
#                             Comp.type))
# 
# 
# ## save anova table
# write_csv(neuron.cluster.aov,
#           file = 'neurons/neuron_seurat_cluster_anova.csv')

#### GLM binomial
### run glm on every cluster
## use with proportion data
# create empty data frame
neuron.cluster.glm = data.frame(matrix(ncol = 7, nrow = 0))

#provide column names
colnames(neuron.cluster.glm) <- c("contrast",
                                  "Estimate",
                                  "Std.error",
                                  "z.value",
                                  "adj.pvalue",
                                  "z.ratio",
                                  "seurat_clusters")
## loop through each cluster
for (i in unique(neuron.cluster.genotype.count$seurat_clusters)) {
  ## run interaction binomial GLM on sex and status
  tmp = glm(cbind(Freq, Others) ~ Sex*Status, 
               data = neuron.cluster.genotype.count %>% 
               filter(seurat_clusters == i) %>% 
               mutate(Others = total.neuron.count - Freq,
                      Sex = as.factor(Sex),
                      Status = as.factor(Status)), 
             family = binomial)
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  
  # save results
  tmp.res.df = data.frame(Estimate = tmp.res$test$coefficients,
                          Std.error = tmp.res$test$sigma,
                          z.value = tmp.res$test$tstat,
                          adj.pvalue = tmp.res$test$pvalues,
                          z.ratio = NA) %>% 
    rownames_to_column('contrast') %>% 
    mutate(seurat_clusters = i)
  
  # graph confidence intervals
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint)
  ci_df$group <- row.names(ci_df)
  
  # graph results
  ci_df %>% 
    filter(group != "(Intercept)") %>% 
  ggplot(aes(x = Estimate, 
             y = group)) +
    geom_point() +
    geom_errorbar(aes(xmin = lwr, 
                      xmax = upr), 
                  width = 0.2) +
    geom_vline(xintercept = 0, 
               linetype = 2) +
    labs(x = "Difference in Means", 
         y = "Comparison") +
    theme_classic() +
    ggtitle(paste0("Neuron cluster ",
                   i,
                   ": binomial GLM cell count"))
  ggsave(paste0('neurons/neuron_cluster/',
                i,
                ' binomial GLM cell count.png'))
  
  ## run all pairwise comparison
  # Tukey needed for all pairwise comparisons
  tmp2 = glm(cbind(Freq, Others) ~ orig.ident, 
            data = neuron.cluster.genotype.count %>% 
              filter(seurat_clusters == i) %>% 
              mutate(Others = total.neuron.count - Freq,
                     Sex = as.factor(Sex),
                     Status = as.factor(Status)), 
            family = binomial)
  
  # use emmeans to get pairwise differences
  emm = emmeans(tmp2,
                          ~ "orig.ident")
  emm.res = pairs(emm)


  # save results
  emm.res.df = data.frame(Estimate = summary(emm.res)$estimate,
                    Std.error = summary(emm.res)$SE,
                    z.ratio = summary(emm.res)$z.ratio,
                    adj.pvalue = summary(emm.res)$p.value,
                    contrast = summary(emm.res)$contrast,
                    z.value = NA) %>%
    mutate(seurat_clusters = i)
  
  ## combine in data frame
  neuron.cluster.glm = neuron.cluster.glm %>% 
    rbind(tmp.res.df) %>% 
    rbind(emm.res.df)
}


## FDR correction for all comparisons
neuron.cluster.glm.fdr = neuron.cluster.glm %>% 
  filter(is.na(z.ratio)) %>% 
  filter(contrast != "(Intercept)") %>% 
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                          method = 'fdr'))

# combine with data frame
neuron.cluster.glm = neuron.cluster.glm %>% 
  left_join(neuron.cluster.glm.fdr)



## create rounded p.value to make it easier to read
neuron.cluster.glm = neuron.cluster.glm %>% 
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

## save glm table
write_csv(neuron.cluster.glm,
          file = 'neurons/neuron_seurat_cluster_glm.csv')

## create p value dataframe
neuron.cluster.glm.p = neuron.cluster.glm %>% 
  dplyr::select(seurat_clusters,
         adj.pvalue.round,
         contrast) %>% 
  pivot_wider(names_from = contrast,
              values_from = adj.pvalue.round) 

#### neuron cluster graphs ####
### Number of reads and genes per cluster
library(gghalves)

# https://romanhaa.github.io/projects/scrnaseq_workflow/#number-of-transcripts-and-expressed-genes
temp_labels <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally()

p1 <- ggplot() +
  geom_half_violin(
    data = mouse.snseq.combined.sct.neurons.recluster@meta.data, aes(integrated_snn_res.0.4, nCount_RNA, fill = integrated_snn_res.0.4),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = mouse.snseq.combined.sct.neurons.recluster@meta.data, aes(integrated_snn_res.0.4, nCount_RNA, fill = integrated_snn_res.0.4),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = integrated_snn_res.0.4, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_color_manual(values = custom_colors$discrete) +
  # scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

p2 <- ggplot() +
  geom_half_violin(
    data = mouse.snseq.combined.sct.neurons.recluster@meta.data, aes(integrated_snn_res.0.4, nFeature_RNA, fill = integrated_snn_res.0.4),
    side = 'l', show.legend = FALSE, trim = FALSE
  ) +
  geom_half_boxplot(
    data = mouse.snseq.combined.sct.neurons.recluster@meta.data, aes(integrated_snn_res.0.4, nFeature_RNA, fill = integrated_snn_res.0.4),
    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = integrated_snn_res.0.4, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_color_manual(values = custom_colors$discrete) +
  # scale_fill_manual(values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(
  'neurons/ncount_nfeature_by_cluster.png',
  p1 + p2 + plot_layout(ncol = 1),
  height = 7, width = 14
)

#### create counts table 
## sample by cluster
table_orig.idents_by_clusters <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(orig.ident, integrated_snn_res.0.4) %>%
  summarize(count = n()) %>%
  spread(integrated_snn_res.0.4, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('orig.ident', 'total_cell_count', everything())) %>%
  arrange(factor(orig.ident, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$orig.ident)))
## cluster by sample
table_clusters_by_orig.idents <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::rename('cluster' = 'integrated_snn_res.0.4') %>%
  group_by(cluster, orig.ident) %>%
  summarize(count = n()) %>%
  spread(orig.ident, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$integrated_snn_res.0.4)))

### Percent of cells per cluster by sample
temp_labels <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
  group_by(orig.ident) %>%
  tally()

p1 <- table_orig.idents_by_clusters %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  # mutate(orig.ident = factor(orig.ident, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = orig.ident, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4)

p2 <- table_clusters_by_orig.idents %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$integrated_snn_res.0.4))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels, aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_fill_manual(name = 'orig.ident', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  'neurons/composition_samples_clusters_by_percent.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      mouse.snseq.combined.sct.neurons.recluster@meta.data$orig.ident %>% unique() %>% length(),
      mouse.snseq.combined.sct.neurons.recluster@meta.data$integrated_snn_res.0.4 %>% unique() %>% length()
    )),
  width = 18, height = 8
)


### Count of cells per cluster by orig.ident
temp_labels <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(orig.ident) %>%
  tally()

p1 <- table_orig.idents_by_clusters %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  # mutate(orig.ident = factor(orig.ident, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = orig.ident, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4)

p2 <- table_clusters_by_orig.idents %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  # mutate(cluster = factor(cluster, levels = levels(mouse.snseq.combined.sct.neurons.recluster@meta.data$integrated_snn_res.0.4))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_fill_manual(name = 'orig.ident', values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  'neurons/composition_orig.idents_clusters_by_number.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      mouse.snseq.combined.sct.neurons.recluster@meta.data$orig.ident %>% unique() %>% length(),
      mouse.snseq.combined.sct.neurons.recluster@meta.data$integrated_snn_res.0.4 %>% unique() %>% length()
    )),
  width = 18, height = 8
)

#### Build neuron hierarchical tree ####
# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

### create tree of cell types
mouse.snseq.combined.sct.neurons.recluster <- BuildClusterTree(
  mouse.snseq.combined.sct.neurons.recluster,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE,
  slot = 'data',
  assay = "SCT"
)

## create tree
neuron.tree <- mouse.snseq.combined.sct.neurons.recluster@tools$BuildClusterTree
# tree$tip.label <- paste0("Cluster ", tree$tip.label)

# # convert tree to tibble
# neuron.tree = as_tibble(neuron.tree)
# 
# ## add cell type data
# neuron.tree = full_join(neuron.tree,
#                         mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
#                    dplyr::select(integrated_snn_res.0.4,
#                                  parent_id.exp.prob) %>% 
#                    distinct(),
#                  by = c("label" = "integrated_snn_res.0.4"))
# 
# # convert back to tree
# neuron.tree = treeio::as.treedata(neuron.tree)

### graph tree of celltypes
ggtree::ggtree(neuron.tree, 
               aes(x, 
                   y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(
    # aes(color = parent_id.exp.prob),
                        shape = 16,
                        size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0),
                           'cm'))
ggsave('neurons/Neuron cluster tree.png',
       height = 6,
       width = 6)











#### Neuron marker specificity score ####
### Prep seurat object for running find markers 
## need to run before FindAllMarkers
# needs to recorrect again after filtering down to neurons? 
mouse.snseq.combined.sct.neurons.recluster = PrepSCTFindMarkers(mouse.snseq.combined.sct.neurons.recluster,
                                                                assay = "SCT")

### calculate specificity score with DEGs for each cluster
neuron.markers.df <- FindAllMarkers(mouse.snseq.combined.sct.neurons.recluster, 
                                    logfc.threshold = 0.1, # increase to increase speed
                                    test.use = 'wilcox',
                                    verbose = TRUE,
                                       min.pct = 0.10, 
                                       assay = 'SCT')
# create gene column
neuron.markers.df = neuron.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
neuron.markers.df = neuron.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))

## save data
write_csv(neuron.markers.df,
          file = 'neurons/neuron_seurat_cluster_markers.csv')


# # filter out poor markers
# neuron.markers.df.reduce = neuron.markers.df %>% 
#   filter(specificity > 0.87) %>% 
#   filter(p_val_adj < 0.05) 

### create heatmap of 6 top markers per cluster
## reduce to top 6 genes per cluster
neuron.markers.df.reduce = neuron.markers.df %>% 
  group_by(cluster) %>% 
  slice_max(specificity,
             n = 6) %>% 
  ungroup()

## create heatmap dataframe
## matrix 
# clusters as cols
# rows as genes
# neuron.markers.df.reduce.heatmap = neuron.markers.df %>% 
#   filter(gene %in% unique(neuron.markers.df.reduce$gene)) %>% 
#   select(cluster,
#          gene,
#          specificity) %>% 
#   # arrange(specificity) %>% 
#   arrange(cluster, 
#           -specificity) %>% 
#   pivot_wider(names_from = 'cluster',
#               values_from = 'specificity',
#               id_cols = 'gene') %>% 
#   column_to_rownames('gene') %>% 
#   as.matrix()

neuron.markers.df.reduce.heatmap = neuron.markers.df %>% 
  filter(gene %in% unique(neuron.markers.df.reduce$gene)) %>% 
  select(cluster,
         gene,
         specificity) %>% 
  left_join(neuron.markers.df.reduce %>% 
              select(cluster,
                     gene,
                     specificity) %>% 
              mutate(keep = cluster)) %>% 
  arrange(keep) %>%  
  pivot_wider(names_from = 'cluster',
              values_from = 'specificity',
              id_cols = 'gene') %>% 
  column_to_rownames('gene') %>% 
  as.matrix()
  
# replace NA with 0
neuron.markers.df.reduce.heatmap[is.na(neuron.markers.df.reduce.heatmap)] <- 0


# graph heatmap
neuron.markers.df.reduce.heatmap %>% 
  pheatmap::pheatmap(na_col = 'white',
                     breaks = c(-1, -0.1, 0.1, 1, 10, 65),
                     col = c('lightblue',
                             'white',
                             'yellow',
                             'orange',
                             'red'),
                     cluster_rows = F,
                     cluster_cols = F,
                     main = 'Neuron cluster specificity genes',
                     filename = 'neurons/neuron markers heatmap.png',
                     width = 10,
                     height = 10)


## possible to run this across nodes?


#### neuron cluster limmatrend ####
### create a loop for each neuron cluster 
# limmatrend
#set idents
Idents(object = mouse.snseq.combined.sct.neurons.recluster) <- "seurat_clusters"

## run across all neuron clusters
# remove cluster 12 since not enough for contrasts (it is postion 9 in list)
for (j in unique(mouse.snseq.combined.sct.neurons.recluster@meta.data$seurat_clusters)[-c(9)]) {
  
  # create new folder
  dir.create(paste0('./neurons/neuron_cluster/limmatrend/',
             j))
  
#subset to neuron clusters
tmp = subset(mouse.snseq.combined.sct.neurons.recluster,
                                             idents = j)

# subset with SCT data
DefaultAssay(tmp) = 'SCT'

### extract gene expression and orig.idents from neuron_cluster
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
tmp.var <- FindVariableFeatures(tmp, 
                                                                assay = 'integrated',
                                                                selection.method = "vst", 
                                                                verbose = F)
# identify top variable genes
neuron_cluster.reduce.group.topgenes.prep <- VariableFeatures(tmp.var,
                                                          assay = 'integrated')

# create dummy
tmp.expression = full_join(full_join(tmp@reductions$umap@cell.embeddings %>% 
                                                                       as.data.frame() %>% 
                                                                       rownames_to_column("Cell.id"),
                                                                     tmp@meta.data %>%
                                                                       rownames_to_column("Cell.id")),
                                                           tmp@assays$SCT@data %>% 
                                                             as.data.frame()  %>% 
                                                             t() %>% as.data.frame() %>% 
                                                             rownames_to_column('Cell.id'))


## create vector of factor
neuron_cluster.reduce.DomvsSub.vector.list.sct.prep = tmp.expression %>% 
  mutate(neuron_cluster.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(neuron_cluster.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
neuron_cluster.reduce.vector.count.sct.prep = GetAssayData(tmp,
                                                       assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  tmp.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% neuron_cluster.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
neuron_cluster.reduce.vector.limma.sct.prep = list(count = neuron_cluster.reduce.vector.count.sct.prep,
                                               condt = neuron_cluster.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_neuron_cluster <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
                               MvsFD = treatMale.Dom - treatFemale.Dom,
                               MvsFS = treatMale.Sub - treatFemale.Sub,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
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
  pdf(file= paste0('./neurons/neuron_cluster/limmatrend/',
                   j,
                   '/',
                   j,
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
  # pdf(file= "./neuron_cluster/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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


###run function  
neuron_cluster.limma.results.sct = run_limmatrend_neuron_cluster(neuron_cluster.reduce.vector.limma.sct.prep)

# save results to dataframe
neuron_cluster.limma.results.sct.df = full_join(neuron_cluster.limma.results.sct$ttDvsS %>% 
                                              rename_with(~paste0(.,"_DvsS")) %>% 
                                              rownames_to_column("Gene"),
                                            neuron_cluster.limma.results.sct$ttDvsSM %>% 
                                              rename_with(~paste0(.,"_DvsS_M")) %>% 
                                              rownames_to_column("Gene")) %>% 
  full_join(neuron_cluster.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron_cluster.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron_cluster.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron_cluster.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
neuron_cluster.limma.results.sct.df = neuron_cluster.limma.results.sct.df %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS)) %>% 
  mutate(Sig_DvsS_M = ifelse(adj.P.Val_DvsS_M <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_M = ifelse(logFC_DvsS_M > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_M = ifelse(Sig_DvsS_M == 'Sig',
                                       Direction.type_DvsS_M,
                                       Sig_DvsS_M)) %>% 
  mutate(Sig_DvsS_F = ifelse(adj.P.Val_DvsS_F <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_F = ifelse(logFC_DvsS_F > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_F = ifelse(Sig_DvsS_F == 'Sig',
                                       Direction.type_DvsS_F,
                                       Sig_DvsS_F)) %>% 
  mutate(Sig_MvsF = ifelse(adj.P.Val_MvsF <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_MvsF = ifelse(logFC_MvsF > 0,
                                      'up',
                                      'down'),
         Sig.direction_MvsF = ifelse(Sig_MvsF == 'Sig',
                                     Direction.type_MvsF,
                                     Sig_MvsF)) %>% 
  mutate(Sig_MvsFD = ifelse(adj.P.Val_MvsFD <= 0.05,
                            'Sig',
                            'Not Sig'),
         Direction.type_MvsFD = ifelse(logFC_MvsFD > 0,
                                       'up',
                                       'down'),
         Sig.direction_MvsFD = ifelse(Sig_MvsFD == 'Sig',
                                      Direction.type_MvsFD,
                                      Sig_MvsFD)) %>% 
  mutate(Sig_MvsFS = ifelse(adj.P.Val_MvsFS <= 0.05,
                            'Sig',
                            'Not Sig'),
         Direction.type_MvsFS = ifelse(logFC_MvsFS > 0,
                                       'up',
                                       'down'),
         Sig.direction_MvsFS = ifelse(Sig_MvsFS == 'Sig',
                                      Direction.type_MvsFS,
                                      Sig_MvsFS)) 

# save data
write_csv(neuron_cluster.limma.results.sct.df,
     file = paste0('./neurons/neuron_cluster/limmatrend/',
     j,
     '/',
     j,
     '_neuropeptides.neuron_cluster.limma.results.sct.df.csv'))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF",
                          "MvsFD",
                          "MvsFS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  neuron_cluster.limma.results.sct.df %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('./neurons/neuron_cluster/limmatrend/',
                j,
                '/',
                j,
                ' limma.neuron_cluster.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
}


#### Red ribbon functions ####
# https://github.com/antpiron/RedRibbon

# ## libraries 
# library(RedRibbon)
# library(tidyverse)

### This melt.matrix function needed to be set up separetly for some reason
melt.matrix <- function (data, ..., na.rm = FALSE, value.name = "value")
{
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

### add value to set as max log scale for colors
#' Render the RRHO map.
#' 
#' You can choose to render the RedRibbon level map using \code{ggplot2}. 
#' 
#' @param self is the RedRibbon object
#' @param n is the number of coordinates to compute on the x and y axis (Default = sqrt(len))
#' @param labels is a list of labels of list a and b. Default is c("a", "b").
#' @param show.quadrants is a flag if set the quadrants lines are drawn
#' @param quadrants is the object returned by 'quadrants()' method
#' @param show.pval is a flag to show the P-values on the plot
#' @param repel.force is the value of the repel force for the p-value ploting (default = 150)
#' @param base_size is the size of the text fields (default = 20)
#' @param .log10 output log10 pval (default = FALSE)
#' @param ... The rest of the parameters
#' 
#' @return A \code{ggplot} object.
#' 
#' @method ggRedRibbon rrho
#' @export
ggRedRibbon.rrho.scale <- function (self, n = NULL, labels = c("a", "b"), show.quadrants=TRUE, quadrants=NULL, 
                                    show.pval=TRUE,
                                    repel.force=150, base_size=20, .log10=FALSE,
                                    new.max.log, # add value
                                    ...)
{
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
    ggplot2::theme(axis.title = ggplot2::element_text(size=base_size,face="bold"),
                   legend.title = ggplot2::element_text(size = base_size * 7 / 10),
                   legend.text = ggplot2::element_text(size = base_size * 1 / 2),
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=base_size* 7 / 10, face="bold"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=base_size * 7 / 10, face="bold", angle=90),
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

#### create function to output redribbon graph and save quadrant data
### just input both gene lists
#' @param celltype is the title of the analysis
#' @param dataset.a is the DEG gene list with one column the gene names and the second column the logFC or p.value results with positive and negative indicating directionality 
#' @param dataset.b is the same as dataset.a
#' @param dataset.a.type is the name given to the 'a' gene list analysis
#' @param dataset.b.type is the same as @param dataset.a.type
#' @param a.variable is the variable name with the values from dataset.a
#' @param b.variable is the same as a.variable
#' @param new.max.log is the value to set the max and min for the log scale in the heatmap graph. Can be left 'NULL' if scale should be determined by the redribbon function. The new.max.log needs to be at least the value of the log scale from the redribbon function. Useful if you want to set the scale to the max value across mutlitple redribbon objects
#' @param file.name is the folder file path where the figure and quadrant data should be saved
RedRibbon.all <- function(celltype = "Analysis.title",
                          dataset.a = gene.list.a,
                          dataset.b = gene.list.b,
                          dataset.a.type = "Name.a",
                          dataset.b.type = "Name.b",
                          a.variable = "Value.name.a",
                          b.variable = "Value.name.b",
                          new.max.log = NULL,
                          file.name = './file.name/')
{
  # create data frame for red ribbon
  # needs to have an id (gene) col and one called 'a' and one called 'b'
  df = dataset.a %>% 
    rename('a' =  a.variable) %>% 
    full_join(dataset.b %>% 
                rename('b' = b.variable))
  
  ## Create RedRibbon object
  rr <- RedRibbon(df, 
                  enrichment_mode="hyper-two-tailed")
  
  ## Run the overlap using evolutionnary algorithm,
  ## computing permutation adjusted p-value for the four quadrants
  quad <- quadrants(rr, 
                    algorithm="ea",
                    permutation=TRUE, 
                    whole=FALSE)
  
  ### compare RRHO2 to Redribbon
  # create list of RRHO outcomes
  # add NA to deal with empty quadrants
  RR.list <- data.frame(Gene = c(df[quad$upup$positions,]$gene,
                                 NA),
                        RRquadrant = 'upup') %>% 
    rbind(data.frame(Gene = c(df[quad$downdown$positions,]$gene,
                              NA),
                     RRquadrant = 'downdown')) %>% 
    rbind(data.frame(Gene = c(df[quad$updown$positions,]$gene,
                              NA),
                     RRquadrant = 'updown')) %>% 
    rbind(data.frame(Gene = c(df[quad$downup$positions,]$gene,
                              NA),
                     RRquadrant = 'downup')) %>% 
    mutate(Sample = celltype) %>% 
    na.omit()
  
  # save file
  write_csv(RR.list,
            file = paste0(file.name,
                          celltype,
                          ' quadrant genes.csv'))
  ## Plots the results
  ggRedRibbon(rr,
              quadrants=quad) + 
    coord_fixed(ratio = 1, 
                clip = "off") +
    xlab(dataset.a) +
    ylab(dataset.b) +
    ggtitle(celltype)
  
  if (is.null(new.max.log))
  {
    ggRedRibbon(rr, 
                quadrants=quad) + 
      coord_fixed(ratio = 1, 
                  clip = "off") +
      xlab(dataset.a.type) +
      ylab(dataset.b.type) +
      ggtitle(celltype)
    ggsave(paste0(file.name,
                  celltype,
                  ' RedRibbon.png'),
           height = 10,
           width = 10)
  } else
  {
    # scaled
    ggRedRibbon.rrho.scale(rr, 
                           quadrants=quad,
                           new.max.log = new.max.log) + 
      coord_fixed(ratio = 1, 
                  clip = "off") +
      xlab(dataset.a.type) +
      ylab(dataset.b.type) +
      ggtitle(celltype)
    ggsave(paste0(file.name,
                  celltype,
                  ' RedRibbon_scaled.png'),
           height = 10,
           width = 10)
  }
  
  
}

#### RedRibbon neuron clusters ####
### loop through all neuron clusters
# remove 12
for (j in c(0:11))
{
  ## load data 
  tmp = read.csv(paste0('./neurons/neuron_cluster/limmatrend/',
                              j,
                              '/',
                              j,
                              '_neuropeptides.neuron_cluster.limma.results.sct.df.csv'))

### compare dom vs sub across sexes
##male data
neuron_cluster.DvsS.M.rrho2 = data.frame(gene = tmp$Gene,
                                 value = -log10(tmp$P.Value_DvsS_M),
                                 direction = tmp$Direction.type_DvsS_M, 
                                 stringsAsFactors = FALSE)
#set positive negative
neuron_cluster.DvsS.M.rrho2 = neuron_cluster.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
neuron_cluster.DvsS.F.rrho2 = data.frame(gene = tmp$Gene,
                                 value = -log10(tmp$P.Value_DvsS_F),
                                 direction = tmp$Direction.type_DvsS_F,
                                 stringsAsFactors = FALSE)
#set positive negative
neuron_cluster.DvsS.F.rrho2 = neuron_cluster.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

#### RedRibbon neuron_cluster clusters ####
## neuron_cluster
RedRibbon.all(celltype = paste0("neuron_cluster_",
                                j),
              dataset.a = neuron_cluster.DvsS.M.rrho2,
              dataset.b = neuron_cluster.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = './neurons/neuron_cluster/limmatrend/RedRibbon/')
}

## create counts per quadrant for each one
for (j in c(0:11))
{
  ## load data 
  tmp = read.csv(paste0('./neurons/neuron_cluster/limmatrend/RedRibbon/neuron_cluster_',
                        j,
                        ' quadrant genes.csv'))
  
  # graph count of genes per quadrant
  tmp %>% 
    group_by(RRquadrant) %>% 
    summarise(Total = n()) %>% 
    as.data.frame() %>% 
    ggplot(aes(x = reorder(RRquadrant,
                           -Total),
               y = Total,
               label = Total)) +
    geom_label() +
    theme_classic() +
    ylab('Number of genes per quadrant') +
    xlab('') +
    ggtitle(paste0(j,
                   ' neuron cluster RR quadrant'))
  ggsave(paste0('./neurons/neuron_cluster/limmatrend/RedRibbon/neuron_cluster_',
         j,
         ' _rr_number_gene_quad.png'))
  
}


#### neuron WGCNA ####
#gene_select parameter:

#   variable: use the genes stored in the Seurat object’s VariableFeatures.
# fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# custom: use genes that are specified in a custom list.

# mouse.snseq.combined.sct.neurons.recluster <- SetupForWGCNA(
#   mouse.snseq.combined.sct.neurons.recluster,
#   gene_select = "variable", # the gene selection approach
#   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# )


# create new variable features list for neuron data
mouse.snseq.combined.sct.neurons.recluster <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.recluster, 
                             selection.method = "vst", 
                             nfeatures = 3000, 
                             verbose = FALSE,
                             assay = 'integrated')

# gets only 2221 genes

## set up object for analysis
mouse.snseq.combined.sct.neurons.recluster = SetupForWGCNA(mouse.snseq.combined.sct.neurons.recluster,
                                                             features = VariableFeatures(mouse.snseq.combined.sct.neurons.recluster),
                                                             wgcna_name = "neurons"
)


# construct metacells  in each group
mouse.snseq.combined.sct.neurons.recluster <- MetacellsByGroups(
  seurat_obj = mouse.snseq.combined.sct.neurons.recluster,
  group.by = c("Genotype",
               "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident', # set the Idents of the metacell seurat object
  slot = "data", #set for counts
  assay = "SCT", #set assay type
  mode = "average", #metacell expression profile determined by average of constituent single cells
  min_cells = 50 #set minimum number of cells per metacell
)

## setup expression matrix
#use metacell expression data
## setup expression matrix
#use metacell expression data
mouse.snseq.combined.sct.neurons.recluster <- SetDatExpr(
  mouse.snseq.combined.sct.neurons.recluster,
  assay = 'SCT', # using SCT assay
  slot = 'data' # using normalized data
)

##select soft-power threshold
mouse.snseq.combined.sct.neurons.recluster <- TestSoftPowers(
  mouse.snseq.combined.sct.neurons.recluster,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.neuron <- PlotSoftPowers(mouse.snseq.combined.sct.neurons.recluster)

#assemble with patchwork
png('./neurons/hdWGCNA/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron,
           ncol=2)
dev.off()

#use soft power threshold 7


### construct co-expression network
# construct co-expression network:
mouse.snseq.combined.sct.neurons.recluster <- ConstructNetwork(mouse.snseq.combined.sct.neurons.recluster, 
                                                                 soft_power=7,
                                                                 # setDatExpr=FALSE,
                                                                 tom_name = 'neuron', # name of the topoligical overlap matrix written to disk
                                                                 overwrite_tom = TRUE
)

#graph dendrogram
png('./neurons/hdWGCNA/WGCNA dendrogram hypomap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(mouse.snseq.combined.sct.neurons.recluster,
               main='neurons hdWGCNA Dendrogram')
dev.off()

#graph dendrogram for poster
# pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
#     width = 10,
#     height = 2.5)
# PlotDendrogram(mouse.snseq.combined.sct.neurons.recluster,
#                main='Neurons hdWGCNA Dendrogram')
# dev.off()
# ended up exporting the right size with plot viewer

##compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
mouse.snseq.combined.sct.neurons.recluster <- ScaleData(mouse.snseq.combined.sct.neurons.recluster,
                                                              features=VariableFeatures(mouse.snseq.combined.sct.neurons.recluster))

# compute all MEs in the full single-cell dataset
mouse.snseq.combined.sct.neurons.recluster <- ModuleEigengenes(
  mouse.snseq.combined.sct.neurons.recluster,
  exclude_grey = TRUE, 
  group.by.vars = 'Genotype'
)

# module eigengenes:
MEs.neurons <- GetMEs(mouse.snseq.combined.sct.neurons.recluster,
                      harmonized=FALSE)


# compute eigengene-based connectivity (kME):
mouse.snseq.combined.sct.neurons.recluster <- ModuleConnectivity(
  mouse.snseq.combined.sct.neurons.recluster)


#### Graph WGCNNA for neurons ####

# plot genes ranked by kME for each module
png('./neurons/hdWGCNA/Plot kMEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotKMEs(mouse.snseq.combined.sct.neurons.recluster,
         ncol=2)
dev.off()

# correlation between modules
png('./neurons/hdWGCNA/Correlogram modules neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleCorrelogram(mouse.snseq.combined.sct.neurons.recluster,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001),
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.neurons.me <- ModuleFeaturePlot(
  mouse.snseq.combined.sct.neurons.recluster,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('./neurons/hdWGCNA/UMAP MEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neurons.me, ncol=2)
dev.off()



# for poster
# stitch together with patchwork
# pdf('./neurons/hdWGCNA/UMAP MEs neurons poster.pdf',
#     width = 10,
#     height = 5)
# wrap_plots(plot_list.neurons.me, ncol=4)
# dev.off()

## create one umap plot
# combine umap data and meta data
# neuron.umap = mouse.snseq.combined.sct.neurons.recluster@reductions$umap@cell.embeddings %>% 
#   as.data.frame() %>% 
#   rownames_to_column('Cell.id') %>% 
#   full_join(mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
#               as.data.frame() %>% 
#               rownames_to_column('Cell.id'))

# add cluster_color for each cell
# neuron.umap.long = neuron.umap %>% 
#   select(c(Cell.id,
#            UMAP_1,
#            UMAP_2,
#            orig.ident,
#            Genotype.id,
#            Cell.type,
#            integrated_snn_res.0.8,
#            turquoise,
#            green,
#            yellow,
#            red,
#            pink,
#            blue,
#            brown,
#            black,
#            grey)) %>% 
#   pivot_longer(cols = c(turquoise,
#                         green,
#                         yellow,
#                         red,
#                         pink,
#                         blue,
#                         brown,
#                         black,
#                         grey),
#                names_to = 'module',
#                values_to = 'ME_score') %>% 
#   group_by(Cell.id) %>% 
#   mutate(max.ME_score = max(ME_score),
#          second.max.ME_score = max(ME_score[ME_score != max(ME_score)]),
#          max.module = ifelse(max.ME_score == ME_score,
#                              module,
#                              NA),
#          second.max.module = ifelse(second.max.ME_score == ME_score,
#                              module,
#                              NA),
#          ratio.max.mes = second.max.ME_score/max.ME_score,
#          keep = case_when(
#            ratio.max.mes < 0 ~ max.module,
#            ratio.max.mes < 0.75 ~ max.module,
#            ratio.max.mes > 1  ~ max.module,
#            TRUE  ~ paste(max.module,
#                          second.max.module,
#                          sep = ':'),
#          )) %>% 
#   filter(!is.na(max.module)) 
# 
# # graph umap
# neuron.umap.long %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = max.module)) +
#   geom_point() +
#   theme_classic()


# get mods from object
mods.neurons <- colnames(MEs.neurons); mods.neurons <- mods.neurons[mods.neurons != 'grey']

# add MEs to Seurat meta-data:
mouse.snseq.combined.sct.neurons.recluster@meta.data <- cbind(mouse.snseq.combined.sct.neurons.recluster@meta.data, 
                                                                MEs.neurons)

##neuropeptides
modules.neurons <- GetModules(mouse.snseq.combined.sct.neurons.recluster)

#heatmap of kme and neuropeptide
neuropeptides.df = data.frame(gene_name = neuropeptides.list %>% 
                                filter(!is.na(Gene.name.nile.tilapia)) %>% 
                                pull(Gene.name.nile.tilapia) %>% 
                                unique(),
                              neuropeptide = TRUE)

#combine modules with neuropeptide
modules.neurons = modules.neurons %>% 
  full_join(neuropeptides.df)

# #convert to long
# modules.long.neurons = modules.neurons %>% 
#   pivot_longer(kME_turquoise:kME_red ,
#                names_to = 'module.kME',
#                values_to = 'kME')
# 
# 
# #graph
# modules.long.neurons %>% 
#   filter(neuropeptide == TRUE) %>% 
#   drop_na() %>%
#   mutate(module.fct=as.integer(module)) %>% 
#   ggplot(aes(x = module.kME,
#              y = reorder(gene_name,
#                          module.fct),
#              label = round(kME,
#                            2))) +
#   geom_tile(aes(fill = kME)) +
#   geom_text() +
#   geom_point(aes(x=-Inf,
#                  color = module),
#              size = 5) +
#   scale_fill_gradientn(colours=c("blue",
#                                  "white",
#                                  "red"), 
#                        limits = c(-0.6, 0.6)) +
#   theme_classic() +
#   ylab('Neuropeptides') +
#   scale_color_manual(values = c("blue"="blue",
#                                 "grey"="grey",
#                                 "turquoise"="turquoise",
#                                 "brown"="brown",
#                                 "red"="red",
#                                 "green"="green",
#                                 "black"="black",
#                                 "yellow"="yellow")) +
#   coord_cartesian(clip = 'off')+ 
#   theme(axis.ticks.y = element_blank())
# ggsave('./neurons/hdWGCNA/Module and neuropeptide heatmap.png',
#        width = 10,
#        height = 10)

## create dotplot of modules ME and clusters
DotPlot(mouse.snseq.combined.sct.neurons.recluster,
        features = c("blue", 
                     "red", 
                     "yellow", 
                     "green",
                     "turquoise",
                     "brown"),
        cols = c('grey',
                 'red'),
        group.by = "integrated_snn_res.0.8",
        col.min = 2,
        dot.min = .5,
        scale = FALSE) +
  RotatedAxis() +
  theme_bw()
ggsave('./neurons/hdWGCNA/ME by cluster dotplot.png',
       width = 10,
       height = 10)

## get dot plot data
# relabel column 
mouse.neuron.wgcna.module.cluster = DotPlot.data(mouse.snseq.combined.sct.neurons.recluster,
                                                   features = c("blue", 
                                                                "red", 
                                                                "yellow", 
                                                                "green",
                                                                "turquoise",
                                                                "brown"),
                                                   cols = c('grey',
                                                            'red'),
                                                   group.by = "integrated_snn_res.0.8",
                                                   col.min = 2,
                                                   dot.min = .5,
                                                   scale = FALSE)


## trait correlation
Idents(mouse.snseq.combined.sct.neurons.recluster) <- mouse.snseq.combined.sct.neurons.recluster$orig.ident

# convert orig.ident to factor
mouse.snseq.combined.sct.neurons.recluster$orig.ident.fct <- as.factor(mouse.snseq.combined.sct.neurons.recluster$orig.ident)

# list of traits to correlate
cur_traits.neurons <- c('orig.ident.fct', 
                        'nCount_SCT', 
                        'nFeature_SCT')

mouse.snseq.combined.sct.neurons.recluster <- ModuleTraitCorrelation(
  mouse.snseq.combined.sct.neurons.recluster,
  traits = cur_traits.neurons,
  features = 'MEs'
)

#plot
png('./neurons/hdWGCNA/Module and trait correlation heatmap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotModuleTraitCorrelation(
  mouse.snseq.combined.sct.neurons.recluster,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 5,
  text_digits = 5,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'white',
  low_color = 'blue',
  plot_max = 0.5,
  combine=TRUE
)
dev.off()

#### Network WGCNA for neurons ####
# visualize network with UMAP
mouse.snseq.combined.sct.neurons.recluster <- RunModuleUMAP(
  mouse.snseq.combined.sct.neurons.recluster,
  n_hubs = 5, #number of hub genes to include for the umap embedding
  n_neighbors=15, #neighbors parameter for umap
  min_dist=0.3, #min distance between points in umap space
  spread=5
)


# compute cell-type marker genes with Seurat:
#set idents to cluster
Idents(mouse.snseq.combined.sct.neurons.recluster) <- mouse.snseq.combined.sct.neurons.recluster$integrated_snn_res.0.8

#calculate marker genes per cluster
# with 3000 variable genes
markers.neuron <- Seurat::FindAllMarkers(
  mouse.snseq.combined.sct.neurons.recluster,
  only.pos = TRUE,
  logfc.threshold=1,
  features = VariableFeatures(mouse.snseq.combined.sct.neurons.recluster)
)

# compute marker gene overlaps
mouse.snseq.combined.sct.neurons.recluster.overlap_df <- OverlapModulesDEGs(
  mouse.snseq.combined.sct.neurons.recluster,
  deg_df = markers.neuron,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cluster
plot_list.neuron.overlap <- OverlapBarPlot(mouse.snseq.combined.sct.neurons.recluster.overlap_df)

# stitch plots with patchwork
png('./neurons/hdWGCNA/Vlnplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron.overlap, 
           ncol=4)
dev.off()

# plot odds ratio of the overlap as a dot plot
png('./neurons/hdWGCNA/Dotplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
OverlapDotPlot(mouse.snseq.combined.sct.neurons.recluster.overlap_df,
               plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cluster markers')
dev.off()

##graph network
#umap
png('./neurons/hdWGCNA/gene_network/Gene network UMAP.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleUMAPPlot(mouse.snseq.combined.sct.neurons.recluster,
               edge.alpha=0.5,
               sample_edges=TRUE,
               keep_grey_edges=FALSE,
               edge_prop=0.075, # taking the top 20% strongest edges in each module
               label_hubs=10 # how many hub genes to plot per module?
)
dev.off()

# for poster
##graph network
#umap
# pdf('./neurons/hdWGCNA/gene_network/Gene network UMAP poster.pdf',
#     width = 10,
#     height = 10)
# ModuleUMAPPlot.size(mouse.snseq.combined.sct.neurons.recluster,
#                edge.alpha= 0.75,
#                sample_edges=TRUE,
#                keep_grey_edges=FALSE,
#                edge_prop=0.075, # taking the top 20% strongest edges in each module
#                label_hubs=0, # how many hub genes to plot per module?
#                vertex.label.cex = 1,
#                dot.size = 8,
#                edge.size = 4
# )
# dev.off()

# module specific network
ModuleNetworkPlot(mouse.snseq.combined.sct.neurons.recluster,
                  outdir = './neurons/hdWGCNA/gene_network/')

# hubgene network
png('./neurons/hdWGCNA/gene_network/Hubgene network.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
HubGeneNetworkPlot(mouse.snseq.combined.sct.neurons.recluster,
                   n_hubs = 5, 
                   n_other=50,
                   edge_prop = 0.75,
                   mods = 'all'
)
dev.off()






#### DMEs WGCNA for neurons ####
## Differential module
# get list of dom and sub cells
group1.neurons.dom <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'dom_mouse_snseq') %>%
  rownames
group2.neurons.sub <- mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'sub_mouse_snseq') %>% 
  rownames

# calculate differential module eigengene
DMEs.neurons <- FindDMEs(
  mouse.snseq.combined.sct.neurons.recluster,
  barcodes1 = group1.neurons.dom,
  barcodes2 = group2.neurons.sub,
  test.use='wilcox',
  harmonized = FALSE
)

# graph DME
PlotDMEsVolcano(
  mouse.snseq.combined.sct.neurons.recluster,
  DMEs.neurons) +
  theme_classic() +
  xlim(-0.5, 0.5)
ggsave('./neurons/hdWGCNA/ME by cluster volcanoplot.png',
       width = 10,
       height = 10)

# for poster
# create mod color list
mod_colors <- DMEs.neurons$module
names(mod_colors) <- as.character(DMEs.neurons$module)
# graph DME
DMEs.neurons %>%
  mutate(anno = ifelse(p_val_adj < 0.05,
                       module,
                       "")) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj)))  +
  geom_rect(aes(xmin = -Inf,
                xmax = Inf,
                ymin = -Inf,
                ymax = -log10(0.05)),
            fill = "grey75",
            alpha = 0.8,
            color = NA) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey75",
             alpha = 0.8) +
  geom_point(aes(fill = module),
             size = 8,
             pch = 21,
             color = "black") +
  geom_text_repel(aes(label = anno),
                  color = "black",
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 8,
                  point.padding = 4,
                  nudge_x = 0.06) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlim(c(-0.45,
         0.45)) +
  scale_fill_manual(values = mod_colors)+
  xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
  ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
  ggtitle('Differential module eigengene analysis') +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('./neurons/hdWGCNA/ME by cluster volcanoplot poster.pdf',
       width = 10,
       height = 5,
       units = "in",
       dpi = 320)

### graph MEs for relevant clusters
for (i in mods.neurons) {
  #create list of clusters with module eigengene above:
  # average scaled expression of 1 and 50% percent expression
  cluster.list = mouse.neuron.wgcna.module.cluster %>% 
    filter(features.plot == i) %>% 
    mutate(keep = ifelse(avg.exp.scaled >= 1 & pct.exp >= 50,
                         "keep",
                         "remove")) %>% 
    filter(keep == "keep") %>%
    mutate(id = as.character(id)) %>% 
    pull(id) 
  
  ##graph ME per cluster across social status
  mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    droplevels() %>%
    ggplot(aes_string(x = 'integrated_snn_res.0.8',
                      y = i,
                      color = 'orig.ident')) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.25,
                                             jitter.height = 0),
               alpha = 0.2)+
    theme_classic() +
    ggtitle(i)
  ggsave(paste('./neurons/hdWGCNA/modules/',
               i,
               ' ME by cluster boxplot.png',
               sep = ''),
         width = 10,
         height = 10)
  
  ## ME across social status
  mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    droplevels() %>%
    ggplot(aes_string(x = 'orig.ident',
                      y = i)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter(width = 0.25,
                height = 0)+
    theme_classic() +
    ggtitle(i)
  ggsave(paste('./neurons/hdWGCNA/modules/',
               i,
               ' ME violinplot.png',
               sep = ''),
         width = 10,
         height = 10)
  
  # cluster bias 
  mouse.snseq.combined.sct.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    mutate(keep = ifelse(get(i) >= 1,
                         "present",
                         "absence")) %>% 
    droplevels() %>%
    select(c(cluster,
             keep,
             orig.ident)) %>% 
    table() %>% 
    as.data.frame() %>% 
    pivot_wider(id_cols = c(cluster,
                            keep),
                names_from = orig.ident,
                values_from = Freq) %>% 
    filter(keep == "present") %>% 
    mutate(dom_mouse_snseq.scaled = 0.783645 * dom_mouse_snseq) %>% 
    ggplot(aes(x = dom_mouse_snseq.scaled,
               y = sub_mouse_snseq,
               label = cluster)) +
    geom_abline(slope = 1,
                intercept = 0) +
    geom_text() +
    theme_classic() +
    ggtitle(paste(i,
                  'sub vs dom scaled counts per cluster'))
  ggsave(paste('./neurons/hdWGCNA/modules/',
               i,
               ' ME count bias.png',
               sep = ''),
         width = 10,
         height = 10)
  
}


#### male neuron WGCNA ####
## subset males
mouse.snseq.combined.sct.neurons.recluster.males = mouse.snseq.combined.sct.neurons.recluster

#set idents
  Idents(object = mouse.snseq.combined.sct.neurons.recluster.males) <- "orig.ident"

#subset to neurons
  mouse.snseq.combined.sct.neurons.recluster.males = subset(mouse.snseq.combined.sct.neurons.recluster.males,
                                          idents = c("Male.Dom",
                                                     "Male.Sub"))
  
#set idents
Idents(object = mouse.snseq.combined.sct.neurons.recluster.males) <- "seurat_clusters"

#gene_select parameter:

#   variable: use the genes stored in the Seurat object’s VariableFeatures.
# fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# custom: use genes that are specified in a custom list.

# mouse.snseq.combined.sct.neurons.recluster <- SetupForWGCNA(
#   mouse.snseq.combined.sct.neurons.recluster,
#   gene_select = "variable", # the gene selection approach
#   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# )


mouse.snseq.combined.sct.neurons.recluster.males = SetupForWGCNA(mouse.snseq.combined.sct.neurons.recluster.males,
                                                           features = VariableFeatures(mouse.snseq.combined.sct.neurons.recluster.males),
                                                           wgcna_name = "neurons.males"
)

# construct metacells  in each group
mouse.snseq.combined.sct.neurons.recluster.males <- MetacellsByGroups(
  seurat_obj = mouse.snseq.combined.sct.neurons.recluster.males,
  group.by = c("orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident', # set the Idents of the metacell seurat object
  slot = "data", #set for counts
  assay = "SCT", #set assay type
  mode = "average", #metacell expression profile determined by average of constituent single cells
  min_cells = 100 #set minimum number of cells per metacell
)


## setup expression matrix
#use metacell expression data

##select soft-power threshold
mouse.snseq.combined.sct.neurons.recluster.males <- TestSoftPowers(
  mouse.snseq.combined.sct.neurons.recluster.males,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.neuron.males <- PlotSoftPowers(mouse.snseq.combined.sct.neurons.recluster.males)

#assemble with patchwork
png('./neurons/hdWGCNA/males/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron.males,
           ncol=2)
dev.off()

#use soft power threshold 9


### construct co-expression network
# construct co-expression network:
mouse.snseq.combined.sct.neurons.recluster.males <- ConstructNetwork(mouse.snseq.combined.sct.neurons.recluster.males, 
                                                               soft_power=9,
                                                               # setDatExpr=FALSE,
                                                               tom_name = 'neuron', # name of the topoligical overlap matrix written to disk
                                                               overwrite_tom = TRUE
)

#graph dendrogram
png('./neurons/hdWGCNA/males/WGCNA dendrogram hypomap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(mouse.snseq.combined.sct.neurons.recluster.males,
               main='neurons males hdWGCNA Dendrogram')
dev.off()

#graph dendrogram for poster
# pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
#     width = 10,
#     height = 2.5)
# PlotDendrogram(mouse.snseq.combined.sct.neurons.recluster,
#                main='Neurons hdWGCNA Dendrogram')
# dev.off()
# ended up exporting the right size with plot viewer

##compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
# mouse.snseq.combined.sct.neurons.recluster <- ScaleData(mouse.snseq.combined.sct.neurons.recluster, 
#                                                               features=VariableFeatures(mouse.snseq.combined.sct.neurons.recluster))

# compute all MEs in the full single-cell dataset

mouse.snseq.combined.sct.neurons.recluster.males <- ModuleEigengenes(
  mouse.snseq.combined.sct.neurons.recluster.males,
  exclude_grey = TRUE
  # , group.by.vars = 'Genotype.id'
)

# module eigengenes:
MEs.neurons.males <- GetMEs(mouse.snseq.combined.sct.neurons.recluster.males,
                      harmonized=FALSE)


# compute eigengene-based connectivity (kME):
mouse.snseq.combined.sct.neurons.recluster.males <- ModuleConnectivity(
  mouse.snseq.combined.sct.neurons.recluster.males)


#### Graph WGCNNA for neurons males ####

# plot genes ranked by kME for each module
png('./neurons/hdWGCNA/males/Plot kMEs neurons males.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotKMEs(mouse.snseq.combined.sct.neurons.recluster.males, 
         ncol=2)
dev.off()

# correlation between modules
png('./neurons/hdWGCNA/males/Correlogram modules neurons males.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleCorrelogram(mouse.snseq.combined.sct.neurons.recluster.males,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001), 
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.neurons.males.me <- ModuleFeaturePlot(
  mouse.snseq.combined.sct.neurons.recluster.males,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('./neurons/hdWGCNA/males/UMAP MEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neurons.males.me, ncol=2)
dev.off()




#### AVP and oxt analysis ####
### graph avp and oxt
VlnPlot(mouse.snseq.combined.sct.neurons.recluster,
        features = c('Avp',
                     'Oxt'),
        assay = 'SCT',
        slot = 'data',
        group.by = 'orig.ident')

### get meta data and expression data for oxt/avp cells
mouse.snseq.combined.sct.neurons.recluster.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                              as.data.frame() %>% 
                                                                              rownames_to_column("Cell.id"),
                                                                            mouse.snseq.combined.sct.neurons.recluster@meta.data %>%
                                                                              rownames_to_column("Cell.id")),
                                                                  mouse.snseq.combined.sct.neurons.recluster@assays$SCT@data %>% 
                                                                    as.data.frame() %>% 
                                                                    filter(rownames(mouse.snseq.combined.sct.neurons.recluster@assays$SCT@data) %in% c('Avp',
                                                                                                                                                       'Oxt')) %>% 
                                                                    t() %>% as.data.frame() %>% 
                                                                    rownames_to_column('Cell.id'))


## assign AVP and oxt cell type cell type
mouse.snseq.combined.sct.neurons.recluster.expression = mouse.snseq.combined.sct.neurons.recluster.expression %>% 
  mutate(AVP.oxt.cell = case_when(Avp > 1 & Oxt > 1 ~ "Both",
                                  Avp > 1 & Oxt < 1 ~ "Avp",
                                  Avp < 1 & Oxt > 1 ~ "Oxt",
                                  Avp < 1 & Oxt < 1 ~ "None",
                                  TRUE ~ 'NA'))


## graph results
mouse.snseq.combined.sct.neurons.recluster.expression %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>% 
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  ggplot(aes(x = orig.ident,
             y = Percent,
             fill = AVP.oxt.cell,
             group = AVP.oxt.cell)) +
  geom_bar(position="dodge",
           stat = 'identity') +
  labs(fill = 'Neuron type') 
ggsave('./neurons/neuropeptides/Percentage neurons expressing AVP and OXT.png',
       height = 10,
       width = 10)

### just avp cells
# graph results
mouse.snseq.combined.sct.neurons.recluster.expression %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>%
  filter(AVP.oxt.cell == c('Avp', 'Both')) %>%  
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  ggplot(aes(x = orig.ident,
             y = Percent,
             fill = AVP.oxt.cell,
             group = AVP.oxt.cell)) +
  geom_bar(position="dodge",
           stat = 'identity') +
  labs(fill = 'AVP neuron type') 
ggsave('./neurons/neuropeptides/Percentage AVP neurons expressing OXT.png',
       height = 10,
       width = 10)


# add color
mouse.snseq.combined.sct.neurons.recluster.expression %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>%
  filter(AVP.oxt.cell == c('Avp', 
                           'Both')) %>% 
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  filter(AVP.oxt.cell == 'Both') %>% 
  filter(!is.na(Percent)) %>% 
  mutate(orig.ident = tolower(orig.ident)) %>% 
  separate(orig.ident,
           c("Sex",
             "Status"),
           remove = F) %>% 
  ggplot(aes(x = reorder(orig.ident,
                         -Percent))) +
  geom_bar_pattern(aes(y=Percent,
                       fill = Status,
                       pattern = Sex), 
                   stat="identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  theme_classic() +
  ylab('Percentage AVP neurons OXT+') +
  xlab('Sample') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46")) +
  scale_pattern_manual(values = c(female = "stripe", 
                                  male = "none")) +
  theme(text = element_text(size = 20)) +
  ylim(0,100)
ggsave('./neurons/neuropeptides/Bar chart overlap AVP neurons expressing OXT.png',
       height = 10,
       width = 15)

# add scale bar
mouse.snseq.combined.sct.neurons.recluster.expression %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(AVP.oxt.cell == c('Avp', 
                           'Both')) %>% 
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  filter(AVP.oxt.cell == 'Both') %>% 
  filter(!is.na(Percent)) %>% 
  mutate(orig.ident = tolower(orig.ident)) %>% 
  ggplot(aes(x = reorder(orig.ident,
                         -Percent))) +
  geom_bar(aes(y=Percent,
                       fill = orig.ident), 
                   stat="identity",
           color = 'black') +
  theme_classic() +
  ylab('Percentage AVP neurons OXT+') +
  xlab('') + 
  scale_fill_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF")) +
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  ylim(0,100)
ggsave('./neurons/neuropeptides/Bar chart overlap AVP neurons expressing OXT color.png',
       height = 10,
       width = 15)




#### all WGCNA ####
#gene_select parameter:

#   variable: use the genes stored in the Seurat object’s VariableFeatures.
# fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# custom: use genes that are specified in a custom list.

# mouse.snseq.combined.sct.neurons.recluster <- SetupForWGCNA(
#   mouse.snseq.combined.sct.neurons.recluster,
#   gene_select = "variable", # the gene selection approach
#   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# )

## set up object for analysis
mouse.snseq.combined.sct = SetupForWGCNA(mouse.snseq.combined.sct,
                                                           features = VariableFeatures(mouse.snseq.combined.sct),
                                                           wgcna_name = "hdWGCNA"
)


# construct metacells  in each group
mouse.snseq.combined.sct <- MetacellsByGroups(
  seurat_obj = mouse.snseq.combined.sct,
  group.by = c("orig.ident",
               "parent_id.broad"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'parent_id.broad', # set the Idents of the metacell seurat object
  slot = "data", #set for counts
  assay = "SCT", #set assay type
  mode = "average", #metacell expression profile determined by average of constituent single cells
  # min_cells = 50 #set minimum number of cells per metacell
)

## setup expression matrix
#use metacell expression data

##select soft-power threshold
mouse.snseq.combined.sct <- TestSoftPowers(
  mouse.snseq.combined.sct,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.all <- PlotSoftPowers(mouse.snseq.combined.sct)

#assemble with patchwork
png('./hdwgcna/all/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.all,
           ncol=2)
dev.off()

#use soft power threshold 5


### construct co-expression network
# construct co-expression network:
mouse.snseq.combined.sct <- ConstructNetwork(mouse.snseq.combined.sct, 
                                                               soft_power=5,
                                                               # setDatExpr=FALSE,
                                                               tom_name = 'all', # name of the topoligical overlap matrix written to disk
                                                               overwrite_tom = TRUE
)

#graph dendrogram
png('./hdwgcna/all/WGCNA dendrogram hypomap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(mouse.snseq.combined.sct,
               main='all hdWGCNA Dendrogram')
dev.off()

#graph dendrogram for poster
# pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
#     width = 10,
#     height = 2.5)
# PlotDendrogram(mouse.snseq.combined.sct.neurons.recluster,
#                main='Neurons hdWGCNA Dendrogram')
# dev.off()
# ended up exporting the right size with plot viewer

##compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
# mouse.snseq.combined.sct.neurons.recluster <- ScaleData(mouse.snseq.combined.sct.neurons.recluster, 
#                                                               features=VariableFeatures(mouse.snseq.combined.sct.neurons.recluster))

# compute all MEs in the full single-cell dataset

mouse.snseq.combined.sct <- ModuleEigengenes(
  mouse.snseq.combined.sct,
  exclude_grey = TRUE
  # , group.by.vars = 'Genotype.id'
)

# module eigengenes:
MEs.all <- GetMEs(mouse.snseq.combined.sct,
                      harmonized=FALSE)


# compute eigengene-based connectivity (kME):
mouse.snseq.combined.sct <- ModuleConnectivity(
  mouse.snseq.combined.sct)


#### Graph WGCNNA for all ####

# plot genes ranked by kME for each module
png('./hdwgcna/all/Plot kMEs all.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotKMEs(mouse.snseq.combined.sct, 
         ncol=2)
dev.off()

# correlation between modules
png('./hdwgcna/all/Correlogram modules all.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleCorrelogram(mouse.snseq.combined.sct,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001), 
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.all.me <- ModuleFeaturePlot(
  mouse.snseq.combined.sct,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('./hdwgcna/all/UMAP MEs all.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.all.me, ncol=2)
dev.off()




# get mods from object
mods.all <- colnames(MEs.all); mods.all <- mods.all[mods.all != 'grey']

# add MEs to Seurat meta-data:
mouse.snseq.combined.sct@meta.data <- cbind(mouse.snseq.combined.sct@meta.data, 
                                                              MEs.all)

#### hdwgcna: GLU ####

### setup GLU wgcna
mouse.snseq.combined.sct = SetupForWGCNA(mouse.snseq.combined.sct,
                                         features = VariableFeatures(mouse.snseq.combined.sct),
                                         metacell_location = "hdWGCNA",
                                         wgcna_name = 'GLU'
)

# get data for each cell type wgcna
mouse.snseq.combined.sct= SetDatExpr(
  mouse.snseq.combined.sct,
  group_name = 'C7-1: GLU',
  use_metacells = TRUE,
  group.by = 'parent_id.broad',
  wgcna_name = "GLU",
  assay = 'SCT',
  slot = "data")


##select soft-power threshold
mouse.snseq.combined.sct <- TestSoftPowers(
  mouse.snseq.combined.sct,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.GLU <- PlotSoftPowers(mouse.snseq.combined.sct)

#assemble with patchwork
png('./hdwgcna/GLU/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.GLU,
           ncol=2)
dev.off()


#use soft power threshold 10


### construct co-expression network
# construct co-expression network:
mouse.snseq.combined.sct <- ConstructNetwork(mouse.snseq.combined.sct, 
                                             soft_power=10,
                                             # setDatExpr=FALSE,
                                             tom_name = 'GLU', # name of the topoligical overlap matrix written to disk
                                             overwrite_tom = TRUE
)

#graph dendrogram
png('./hdwgcna/GLU/WGCNA dendrogram hypomap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(mouse.snseq.combined.sct,
               main='GLU hdWGCNA Dendrogram')
dev.off()

# compute  MEs in the full single-cell dataset
mouse.snseq.combined.sct <- ModuleEigengenes(
  mouse.snseq.combined.sct,
  exclude_grey = TRUE,
  wgcna_name = 'GLU'
)

# module eigengenes:
MEs.GLU <- GetMEs(mouse.snseq.combined.sct,
                  harmonized=FALSE,
                  wgcna_name = 'GLU')


# compute eigengene-based connectivity (kME):
mouse.snseq.combined.sct <- ModuleConnectivity(
  mouse.snseq.combined.sct,
  wgcna_name = 'GLU')



# make a featureplot of hMEs for each module
plot_list.GLU.me <- ModuleFeaturePlot(
  mouse.snseq.combined.sct,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('./hdwgcna/GLU/UMAP MEs GLU.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.GLU.me, ncol=2)
dev.off()



#### hdwgcna: GABA ####

### setup GABA wgcna
mouse.snseq.combined.sct = SetupForWGCNA(mouse.snseq.combined.sct,
                                         features = VariableFeatures(mouse.snseq.combined.sct),
                                         metacell_location = "hdWGCNA",
                                         wgcna_name = 'GABA'
)

# get data for each cell type wgcna
mouse.snseq.combined.sct= SetDatExpr(
  mouse.snseq.combined.sct,
  group_name = 'C7-2: GABA',
  use_metacells = TRUE,
  group.by = 'parent_id.broad',
  wgcna_name = "GABA",
  assay = 'SCT',
  slot = "data")


##select soft-power threshold
mouse.snseq.combined.sct <- TestSoftPowers(
  mouse.snseq.combined.sct,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.GABA <- PlotSoftPowers(mouse.snseq.combined.sct)

#assemble with patchwork
png('./hdwgcna/GABA/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.GABA,
           ncol=2)
dev.off()


#use soft power threshold 9


### construct co-expression network
# construct co-expression network:
mouse.snseq.combined.sct <- ConstructNetwork(mouse.snseq.combined.sct, 
                                             soft_power=9,
                                             # setDatExpr=FALSE,
                                             tom_name = 'GABA', # name of the topoligical overlap matrix written to disk
                                             overwrite_tom = TRUE
)

#graph dendrogram
png('./hdwgcna/GABA/WGCNA dendrogram hypomap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(mouse.snseq.combined.sct,
               main='GABA hdWGCNA Dendrogram')
dev.off()

# compute  MEs in the full single-cell dataset
mouse.snseq.combined.sct <- ModuleEigengenes(
  mouse.snseq.combined.sct,
  exclude_grey = TRUE,
  wgcna_name = 'GABA'
)

# module eigengenes:
MEs.GABA <- GetMEs(mouse.snseq.combined.sct,
                  harmonized=FALSE,
                  wgcna_name = 'GABA')


# compute eigengene-based connectivity (kME):
mouse.snseq.combined.sct <- ModuleConnectivity(
  mouse.snseq.combined.sct,
  wgcna_name = 'GABA')




# make a featureplot of hMEs for each module
plot_list.GABA.me <- ModuleFeaturePlot(
  mouse.snseq.combined.sct,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('./hdwgcna/GABA/UMAP MEs GABA.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.GABA.me, ncol=2)
dev.off()



#### save point ####
# save(mouse.snseq.combined.sct,
#      file = "./hdwgcna/mouse.snseq.combined.sct.RData")
# load('./hdwgcna/mouse.snseq.combined.sct.RData')
