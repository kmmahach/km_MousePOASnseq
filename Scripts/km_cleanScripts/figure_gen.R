#### Set up environment ####
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# load functions 
source("./functions/DGE_fun.R")
source("./functions/QC_filtering_fun.R")

load_packages(c("tidyverse", "limma", "tidytext", "ggrepel", "ggplot2", "ggpol", "scales", 
                "svglite", "viridis", "patchwork", "Seurat", "RRHO2", "RedRibbon"),
              out_prefix = "figures")

#### Figure 2 ####

## cell type annotation based on ScType & HypoMap predictions
setwd(paste0(root.dir, "/HypoMap"))
load("./data/integrated_seurat_withHypoMap_predictions.rda")

# split OPCs and oligodendrocytes and remove ParsTuber
int.ldfs@meta.data %>%
  mutate(Cell_Type = case_when(
    
    parent_id.broad.prob == "C7-4: Oligo+Precursor" & 
      sctype.integrated == "Oligodendrocyte precursor cells" ~ "OPCs",
    
    parent_id.broad.prob == "C7-4: Oligo+Precursor" & 
      sctype.integrated != "Oligodendrocyte precursor cells" ~ "Oligodendrocytes",
    
    parent_id.broad.prob == "C7-1: GLU" ~ "Glutamatergic neurons",
    
    parent_id.broad.prob == "C7-2: GABA" ~ "GABAergic neurons",
    
    parent_id.broad.prob == "C7-3: Astro-Ependymal" & 
      sctype.integrated != "Radial glial cells" ~ "Astrocytes",
    
    parent_id.broad.prob == "C7-5: Immune" ~ "Microglia",
    
    parent_id.broad.prob == "C7-7: Vascular" ~ "Epithelial Cells",
    
    TRUE ~ "Unknown" 
    
  )) %>%
  pull(Cell_Type) -> int.ldfs@meta.data$Cell_Type 

# UMAP - figure 2A
DimPlot(int.ldfs,
        group.by = "Cell_Type",
        cols = c('#7570b3',
                 '#a6761d',
                 '#d95f02',
                 '#1b9e77',
                 '#66a61e',
                 '#e7298a',
                 '#e6ab02',
                 'grey')) +
  ggtitle("Cell Type Annotation") 

ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/assigned.cell_types.UMAP.svg'),
       height = 6, width = 8)  

# Neuron clusters
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

# UMAP - neurons only - figure 2B
MSCneurons.reclust %>%
  DimPlot() +
  ggtitle('Neuron Clusters') +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()

ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/neurons.recluster.clusters.svg'),
       height = 6, width = 6)


# Differential expression of neuron cluster marker genes
setwd(paste0(root.dir, "/DGE_CellTypes/")) 

# data 
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")
load("./neurons/all_neurons/limma_perm/rrho_results.rda")
neuron_clusterMarkers <- read_csv('neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')

# color by enrichment quadrant from RRHO to indicate direction
femaleDEGs <- limma_list$Fdom_vs_sub_all
maleDEGs <- limma_list$Mdom_vs_sub_all

sig_maleDEGs <- maleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

sig_femaleDEGs <- femaleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

conc_genes = intersect(sig_femaleDEGs$gene,
                       sig_maleDEGs$gene)

result2 <- rrho_results %>%
  .[str_detect(names(.), "^cluster\\d+$")] %>%
  imap_dfr(~ {
    cluster_name <- .y
    df <- .x$df
    map_dfr(names(.x$RedRibbon.quads), function(q) {
      pos <- .x$RedRibbon.quads[[q]]$positions
      
      df[pos, ] %>%
        transmute(
          gene,
          cluster = str_extract(cluster_name, "\\d+"),
          quadrant = q,
          avg.logFC = (abs(a) + abs(b)) / 2
        ) 
    })
  })

neuron_clusterMarkers_hiSpec <- neuron_clusterMarkers %>% 
  filter(specificity > 0.1)

result2 <- result2 %>% 
  filter(gene %in% neuron_clusterMarkers_hiSpec$gene &
           avg.logFC > 1) %>% 
  complete(gene, cluster,
           fill = list(avg.logFC = NA))


result3 <- result2 %>%
  filter(!is.na(quadrant)) %>% 
  count(gene, cluster) %>%
  group_by(gene) %>%
  summarise(total_n = sum(n)) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, total_n)) %>% 
  left_join(result2)


# heatmap of cluster markers by avg DE - figure 2C
result3 %>% 
  filter(is.na(quadrant)) %>% 
  ggplot(aes(x = factor(cluster),
             y = reorder(gene, total_n),
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  scale_fill_viridis(na.value = "gray95") +
  geom_tile(data = subset(result3, quadrant == "upup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#7e38b7",
                      na.value = "gray95",
                      limits = c(0, 2),
                      name = bquote(plain(" \u2191") * bold("Dom") *
                                      "\u2640  \u2191" * bold("Dom ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "downdown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#408D8E",
                      na.value = "gray95",
                      limits = c(0, 2),
                      name = bquote(plain(" \u2191") * bold("Dom") *
                                      "\u2640  \u2191" * bold("Sub ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "updown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#f94449",
                      na.value = "gray95",
                      limits = c(0, 2),
                      name = bquote(plain(" \u2191") * bold("Sub") *
                                        "\u2640  \u2191" * bold("Dom ") * "\u2642  " *
                                        plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "downup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#ff7d00",
                      na.value = "gray95",
                      limits = c(0, 2),
                      name = bquote(plain(" \u2191") * bold("Sub") *
                                      "\u2640   \u2191" * bold("Sub ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(drop_na(result3), avg.logFC, quadrant),
            aes(x = 1,
                y = 1,
                fill = quadrant),
            alpha = 0) +
  scale_fill_manual(name = "effect of status across sex",
                    values = c("upup" = "#7e38b7",
                               "downdown" = "#408D8E",
                               "updown" = "#f94449",
                               "downup" = "#ff7d00")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()

ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/cluster_marker.quads.svg"),
       device = "svg", width = 9.25, height = 16)
ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/cluster_marker.quads.png"),
       width = 9.25, height = 16)

# cluster composition by group + number of cells - figure 2D
# create counts table - sample by cluster
table_orig.idents_by_clusters <- MSCneurons.reclust@meta.data %>%
  group_by(orig.ident, integrated_snn_res.0.4) %>%
  summarize(count = n()) %>%
  spread(integrated_snn_res.0.4, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('orig.ident', 'total_cell_count', everything())) %>%
  arrange(factor(orig.ident, levels = levels(MSCneurons.reclust@meta.data$orig.ident))) %>% 
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident))

# cluster by sample
table_clusters_by_orig.idents <- MSCneurons.reclust@meta.data %>%
  dplyr::rename('cluster' = 'integrated_snn_res.0.4') %>%
  group_by(cluster, orig.ident) %>%
  summarize(count = n()) %>% 
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  spread(orig.ident, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(MSCneurons.reclust@meta.data$integrated_snn_res.0.4)))

# cells per cluster by orig.ident
MSCneurons.reclust@meta.data %>% 
  group_by(orig.ident) %>%
  tally() %>% 
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) -> temp_labels

p1 <- table_orig.idents_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = temp_labels, 
            aes(x = orig.ident, 
                y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Cluster") +
  theme(legend.position = 'left',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1),
        plot.margin = margin(t = 20, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = 'pt'))

MSCneurons.reclust@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4) -> temp_labels

p4 <- table_clusters_by_orig.idents %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Group") +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

p5 <- p4 + 
  scale_fill_manual(values = adjustcolor(c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"),
                                         alpha.f = 0.75)) +
  labs(x = "Neuron Cluster",
       y = "Total Nuclei") +
  theme(axis.title = element_text(size = 16))

p6 <- p1 +
  scale_fill_viridis_d(option = "plasma",
                       alpha = 0.75)


ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/composition_neuron.clusters_byPct_Num.svg'),
       p6 + p5 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 18, height = 8)


# RRHO heatmap sex + across status 
# concordance/discordance of gene expression in all neurons - figure 2E
setwd(paste0(root.dir, "/DGE_CellTypes/")) 
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")

group.by = c("Male", "Female")
compare.across = c("Dom", "Sub")

# Reformat for RRHO
suffixes <- sub(".*sub_?", "", names(limma_list))
limma_list <- split(limma_list, suffixes)
limma_list <- lapply(limma_list, function(sublist) {
  newnames <- vapply(names(sublist), function(nm) {
    if (grepl("^F", nm)) {
      "ttFemale_Dom_vs_Female_Sub"
    } else if (grepl("^M", nm)) {
      "ttMale_Dom_vs_Male_Sub"
    } else {
      nm  # fallback if it doesn't match
    }
  }, character(1))
  
  names(sublist) <- newnames
  sublist
})

All_Neurons = limma_list$all
comparisons = names(All_Neurons)

  dataset <- All_Neurons
  comps <- comparisons
  
  # pull comparison dfs
  left_df_raw  <- All_Neurons[[comps[1]]] %>% 
    column_to_rownames(var = "gene")
  right_df_raw <- dataset[[comps[2]]] %>%
    column_to_rownames(var = "gene")
  
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
                         value = left_df_raw$logFC * left_sign,
                         stringsAsFactors = FALSE)
  
  right_df <- data.frame(gene = rownames(right_df_raw),
                         value = right_df_raw$logFC * right_sign,
                         stringsAsFactors = FALSE)

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
  
  # clean_names <- sub("^tt", "", comps[grepl("^tt", comps)])
  
  # RedRibbon plot
  prr <- ggRedRibbon(rr, 
                    quadrants = RedRibbon.quads) +
    coord_fixed(ratio = 1,
                clip = "off") +
    xlab(gsub("_", " ", comps[1])) +    
    ylab(gsub("_", " ", comps[2])) +
    ggtitle("All Neurons") +
    theme(plot.title = element_text(size = 25,
                                    face = "bold",
                                    hjust = 0.5),
          plot.margin = unit(c(0,0.5,0,0.5), units = "cm")) # top, right, bottom, left
  
  prr$scales$scales[[2]]$labels = c("higher in Sub", "higher in Dom")
  prr$scales$scales[[2]]$breaks = c(85, 425)
  prr$scales$scales[[2]]$name = "Females"
  prr$scales$scales[[3]]$labels = c("higher in Sub", "higher in Dom")
  prr$scales$scales[[3]]$breaks = c(150, 500)
  prr$scales$scales[[3]]$name = "Males"
    
    ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/AllNeurons_RedRibbon_concordance_by_sex.svg"),
           height = 7.5, width = 7.5, dpi = 300)


#### Most common DEGs across clusters ####
setwd(paste0(root.dir, "/DGE_CellTypes/")) 

# data 
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")
load("./neurons/all_neurons/limma_perm/rrho_results.rda")

subset_by_sex <- function(limma_list,
                          logFC_threshold = 0.2,
                          min_pval = 0.05) 
{
  sublist <- lapply(limma_list, \(lst) {
    lst <- subset(lst, abs(lst$logFC) > logFC_threshold &
                    lst$P.Value < min_pval)
  }
  )
  
  Females <- subset(sublist, grepl("^[F]", names(sublist))==TRUE)  
  Males <- subset(sublist, grepl("^[M]", names(sublist))==TRUE)  
  
  newlist <- as_named_list(Females, Males)
  
  np_DGE <- lapply(newlist, \(df) {
    
    results = bind_rows(df, .id = "contrast") %>% 
      mutate(contrast = sub(".*_", "", contrast))
    
    newdf <- list( results = results,
                   num_clust = results %>%
                     distinct(gene, contrast) %>%
                     count(gene, name = "n_clust"),      
                   num_genes = results %>% 
                     distinct(contrast, gene) %>% 
                     count(contrast, name = "n_genes") )
    return(newdf)
  }
  ) 
  return(np_DGE)
}


DGE_list <- subset_by_sex(limma_list)


# DE in > 6 clusters/subtypes, get top genes DE in the most clusters
topDEGs <- lapply(DGE_list, \(lst) {
  
  subset_genes <- subset(lst$num_clust, n_clust > 6)
  top_genes <- subset(lst$results, gene %in% subset_genes$gene) %>% 
    left_join(lst$num_clust,
              by = "gene") %>% 
    group_by(contrast) %>%
    slice_max(order_by = n_clust, n = 7)
  
  results <- subset(lst$results, gene %in% top_genes$gene)
  
  lst.df <- list( results = results,
                  num_clust = results %>%
                    distinct(gene, contrast) %>%
                    count(gene, name = "n_clust"),      
                  num_genes = results %>% 
                    distinct(contrast, gene) %>% 
                    count(contrast, name = "n_genes") )


  return(lst.df)
  }
)

# make df with contrast(cluster), sex, status, gene, logFC, num_genes, and n_clust
topDEGs <- lapply(topDEGs, \(lst) {
  df = lst$results %>% 
    left_join(lst$num_clust,
              by = "gene") %>% 
    left_join(lst$num_genes,
              by = "contrast")
})

topDEGs.test <- topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) 

topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) %>% 
  mutate(contrast = reorder_within(contrast, -n_genes, Sex),   
         gene = reorder_within(gene, n_clust, Sex)) %>%
  ggplot(aes(x = contrast, 
             y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey90") +
  labs(x = "neuronal subgroups", 
       y = "top DEGs") +
  scale_fill_gradient2(name = bquote(~italic(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(topDEGs.test, Status, logFC),
            aes(x = 1, 
                y = 1, 
                fill = Status), 
            alpha = 0) +
  scale_fill_manual(name = paste0(sprintf("\u2191 \u2191"), " expr in"),
                    values = c("Dom" = muted("blue"), "Sub" = muted("red"))) +
  guides(fill = guide_legend(override.aes = list(values = c(muted("blue"),
                                                            muted("red")),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  facet_wrap(~ Sex, scales = "free") +
  scale_x_reordered() +
  scale_y_reordered() +
  ggtitle("DEGs across clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.65,
                                   face = "italic", size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15, face = "italic"),
        legend.text = element_text(size = 12, face = "bold"))

ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/topDEGs_byCluster.svg"),
       device = "svg", width = 12, height = 6) 

#### figure out where to put this figure (if keeping it)
# overlap with cluster marker genes 
neuron_clusterMarkers <- read_csv('neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')

sig_maleDEGs <- maleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

sig_femaleDEGs <- femaleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

conc_genes = intersect(sig_femaleDEGs$gene,
                       sig_maleDEGs$gene)

result <- limma_list %>%
  .[str_detect(names(.), "_cluster\\d+$")] %>%
  imap_dfr(~ {
    .x %>%
      filter(gene %in% conc_genes) %>%
      select(gene, logFC, P.Value) %>%
      mutate(
        sex = str_to_lower(str_sub(.y, 1, 1)),  
        cluster = str_extract(.y, "(?<=_cluster)\\d+")
      )
  }) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(logFC, P.Value),
    names_glue = "{.value}.{sex}"
  ) %>%
  filter(!is.na(logFC.f) & !is.na(logFC.m)) %>%
  group_by(cluster) %>%
  filter(all(!is.na(logFC.f)) & all(!is.na(logFC.m))) %>%
  ungroup() %>% 
  mutate(
    logFC.f   = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, logFC.f, NA),
    logFC.m   = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, logFC.m, NA),
    P.Value.f = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, P.Value.f, NA),
    P.Value.m = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, P.Value.m, NA),
    avg.logFC = (abs(logFC.f) + abs(logFC.m)) / 2 
  ) %>% 
  filter(!(is.na(logFC.f) | is.na(logFC.m)))

neuron_clusterMarkers_conc <- neuron_clusterMarkers %>% 
  filter(gene %in% result$gene) 

result.test = left_join(result, neuron_clusterMarkers_conc[, c("gene", "specificity")]) 

result.test %>% 
  # slice_max(avg.logFC, n = 6) %>% 
  dplyr::select(gene, cluster, logFC.f, logFC.m, avg.logFC, specificity) %>%
  complete(gene, cluster,
           fill = list(avg.logFC = NA)) %>%
  group_by(cluster) %>% 
  mutate(avg.logFC = ifelse(logFC.f * logFC.m < 0, -1 * avg.logFC, avg.logFC),
         direction = ifelse(avg.logFC > 0, "concordant", "discordant")) -> result.test

result.test %>% 
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  scale_fill_gradient2(low = "#408D8E",
                       high = "#7e38b7",
                       na.value = "gray95",
                       name = bquote(avg~abs(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(drop_na(result.test), avg.logFC, direction),
            aes(x = 1,
                y = 1,
                fill = direction),
            alpha = 0) +
  scale_fill_manual(name = "effect of status across sex",
                    values = c("concordant" = "#7e38b7", "discordant" = "#408D8E")) +
  guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
                                                            "#408D8E"),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/cluster_marker.DEGs.svg"),
       device = "svg", width = 9.25, height = 16)



#### Oxytocin and Vasopressin DE by cluster #### 
setwd(paste0(root.dir, "/DGE_CellTypes/")) 
## if starting here:
load("./neurons/all_neurons/limma_perm/norm_counts.rda")
metadata <- read_csv("./neurons/all_neurons/limma_perm/metadata.csv")

norm_counts <- data.frame(v.dl$E) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(gene %in% c("Oxt", "Avp")) %>% 
  pivot_longer(cols = !gene, 
               values_to = "norm_expr",
               names_to = "pb_colname") %>% 
  left_join(metadata, by = "pb_colname") %>% 
  mutate(group = paste(Sex, Status, sep = " "),
         id = paste(Sex, Status, indiv_genotype, sep = "_"),
         seurat_clusters = paste("cluster", seurat_clusters, sep = " ")) %>% 
  group_by(gene, group, id, seurat_clusters) %>% 
  summarise(norm_expr = mean(norm_expr))

norm_counts %>% 
  filter(gene == "Oxt") %>% 
  ggplot(aes(group, norm_expr, 
             color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 30), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15)) + 
  ggtitle("OXT expr in neurons")

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/oxtExpr_byCluster.svg"),
       device = "svg", width = 10, height = 7)


norm_counts %>% 
  filter(gene == "Avp") %>% 
  ggplot(aes(group, norm_expr, 
             color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 30), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15)) + 
  ggtitle("AVP expr in neurons")

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/avpExpr_byCluster.svg"),
       device = "svg", width = 10, height = 7)

#### Supplemental Figure 1 ####
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/"))
load("./data/rawdata_withGeneByCell_counts.rda")
load("./data/souporcell_output.rda")

l.dfs <- make.seurat.obj(l.df2) 
# add metadata
l.dfs <- map2(l.dfs, l.clust, ~ { AddMetaData(.x, metadata = .y) })

# label mitochondrial genes
lapply(l.dfs, \(x) {
  PercentageFeatureSet(x, pattern = "^mt-",
                       col.name = "percent.mt") 
 } 
) -> l.dfs

# choose filter cutoffs 
filter.categories <- c("nFeature_min", "nFeature_max", "mito.max")
  fdf <- c(700, 1400, 5)
  fsf <- c(700, 1500, 5)
  mdf <- c(700, 3000, 5)
  msf <- c(700, 2500, 5)

filters <- matrix(rbind(fdf,fsf,mdf,msf),
                  nrow = length(l.dfs),
                  ncol = length(filter.categories),
                  dimnames = list(names(l.dfs),
                                  filter.categories))

# plot filtered QC metrics - supp fig 1
plot.filters(l.dfs, filters, 
             outdir = paste0(root.dir, "/manuscriptFigures/supplemental/"))


#### Supplemental Figure 4 ####
setwd(paste0(root.dir, "/HypoMap"))
load("./data/integrated_seurat_withHypoMap_predictions.rda")

# exc/inh neuron clusters - supp fig 4D
MSCneurons.reclust %>%
  DimPlot(group.by = 'parent_id.broad.prob') +
  ggtitle('Excitatory vs. Inhibitory Neurons') +
  scale_color_manual(values = c('#1b9e77',
                                '#d95f02')) +
  theme_classic()

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster.parent_id.broad.svg'),
       height = 5, width = 6)

# colored by sample pool (orig.ident)
# randomize overlay of points/make sure all are relatively visible
MSCneurons.reclust@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Cell_ID") %>% 
  left_join(data.frame(MSCneurons.reclust@meta.data) %>% 
              dplyr::select(orig.ident, Cell_ID),
            by = "Cell_ID") %>% 
  mutate(orig.ident = sub("\\.data$", "", orig.ident)) -> umap

umap = umap[sample(nrow(umap)),]
rownames(umap) = seq(1:nrow(umap))


ggplot(data = subset(umap, !orig.ident == "female.sub"),
       aes(x = UMAP_1,
           y = UMAP_2,
           color = orig.ident)) +
  geom_point(size = 0.3) +
  # scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  geom_point(data = subset(umap, orig.ident == "female.sub"),
             size = 0.3) +
  theme_classic() +
  labs(x = "UMAP_1",
       y = "UMAP_2") +
  guides(color = guide_legend(override.aes = 
                                list(size = 3))) +
  theme(legend.text = element_text(face = "bold"),
        # axis.text = element_text(size = 10),
        # axis.title = element_text(size = 11),
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm")
  ) + NoAxes()

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster.orig.ident.svg'),
       height = 5.5, width = 7)


#### Supplemental Figure 5 ####
setwd(paste0(root.dir, "/HypoMap"))
load("./data/integrated_seurat_withHypoMap_predictions.rda")

# specific sub-types - supp fig 5A
MSCneurons.reclust %>%
  DimPlot(group.by = 'predicted.prob') +
  ggtitle('Neuron Subtypes') +
  theme_classic()

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster_predicted.prob.svg'),
       height = 5, width = 11)


# volcano plots 
setwd(paste0(root.dir, "/DGE_CellTypes"))
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")

sex = c("Male", "Female")
status = c("Dom", "Sub")

femaleDEGs <- limma_list$Fdom_vs_sub_all
maleDEGs <- limma_list$Mdom_vs_sub_all

up <- "Dom"
down <- "Sub"

# males - all neurons - supp fig 5B
dcm <- maleDEGs %>% 
  mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
  mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
  mutate(rank = logFC * log10) %>% 
  mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))

max.x <- round(max(dcm$logFC) * 2) / 2
min.x <- round(min(dcm$logFC) * 2) / 2 
breaks = seq(min.x, max.x, 0.5)

upgene_labels <- dcm %>% slice_max(order_by = rank, n = 10)
downgene_labels <- dcm %>% slice_max(order_by = -rank, n = 10)

dcm_jitters <- subset(dcm, dcm$log10 > 3 & 
                        !gene %in% c(upgene_labels$gene, downgene_labels$gene))

dcm_jitters_up <- subset(dcm, dcm$gene %in% upgene_labels)
dcm_jitters_down <- subset(dcm, dcm$gene %in% downgene_labels)
pos_jitter <- position_jitter(width = 0.4, height = 0.1, seed = 123)

ggplot() +
  geom_point(data = subset(dcm, !gene %in% c(dcm_jitters_up$gene,
                                             dcm_jitters_down$gene,
                                             dcm_jitters$gene)),
             aes(x = logFC, y = log10, color = diffexpressed), 
             alpha = 0.25, size = 3.5) +
  geom_jitter(data = dcm_jitters, 
              aes(x = logFC, y = log10, color = diffexpressed), 
              height = 0.4, width = 0.1, 
              alpha = 0.25, size = 3.5) +
  geom_jitter(data = upgene_labels,
              aes(x = logFC, y = log10, color = diffexpressed), 
              position = pos_jitter,
              alpha = 0.25, size = 3) +
  geom_jitter(data = downgene_labels,
              aes(x = logFC, y = log10, color = diffexpressed), 
              position = pos_jitter,
              alpha = 0.25, size = 3) +
  geom_text_repel(data = upgene_labels,
                  aes(x = logFC, y = log10, label = gene),
                  size = 5, color = "black", hjust = 0.2, vjust = 0.7,
                  max.overlaps = Inf) +
  geom_text_repel(data = downgene_labels,
                  aes(x = logFC, y = log10, label = gene),
                  size = 5, color = "black", vjust = 0.5, hjust = -0.5,
                  max.overlaps = Inf) +
  geom_hline(yintercept = 1.301, lty = 4, col = "grey", lwd = 0.8) +
  labs(x = "log2 Fold Change",
       y = bquote(~-Log[10]~italic(eFDR))) +
  annotate(geom = "text", x = 1.5, y = 0.5, 
           label = paste0("\u2191\u2191", up), 
           fontface = "bold", color = "black", size = 8) +
  annotate(geom = "text", x = -1.5, y = 0.5, 
           label = paste0("\u2191\u2191", down), 
           fontface = "bold", color = "black", size = 8) +
  scale_color_manual(values = c("UP" =  "#ff7d00",  "DOWN" = "#7e38b7", "NO" = "grey")) +
  scale_x_continuous(limits = c(min.x, max.x), breaks = breaks) +
  ylim(0, 4) +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        axis.text.y = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(color = "#3C3C3C", size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 25)) +
  ggtitle("Males only")


ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/onlyMale_all_volcano_plot.svg"),
       width = 27, height = 22, units = "cm", dpi = 600)

# females - all neurons - supp fig 5C
dcf <- femaleDEGs %>% 
  mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
  mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
  mutate(rank = logFC * log10) %>% 
  mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))

max.x <- round(max(dcf$logFC) * 2) / 2
min.x <- round(min(dcf$logFC) * 2) / 2 
breaks = seq(min.x - 1, max.x, 0.5)

upgene_labels <- dcf %>% slice_max(order_by = rank, n = 10)
downgene_labels <- dcf %>% slice_max(order_by = -rank, n = 10)

dcf_jitters <- subset(dcf, dcf$log10 > 3 & 
                        !gene %in% c(upgene_labels$gene, downgene_labels$gene))

dcf_jitters_up <- subset(dcf, dcf$gene %in% upgene_labels)
dcf_jitters_down <- subset(dcf, dcf$gene %in% downgene_labels)
pos_jitter <- position_jitter(width = 0.4, height = 0.1, seed = 123)

ggplot() +
  geom_point(data = subset(dcf, !gene %in% c(dcf_jitters_up$gene,
                                             dcf_jitters_down$gene,
                                             dcf_jitters$gene)),
             aes(x = logFC, y = -log10(P.Value), color = diffexpressed), 
             alpha = 0.25, size = 3.5) +
  geom_jitter(data = dcf_jitters, 
              aes(x = logFC, y = log10, color = diffexpressed), 
              height = 0.4, width = 0.1, 
              alpha = 0.25, size = 3.5) +
  geom_jitter(data = upgene_labels,
              aes(x = logFC, y = log10, color = diffexpressed), 
              position = pos_jitter,
              alpha = 0.25, size = 3) +
  geom_jitter(data = downgene_labels,
              aes(x = logFC, y = log10, color = diffexpressed), 
              position = pos_jitter,
              alpha = 0.25, size = 3) +
  geom_text_repel(data = upgene_labels,
                  aes(x = logFC, y = -log10(P.Value), label = gene),
                  size = 5, color = "black", hjust = -1, vjust = 0.5,
                  max.overlaps = Inf) +
  geom_text_repel(data = downgene_labels,
                  aes(x = logFC, y = -log10(P.Value), label = gene),
                  size = 5, color = "black", vjust = -1.2, hjust = -0.5,
                  max.overlaps = Inf) +
  geom_hline(yintercept = 1.301, lty = 4, col = "grey", lwd = 0.8) +
  labs(x = "log2 Fold Change",
       y = bquote(~-Log[10]~italic(eFDR))) +
  annotate(geom = "text", x = 1.5, y = 0.5, 
           label = paste0("\u2191\u2191", up), 
           fontface = "bold", color = "black", size = 8) +
  annotate(geom = "text", x = -1.5, y = 0.5, 
           label = paste0("\u2191\u2191", down), 
           fontface = "bold", color = "black", size = 8) +
  scale_color_manual(values = c("UP" = "#f94449", "DOWN" = "#408D8E", "NO" = "grey")) +
  scale_x_continuous(limits = c(min.x -1, max.x), breaks = breaks) +
  ylim(0, 4) +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, size = 15),
        axis.text.y = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(color = "#3C3C3C", size = 20),
        axis.title = element_text(size = 20),   
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size = 25)) +
  ggtitle("Females only")


ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/onlyFem_all_volcano_plot.svg"),
       width = 27, height = 22, units = "cm", dpi = 600)


#### Supplemental Figure 6 ####
setwd(paste0(root.dir, "/DGE_CellTypes/")) 

# data for DEG x clusters heatmap - figure 6C 
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")
load("./neurons/all_neurons/limma_perm/rrho_results.rda")

subset_by_sex <- function(limma_list,
                          logFC_threshold = 0.2,
                          min_pval = 0.05) 
{
  sublist <- lapply(limma_list, \(lst) {
    lst <- subset(lst, abs(lst$logFC) > logFC_threshold &
                    lst$P.Value < min_pval)
  }
  )
  
  Females <- subset(sublist, grepl("^[F]", names(sublist))==TRUE)  
  Males <- subset(sublist, grepl("^[M]", names(sublist))==TRUE)  
  
  newlist <- as_named_list(Females, Males)
  
  np_DGE <- lapply(newlist, \(df) {
    
    results = bind_rows(df, .id = "contrast") %>% 
      mutate(contrast = sub(".*_", "", contrast))
    
    newdf <- list( results = results,
                   num_clust = results %>%
                     distinct(gene, contrast) %>%
                     count(gene, name = "n_clust"),      
                   num_genes = results %>% 
                     distinct(contrast, gene) %>% 
                     count(contrast, name = "n_genes") )
    return(newdf)
  }
  ) 
  return(np_DGE)
}


DGE_list <- subset_by_sex(limma_list)

# DE in > 6 clusters/subtypes, get top genes DE in the most clusters
topDEGs <- lapply(DGE_list, \(lst) {
  
  subset_genes <- subset(lst$num_clust, n_clust > 6)
  top_genes <- subset(lst$results, gene %in% subset_genes$gene) %>% 
    left_join(lst$num_clust,
              by = "gene") %>% 
    group_by(contrast) %>%
    slice_max(order_by = n_clust, n = 7)
  
  results <- subset(lst$results, gene %in% top_genes$gene)
  
  lst.df <- list( results = results,
                  num_clust = results %>%
                    distinct(gene, contrast) %>%
                    count(gene, name = "n_clust"),      
                  num_genes = results %>% 
                    distinct(contrast, gene) %>% 
                    count(contrast, name = "n_genes") )
  
  
  return(lst.df)
}
)

# make df with contrast(cluster), sex, status, gene, logFC, num_genes, and n_clust
topDEGs <- lapply(topDEGs, \(lst) {
  df = lst$results %>% 
    left_join(lst$num_clust,
              by = "gene") %>% 
    left_join(lst$num_genes,
              by = "contrast")
})

topDEGs.test <- topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) 

topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) %>% 
  mutate(contrast = reorder_within(contrast, -n_genes, Sex),   
         gene = reorder_within(gene, n_clust, Sex)) %>%
  ggplot(aes(x = contrast, 
             y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey90") +
  labs(x = "neuronal subgroups", 
       y = "top DEGs") +
  scale_fill_gradient2(name = bquote(~italic(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(topDEGs.test, Status, logFC),
            aes(x = 1, 
                y = 1, 
                fill = Status), 
            alpha = 0) +
  scale_fill_manual(name = paste0(sprintf("\u2191 \u2191"), " expr in"),
                    values = c("Dom" = muted("blue"), "Sub" = muted("red"))) +
  guides(fill = guide_legend(override.aes = list(values = c(muted("blue"),
                                                            muted("red")),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  facet_wrap(~ Sex, scales = "free") +
  scale_x_reordered() +
  scale_y_reordered() +
  ggtitle("DEGs across clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.65,
                                   face = "italic", size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15, face = "italic"),
        legend.text = element_text(size = 12, face = "bold"))

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/topDEGs_byCluster.svg"),
       width = 12, height = 6) 


# overlap with cluster marker genes 
neuron_clusterMarkers <- read_csv('neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')

sig_maleDEGs <- maleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

sig_femaleDEGs <- femaleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

conc_genes = intersect(sig_femaleDEGs$gene,
                       sig_maleDEGs$gene)

result <- limma_list %>%
  .[str_detect(names(.), "_cluster\\d+$")] %>%
  imap_dfr(~ {
    .x %>%
      filter(gene %in% conc_genes) %>%
      select(gene, logFC, P.Value) %>%
      mutate(
        sex = str_to_lower(str_sub(.y, 1, 1)),  
        cluster = str_extract(.y, "(?<=_cluster)\\d+")
      )
  }) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(logFC, P.Value),
    names_glue = "{.value}.{sex}"
  ) %>%
  filter(!is.na(logFC.f) & !is.na(logFC.m)) %>%
  group_by(cluster) %>%
  filter(all(!is.na(logFC.f)) & all(!is.na(logFC.m))) %>%
  ungroup() %>% 
  mutate(
    logFC.f   = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, logFC.f, NA),
    logFC.m   = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, logFC.m, NA),
    P.Value.f = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, P.Value.f, NA),
    P.Value.m = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, P.Value.m, NA),
    avg.logFC = (abs(logFC.f) + abs(logFC.m)) / 2 
  ) %>% 
  filter(!(is.na(logFC.f) | is.na(logFC.m)))

neuron_clusterMarkers_conc <- neuron_clusterMarkers %>% 
  filter(gene %in% result$gene) 

result.test = left_join(result, neuron_clusterMarkers_conc[, c("gene", "specificity")]) 

result.test %>% 
  # slice_max(avg.logFC, n = 6) %>% 
  dplyr::select(gene, cluster, logFC.f, logFC.m, avg.logFC, specificity) %>%
  complete(gene, cluster,
           fill = list(avg.logFC = NA)) %>%
  group_by(cluster) %>% 
  mutate(avg.logFC = ifelse(logFC.f * logFC.m < 0, -1 * avg.logFC, avg.logFC),
         direction = ifelse(avg.logFC > 0, "concordant", "discordant")) -> result.test

result.test %>% 
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  scale_fill_gradient2(low = "#408D8E",
                       high = "#7e38b7",
                       na.value = "gray95",
                       name = bquote(avg~abs(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(drop_na(result.test), avg.logFC, direction),
            aes(x = 1,
                y = 1,
                fill = direction),
            alpha = 0) +
  scale_fill_manual(name = "effect of status across sex",
                    values = c("concordant" = "#7e38b7", "discordant" = "#408D8E")) +
  guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
                                                            "#408D8E"),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()

ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/cluster_marker.DEGs.svg"),
       device = "svg", width = 9.25, height = 16)


# color by enrichment quadrant from RRHO to indicate direction
result2 <- rrho_results %>%
  # keep only elements that look like clusters
  .[str_detect(names(.), "^cluster\\d+$")] %>%
  imap_dfr(~ {
    cluster_name <- .y
    df <- .x$df
    
    # iterate over quadrants
    map_dfr(names(.x$RedRibbon.quads), function(q) {
      pos <- .x$RedRibbon.quads[[q]]$positions
      
      df[pos, ] %>%
        transmute(
          gene,
          cluster = str_extract(cluster_name, "\\d+"),
          quadrant = q,
          avg.logFC = (abs(a) + abs(b)) / 2
        ) 
    })
  })

neuron_clusterMarkers_hiSpec <- neuron_clusterMarkers %>% 
  filter(specificity > 0.1)

result2 <- result2 %>% 
  filter(gene %in% neuron_clusterMarkers_hiSpec$gene &
           avg.logFC > 1) %>% 
  complete(gene, cluster,
           fill = list(avg.logFC = NA))


result2 %>% 
  filter(is.na(quadrant)) %>% 
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  scale_fill_viridis(na.value = "gray95") +
  geom_tile(data = subset(result2, quadrant == "upup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#7e38b7",
                      na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'upup' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "downdown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#408D8E",
                      na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'downdown' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "updown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#f94449",
                      na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'updown' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "downup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#ff7d00",
                      na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'downup' quadrant") +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  # ggnewscale::new_scale_fill() +
  # geom_tile(data = distinct(drop_na(result2), avg.logFC, quadrant),
  #           aes(x = 1,
  #               y = 1,
  #               fill = quadrant),
  #           alpha = 0) +
  # scale_fill_manual(name = "effect of status across sex",
  #                   values = c("upup" = "#7e38b7",
  #                              "downdown" = "#408D8E",
  #                              "updown" = "#f94449",
  #                              "downup" = "#ff7d00")) +
  # guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
  #                                                           "#408D8E",
  #                                                           "#f94449",
  #                                                           "#ff7d00")),
  #                                                alpha = c(1,1),
  #                                                shape = c(1,1))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()

ggsave(paste0(root.dir, "/manuscriptFigures/main_figs/cluster_marker.quads.svg"),
       device = "svg", width = 9.25, height = 16)
