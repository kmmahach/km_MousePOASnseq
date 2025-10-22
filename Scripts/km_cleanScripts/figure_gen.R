#### Set up environment ####
net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# load functions 
source("./functions/DGE_fun.R")
source("./functions/QC_filtering_fun.R")

load_packages(c("tidyverse", "limma", "tidytext", "ggrepel", "ggplot2", "cowplot", "ggpol", "scales", "gghalves",
                "svglite", "viridis", "patchwork", "Seurat", "RRHO2", "RedRibbon", "dittoSeq", "ggraph", "ggtext",
                "clusterProfiler", "enrichplot", "HDtest", "PEtests", "doParallel"),
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

# color by direction
overlapping_genes <- limma_list %>%
  .[str_detect(names(.), "_cluster\\d+$")] %>%
  imap_dfr(~ {
    .x %>%
      mutate(
        sex = str_to_upper(str_sub(.y, 1, 1)),        # "M" or "F"
        cluster = str_extract(.y, "(?<=_cluster)\\d+"),
        direction = case_when(
          P.Value < 0.05 & logFC >  0.2  ~ "dom",
          P.Value < 0.05 & logFC < -0.2  ~ "sub",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(direction)) %>%
      select(gene, sex, cluster, direction, logFC)
  })


result <- overlapping_genes %>%
  select(gene, sex, cluster, logFC) %>%
  pivot_wider(
    names_from = sex,
    values_from = logFC,
    names_prefix = "logFC."
  ) %>%
  mutate(
    both_present = !is.na(logFC.M) & !is.na(logFC.F),
    dir_M = ifelse(both_present, ifelse(logFC.M > 0, "up", "down"), NA),
    dir_F = ifelse(both_present, ifelse(logFC.F > 0, "up", "down"), NA),
    quadrant = ifelse(both_present, paste(dir_M, dir_F, sep = ""), NA),
    avg.logFC = ifelse(both_present, 
                       rowMeans(abs(select(., logFC.M, logFC.F)), na.rm = TRUE),
                       NA)
  ) %>%
  # filter(both_present) %>% 
  select(gene, cluster, quadrant, avg.logFC)


neuron_clusterMarkers_hiSpec <- neuron_clusterMarkers %>%
  group_by(cluster) %>%
  mutate(top10 = quantile(specificity, probs = 0.95)) %>%
  ungroup() %>%
  filter(specificity >= top10)


result2 <- result %>% 
  filter(gene %in% neuron_clusterMarkers_hiSpec$gene) %>% 
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

# heatmap of cluster markers by avg DE - figure 2D
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
                      high = "#4b306e",
                      na.value = "gray95",
                      limits = c(0.5, 2.2),
                      name = bquote(plain(" \u2191") * bold("Dom") *
                                      "\u2640  \u2191" * bold("Dom ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "downdown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#ffa600",
                      na.value = "gray95",
                      limits = c(0.5, 2.2),
                      name = bquote(plain(" \u2191") * bold("Sub") *
                                      "\u2640  \u2191" * bold("Sub ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "updown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#005249",
                      na.value = "gray95",
                      limits = c(0.5, 2.2),
                      name = bquote(plain(" \u2191") * bold("Dom") *
                                      "\u2640  \u2191" * bold("Sub ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result3, quadrant == "downup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                      high = "#8a124b",
                      na.value = "gray95",
                      limits = c(0.5, 2.2),
                      name = bquote(plain(" \u2191") * bold("Sub") *
                                      "\u2640   \u2191" * bold("Dom ") * "\u2642  " *
                                      plain("avg. ") * italic(log[2]*FC))) +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(drop_na(result3), avg.logFC, quadrant),
            aes(x = 1,
                y = 1,
                fill = quadrant),
            alpha = 0) +
  scale_fill_manual(name = "effect of status across sex",
                    values = c("upup" = "#4b306e",
                               "downdown" = "#ffa600",
                               "updown" = 	"#005249",
                               "downup" = "#8a124b")) +
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

# cluster composition by group + number of cells - figure 2E
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


#### Figure 3 ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# fig 3A
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda") 
neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"
# subset with RNA counts
DefaultAssay(int.ldfs) = "RNA"
# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

subset_by_sex <- function(limma_list,
                          subset_genes = top25) {
  
  genes <- subset_genes$gene
  
  sublist <- lapply(limma_list, \(lst) {
    lst <- subset(lst, lst$gene %in% genes)
  }
  )
  
  Females <- subset(sublist, grepl("^[F]", names(sublist))==TRUE)  
  Males <- subset(sublist, grepl("^[M]", names(sublist))==TRUE)  
  
  newlist <- as_named_list(Females, Males)
  
  np_DGE <- lapply(newlist, \(x) {
    x <- bind_rows(x, .id = "contrast") %>% 
      left_join(subset_genes, 
                by = 'gene') %>% 
      mutate(logFC = ifelse(P.Value < 0.05, logFC, 0)) %>% 
      mutate(contrast = sub(".*_", "", contrast))
  }
  ) 
  return(np_DGE)
}


# genes
neuropeptides %>% 
  rename(gene = Gene.name) %>% 
  mutate(gene = str_to_title(gene)) %>% 
  filter(gene %in% rownames(int.ldfs)) -> neuropeptides.genes

unlist(as.vector(neuropeptides.genes$gene)) -> neuropeptides.genes

# subset Seurat object by neuropeptide.genes
sub.MSCneurons <- subset_by_gene(int.ldfs,
                                 neuropeptides.genes,
                                 slot = "counts",
                                 min_count = 2)

# get presence/absence (1/0) for neuropeptide.genes
umap = data.frame(int.ldfs@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(int.ldfs@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

sum_neur <- select_neur %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

# get top 25 NPs 
sum_neur %>% 
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "cell_count"
  ) %>% 
  dplyr::select(gene, cell_count) %>% 
  filter(cell_count >= quantile(cell_count, 0.75)) -> top25

np_DGE <- subset_by_sex(limma_list) %>% 
  bind_rows(., .id = "sex") %>% 
  mutate(group = ifelse(logFC > 0, "Dom", "Sub"))

ggplot(np_DGE, aes(x = reorder(gene, -cell_count), 
                   y = contrast,
                   fill = logFC)) +
  labs(x = "",
       y = "neuronal subgroup") +
  geom_tile(color = "grey90") + 
  scale_fill_gradient2(name = bquote(~italic(log[2]~FC)),
                       low = "#4b306e",
                       high = "#8a3324") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(np_DGE, group, logFC),
            aes(x = 1, y = 1, 
                fill = group), 
            alpha = 0) +
  scale_fill_manual(name = paste0(sprintf("\u2191 \u2191"), " expr in"),
                    values = c("Dom" = "#8a3324", 
                               "Sub" = "#4b306e")) +
  guides(fill = guide_legend(override.aes = 
                               list(values = c("#4b306e",
                                               "#8a3324"),
                                    alpha = c(1,1),
                                    shape = c(1,1)))) +
  theme_classic() +
  facet_wrap(~ sex) +
  ggtitle("Differential Expression of Neuropeptides in Neuronal Nuclei") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.75,
                                   face = "italic",
                                   size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15,
                                    face = "italic"),
        legend.text = element_text(size = 12,
                                   face = "bold"))


ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/DEGneuropeptides_byCluster.svg'),
       width = 10, height = 5)


# fig 3B
# load data
load(paste0(root.dir, "/Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda"))
DefaultAssay(MSCneurons.reclust) = 'RNA'
neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))

neuropeptides %>% 
  mutate(Gene.name = str_to_title(Gene.name)) %>% 
  filter(Gene.name %in% rownames(MSCneurons.reclust)) %>% 
  pull(Gene.name) -> np.genes


sub.MSCneurons = subset_by_gene(MSCneurons.reclust,
                                np.genes, 
                                "counts",
                                min_count = 2)


# get presence/absence (1/0) for neuropeptide.genes
umap = data.frame(MSCneurons.reclust@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(MSCneurons.reclust@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

sum_neur_indiv <- select_neur %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(total.neuron.count = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(Freq = ~ sum(.x),
                        Percent = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")


# format data
sum_neur_indiv %>%
  pivot_longer(
    cols = matches("^(Percent|Freq)_"),
    names_to = c(".value", "neuropeptide"),
    names_pattern = "(Percent|Freq)_(.*)"
  ) %>%
  dplyr::select(orig.ident,
                indiv_genotype,
                neuropeptide,
                Freq,
                total.neuron.count,
                Percent) %>%
  mutate(Sex =
           ifelse(grepl("female", orig.ident),
                  "Female",
                  "Male")) %>%
  mutate(Status =
           ifelse(grepl("dom", orig.ident),
                  "Dom",
                  "Sub")) -> dat.for.stats



dat.for.stats %>% 
  group_by(Sex,
           Status,
           neuropeptide) %>% 
  summarise(percent.avg = mean(Percent)) -> avg.NPs_per_neuron


neuron.np.pairwise <- read.csv("neurons/neuropeptides/stats/neuron_np_pairwise.csv")

neuron.np.pairwise %>% 
  filter(contrast == 'male.dom.data - male.sub.data' & adj.pval <= 0.05) %>% 
  mutate(Male.direction = ifelse(Estimate > 0,
                                 "Male.Dom.Bias",
                                 "Male.Sub.Bias")) %>% 
  dplyr::select(Male.direction,
                neuron.np) %>% 
  dplyr::rename(neuropeptide = neuron.np) -> male.np.results

# plot direction of bias
# males
avg.NPs_per_neuron %>% 
  filter(Sex == 'Male') %>% 
  pivot_wider(id_cols = c('neuropeptide'),
              names_from = 'Status',
              values_from = 'percent.avg') %>% 
  full_join(male.np.results, by = "neuropeptide") %>% 
  mutate(label = ifelse(is.na(Male.direction) & neuropeptide != "Scg2", NA, neuropeptide),
         Male.direction = ifelse(is.na(Male.direction), 'No Bias', Male.direction)) %>% 
  ggplot(aes(x = Sub/100,
             y = Dom/100,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0, 
              linetype = "dotted", 
              linewidth = 1.5) +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(0, 1)) +
  geom_point(aes(color = Male.direction),
             size = 5) +
  geom_text_repel(size = 8, max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = c("#ff7d00", 
                                "#7e38b7",
                                "grey"),
                     name = "Status Bias",
                     labels = c("Male.Dom.Bias" = "Dominant \u2642",
                                "Male.Sub.Bias" = "Subordinate \u2642",
                                "No Bias" = "No Bias")) +
  scale_size(guide = 'none') +
  theme(legend.position = c(0.8, 0.2)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 23.5),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25)) +
  xlab('Subordinate Male Neurons (%)') +
  ylab('Dominant Male Neurons (%)') +
  ggtitle('Male Neuropeptide Expression') 


ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/male_countsGLM_statusBias.svg'),
       height = 8, width = 8.5)



# females
neuron.np.pairwise %>% 
  filter(contrast == 'female.dom.data - female.sub.data' & adj.pval <= 0.05) %>% 
  mutate(Female.direction = ifelse(Estimate > 0,
                                   "Female.Dom.Bias",
                                   "Female.Sub.Bias")) %>% 
  dplyr::select(Female.direction,
                neuron.np) %>% 
  dplyr::rename(neuropeptide = neuron.np) -> female.np.results

avg.NPs_per_neuron %>% 
  filter(Sex == 'Female') %>% 
  pivot_wider(id_cols = c('neuropeptide'),
              names_from = 'Status',
              values_from = 'percent.avg') %>% 
  full_join((female.np.results %>% 
               add_row(Female.direction = "Female.Dom.Bias", neuropeptide = NA)), 
            by = "neuropeptide") %>% 
  mutate(label = ifelse(is.na(Female.direction), NA, neuropeptide),
         Female.direction = ifelse(is.na(Female.direction), 'No Bias', Female.direction)) %>% 
  ggplot(aes(x = Sub/100,
             y = Dom/100,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0, 
              linetype = "dotted", 
              linewidth = 1.5) +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(0, 1)) +
  geom_point(aes(color = Female.direction),
             size = 5) +
  geom_text_repel(size = 8, max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "grey"),
                     name = "Status Bias",
                     labels = c("Female.Dom.Bias" = paste0("Dominant \u2640"),
                                "Female.Sub.Bias" = "Subordinate \u2640",
                                "No Bias" = "No Bias")) +
  scale_size(guide = 'none') +
  theme(legend.position = c(0.8, 0.2)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 23.5),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25)) +
  xlab('Subordinate Female Neurons (%)') +
  ylab('Dominant Female Neurons (%)') +
  ggtitle('Female Neuropeptide Expression') 

ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/female_countsGLM_statusBias.svg'),
       height = 8, width = 8.5)

# distributions
select_neur[,c(colnames(select_neur) %in% c(np.genes, 
                                            "orig.ident", 
                                            "indiv_genotype"))] %>% 
  mutate(count.expr = rowSums(across(where(is.numeric)))) %>% 
  dplyr::select(c("orig.ident", "indiv_genotype", "count.expr")) -> np.coexpr.counts

np.coexpr.counts %>% 
  mutate(count.expr = as.character(count.expr)) %>% 
  group_by(orig.ident, indiv_genotype, count.expr) %>% 
  table() %>% as.data.frame() -> np.coexpr.summary

np.coexpr.summary %>% 
  group_by(orig.ident, indiv_genotype) %>% 
  mutate(total.neur = Reduce(`+`, Freq),
         Percent = 100*(Freq/total.neur),
         Cum.Percent = cumsum(Percent)) %>% 
  ungroup() -> np.coexpr.summary



# percent
np.coexpr.summary %>%
  group_by(orig.ident, count.expr) %>% 
  mutate(avg.Percent = mean(Percent)) %>% 
  ggplot() +
  geom_point(aes(x = count.expr,
                 y = avg.Percent/100,
                 group = orig.ident,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3,
             shape = 3) +
  geom_smooth(aes(x = count.expr,
                  y = avg.Percent/100,
                  group = orig.ident,
                  color = orig.ident,
                  fill = orig.ident),
              se = FALSE,
              span = 0.4) + 
  geom_point(data = np.coexpr.summary,
             aes(x = count.expr,
                 y = Percent/100,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3, alpha = 0.5) +
  geom_line(data = np.coexpr.summary,
            aes(x = count.expr,
                y = Percent/100,
                color = orig.ident)) +
  labs(y = "Percent of Neurons (exclusive)",
       x = "Number of Co-expressed Neuropeptide Genes") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "#ff7d00", 
                                "#7e38b7"),
                     name = "Group",
                     labels = c("Female Dom",
                                "Female Sub",
                                "Male Dom",
                                "Male Sub")) +
  scale_fill_manual(values = c("#f94449", 
                                "#408D8E",
                                "#ff7d00", 
                                "#7e38b7"),
                    name = "Group",
                    labels = c("Female Dom",
                               "Female Sub",
                               "Male Dom",
                               "Male Sub"))


ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/coexprNPs_distribution.svg'),
       height = 6, width = 10)



# create cumulative freq plot
np.coexpr.summary %>% 
  group_by(orig.ident, count.expr) %>% 
  mutate(avg.Cum.Percent = mean(Cum.Percent)) %>% 
  ggplot() +
  geom_point(aes(x = count.expr,
                 y = avg.Cum.Percent/100,
                 group = orig.ident,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3,
             shape = 3) +
  geom_smooth(aes(x = count.expr,
                  y = avg.Cum.Percent/100,
                  group = orig.ident,
                  color = orig.ident,
                  fill = orig.ident),
              se = FALSE,
              span = 0.5) + 
  geom_point(data = np.coexpr.summary,
             aes(x = count.expr,
                 y = Cum.Percent/100,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3, alpha = 0.5) +
  geom_line(data = np.coexpr.summary,
            aes(x = count.expr,
                y = Cum.Percent/100,
                color = orig.ident)) +
  labs(y = "Percent of Neurons (cumulative)",
       x = "Number of Co-expressed Neuropeptide Genes") +
  theme_classic() + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "#ff7d00", 
                                "#7e38b7"),
                     name = "Group",
                     labels = c("Female Dom",
                                "Female Sub",
                                "Male Dom",
                                "Male Sub")) +
  scale_fill_manual(values = c("#f94449", 
                               "#408D8E",
                               "#ff7d00", 
                               "#7e38b7"),
                    name = "Group",
                    labels = c("Female Dom",
                               "Female Sub",
                               "Male Dom",
                               "Male Sub"))

ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/coexprNPs_distribution_cumulative.svg'),
       height = 6, width = 10)


# fig 3C
load("./data/integrated_seurat_onlyNeurons.rda")
DefaultAssay(MSCneurons.reclust) = 'RNA'

neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))

neuropeptides %>% 
  mutate(Gene.name = str_to_title(Gene.name)) %>% 
  filter(Gene.name %in% rownames(MSCneurons.reclust)) %>% 
  pull(Gene.name) -> np.genes


sub.MSCneurons = subset_by_gene(MSCneurons.reclust,
                                np.genes, 
                                "counts",
                                min_count = 2)

# get presence/absence (1/0) for neuropeptide.genes
umap = data.frame(MSCneurons.reclust@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(MSCneurons.reclust@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

select_neur %>% 
  rowwise() %>% 
  filter(sum(c_across(intersect(np.genes, colnames(.)))) > 0) -> filtered_neur

select_neur %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  pivot_longer(cols = starts_with("pct_"),
               names_to = "gene") %>% 
  filter(value > 1) %>% 
  mutate(gene = sub("pct_", "", gene)) %>% pull(gene) -> topNP 


filtered_neur %>% 
  column_to_rownames(var = "Cell_ID") %>% 
  dplyr::select(all_of(topNP)) -> filtered_neur_topNP


# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes.all <- (scale01(apply(t(filtered_neur_topNP), 1, mean)) + 1.0) 


# load in plot object
setwd(paste0(root.dir, "/manuscriptFigures"))
load("./plotlist.rda") 


ggsave(file = "./main_figs/coexpression_networks.svg",
       wrap_plots(plotlist, ncol = 2),
       width = 8.5, height = 8)

#### Figure 4 ####
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

# subset with norm RNA data
DefaultAssay(MSCneurons.reclust) = "RNA"
MSCneurons.reclust = NormalizeData(MSCneurons.reclust)


# box plots of norm counts in all neurons
np.sub <- subset_by_gene(MSCneurons.reclust,
                         subset_genes = c("Oxt", "Avp", "Scg2"),
                         min_count = 2)


np.expr <- lapply(np.sub, \(obj) {
  df = as.data.frame(obj@assays$RNA@data) %>% 
    filter(rownames(.) == obj@project.name) %>% 
    pivot_longer(col = everything(),
                 names_to = "Cell_ID",
                 values_to = "expr") 
  
  meta = as.data.frame(obj@meta.data) %>% 
    filter(Cell_ID %in% df$Cell_ID) %>% 
    select(orig.ident, indiv_genotype, Cell_ID)
  
  expr.df = left_join(df, meta, 
                      by = "Cell_ID")
  
  expr.df %>% 
    group_by(orig.ident, indiv_genotype) %>% 
    summarise(avg_expr = mean(expr))
  
}) %>% 
  bind_rows(.id = 'gene') %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))


np.expr %>% 
  ggplot(aes(orig.ident, avg_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 1.75, size = .5, errorbar.draw = TRUE,
                 position = position_dodge(2)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"),
                     name = "Group") +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"),
                    name = "Group") +
  facet_wrap(~ gene, nrow = 1) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Counts", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 12)) + 
  ggtitle("Expression in Neuronal Nuclei with >1 read") 

ggsave("manuscriptFigures/main_figs/top3_NPs_NormCounts_RNAdata.svg",
       width = 7, height = 3)


DotPlot(MSCneurons.reclust,
        features = c("Oxt", "Avp", "Scg2"),
        group.by = 'orig.ident',
        min_count = 2,
        col.min = 0,
        scale.min = 0,
        scale.max = 100,
        scale = FALSE) +
  scale_y_discrete(labels = sub("\\.data", "",
                                unique(MSCneurons.reclust$orig.ident))) +
  theme_classic() +
  labs(y = "", x = "") +
  ggtitle('Expression in Neurons') +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 12,
                                 face = "bold"),
        axis.text.y = element_text(size = 12))

ggsave('manuscriptFigures/main_figs/top3_NPs_dotplot.svg',
       height = 3.5, width = 5.75)


# percent of neurons expressing each combination:
# avp only, oxt only, scg2 only; avp + oxt, avp + scg2, oxt + scg2; oxt + avp + scg2
top3np <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp", "Scg2")) %>% 
  t() %>% as.data.frame() %>% 
  mutate(Oxt = ifelse(Oxt > 1, 1, 0),
         Avp = ifelse(Avp > 1, 1, 0),
         Scg2 = ifelse(Scg2 > 1, 1, 0),
         Oxt_Avp = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 0, 1, 0),
         Oxt_Scg2 = ifelse(Oxt == 1 & Avp == 0 & Scg2 == 1, 1, 0),
         Avp_Scg2 = ifelse(Oxt == 0 & Avp == 1 & Scg2 == 1, 1, 0),
         Oxt_Avp_Scg2 = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 1, 1, 0)) %>% 
  rownames_to_column(var = "Cell_ID") 

meta = as.data.frame(MSCneurons.reclust@meta.data) %>% 
  filter(Cell_ID %in% top3np$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, Cell_ID) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))


top3np.excl <- left_join(top3np, meta,
                         by = "Cell_ID") %>% 
  mutate(Oxt = ifelse(Oxt == 1 &
                        rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Oxt),
         Avp = ifelse(Avp == 1 &
                        rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Avp),
         Scg2 = ifelse(Scg2 == 1 &
                         rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Scg2))


# plot pcts for each NP
sum_np_group <- top3np.excl %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n(),
            across(all_of(c(2:8)),
                   list(sum = ~ sum(.x)),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

# split by oxt/avp/scg2
# use for combined legend before saving
labels = c("sum_Scg2"="Scg2 only",
           "sum_Oxt"="Oxt only",
           "sum_Avp"="Avp only",
           "sum_Oxt_Avp"="Oxt + Avp",
           "sum_Oxt_Scg2"="Oxt + Scg2",
           "sum_Avp_Scg2"="Avp + Scg2",
           "sum_Oxt_Avp_Scg2"="Oxt + Avp + Scg2")

# oxytocin 
sum_oxt_group = sum_np_group %>% 
  mutate(across(!contains(c("Oxt", "cell_count", "orig.ident")), ~ 0)) %>% 
  rowwise() %>% 
  mutate(cell_count = sum(c_across(starts_with("sum")), 
                          na.rm = TRUE))

oxt.temp_labels = sum_oxt_group[,c("orig.ident", "cell_count")]

sum_oxt_group %>% 
  select(!"cell_count") %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = oxt.temp_labels,
            aes(x = orig.ident,
                y = Inf,
                label = paste0('n = ', format(cell_count, big.mark = ',', trim = TRUE)),
                vjust = -0.5),
            color = 'black',
            size = 5) +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  scale_fill_viridis_d(option = "plasma", labels = labels, end = 0.9) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Co-expression") +
  ggtitle("Oxytocin Neurons") +
  theme(legend.position = 'left',
        plot.title = element_text(hjust = 0.5, vjust = 5),
        text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1,
                                   size = 25),
        plot.margin = margin(t = 1, 
                             r = 1, 
                             b = 1, 
                             l = 1, 
                             unit = 'cm')) -> p1


# vasopressin
sum_avp_group = sum_np_group %>% 
  mutate(across(!contains(c("Avp", "cell_count", "orig.ident")), ~ 0)) %>% 
  rowwise() %>% 
  mutate(cell_count = sum(c_across(starts_with("sum")), 
                          na.rm = TRUE))

avp.temp_labels = sum_avp_group[,c("orig.ident", "cell_count")]

sum_avp_group %>% 
  select(!"cell_count") %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = avp.temp_labels,
            aes(x = orig.ident,
                y = Inf,
                label = paste0('n = ', format(cell_count, big.mark = ',', trim = TRUE)),
                vjust = -0.5),
            color = 'black',
            size = 5) +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  scale_fill_viridis_d(option = "plasma", labels = labels, end = 0.9) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Co-expression") +
  ggtitle("Vasopressin Neurons") +
  theme(legend.position = 'left',
        plot.title = element_text(hjust = 0.5, vjust = 5),
        text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1,
                                   size = 25),
        plot.margin = margin(t = 1, 
                             r = 1, 
                             b = 1, 
                             l = 1, 
                             unit = 'cm')) -> p2


# Secretogranin2
sum_scg2_group = sum_np_group %>% 
  mutate(across(!contains(c("Scg2", "cell_count", "orig.ident")), ~ 0)) %>% 
  rowwise() %>% 
  mutate(cell_count = sum(c_across(starts_with("sum")), 
                          na.rm = TRUE))

scg2.temp_labels = sum_scg2_group[,c("orig.ident", "cell_count")]

sum_scg2_group %>% 
  select(!"cell_count") %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = scg2.temp_labels,
            aes(x = orig.ident,
                y = Inf,
                label = paste0('n = ', format(cell_count, big.mark = ',', trim = TRUE)), 
                vjust = -0.5),
            color = 'black',
            size = 5) +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  scale_fill_viridis_d(option = "plasma", labels = labels, end = 0.9) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Co-expression") +
  ggtitle("Secretogranin2 Neurons") +
  theme(legend.position = 'left',
        plot.title = element_text(hjust = 0.5, vjust = 5),
        text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1,
                                   size = 25),
        plot.margin = margin(t = 1, 
                             r = 1, 
                             b = 1, 
                             l = 1, 
                             unit = 'cm')) -> p3


ggsave(paste0(root.dir, '/manuscriptFigures/main_figs/coexpression_neuronComposition.png'),
       p1 + p2 + p3 +
         plot_layout(ncol = 3, axes = "collect", guides = "collect"),
       width = 26, height = 8)

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
plot.filters(l.dfs, filters, ".svg",
             outdir = paste0(root.dir, "/manuscriptFigures/supplemental/"))


#### Supplemental Figure 2 ####
load("./data/integrated_seurat_withHypoMap.rda")

# across samples cluster
DimPlot(int.ldfs,
        reduction = "umap",
        split.by = 'orig.ident',
        pt.size = 1, 
        cols = c('#d95f02',
                 '#f8766d',
                 '#7570b3',
                 '#e7298a',
                 '#1b9e77',
                 '#00bdd0',
                 '#00b0f6',
                 '#a6761d',
                 '#ac88ff',
                 '#a3a500',
                 '#66a61e',
                 '#e76bf3',
                 '#7997ff',
                 '#9acd32',
                 '#ff61c9',
                 '#00bc59',
                 '#e6ab02',
                 '#00008b',
                 '#341539',
                 'grey',
                 '#bb9d00',
                 '#ed8141')) +
  theme(text = element_text(size = 18)) -> d1

d1$data$orig.ident = sub("\\.data$", ".clusters", d1$data$orig.ident)

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/domVsub.clusters.dimplot.comparison.all.svg"),
       d1, width = 20, height = 5)

# across samples genotype
DimPlot(int.ldfs,
        reduction = "umap",
        split.by = 'orig.ident',
        group.by = 'indiv_genotype',
        pt.size = 1) +
  theme(text = element_text(size = 18)) +
  ggtitle(NULL) -> d2

d2$data$orig.ident = sub("\\.data$", ".indiv_genotypes", d2$data$orig.ident)

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/domVsub.genotype.dimplot.comparison.all.svg"),
       d2, width = 20, height = 5)

# across samples nfeature
FeaturePlot(int.ldfs,
            reduction = "umap",
            split.by = 'orig.ident',
            features = 'nFeature_SCT',
            pt.size = 1,
            combine = TRUE) +
  plot_layout(axes = "collect",
              guides = "collect") &
  labs(color = "nFeature_SCT") &
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title.y.right = element_blank()) -> d3

for (i in 1:length(d3)) {
  d3[[i]]$labels$title <- sub("\\.data$", ".feature_count", d3[[i]]$labels$title)
}

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/domVsub.nfeature.dimplot.comparison.all.svg"),
       d3, width = 20, height = 5)


#### Supplemental Figure 3 ####
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# supp fig 3A
DimPlot(int.ldfs,
        group.by = "parent_id.broad.prob",
        cols = c('#1b9e77',
                 '#d95f02',
                 '#7570b3',
                 '#e7298a',
                 '#66a61e',
                 '#e6ab02',
                 '#a6761d',
                 'grey')) +
  NoAxes() 

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/predicted.broad.UMAP.svg'),
       height = 5, width = 6)  

# supp fig 3B
DimPlot(int.ldfs, 
        reduction = "umap", 
        group.by = 'sctype.integrated', 
        cols = c('#7570b3',
                 '#a6761d',
                 '#1b9e77',
                 '#66a61e',
                 '#e6ab02',
                 '#d95f02',
                 '#e7298a',
                 '#341539',
                 'grey')) + 
  NoAxes()


ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/sctype.dimplot.colormatch.svg'),
       height = 5, width = 6)  

# supp fig 3C
DimPlot(int.ldfs,
        reduction = "umap",
        cols = c('#d95f02',
                 '#f8766d',
                 '#7570b3',
                 '#e7298a',
                 '#1b9e77',
                 '#00bdd0',
                 '#00b0f6',
                 '#a6761d',
                 '#ac88ff',
                 '#a3a500',
                 '#66a61e',
                 '#e76bf3',
                 '#7997ff',
                 '#9acd32',
                 '#ff61c9',
                 '#00bc59',
                 '#e6ab02',
                 '#00008b',
                 '#341539',
                 'grey',
                 '#bb9d00',
                 '#ed8141')) +
  ggtitle("integrated.clusters") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/clusters.dimplot.integrated.svg"),
       height = 5, width = 5.75)

# supp fig 3D
DimPlot(int.ldfs,
        reduction = "umap",
        group.by = "orig.ident") +
  NoAxes() -> d0

d0$data$orig.ident = sub("\\.data$", "", d0$data$orig.ident)

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/dimplot.ident.integrated.svg"),
       d0, height = 5, width = 6)

#### Supplemental Figure 4 ####
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

### Number of reads and genes per cluster - supp fig 4A
MSCneurons.reclust@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() -> temp_labels

p1 <- ggplot() +
  geom_half_violin(data = MSCneurons.reclust@meta.data, 
                   aes(integrated_snn_res.0.4, nCount_RNA, fill = integrated_snn_res.0.4),
                   side = 'l', show.legend = FALSE, trim = FALSE) +
  geom_half_boxplot(data = MSCneurons.reclust@meta.data, 
                    aes(integrated_snn_res.0.4, nCount_RNA, fill = integrated_snn_res.0.4),
                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE) +
  geom_text(data = temp_labels, 
            aes(x = integrated_snn_res.0.4, 
                y = -Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1), 
            color = 'black', 
            size = 2.8) +
  scale_y_continuous(name = 'Number of transcripts', 
                     labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("nCount_RNA") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot() +
  geom_half_violin(data = MSCneurons.reclust@meta.data, 
                   aes(integrated_snn_res.0.4, nFeature_RNA, fill = integrated_snn_res.0.4),
                   side = 'l', show.legend = FALSE, trim = FALSE) +
  geom_half_boxplot(data = MSCneurons.reclust@meta.data, 
                    aes(integrated_snn_res.0.4, nFeature_RNA, fill = integrated_snn_res.0.4),
                    side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE) +
  geom_text(data = temp_labels, 
            aes(x = integrated_snn_res.0.4, 
                y = -Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  scale_y_continuous(name = 'Number of expressed genes', 
                     labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("nFeature_RNA") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/nCount_nFeature_byCluster.svg'), 
       p1 + p2 + plot_layout(ncol = 1),
       height = 8, width = 14)



# hierarchical cluster tree - supp fig 4B
# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

BuildClusterTree(MSCneurons.reclust,
                 dims = 1:15,
                 reorder = FALSE,
                 reorder.numeric = FALSE,
                 slot = 'data',
                 assay = "SCT") -> MSCneurons.reclust

# create tree
neuron.tree <- MSCneurons.reclust@tools$BuildClusterTree

# graph tree of celltypes
ggtree::ggtree(neuron.tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(shape = 16, 
                        size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neuron.cluster.tree.svg'),
       height = 6, width = 6)

# scaled counts of neurons in each cluster (as percent) - supp fig 4C
# subset by select_clusters
Idents(MSCneurons.reclust) <- "seurat_clusters"
select_clusters <- unique(MSCneurons.reclust$seurat_clusters)

sub.MSCneurons <- subset_by_ident(MSCneurons.reclust,
                                  select_clusters,
                                  cluster = TRUE)

# get presence/absence (1/0) for clusters
umap = data.frame(MSCneurons.reclust@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(MSCneurons.reclust@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

# counts & percentages
sum_neur_indiv <- select_neur %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(total.neuron.count = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(Freq = ~ sum(.x),
                        Percent = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

# percent of neurons expr select genes by individual
sum_neur_indiv %>%
  pivot_longer(
    cols = starts_with("Percent_"),
    names_to = "seurat_clusters",
    names_prefix = "Percent_",
    values_to = "Percent"
  ) %>% 
  dplyr::select(orig.ident,
                indiv_genotype,
                seurat_clusters,
                total.neuron.count,
                Percent) %>% 
  mutate(seurat_clusters = str_remove(seurat_clusters, 
                                      "cluster_") %>% 
           as.numeric()) %>% 
  ggplot(aes(x = reorder(as.integer(seurat_clusters),
                         as.numeric(seurat_clusters)),
             y = Percent/100,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1.5) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 4)+
  facet_grid(~seurat_clusters,
             scales = "free_x") +
  theme_classic() +
  xlab("Neuron Cluster") +
  ylab("Percent of Neurons") +
  labs(color = "Group",
       shape = "Genotype") +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 0.5)) +
  theme(legend.position = c(0.93, 0.75),
        legend.key.size = unit(0.75, "cm"),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 23),
        strip.text.x.top = element_blank(),
        axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 25)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"),
                     # labels = c("Female Dom",
                     #            "Female Sub",
                     #            "Male Dom",
                     #            "Male Sub")
                     labels = sub("\\.data", "", unique(sum_neur_indiv$orig.ident)))


ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/clustersVorig.ident.neurons.recluster.genotype.scaled.svg'),
       height = 7, width = 20)



# exc/inh neuron clusters - supp fig 4D
MSCneurons.reclust %>%
  DimPlot(group.by = 'parent_id.broad.prob') +
  ggtitle('Excitatory vs. Inhibitory Neurons') +
  scale_color_manual(values = c('#1b9e77',
                                '#d95f02')) +
  theme_classic()

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster.parent_id.broad.svg'),
       height = 5, width = 6)

# colored by sample pool (orig.ident) - supp fig 4E
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
setwd(paste0(root.dir, "/DGE_CellTypes/")) 

# cluster marker specificity scores - supp fig 5A
neuron.markers.df <- read_csv('neurons/cluster_stats/neuron_clusterMarkers.csv')
load(paste0(root.dir, "/Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda"))

# add weights for neuron proportions to specificity
MSCneurons.reclust@meta.data %>% 
  group_by(seurat_clusters) %>% 
  summarise(n_cells = n()) %>% 
  mutate(prop = n_cells / sum(n_cells),
         seurat_clusters = as.numeric(seurat_clusters) - 1) %>% 
  rename(cluster = seurat_clusters) -> neuron.props

neuron.markers.df %>% 
  left_join(neuron.props,
            by = "cluster") %>% 
  mutate(specificity = ifelse(specificity < 0, 
                              specificity * prop * -1,
                              specificity * prop)) -> markers.weighted.spec

top6_genes = markers.weighted.spec %>% 
  group_by(cluster) %>% 
  slice_max(specificity, 
            n = 6) %>% 
  ungroup() 

top6_genes %>% 
  group_by(gene) %>%
  slice_max(specificity, n = 1) %>%
  dplyr::select(gene, cluster) %>% 
  rename(top_clust = cluster) %>% 
  distinct() -> markers.weighted.summary


markers.weighted.spec %>%
  dplyr::select(cluster,gene,specificity) %>% 
  filter(gene %in% top6_genes$gene) %>% 
  complete(gene, cluster, 
           fill = list(specificity = NA)) %>% 
  left_join(markers.weighted.summary, by = "gene") -> markers.df


markers.df %>% 
  mutate(specificity = ifelse(specificity < 0.02, NA, specificity)) %>%
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = specificity)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  labs(x = "Neuron clusters", y = "Marker genes") +
  scale_fill_viridis(option = "plasma", 
                     na.value = "gray95",
                     name = "Specificity score",
                     direction = -1, 
                     begin = 0.0, end = 0.9) +
  facet_grid(top_clust ~ ., 
             scales = "free_y", 
             space = "free_y") +
  theme_classic() +
  ggtitle("Cluster Marker Genes") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0.5, size = 25),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 11),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 17),
        legend.title = element_text(size = 20, face = "italic"),
        legend.text = element_text(size = 17, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/cluster_markers_spec.svg"),
       width = 7.5, height = 15)


# supp fig 5B - GABA/GLU split by cluster
MSCneurons.reclust@meta.data %>% 
  group_by(seurat_clusters, parent_id.broad) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         cluster = as.numeric(cluster) - 1) -> neuron.props


# cells per cluster by orig.ident
MSCneurons.reclust@meta.data %>% 
  group_by(orig.ident, seurat_clusters, parent_id.broad) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) -> neuron.props.ident


MSCneurons.reclust@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4) -> temp_labels

p2 <- neuron.props.ident %>% 
  ggplot(aes(cluster, n_cells)) +
  geom_bar(aes(fill = parent_id.broad), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_fill_discrete(type =  adjustcolor(c('#1b9e77', '#d95f02'), 
                                          alpha.f = 0.75)) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Subtype",
       x = "Neuron Cluster",
       y = "Number of Nuclei") +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, vjust = -2),
    axis.title.x = element_text(size = 22),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )



MSCneurons.reclust@meta.data %>%
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

MSCneurons.reclust@meta.data %>% 
  group_by(orig.ident, parent_id.broad) %>% 
  summarise(n_cells = n()) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) -> neuron.props.group


p1 <- neuron.props.group %>%
  ggplot(aes(orig.ident, prop)) +
  geom_bar(aes(fill = parent_id.broad), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = temp_labels,
            aes(x = orig.ident, 
                y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  scale_fill_discrete(type =  adjustcolor(c('#1b9e77', '#d95f02'), 
                                          alpha.f = 0.9)) +
  scale_y_continuous(name = 'Percent of Nuclei', 
                     labels = scales::percent_format(), 
                     expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "parent_id.broad") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt'))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/composition_GABAvsGLU_ClustersByNum.svg'),
       p1 + p2 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 16, height = 5.5)

# data for DEG x clusters heatmap - supp fig 5C 
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

# get top 1% of genes DE in the most clusters
topDEGs <- lapply(DGE_list, \(lst) {
  
  subset_genes <- subset(lst$num_clust, 
                         n_clust > quantile(n_clust, probs = 0.99))
  
  top_genes <- subset(lst$results, gene %in% subset_genes$gene) %>% 
    left_join(lst$num_clust,
              by = "gene") %>% 
    group_by(contrast) 
  
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


# supp fig 5D
# combine with spatial data
data.spatial = read.csv(paste0(root.dir, '/HypoMap/data/metadata_with_predictedRegions.csv'))
pred.reg = data.spatial[, c("Cell_ID", "Region_summarized")] %>% 
  column_to_rownames(var = "Cell_ID")

MSCneurons.reclust = AddMetaData(MSCneurons.reclust,
                                 pred.reg)

## add unknown 
### use prediction probability to assign cell type
predicted.prob.df = MSCneurons.reclust@meta.data

# celltype call
predicted.prob.df <- predicted.prob.df %>%
  mutate(predicted.prob = ifelse(prediction_probability >= 0.25,
                                 predicted, 'Unknown'))


# cells per cluster by group
predicted.prob.df %>%
  filter(!(predicted.prob == "Unknown")) %>% 
  group_by(orig.ident, Region_summarized) %>% 
  summarise(n_cells = n()) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  mutate(Region_summarized = as.factor(ifelse(is.na(Region_summarized),
                                              "Unknown", Region_summarized))) -> region.props.group

region.props.group$Region_summarized <- fct_relevel(region.props.group$Region_summarized,
                                                    "Unknown", after = Inf)

MSCneurons.reclust@meta.data %>%
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- region.props.group %>%
  # arrange(Region_summarized == "Unknown") %>% 
  ggplot(aes(orig.ident, prop)) +
  geom_bar(aes(color = Region_summarized,
               fill = Region_summarized), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = temp_labels,
            aes(x = orig.ident, 
                y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  scale_color_viridis_d(option = "viridis",
                        direction = -1) +
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +
  scale_y_continuous(name = 'Percent of Nuclei', 
                     labels = scales::percent_format(), 
                     expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "parent_id.broad") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt'))

MSCneurons.reclust@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4) -> temp_labels


# by ident
predicted.prob.df %>%
  filter(!(predicted.prob == "Unknown")) %>% 
  group_by(orig.ident, seurat_clusters, Region_summarized) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  mutate(Region_summarized = as.factor(ifelse(is.na(Region_summarized),
                                              "Unknown", Region_summarized))) -> region.props.ident

region.props.ident$Region_summarized <- fct_relevel(region.props.ident$Region_summarized,
                                                    "Unknown", after = Inf)

p2 <- region.props.ident %>% 
  ggplot(aes(cluster, n_cells)) +
  geom_bar(aes(color = Region_summarized,
               fill = Region_summarized), 
           position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_color_viridis_d(option = "viridis",
                        direction = -1) +
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(color = "Predicted brain region",
       fill = "Predicted brain region",
       x = "Neuron Cluster",
       y = "Number of Nuclei") +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, vjust = -2),
    axis.title.x = element_text(size = 22),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/composition_predictedRegion_ClustersByNum.svg'),
       p1 + p2 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 20, height = 5.5)



# # overlap with cluster marker genes 
# neuron_clusterMarkers <- read_csv('neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')
# femaleDEGs <- limma_list$Fdom_vs_sub_all
# maleDEGs <- limma_list$Mdom_vs_sub_all
# 
# 
# sig_maleDEGs <- maleDEGs %>% 
#   filter(abs(logFC) > 0.2, P.Value < 0.05)
# 
# sig_femaleDEGs <- femaleDEGs %>% 
#   filter(abs(logFC) > 0.2, P.Value < 0.05)
# 
# conc_genes = intersect(sig_femaleDEGs$gene,
#                        sig_maleDEGs$gene)
# 
# result <- limma_list %>%
#   .[str_detect(names(.), "_cluster\\d+$")] %>%
#   imap_dfr(~ {
#     .x %>%
#       filter(gene %in% conc_genes) %>%
#       select(gene, logFC, P.Value) %>%
#       mutate(
#         sex = str_to_lower(str_sub(.y, 1, 1)),  
#         cluster = str_extract(.y, "(?<=_cluster)\\d+")
#       )
#   }) %>%
#   pivot_wider(
#     names_from = sex,
#     values_from = c(logFC, P.Value),
#     names_glue = "{.value}.{sex}"
#   ) %>%
#   filter(!is.na(logFC.f) & !is.na(logFC.m)) %>%
#   group_by(cluster) %>%
#   filter(all(!is.na(logFC.f)) & all(!is.na(logFC.m))) %>%
#   ungroup() %>% 
#   mutate(
#     logFC.f   = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, logFC.f, NA),
#     logFC.m   = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, logFC.m, NA),
#     P.Value.f = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, P.Value.f, NA),
#     P.Value.m = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, P.Value.m, NA),
#     avg.logFC = (abs(logFC.f) + abs(logFC.m)) / 2 
#   ) %>% 
#   filter(!(is.na(logFC.f) | is.na(logFC.m)))
# 
# neuron_clusterMarkers_conc <- neuron_clusterMarkers %>% 
#   filter(gene %in% result$gene) 
# 
# result.test = left_join(result, neuron_clusterMarkers_conc[, c("gene", "specificity")]) 
# 
# result.test %>% 
#   # slice_max(avg.logFC, n = 6) %>% 
#   dplyr::select(gene, cluster, logFC.f, logFC.m, avg.logFC, specificity) %>%
#   complete(gene, cluster,
#            fill = list(avg.logFC = NA)) %>%
#   group_by(cluster) %>% 
#   mutate(avg.logFC = ifelse(logFC.f * logFC.m < 0, -1 * avg.logFC, avg.logFC),
#          direction = ifelse(avg.logFC > 0, "concordant", "discordant")) -> result.test
# 
# result.test %>% 
#   ggplot(aes(x = factor(cluster),
#              y = gene,
#              fill = avg.logFC)) +
#   geom_tile(color = "white",
#             lwd = 0.8) +
#   labs(x = "Neuron clusters", y = "cluster marker genes") +
#   scale_fill_gradient2(low = "#408D8E",
#                        high = "#7e38b7",
#                        na.value = "gray95",
#                        name = bquote(avg~abs(log[2]~FC))) +
#   ggnewscale::new_scale_fill() +
#   geom_tile(data = distinct(drop_na(result.test), avg.logFC, direction),
#             aes(x = 1,
#                 y = 1,
#                 fill = direction),
#             alpha = 0) +
#   scale_fill_manual(name = "effect of status across sex",
#                     values = c("concordant" = "#7e38b7", "discordant" = "#408D8E")) +
#   guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
#                                                             "#408D8E"),
#                                                  alpha = c(1,1),
#                                                  shape = c(1,1)))) +
#   theme_classic() +
#   ggtitle("Differential expression of cluster markers") +
#   theme(strip.text.y = element_blank(),
#         panel.spacing = unit(1, "mm"),
#         plot.title = element_text(hjust = 0, size = 22),
#         axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
#         axis.title.x = element_text(size = 21),
#         axis.title.y = element_text(size = 21),
#         axis.text.y = element_text(size = 12),
#         axis.title = element_text(size = 15),
#         legend.title = element_text(size = 18, face = "italic"),
#         legend.text = element_text(size = 15, face = "bold"),
#         plot.margin = unit(c(1, 1, 1, 1), "cm")) +
#   coord_fixed()
# 
# ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/cluster_marker.DEGs.svg"),
#        width = 9.25, height = 16)


#### Supplemental Figure 6 ####
# supp fig 6A (UMAP)
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

# specific sub-types - supp fig 6A
MSCneurons.reclust %>%
  DimPlot(group.by = 'predicted.prob') +
  ggtitle('Neuron Subtypes') +
  theme_classic() +
  scale_color_viridis_d(option = "turbo",
                        direction = -1)

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster_predicted.prob.svg'),
       height = 5, width = 11)
f6B

# percent of neurons in each subclass
predicted.prob.df = MSCneurons.reclust@meta.data 

predicted.prob.df %>% 
  group_by(predicted.prob,
           orig.ident,
           indiv_genotype) %>% 
  summarise(total_cells = n()) %>% 
  ungroup() %>%
  group_by(orig.ident, 
           indiv_genotype) %>% 
  mutate(n_cells = sum(total_cells)) %>%  
  ungroup() %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident)) -> df


df %>% 
  group_by(orig.ident, 
           predicted.prob) %>% 
  summarise(group.mean = mean(total_cells/n_cells),
            total_cells = sum(total_cells)) %>% 
  ungroup() %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident)) -> df.bar


ggplot() +
  scale_y_continuous(labels = scales::percent, expand = c(0,0.001)) +
  geom_bar(aes(x = reorder(predicted.prob, -total_cells),
               y = group.mean,
               fill = orig.ident), df.bar,
           stat = "identity",
           position = "dodge") + 
  geom_point(position = position_dodge(0.9),
             aes(x = reorder(predicted.prob, -total_cells),
                 y = total_cells/n_cells,
                 shape = indiv_genotype,
                 group = orig.ident), df,
             size = 2, alpha = 0.75) +
  theme_classic() +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  labs(x = "",
       y = "percent of neurons") +
  theme(axis.text.x = element_text(hjust = 1, 
                                   vjust = 1, 
                                   angle = 45),
        plot.margin = unit(c(1,1,1,1), "cm")) 

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/pct_cells_bySubtype.svg"),
       width = 19, height = 7)

# # color some of the interesting ones to highlight 
# # percent of neurons in each subclass
# predicted.prob.df %>% 
#   group_by(predicted.prob,
#            orig.ident,
#            indiv_genotype) %>% 
#   summarise(total_cells = n()) %>% 
#   ungroup() %>%
#   group_by(orig.ident, 
#            indiv_genotype) %>% 
#   mutate(n_cells = sum(total_cells)) %>%  
#   ungroup() %>% 
#   mutate(orig.ident = sub("\\.data", "", orig.ident)) %>% 
#   mutate(color = ifelse(grepl("Lhx8|Trh|Chat|Meis2|Caprin2|Pmch", .$predicted.prob) == TRUE,
#                         .$orig.ident, NA)) -> df
# 
# 
# df %>% 
#   group_by(orig.ident, 
#            predicted.prob) %>% 
#   summarise(group.mean = mean(total_cells/n_cells),
#             total_cells = sum(total_cells)) %>% 
#   ungroup() %>% 
#   mutate(orig.ident = sub("\\.data", "", orig.ident)) %>% 
#   mutate(color = ifelse(grepl("Lhx8|Trh|Chat|Meis2|Caprin2|Pmch", .$predicted.prob) == TRUE,
#                         .$orig.ident, NA)) -> df.bar
# 
# 
# ggplot() +
#   scale_y_continuous(labels = scales::percent, expand = c(0,0.001)) +
#   geom_bar(aes(x = reorder(predicted.prob, -total_cells),
#                y = group.mean,
#                group = orig.ident,
#                fill = color), df.bar,
#            stat = "identity",
#            position = "dodge") + 
#   geom_point(position = position_dodge(0.9),
#              aes(x = reorder(predicted.prob, -total_cells),
#                  y = total_cells/n_cells,
#                  shape = indiv_genotype,
#                  group = orig.ident), df,
#              size = 2, alpha = 0.75) +
#   theme_classic() +
#   scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
#   labs(x = "",
#        y = "percent of neurons") +
#   theme(axis.text.x = element_text(hjust = 1, 
#                                    vjust = 1, 
#                                    angle = 45),
#         plot.margin = unit(c(1,1,1,1), "cm")) 

# supp fig 6C - percent sub-class by cluster
predicted.prob.df %>%
  filter(!(predicted.prob == "Unknown")) %>% 
  group_by(orig.ident, seurat_clusters, predicted.prob) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  mutate(predicted.prob = as.factor(ifelse(is.na(predicted.prob),
                                           "Unknown", predicted.prob))) -> predicted.prob.ident

predicted.prob.ident$predicted.prob <- fct_relevel(predicted.prob.ident$predicted.prob,
                                                   "Unknown", after = Inf)

MSCneurons.reclust@meta.data %>%
mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- predicted.prob.ident %>%
  ggplot(aes(orig.ident, prop)) +
  geom_bar(aes(color = predicted.prob,
               fill = predicted.prob), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = temp_labels,
            aes(x = orig.ident, 
                y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  scale_color_viridis_d(option = "turbo",
                        direction = -1) +
  scale_fill_viridis_d(option = "turbo",
                       direction = -1) +
  scale_y_continuous(name = 'Percent of Nuclei', 
                     labels = scales::percent_format(), 
                     expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "predicted.prob") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt'))
p1

MSCneurons.reclust@meta.data %>%
  group_by(integrated_snn_res.0.4) %>%
  tally() %>%
  dplyr::rename('cluster' = integrated_snn_res.0.4) -> temp_labels


p2 <- predicted.prob.ident %>% 
  ggplot(aes(cluster, n_cells)) +
  geom_bar(aes(color = predicted.prob,
               fill = predicted.prob), 
           position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_color_viridis_d(option = "turbo",
                        direction = -1) +
  scale_fill_viridis_d(option = "turbo",
                       direction = -1) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(color = "Neuron Sub-class",
       fill = "Neuron Sub-class",
       x = "Neuron Cluster",
       y = "Number of Nuclei") +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, vjust = -2),
    axis.title.x = element_text(size = 22),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )
p2

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/composition_predictedSubtype_ClustersByNum.svg'),
       p1 + p2 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 15, height = 5.5)


# # supp fig 6C alternative - predicted brain regions 
# # combine with spatial data
# data.spatial = read.csv(paste0(root.dir, '/HypoMap/data/metadata_with_predictedRegions.csv'))
# pred.reg = data.spatial[, c("Cell_ID", "Region_summarized")] %>% 
#   column_to_rownames(var = "Cell_ID") %>% 
#   mutate(Region_summarized = ifelse(is.na(Region_summarized), "Unknown", Region_summarized))
# 
# pred.reg$Region_summarized = as.factor(pred.reg$Region_summarized)
# MSCneurons.reclust = AddMetaData(MSCneurons.reclust,
#                                  pred.reg)
# 
# MSCneurons.reclust$Region_summarized = fct_relevel(MSCneurons.reclust$Region_summarized, 
#                                                    "Unknown", after = Inf)
# 
# # left - UMAP
# MSCneurons.reclust %>% 
#   DimPlot(group.by = 'Region_summarized') +
#   ggtitle('Predicted brain region') +
#   scale_color_viridis_d(option = "viridis",
#                         direction = -1) +
#   theme_classic()
# 
# ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/neurons.recluster.Region_summarized.svg'),
#        height = 5, width = 9)
# 
# # right - percent of neurons in each brain region by orig.ident + indiv
# predicted.prob.df %>% 
#   group_by(Region_summarized,
#            orig.ident,
#            indiv_genotype) %>% 
#   summarise(total_cells = n()) %>% 
#   ungroup() %>%
#   group_by(orig.ident, 
#            indiv_genotype) %>% 
#   mutate(n_cells = sum(total_cells)) %>%  
#   ungroup() %>% 
#   mutate(Region_summarized = ifelse(is.na(Region_summarized), "Unknown", 
#                                     Region_summarized),
#          orig.ident = sub("\\.data", "", orig.ident)) -> df
# 
# df %>% 
#   group_by(orig.ident, 
#            Region_summarized) %>% 
#   summarise(group.mean = mean(total_cells/n_cells),
#             total_cells = sum(total_cells)) %>% 
#   ungroup() %>% 
#   mutate(orig.ident = sub("\\.data", "", orig.ident)) -> df.bar
# 
# 
# ggplot() +
#   scale_x_continuous(labels = scales::percent) +
#   geom_bar(aes(x = group.mean,
#                y = reorder(Region_summarized, -total_cells),
#                fill = orig.ident), df.bar,
#            stat = "identity",
#            position = "dodge") + 
#   geom_point(position = position_dodge(0.9),
#              aes(x = total_cells/n_cells,
#                  y = reorder(Region_summarized, -total_cells),
#                  shape = indiv_genotype,
#                  group = orig.ident), df,
#              size = 2, alpha = 0.75) +
#   theme_classic() +
#   scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
#   labs(x = "percent of neurons",
#        y = "predicted region") +
#   theme_classic() 
# 
# ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/pct_cells_byRegion.svg"),
#        width = 8, height = 8)



#### Supplemental Figure 7 ####
# volcano plots 
setwd(paste0(root.dir, "/DGE_CellTypes"))
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")

sex = c("Male", "Female")
status = c("Dom", "Sub")

femaleDEGs <- limma_list$Fdom_vs_sub_all
maleDEGs <- limma_list$Mdom_vs_sub_all

up <- "Dom"
down <- "Sub"

# males - all neurons - supp fig 7A
dcm <- maleDEGs %>% 
  mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
  mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
  mutate(rank = logFC * log10) %>% 
  mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))

max.x <- 0.5 + round(max(dcm$logFC) * 2) / 2
min.x <- -0.5 + round(min(dcm$logFC) * 2) / 2 
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
                                             dcm_jitters$gene) & 
                             P.Value > 0.05),
             aes(x = logFC, y = log10, color = diffexpressed), 
             alpha = 0.25, size = 3.5) +
  geom_point(data = subset(dcm, !gene %in% c(dcm_jitters_up$gene,
                                             dcm_jitters_down$gene,
                                             dcm_jitters$gene) & 
                             P.Value < 0.05),
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
  ylim(0, 4.5) +
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
       width = 25, height = 21, units = "cm", dpi = 600)

# females - all neurons - supp fig 7B
dcf <- femaleDEGs %>% 
  mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
  mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
  mutate(rank = logFC * log10) %>% 
  mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))

max.x <- 0.5 + round(max(dcf$logFC) * 2) / 2
min.x <- -1.5 + round(min(dcf$logFC) * 2) / 2 
breaks = seq(min.x, max.x, 0.5)

upgene_labels <- dcf %>% slice_max(order_by = rank, n = 10)
downgene_labels <- dcf %>% slice_max(order_by = -rank, n = 10)

dcf_jitters <- subset(dcf, dcf$log10 > 3 & 
                        !gene %in% c(upgene_labels$gene, downgene_labels$gene))

dcf_jitters_up <- subset(dcf, dcf$gene %in% upgene_labels)
dcf_jitters_down <- subset(dcf, dcf$gene %in% downgene_labels)
pos_jitter <- position_jitter(width = 0.4, height = 0.5, seed = 123)

ggplot() +
  geom_point(data = subset(dcf, !gene %in% c(dcf_jitters_up$gene,
                                             dcf_jitters_down$gene,
                                             dcf_jitters$gene) & 
                             P.Value > 0.05),
             aes(x = logFC, y = log10, color = diffexpressed), 
             alpha = 0.25, size = 3.5) +
  geom_point(data = subset(dcf, !gene %in% c(dcf_jitters_up$gene,
                                             dcf_jitters_down$gene,
                                             dcf_jitters$gene) & 
                             P.Value < 0.05),
             aes(x = logFC, y = log10, color = diffexpressed), 
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
  scale_x_continuous(limits = c(min.x, max.x), breaks = breaks) +
  ylim(0, 4.5) +
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
       width = 25, height = 21, units = "cm", dpi = 600)



# number of DEGs in either direction for each cluster/subtype
DEG_summary <- limma_list %>%
  imap_dfr(~ {.x %>% 
      mutate(contrast_full = .y,
             contrast = str_extract(.y, "(?<=_)[^_]+$"),
             sex = str_to_upper(str_sub(.y, 1, 1)),
             direction = case_when(P.Value < 0.05 & logFC >  0.2  ~ "up",
                                   P.Value < 0.05 & logFC < -0.2  ~ "down",
                                   TRUE ~ NA_character_)) %>%
      filter(!is.na(direction)) %>%
      select(gene, sex, contrast, direction)}) %>%
  group_by(contrast) %>%
  summarise(up_Mdom = sum(sex == "M" & direction == "up"),
            up_Fdom   = sum(sex == "F" & direction == "up"),
            down_Mdom = sum(sex == "M" & direction == "down"),
            down_Fdom = sum(sex == "F" & direction == "down"),
            up_overlap = length(intersect(gene[sex == "M" & direction == "up"],
                                          gene[sex == "F" & direction == "up"])),
            down_overlap = length(intersect(gene[sex == "M" & direction == "down"],
                                            gene[sex == "F" & direction == "down"]))) %>% 
  pivot_longer(cols = -contrast, 
               names_to = "stat", 
               values_to = "n") %>%
  pivot_wider(names_from = contrast, 
              values_from = n) %>% 
  arrange(desc(stat))


write_csv(DEG_summary, file = paste0(root.dir, "/manuscriptFigures/supplemental/DEG_overlap_summary.csv"))

#### Supplemental Figure 8 ####
# GO analysis on RRHO results
setwd(paste0(root.dir, "/DGE_CellTypes/neurons/all_neurons/"))
load("./RRHO/rrho_results.rda")


get.GOtop15.RRHO(rrho_results,
                 quadrants_to_check = c("upup", "downdown"),
                 filetype = ".pdf",
                 analysis = "RedRibbon",
                 img.params = list(width = 14, height = 7, paper = "special"),
                 outdir = paste0(root.dir, "/manuscriptFigures/supplemental"))


#### Supplemental Figure 9 ####
setwd(paste0(root.dir, "/DGE_CellTypes"))
load("./neurons/neuroendocrine_genes/limma_trend/limma_results.rda")

# how to determine max log scale before graphing? 
max.log.scale = 115;
# make sure these match upper/lower case
sex = c("Female", "Male");
status = c("Dom", "Sub");

rrho_results <- get.RRHO(limma_results,
                         group.by = sex,
                         compare.across = status,
                         new.max.log = max.log.scale,
                         filetype = ".svg",
                         analysis = "RedRibbon",
                         img.params = c(width = 10, height = 10),
                         outdir = paste0(root.dir, "/manuscriptFigures/supplemental/neuroendocrine_genes"))

# supp fig 9F
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")
DefaultAssay(MSCneurons.reclust) = "RNA"
MSCneurons.reclust = NormalizeData(MSCneurons.reclust)

# subset genes
subset_genes <- c("Esr1",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2")

DotPlot(MSCneurons.reclust,
        assay = 'RNA',
        features = subset_genes,
        group.by = 'orig.ident',
        min_count = 2,
        col.min = 0,
        scale = FALSE) +
  scale_y_discrete(labels = sub("\\.data", "", 
                                unique(MSCneurons.reclust$orig.ident))) +
  theme_classic() +
  ggtitle('Expression in Neurons') +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(.5,.5,2,.5), "cm"))

ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/neur.endoGenes_neurons_unscaledExpr.svg"),
       height = 4, width = 5.25)

#### Supplemental Figure 10 ####
# load data
load(paste0(root.dir, "/Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda"))
DefaultAssay(MSCneurons.reclust) = 'RNA'
MSCneurons.reclust = NormalizeData(MSCneurons.reclust)
neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))
load(paste0(net.dir, "/seurat/mouse.snseq.combined.sct.RData"))

neuropeptides %>% 
  mutate(Gene.name = str_to_title(Gene.name)) %>% 
  filter(Gene.name %in% rownames(MSCneurons.reclust)) %>% 
  pull(Gene.name) -> np.genes


sub.MSCneurons = subset_by_gene(MSCneurons.reclust,
                                np.genes, 
                                "counts",
                                min_count = 2)


# get presence/absence (1/0) for neuropeptide.genes
umap = data.frame(MSCneurons.reclust@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(MSCneurons.reclust@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

select_neur %>% 
  rowwise() %>% 
  filter(sum(c_across(intersect(np.genes, colnames(.)))) > 0) -> filtered_neur

select_neur %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  pivot_longer(cols = starts_with("pct_"),
               names_to = "gene") %>% 
  filter(value > 1) %>% 
  mutate(gene = sub("pct_", "", gene)) %>% pull(gene) -> topNP 



sapply(unique(MSCneurons.reclust$orig.ident), 
       \(ident, df, genes, seurat_obj) {
         
         cells = df %>% 
           filter(orig.ident == ident) %>% 
           pull("Cell_ID")
         
         cov.mat = t(as.matrix(seurat_obj@assays$RNA@counts[genes, cells] %>% 
                                 data.frame() %>% 
                                 mutate(across(where(is.numeric), 
                                               function(x) ifelse(x > 1, 1, 0)))))
         
       }, 
       df = filtered_neur, genes = topNP, seurat_obj = MSCneurons.reclust) -> cov.matList



### run hd test
males.hd = testCov(cov.matList[["male.dom.data"]],
                   cov.matList[["male.sub.data"]], 
                   method = "HD",
                   J = 1000, 
                   alpha = 0.05, 
                   n.core = 12)

females.hd = testCov(cov.matList[["female.dom.data"]],
                     cov.matList[["female.sub.data"]], 
                     method = "HD",
                     J = 1000, 
                     alpha = 0.05, 
                     n.core = 12)

doms.hd = testCov(cov.matList[["male.dom.data"]],
                  cov.matList[["female.dom.data"]], 
                  method = "HD",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)

subs.hd = testCov(cov.matList[["male.sub.data"]],
                  cov.matList[["female.sub.data"]], 
                  method = "HD",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)

## run clx test
males.CLX = testCov(cov.matList[["male.dom.data"]],
                   cov.matList[["male.sub.data"]], 
                   method = "CLX",
                   J = 1000, 
                   alpha = 0.05, 
                   n.core = 12)

females.CLX = testCov(cov.matList[["female.dom.data"]],
                     cov.matList[["female.sub.data"]], 
                     method = "CLX",
                     J = 1000, 
                     alpha = 0.05, 
                     n.core = 12)

doms.CLX = testCov(cov.matList[["male.dom.data"]],
                  cov.matList[["female.dom.data"]], 
                  method = "CLX",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)

subs.CLX = testCov(cov.matList[["male.sub.data"]],
                  cov.matList[["female.sub.data"]], 
                  method = "CLX",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)


## run LC test
males.LC = covtest.lc(cov.matList[["male.dom.data"]],
                      cov.matList[["male.sub.data"]])

females.LC = covtest.lc(cov.matList[["female.dom.data"]],
                        cov.matList[["female.sub.data"]])

doms.LC = covtest.lc(cov.matList[["male.dom.data"]],
                     cov.matList[["female.dom.data"]])

subs.LC = covtest.lc(cov.matList[["male.sub.data"]],
                     cov.matList[["female.sub.data"]])


#### create dataframe
results.testcov = data.frame(comparison = c('Males','Females','Doms','Subs'),
                             HD.p.val = c(males.hd[[2]],
                                          females.hd[[2]],
                                          doms.hd[[2]],
                                          subs.hd[[2]]),
                             HD.stat = c(males.hd[[1]],
                                         females.hd[[1]],
                                         doms.hd[[1]],
                                         subs.hd[[1]]),
                             CLX.p.val = c(males.CLX[[2]],
                                           females.CLX[[2]],
                                           doms.CLX[[2]],
                                           subs.CLX[[2]]),
                             CLX.stat = c(males.CLX[[1]],
                                          females.CLX[[1]],
                                          doms.CLX[[1]],
                                          subs.CLX[[1]]),
                             LC.p.val = c(males.LC$pval,
                                          females.LC$pval,
                                          doms.LC$pval,
                                          subs.LC$pval),
                             LC.stat = c(males.LC$stat,
                                         females.LC$stat,
                                         doms.LC$stat,
                                         subs.LC$stat))
results.testcov
  write_csv(results.testcov, paste0(root.dir, "/manuscriptFigures/supplemental/network.stats.topNPexpr.csv"))


## all expressed NP genes
select_neur %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  pivot_longer(cols = starts_with("pct_"),
               names_to = "gene") %>% 
  mutate(gene = sub("pct_", "", gene)) %>% pull(gene) -> allNPexpr 


sapply(names(sub.MSCneurons), \(NP) {
  
  if (sum(dim(table(sub.MSCneurons[[NP]]@meta.data$orig.ident, 
                    sub.MSCneurons[[NP]]@meta.data$indiv_genotype))) < 6) {
    return(NA)
  } else {
    return(NP)
  }
  
}) -> allNPs

allNPs = subset(allNPs, !is.na(allNPs))

sapply(unique(MSCneurons.reclust$orig.ident), 
       \(ident, df, genes, seurat_obj) {
         
         cells = df %>% 
           filter(orig.ident == ident) %>% 
           pull("Cell_ID")
         
         cov.mat = t(as.matrix(seurat_obj@assays$RNA@counts[genes, cells] %>% 
                                 data.frame() %>% 
                                 mutate(across(where(is.numeric), 
                                               function(x) ifelse(x > 1, 1, 0)))))
         
       }, 
       df = filtered_neur, genes = allNPs, seurat_obj = MSCneurons.reclust) -> cov.matList.allNP

### run hd test
males.hd = testCov(cov.matList.allNP[["male.dom.data"]],
                   cov.matList.allNP[["male.sub.data"]], 
                   method = "HD",
                   J = 1000, 
                   alpha = 0.05, 
                   n.core = 12)

females.hd = testCov(cov.matList.allNP[["female.dom.data"]],
                     cov.matList.allNP[["female.sub.data"]], 
                     method = "HD",
                     J = 1000, 
                     alpha = 0.05, 
                     n.core = 12)

doms.hd = testCov(cov.matList.allNP[["male.dom.data"]],
                  cov.matList.allNP[["female.dom.data"]], 
                  method = "HD",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)

subs.hd = testCov(cov.matList.allNP[["male.sub.data"]],
                  cov.matList.allNP[["female.sub.data"]], 
                  method = "HD",
                  J = 1000, 
                  alpha = 0.05, 
                  n.core = 12)

## run clx test
males.CLX = testCov(cov.matList.allNP[["male.dom.data"]],
                    cov.matList.allNP[["male.sub.data"]], 
                    method = "CLX",
                    J = 1000, 
                    alpha = 0.05, 
                    n.core = 12)

females.CLX = testCov(cov.matList.allNP[["female.dom.data"]],
                      cov.matList.allNP[["female.sub.data"]], 
                      method = "CLX",
                      J = 1000, 
                      alpha = 0.05, 
                      n.core = 12)

doms.CLX = testCov(cov.matList.allNP[["male.dom.data"]],
                   cov.matList.allNP[["female.dom.data"]], 
                   method = "CLX",
                   J = 1000, 
                   alpha = 0.05, 
                   n.core = 12)

subs.CLX = testCov(cov.matList.allNP[["male.sub.data"]],
                   cov.matList.allNP[["female.sub.data"]], 
                   method = "CLX",
                   J = 1000, 
                   alpha = 0.05, 
                   n.core = 12)


## run LC test
males.LC = covtest.lc(cov.matList.allNP[["male.dom.data"]],
                      cov.matList.allNP[["male.sub.data"]])

females.LC = covtest.lc(cov.matList.allNP[["female.dom.data"]],
                        cov.matList.allNP[["female.sub.data"]])

doms.LC = covtest.lc(cov.matList.allNP[["male.dom.data"]],
                     cov.matList.allNP[["female.dom.data"]])

subs.LC = covtest.lc(cov.matList.allNP[["male.sub.data"]],
                     cov.matList.allNP[["female.sub.data"]])


#### create dataframe
results.testcov = data.frame(comparison = c('Males','Females','Doms','Subs'),
                             HD.p.val = c(males.hd[[2]],
                                          females.hd[[2]],
                                          doms.hd[[2]],
                                          subs.hd[[2]]),
                             HD.stat = c(males.hd[[1]],
                                         females.hd[[1]],
                                         doms.hd[[1]],
                                         subs.hd[[1]]),
                             CLX.p.val = c(males.CLX[[2]],
                                           females.CLX[[2]],
                                           doms.CLX[[2]],
                                           subs.CLX[[2]]),
                             CLX.stat = c(males.CLX[[1]],
                                          females.CLX[[1]],
                                          doms.CLX[[1]],
                                          subs.CLX[[1]]),
                             LC.p.val = c(males.LC$pval,
                                          females.LC$pval,
                                          doms.LC$pval,
                                          subs.LC$pval),
                             LC.stat = c(males.LC$stat,
                                         females.LC$stat,
                                         doms.LC$stat,
                                         subs.LC$stat))
results.testcov
  write_csv(results.testcov, paste0(root.dir, "/manuscriptFigures/supplemental/network.stats.allNPexpr.csv"))

# from Neuropeptide_network_snseq_mouse_IMC.R line 3640-3905 (end script), 'testing covariance matrices'

#### Supplemental Figure 11 #### 
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

# percent of neurons expressing each combination:
# avp only, oxt only, scg2 only; avp + oxt, avp + scg2, oxt + scg2; oxt + avp + scg2
top3np <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp", "Scg2")) %>% 
  t() %>% as.data.frame() %>% 
  mutate(Oxt = ifelse(Oxt > 1, 1, 0),
         Avp = ifelse(Avp > 1, 1, 0),
         Scg2 = ifelse(Scg2 > 1, 1, 0),
         Oxt_Avp = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 0, 1, 0),
         Oxt_Scg2 = ifelse(Oxt == 1 & Avp == 0 & Scg2 == 1, 1, 0),
         Avp_Scg2 = ifelse(Oxt == 0 & Avp == 1 & Scg2 == 1, 1, 0),
         Oxt_Avp_Scg2 = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 1, 1, 0)) %>% 
  rownames_to_column(var = "Cell_ID") 

meta = as.data.frame(MSCneurons.reclust@meta.data) %>% 
  filter(Cell_ID %in% top3np$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, Cell_ID) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))


top3np.excl <- left_join(top3np, meta,
                         by = "Cell_ID") %>% 
  mutate(Oxt = ifelse(Oxt == 1 &
                        rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Oxt),
         Avp = ifelse(Avp == 1 &
                        rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Avp),
         Scg2 = ifelse(Scg2 == 1 &
                         rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Scg2))

## supp fig 11A
# as stacked barplot
labels = c("sum_Scg2"="Scg2 only",
           "sum_Oxt"="Oxt only",
           "sum_Avp"="Avp only",
           "sum_Oxt_Avp"="Oxt + Avp",
           "sum_Oxt_Scg2"="Oxt + Scg2",
           "sum_Avp_Scg2"="Avp + Scg2",
           "sum_Oxt_Avp_Scg2"="Oxt + Avp + Scg2")

# plot pcts for each NP
sum_np_group <- top3np.excl %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n(),
            across(all_of(c(2:8)),
                   list(sum = ~ sum(.x)),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

temp_labels = sum_np_group[,c("orig.ident", "cell_count")]

sum_np_group %>% 
  rowwise() %>% 
  mutate(cell_count = sum(c_across(starts_with("sum")), 
                          na.rm = TRUE)) %>% 
  dplyr::select(-c('cell_count'),
                -c(starts_with('pct'))) %>% 
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill',
           stat = 'identity') +
  geom_text(data = temp_labels,
            aes(x = orig.ident,
                y = Inf,
                label = paste0('n = ', format(cell_count, big.mark = ',', trim = TRUE)),
                vjust = -1),
            color = 'black',
            size = 2.8) +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  scale_fill_viridis_d(option = "plasma", labels = labels, end = 0.9) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Co-expression") +
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
                             unit = 'pt')) -> p1


# as pct of each NP type - by indiv
pct_oxt_indiv <- top3np.excl %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(across(all_of(contains("Oxt")),
                   list(sum = ~ sum(.x)),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  mutate(cell_count = rowSums(select(., contains("Oxt")))) %>% 
  ungroup() %>% 
  mutate(across(all_of(contains("Oxt")),
                list(pct = ~ (.x / cell_count) * 100),
                .names = "{fn}_{.col}")) %>% 
  select(c("orig.ident", "indiv_genotype", starts_with("pct"))) %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "coexpr",
               values_to = "pct_ofOxt_neurons")

pct_avp_indiv <- top3np.excl %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(across(all_of(contains("Avp")),
                   list(sum = ~ sum(.x)),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  mutate(cell_count = rowSums(select(., contains("Avp")))) %>% 
  ungroup() %>% 
  mutate(across(all_of(contains("Avp")),
                list(pct = ~ (.x / cell_count) * 100),
                .names = "{fn}_{.col}")) %>% 
  select(c("orig.ident", "indiv_genotype", starts_with("pct"))) %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "coexpr",
               values_to = "pct_ofAvp_neurons")


pct_scg2_indiv <- top3np.excl %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(across(all_of(contains("Scg2")),
                   list(sum = ~ sum(.x)),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  mutate(cell_count = rowSums(select(., contains("Scg2")))) %>% 
  ungroup() %>% 
  mutate(across(all_of(contains("Scg2")),
                list(pct = ~ (.x / cell_count) * 100),
                .names = "{fn}_{.col}")) %>% 
  select(c("orig.ident", "indiv_genotype", starts_with("pct"))) %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "coexpr",
               values_to = "pct_ofScg2_neurons")

np_list = list(pct_avp_indiv, pct_oxt_indiv, pct_scg2_indiv)

pct_np_indiv <- np_list %>% 
  reduce(full_join, by = c("orig.ident", 
                           "indiv_genotype", 
                           "coexpr")) %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "n_type",
               values_to = "percentage",
               values_drop_na = TRUE) %>% 
  mutate(coexpr = sub("pct_sum_", "", coexpr),
         n_type = sub("pct_of", "", n_type))

labels = c("Scg2"="Scg2 only",
           "Oxt"="Oxt only",
           "Avp"="Avp only",
           "Oxt_Avp"="Oxt + Avp",
           "Oxt_Scg2"="Oxt + Scg2",
           "Avp_Scg2"="Avp + Scg2",
           "Oxt_Avp_Scg2"="Oxt + Avp + Scg2")

pct_np_indiv %>% 
  ggplot(aes(x = factor(coexpr, levels = c("Avp", "Oxt", "Scg2", 
                                           "Oxt_Avp", "Avp_Scg2", "Oxt_Scg2", 
                                           "Oxt_Avp_Scg2")),
             y = percentage/100,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0),
                     limits = c(0, 1)) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 16)) +
  scale_x_discrete(labels = labels) +
  facet_wrap(~ n_type, scales = "free_x") +
  ggtitle('Coexpression in neurons expressing top 3 neuropeptides') -> p2 

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/top3_NPs_pctNeur_indiv.svg'),
       p1 + p2 + 
         plot_layout(ncol = 2, axes = "collect", guides = "collect",
                     widths = c(1, 3)),
       height = 7, width = 15)

# supp figure 11B
top3np <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp", "Scg2")) %>% 
  t() %>% as.data.frame() %>% 
  mutate(Oxt = ifelse(Oxt > 1, 1, 0),
         Avp = ifelse(Avp > 1, 1, 0),
         Scg2 = ifelse(Scg2 > 1, 1, 0)) %>% 
  rownames_to_column(var = "Cell_ID") 

meta = as.data.frame(MSCneurons.reclust@meta.data) %>% 
  filter(Cell_ID %in% top3np$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, Cell_ID, seurat_clusters) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))

top3np = left_join(meta, top3np, 
                   by = "Cell_ID")


# counts & percentages
sum_np_indiv <- top3np %>%
  group_by(orig.ident,
           indiv_genotype,
           seurat_clusters) %>%
  filter(seurat_clusters %in% c(0:9)) %>% 
  mutate(seurat_clusters = paste("cluster", seurat_clusters, sep = " ")) %>% 
  summarise(cell_count = n(),
            across(all_of(c("Avp", "Oxt", "Scg2")),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

# percent of all neurons expr select genes by individual
sum_np_indiv %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "gene",
    names_prefix = "pct_",
    values_to = "percentage"
  ) -> sum.stats

sum.stats %>% 
  group_by(orig.ident, gene, seurat_clusters) %>% 
  summarise(percentage = mean(percentage)) -> avgs

avgs %>% 
  ggplot(aes(x = factor(gene, levels = c("Oxt", "Avp", "Scg2")),
             y = percentage/100)) +
  geom_bar(aes(fill = orig.ident), avgs,
           stat = 'identity',
           position = "dodge") +
  geom_point(position = position_dodge(0.9),
             aes(shape = indiv_genotype,
                 group = orig.ident), sum.stats,
             size = 2, alpha = 0.75) +
  theme_classic() +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  xlab('') +
  scale_y_continuous(name = 'Percent of Nuclei',
                     labels = scales::percent_format(), expand = c(0.01,0),
                     limits = c(0, 1)) +
  facet_wrap( ~ seurat_clusters, nrow = 2) +
  ggtitle('Percent of neurons expressing top 3 neuropeptides')

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/top3_NPs_pctNeur_byCluster.svg'),
       height = 3.5, width = 8)


# supp fig 11 C
## if starting here:
# load("./neurons/all_neurons/limma_perm/norm_counts.rda")
# metadata <- read_csv("./neurons/all_neurons/limma_perm/metadata.csv")
# 
# norm_counts <- data.frame(v.dl$E) %>% 
#   rownames_to_column(var = "gene") %>% 
#   filter(gene %in% c("Oxt", "Avp", "Scg2")) %>% 
#   pivot_longer(cols = !gene, 
#                values_to = "norm_expr",
#                names_to = "pb_colname") %>% 
#   left_join(metadata, by = "pb_colname") %>% 
#   mutate(group = paste(Sex, Status, sep = " "),
#          id = paste(Sex, Status, indiv_genotype, sep = "_"),
#          seurat_clusters = paste("cluster", seurat_clusters, sep = " ")) %>% 
#   group_by(gene, group, id, seurat_clusters) %>% 
#   summarise(norm_expr = mean(norm_expr))

# Seurat-normalized instead of limma
DefaultAssay(MSCneurons.reclust) = "RNA"
MSCneurons.reclust = NormalizeData(MSCneurons.reclust)

np.sub <- subset_by_gene(MSCneurons.reclust,
                         subset_genes = c("Oxt", "Avp", "Scg2"),
                         slot = "counts",
                         min_count = 2)

norm_countsOXT <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Oxt") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Oxt_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(seurat_clusters %in% c(0:9) & Cell_ID %in% np.sub$Oxt@meta.data$Cell_ID)

oxt.exp <- norm_countsOXT %>% 
  group_by(orig.ident, seurat_clusters, indiv_genotype) %>% 
  summarise(avg.Oxt_expr = mean(Oxt_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident),
         seurat_clusters = paste("cluster", seurat_clusters, sep = " "))
  

oxt.exp %>% 
  ggplot(aes(orig.ident, avg.Oxt_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .4, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(limits = c(-0.2, 4.2)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 28), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("OXT") -> p1


norm_countsAVP <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Avp") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Avp_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(seurat_clusters %in% c(0:9) & Cell_ID %in% np.sub$Avp@meta.data$Cell_ID)


avp.exp <- norm_countsAVP %>% 
  group_by(orig.ident, seurat_clusters, indiv_genotype) %>% 
  summarise(avg.Avp_expr = mean(Avp_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident),
         seurat_clusters = paste("cluster", seurat_clusters, sep = " "))


avp.exp %>% 
  ggplot(aes(orig.ident, avg.Avp_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .4, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(limits = c(-0.2, 4.2)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 28), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("AVP") -> p2


norm_countsSCG2 <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Scg2") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Scg2_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(seurat_clusters %in% c(0:9) & Cell_ID %in% np.sub$Scg2@meta.data$Cell_ID)

scg2.exp <- norm_countsSCG2 %>% 
  group_by(orig.ident, seurat_clusters, indiv_genotype) %>% 
  summarise(avg.Scg2_expr = mean(Scg2_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident),
         seurat_clusters = paste("cluster", seurat_clusters, sep = " "))


scg2.exp %>% 
  ggplot(aes(orig.ident, avg.Scg2_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .4, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(limits = c(-0.2, 4.2)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 28), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("SCG2") -> p3


ggsave(paste0(root.dir, "/manuscriptFigures/supplemental/top3NPs_byCluster.svg"),
       p1 + p2 + p3 +
         plot_layout(ncol = 3, axes = "collect", guides = "collect"),
       width = 23, height = 5.5)


#### Supplemental Figure 12 #### 
setwd(root.dir)
load("./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda")

# Seurat-normalized instead of limma
DefaultAssay(MSCneurons.reclust) = "RNA"
MSCneurons.reclust = NormalizeData(MSCneurons.reclust)

# percent of neurons in each subclass
predicted.prob.df = MSCneurons.reclust@meta.data 

predicted.prob.df %>% 
  group_by(predicted.prob) %>% 
  summarise(cell_count = n()) %>%  
  mutate(total_cells = sum(cell_count),
         pct_cells = cell_count/total_cells) %>% 
  filter(pct_cells > 0.01) %>% 
  pull(predicted.prob) -> top_neuron_subclasses


np.sub <- subset_by_gene(MSCneurons.reclust,
                         subset_genes = c("Oxt", "Avp", "Scg2"),
                         slot = "counts",
                         min_count = 2)

# Oxytocin
norm_countsOXT <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Oxt") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Oxt_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(Cell_ID %in% np.sub$Oxt@meta.data$Cell_ID)

oxt.exp <- norm_countsOXT %>% 
  group_by(orig.ident, predicted.prob, indiv_genotype) %>% 
  summarise(avg.Oxt_expr = mean(Oxt_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))

# oxt.exp.filtered <- oxt.exp %>% 
#   group_by(predicted.prob) %>% 
#   count(predicted.prob) %>% 
#   filter(n == 12) %>% 
#   ungroup()

### pct of neurons 
top3np <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp", "Scg2")) %>% 
  t() %>% as.data.frame() %>% 
  mutate(Oxt = ifelse(Oxt > 1, 1, 0),
         Avp = ifelse(Avp > 1, 1, 0),
         Scg2 = ifelse(Scg2 > 1, 1, 0)) %>% 
  rownames_to_column(var = "Cell_ID") 

meta = as.data.frame(MSCneurons.reclust@meta.data) %>% 
  filter(Cell_ID %in% top3np$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, Cell_ID, predicted.prob) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))

top3np = left_join(meta, top3np, 
                   by = "Cell_ID")


# counts & percentages
sum_np_indiv <- top3np %>%
  group_by(orig.ident,
           indiv_genotype,
           predicted.prob) %>%
  filter(predicted.prob %in% top_neuron_subclasses) %>% 
  summarise(cell_count = n(),
            across(all_of(c("Avp", "Oxt", "Scg2")),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")


# percent of all neurons expr select genes by individual
sum_np_indiv %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "gene",
    names_prefix = "pct_",
    values_to = "percentage"
  ) -> sum.stats

sum.stats %>% 
  group_by(orig.ident, gene, predicted.prob) %>% 
  summarise(percentage = mean(percentage),
            total_cells = sum(cell_count)) %>% 
  mutate(top = percentage ) -> avgs

facet_order <- avgs %>%
  group_by(predicted.prob) %>%
  summarise(total_cells = mean(total_cells)) %>% 
  arrange(desc(total_cells)) %>% 
  pull(predicted.prob)


avgs$predicted.prob <- factor(avgs$predicted.prob, levels = facet_order)
sum.stats$predicted.prob <- factor(sum.stats$predicted.prob, levels = facet_order)

facet_labels <- avgs %>%
  distinct(predicted.prob) %>%
  mutate(label_html = str_replace(predicted.prob,
                                  "^(.*?):(.*)$", 
                                  "<i>\\1</i>:<b>\\2</b>")) %>%
  { setNames(.$label_html, .$predicted.prob) }


avgs %>%
  filter(gene == "Oxt") %>% 
  mutate(ymin_gray = 0,
         ymax_gray = pmin(percentage, 50),
         ymin_color = 50,
         ymax_color = ifelse(percentage > 50, percentage, 50)) -> oxt.avgs


ggplot() +
  geom_rect(data = oxt.avgs,
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_gray,
                ymax = ymax_gray),
            fill = "gray80") +
  geom_rect(data = subset(oxt.avgs, percentage > 50),
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_color,
                ymax = ymax_color,
                fill = orig.ident)) +
  geom_jitter(data = subset(sum.stats, gene == "Oxt"),
              aes(x = orig.ident, y = percentage, shape = indiv_genotype),
              width = 0.15, size = 4, alpha = 0.75) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0.01, 0),
                     limits = c(-1, 101)) +
  facet_wrap(~ predicted.prob, nrow = 3,
             labeller = labeller(predicted.prob = facet_labels)) +
  theme_classic() +
  labs(x = NULL,
       y = "Percent of Nuclei",
       title = "Neurons in each sub-class expressing >1 OXT transcript") +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/OXT_pctNeur_bySubclass.svg'),
       height = 10, width = 20)


  
# norm counts
oxt.exp %>% 
  mutate(predicted.prob = factor(predicted.prob, levels = facet_order)) %>% 
  filter(predicted.prob %in% top_neuron_subclasses) %>% 
  ggplot(aes(orig.ident, avg.Oxt_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = 1.2, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_y_continuous(limits = c(0, 5.5)) +
  facet_wrap(~ predicted.prob, nrow = 3,
             labeller = labeller(predicted.prob = facet_labels)) +
  labs(x = NULL,
       y = "Normalized Expression",
       title = "Expression in neurons with >1 OXT transcript") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 20),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/OXTexpr_bySubclass.svg'),
       height = 10, width = 18)



## vasopressin
norm_countsAVP <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Avp") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Avp_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(Cell_ID %in% np.sub$Avp@meta.data$Cell_ID)

avp.exp <- norm_countsAVP %>% 
  group_by(orig.ident, predicted.prob, indiv_genotype) %>% 
  summarise(avg.Avp_expr = mean(Avp_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))

# avp.exp.filtered <- avp.exp %>% 
#   group_by(predicted.prob) %>% 
#   count(predicted.prob) %>% 
#   filter(n == 12) %>% 
#   ungroup()

avgs %>%
  filter(gene == "Avp") %>% 
  mutate(ymin_gray = 0,
         ymax_gray = pmin(percentage, 10),
         ymin_color = 10,
         ymax_color = ifelse(percentage > 10, percentage, 10)) -> avp.avgs


ggplot() +
  geom_rect(data = avp.avgs,
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_gray,
                ymax = ymax_gray),
            fill = "gray80") +
  geom_rect(data = subset(avp.avgs, percentage > 10),
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_color,
                ymax = ymax_color,
                fill = orig.ident)) +
  geom_jitter(data = subset(sum.stats, gene == "Avp"),
              aes(x = orig.ident, y = percentage, shape = indiv_genotype),
              width = 0.15, size = 4, alpha = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0.01, 0),
                     limits = c(-1, 101)) +
  facet_wrap(~ predicted.prob, nrow = 3,
             labeller = labeller(predicted.prob = facet_labels)) +
  theme_classic() +
  labs(x = NULL,
       y = "Percent of Nuclei",
       title = "Neurons in each sub-class expressing >1 AVP transcript") +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))


ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/AVP_pctNeur_bySubclass.svg'),
       height = 10, width = 20)


# norm counts
avp.exp %>% 
  mutate(predicted.prob = factor(predicted.prob, levels = facet_order)) %>% 
  filter(predicted.prob %in% top_neuron_subclasses) %>% 
  ggplot(aes(orig.ident, avg.Avp_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = 1.2, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_y_continuous(limits = c(0, 5.5)) +
  facet_wrap(~ predicted.prob, nrow = 3,
             labeller = labeller(predicted.prob = facet_labels)) +
  labs(x = NULL,
       y = "Normalized Expression",
       title = "Expression in neurons with >1 AVP transcript") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 20),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/AVPexpr_bySubclass.svg'),
       height = 10, width = 18)


## secretogranin2
norm_countsSCG2 <- data.frame(MSCneurons.reclust@assays$RNA@data) %>% 
  filter(rownames(.) == "Scg2") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Scg2_expr") %>% 
  mutate(Cell_ID = sub("\\.", "-", Cell_ID)) %>% 
  left_join(as.data.frame(MSCneurons.reclust@meta.data),
            by = "Cell_ID") %>% 
  filter(Cell_ID %in% np.sub$Scg2@meta.data$Cell_ID)

scg2.exp <- norm_countsSCG2 %>% 
  group_by(orig.ident, predicted.prob, indiv_genotype) %>% 
  summarise(avg.Scg2_expr = mean(Scg2_expr)) %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident))

# scg2.exp.filtered <- scg2.exp %>% 
#   group_by(predicted.prob) %>% 
#   count(predicted.prob) %>% 
#   filter(n == 12) %>% 
#   ungroup()


avgs %>%
  filter(gene == "Scg2") %>% 
  mutate(ymin_gray = 0,
         ymax_gray = pmin(percentage, 10),
         ymin_color = 10,
         ymax_color = ifelse(percentage > 10, percentage, 10)) -> scg2.avgs


ggplot() +
  geom_rect(data = scg2.avgs,
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_gray,
                ymax = ymax_gray),
            fill = "gray80") +
  geom_rect(data = subset(scg2.avgs, percentage > 10),
            aes(xmin = as.numeric(factor(orig.ident)) - 0.45,
                xmax = as.numeric(factor(orig.ident)) + 0.45,
                ymin = ymin_color,
                ymax = ymax_color,
                fill = orig.ident)) +
  geom_jitter(data = subset(sum.stats, gene == "Scg2"),
              aes(x = orig.ident, y = percentage, shape = indiv_genotype),
              width = 0.15, size = 4, alpha = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0.01, 0),
                     limits = c(-1, 101)) +
  facet_wrap(~ predicted.prob, nrow = 3, 
             labeller = labeller(predicted.prob = facet_labels)) +
  theme_classic() +
  labs(x = NULL,
       y = "Percent of Nuclei",
       title = "Neurons in each sub-class expressing >1 SCG2 transcript") +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))


ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/SCG2_pctNeur_bySubclass.svg'),
       height = 10, width = 20)


scg2.exp %>% 
  mutate(predicted.prob = factor(predicted.prob, levels = facet_order)) %>% 
  filter(predicted.prob %in% top_neuron_subclasses) %>% 
  ggplot(aes(orig.ident, avg.Scg2_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = 1.2, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_y_continuous(limits = c(0, 5.5)) +
  facet_wrap(~ predicted.prob, nrow = 3,
             labeller = labeller(predicted.prob = facet_labels)) +
  labs(x = NULL,
       y = "Normalized Expression",
       title = "Expression in neurons with >1 SCG2 transcript") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_markdown(),
        text = element_text(size = 20),
        axis.text.x = element_text(hjust = 0.75, 
                                   vjust = 0.75, 
                                   angle = 45))

ggsave(paste0(root.dir, '/manuscriptFigures/supplemental/SCG2expr_bySubclass.svg'),
       height = 10, width = 18)



# # percent of neurons expressing each combination:
# # avp only, oxt only, scg2 only; avp + oxt, avp + scg2, oxt + scg2; oxt + avp + scg2
# top3np <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
#   filter(rownames(.) %in% c("Oxt", "Avp", "Scg2")) %>% 
#   t() %>% as.data.frame() %>% 
#   mutate(Oxt = ifelse(Oxt > 1, 1, 0),
#          Avp = ifelse(Avp > 1, 1, 0),
#          Scg2 = ifelse(Scg2 > 1, 1, 0),
#          Oxt_Avp = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 0, 1, 0),
#          Oxt_Scg2 = ifelse(Oxt == 1 & Avp == 0 & Scg2 == 1, 1, 0),
#          Avp_Scg2 = ifelse(Oxt == 0 & Avp == 1 & Scg2 == 1, 1, 0),
#          Oxt_Avp_Scg2 = ifelse(Oxt == 1 & Avp == 1 & Scg2 == 1, 1, 0)) %>% 
#   rownames_to_column(var = "Cell_ID") 
# 
# meta = as.data.frame(MSCneurons.reclust@meta.data) %>% 
#   filter(Cell_ID %in% top3np$Cell_ID) %>% 
#   select(orig.ident, indiv_genotype, Cell_ID, predicted.prob) %>% 
#   mutate(orig.ident = sub("\\.data", "", orig.ident))
# 
# 
# top3np.excl <- left_join(top3np, meta,
#                          by = "Cell_ID") %>% 
#   mutate(Oxt = ifelse(Oxt == 1 &
#                         rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Oxt),
#          Avp = ifelse(Avp == 1 &
#                         rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Avp),
#          Scg2 = ifelse(Scg2 == 1 &
#                          rowSums(select(., Oxt, Avp, Scg2), na.rm = TRUE) > 1, 0, Scg2))
# 
# 
# # counts & percentages
# sum_np_indiv <- top3np.excl %>%
#   group_by(orig.ident,
#            indiv_genotype,
#            predicted.prob) %>%
#   summarise(cell_count = n(),
#             across(all_of(c(2:8)),
#                    list(sum = ~ sum(.x),
#                         pct = ~ (sum(.x) / n()) * 100),
#                    .names = "{fn}_{.col}"),
#             .groups = "drop")
# 
# # percent of all neurons expr select genes by individual
# sum_np_indiv %>%
#   pivot_longer(
#     cols = starts_with("pct_"),
#     names_to = "gene",
#     names_prefix = "pct_",
#     values_to = "percentage"
#   ) -> sum.stats
# 
# sum.stats %>% 
#   group_by(orig.ident, gene, predicted.prob) %>% 
#   summarise(percentage = mean(percentage)) -> avgs
# 
# avgs %>% 
#   ggplot(aes(x = factor(gene, levels = c("Oxt", "Avp", "Scg2", 
#                                          "Oxt_Avp", "Oxt_Scg2", "Avp_Scg2", 
#                                          "Oxt_Avp_Scg2")),
#              y = percentage)) +
#   geom_bar(aes(fill = orig.ident), avgs,
#            stat = 'identity',
#            position = "dodge") +
#   geom_point(position = position_dodge(0.9),
#              aes(shape = indiv_genotype,
#                  group = orig.ident), sum.stats,
#              size = 2, alpha = 0.75) +
#   theme_classic() +
#   scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
#   xlab('') + ylab('% of neurons') +
#   # scale_y_continuous(limits = c(0,75)) +
#   facet_wrap( ~ predicted.prob) +
#   ggtitle('Percent of neurons expressing top 3 neuropeptides')
# 
# # ggsave('manuscriptFigures/supplemental/top3_NPs_pctNeur.svg',
# #        height = 3.5, width = 9)




#### Decide whether to include these #### 
# dot plot with cluster assignments!!
sub.MSCneurons.reclust <- subset(MSCneurons.reclust, 
                                 seurat_clusters %in% c(0:9))

dittoSeq::dittoDotPlot(sub.MSCneurons.reclust,
                       assay = 'RNA',
                       # features = c("Oxt"),
                       vars = c("Oxt", "Avp", "Scg2"),
                       group.by = 'seurat_clusters',
                       split.by = 'orig.ident',
                       min = 1,
                       min.color = "lightgrey",
                       max.color = "#0000FF",
                       # cols = c("red", "blue", "green", "purple"),
                       # min_count = 2,
                       # col.min = 0,
                       scale = FALSE) +
  theme_classic() +
  ggtitle('Expression in Neurons') +
  theme(plot.margin = unit(c(.5,.5,2,.5), "cm")) 

# all clusters
dittoSeq::dittoDotPlot(MSCneurons.reclust,
                       assay = 'RNA',
                       vars = c("Oxt", "Avp", "Scg2"),
                       group.by = 'orig.ident',
                       # split.by = 'orig.ident',
                       min = 1,
                       max.color = "#0000FF",
                       min.percent = 0,
                       max.percent = 1,
                       summary.fxn.color = function(x) {mean(x[x>1])},
                       summary.fxn.size = function(x) {mean(x>1)},
                       scale = FALSE) +
  theme_classic() +
  ggtitle('Expression in Neurons') +
  theme(plot.margin = unit(c(.5,.5,2,.5), "cm")) 

# compare to regular dotplot Seurat function
DotPlot(MSCneurons.reclust,
        features = c("Oxt", "Avp", "Scg2"),
        group.by = 'orig.ident',
        min_count = 2,
        col.min = 0,
        scale.min = 0,
        scale.max = 100,
        scale = FALSE) +
  scale_y_discrete(labels = sub("\\.data", "",
                                unique(MSCneurons.reclust$orig.ident))) +
  theme_classic() +
  ggtitle('Expression in Neurons') +
  theme(plot.margin = unit(c(.5,.5,2,.5), "cm"))

ggsave('manuscriptFigures/supplemental/top3_NPs_dotplot_clusters.png',
       height = 4, width = 5.5)


