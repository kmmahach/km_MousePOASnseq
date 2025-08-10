# R-4.3.1, Seurat v.4.4.0
# Reclustering with neuronal nuclei only (exclude all glia and blood cells)

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
set.seed(12345)
compression = "xz" # slower, but usually smallest compression
setwd(root.dir) 

# functions
source("./Scripts/km_cleanScripts/functions/network_graph_fun.R")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "Seurat", "limma", "edgeR", "RedRibbon"), 
              out_prefix = "03", folder = "./Scripts/km_cleanScripts")

#### load data ####
load("./HypoMap/data/integrated_seurat_withHypoMap_predictions.rda")
  
# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# need to set to integrated for clustering
DefaultAssay(int.ldfs) = 'integrated'

# recluster
int.ldfs %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.4) -> MSCneurons.reclust

rm(int.ldfs) # save space in env

  save(MSCneurons.reclust, 
       file = "./Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda", 
       compress = compression)
  
#### graph re-clustered UMAP ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# idents to new clusters
Idents(object = MSCneurons.reclust) <- "seurat_clusters"

# UMAPs - all neuron clusters
MSCneurons.reclust %>%
  DimPlot(label = TRUE) +
  theme_classic()

ggsave('./neurons/UMAP/neurons.recluster.clusters.png',
       height = 10, width = 12)

# broad cell types
MSCneurons.reclust %>%
  DimPlot(group.by = 'parent_id.broad.prob')+
  theme_classic()

ggsave('./neurons/UMAP/neurons.recluster.parent_id.broad.png',
       height = 10, width = 12)

# sex/social status 
MSCneurons.reclust %>%
  DimPlot(group.by = 'orig.ident',
          size = 3) +
  theme_classic() +
  scale_color_manual(values = c("red",
                                'cyan',
                                'black',
                                'green'))

ggsave('./neurons/UMAP/neurons.recluster.orig.ident.png',
       height = 10, width = 12)

# social status
MSCneurons.reclust %>%
  DimPlot(group.by = 'orig.ident',
          size = 3)+
  theme_classic() +
  scale_color_manual(values = c("black",
                                'red',
                                'black',
                                'red'))

ggsave('./neurons/UMAP/neurons.recluster.status.png',
       height = 10, width = 10)

# indiv_genotype (predicted by souporcell) for each group
# dom male
subset(MSCneurons.reclust, orig.ident == "male.dom.data") %>%
  DimPlot(group.by = 'indiv_genotype')+
  theme_classic() +
  ggtitle('Male Dom - indiv_genotype')

ggsave('./neurons/UMAP/neurons.recluster.maleDom_genotypes.png',
       height = 10, width = 12)

# dom female
subset(MSCneurons.reclust, orig.ident == "female.dom.data") %>%
  DimPlot(group.by = 'indiv_genotype')+
  theme_classic() +
  ggtitle('Female Dom - indiv_genotype')

ggsave('./neurons/UMAP/neurons.recluster.femaleDom_genotypes.png',
       height = 10, width = 12)

# sub male
subset(MSCneurons.reclust, orig.ident == "male.sub.data") %>%
  DimPlot(group.by = 'indiv_genotype')+
  theme_classic() + 
  ggtitle('Male Sub - indiv_genotype')

ggsave('./neurons/UMAP/neurons.recluster.maleSub_genotypes.png',
       height = 10, width = 12)

# sub female
subset(MSCneurons.reclust, orig.ident == "female.sub.data") %>%
  DimPlot(group.by = 'indiv_genotype')+
  theme_classic() +
  ggtitle('Female Sub - indiv_genotype')

ggsave('./neurons/UMAP/neurons.recluster.femaleSub_genotypes.png',
       height = 10, width = 12)

# dom candidate genes?
DotPlot(MSCneurons.reclust,
        assay = 'SCT', features = c("Lhx8", "Crym", "Ache", "Myb", "Mog")) +
  theme_classic()

ggsave('./neurons/dom.genes.dotplot.png',
       height = 5, width = 5)

# scaled counts of neurons in each cluster (as percent)
# scale
MSCneurons.reclust@meta.data %>%
  dplyr::select(c(orig.ident, indiv_genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.neuron.count = Freq) -> neuron.genotype.scale

# cluster counts
MSCneurons.reclust@meta.data %>%
  dplyr::select(c(orig.ident, seurat_clusters, indiv_genotype)) %>% 
  table() %>%
  as.data.frame() -> neuron.cluster.genotype.count

# combine counts and scale
neuron.cluster.genotype.count %>% 
  full_join(neuron.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.neuron.count) -> neuron.cluster.genotype.count

# graph
neuron.cluster.genotype.count %>%
  ggplot(aes(x = seurat_clusters,
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 3)+
  facet_grid(~seurat_clusters,
             scales = "free_x") +
  theme_classic() +
  xlab("Neuron Cluster") +
  ylab("Percent of Neurons") +
  labs(color = "Group",
       shape = "Genotype") +
  theme(legend.position = c(0.9, 0.75),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 23),
        strip.text.x.top = element_blank(),
        axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 25),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 25)) +
  scale_color_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF"),
                     labels = c("Female Dom",
                                "Female Sub",
                                "Male Dom",
                                "Male Sub"))


ggsave('./neurons/clustersVorig.ident.neurons.recluster.genotype.scaled.png',
       height = 8, width = 12)

#### Prep for stats ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

#### cluster names
MSCneurons.reclust@meta.data %>%
  dplyr::select(parent_id.exp.prob,
         parent_id.broad.prob,
         seurat_clusters)  %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/total) -> neuron.cluster.names

## save neuron cluster names
write_csv(neuron.cluster.names, 'neurons/stats/neuron_clusterNames.csv')

## add variable for status and sex
neuron.cluster.genotype.count %>% 
  mutate(Sex = ifelse(grepl("female", orig.ident), "female", "male")) %>% 
  mutate(Status = ifelse(grepl("dom", orig.ident), "dom", "sub")) -> neuron.cluster.genotype.count

#### ANOVA ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

### run anova on every cluster
# create empty data frame
neuron.cluster.aov = data.frame(matrix(ncol = 6, nrow = 0))

# provide column names
colnames(neuron.cluster.aov) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "seurat_clusters")

## loop through each cluster
for (i in unique(neuron.cluster.genotype.count$seurat_clusters)) {
  ## run two way anova on sex and status
  tmp <- aov(Percent ~ Sex*Status,
             data = neuron.cluster.genotype.count %>%
                    filter(seurat_clusters == i))

  # get results
  tmp2 = summary(tmp)[[1]] %>%
    as.data.frame() %>%
    mutate(seurat_clusters = i)

  # combine in data frame
  neuron.cluster.aov = neuron.cluster.aov %>%
    rbind(tmp2)
}

# drop residuals
neuron.cluster.aov <- na.omit(neuron.cluster.aov)

## adjust pvalue (FDR)
neuron.cluster.aov %>%
  mutate(p_adj = p.adjust(`Pr(>F)`,
                          method = 'fdr',
                          n = nrow(neuron.cluster.aov))) %>% 
  rownames_to_column('Type') -> neuron.cluster.aov

# clean up
neuron.cluster.aov %>%
  mutate(Comp.type = case_when(grepl("Sex", Type) ~ "Sex",
                               grepl("Status", Type) ~ "Status")) %>%
  mutate(Comp.type = ifelse(grepl("Sex:Status", Type), 
                            "Sex:Status", Comp.type)) -> neuron.cluster.aov


## save anova table
write_csv(neuron.cluster.aov, 'neurons/stats/neuron_clusterANOVA.csv')

#### GLM binomial ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

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
    ggtitle(paste0("Neuron cluster ", i, ": binomial GLM cell count"))
  
  ggsave(paste0('neurons/stats/neuron_clusterGLM/graphs/cluster', i,'_binomial.GLM.cell.count.png'),
         width = 10, height = 6)
  
  message(paste0("saving cell count differences for cluster ", i))
  
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
write_csv(neuron.cluster.glm, 'neurons/stats/neuron_clusterGLM.csv')

## create p value dataframe
neuron.cluster.glm.p = neuron.cluster.glm %>% 
  dplyr::select(seurat_clusters,
                adj.pvalue.round,
                contrast) %>% 
  pivot_wider(names_from = contrast,
              values_from = adj.pvalue.round) 

#### Neuron cluster graphs ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

### Number of reads and genes per cluster

# https://romanhaa.github.io/projects/scrnaseq_workflow/#number-of-transcripts-and-expressed-genes
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

ggsave('neurons/nCount_nFeature_byCluster.png', 
       p1 + p2 + plot_layout(ncol = 1),
       height = 8, width = 14)

#### create counts table 
## sample by cluster
table_orig.idents_by_clusters <- MSCneurons.reclust@meta.data %>%
  group_by(orig.ident, integrated_snn_res.0.4) %>%
  summarize(count = n()) %>%
  spread(integrated_snn_res.0.4, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('orig.ident', 'total_cell_count', everything())) %>%
  arrange(factor(orig.ident, levels = levels(MSCneurons.reclust@meta.data$orig.ident)))

## cluster by sample
table_clusters_by_orig.idents <- MSCneurons.reclust@meta.data %>%
  dplyr::rename('cluster' = 'integrated_snn_res.0.4') %>%
  group_by(cluster, orig.ident) %>%
  summarize(count = n()) %>%
  spread(orig.ident, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(MSCneurons.reclust@meta.data$integrated_snn_res.0.4)))

### Percent of cells per cluster by sample
MSCneurons.reclust@meta.data %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- table_orig.idents_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  # mutate(orig.ident = factor(orig.ident, levels = levels(MSCneurons.reclust@meta.data$orig.ident))) %>%
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
  # scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
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

p2 <- table_clusters_by_orig.idents %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, 
                          levels = levels(MSCneurons.reclust@meta.data$integrated_snn_res.0.4))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), 
           position = 'fill', 
           stat = 'identity') +
  geom_text(data = temp_labels, 
            aes(x = cluster, y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1), 
            color = 'black', 
            size = 2.8) +
  # scale_fill_manual(name = 'orig.ident', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(),
                     expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 20, 
                             r = 0, 
                             b = 0, 
                             l = 10, 
                             unit = 'pt'))

ggsave('neurons/sample_composition_ClustersByPercent.png', 
       p1 + p2 + 
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
         unique() %>% length(), 
         MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
           unique() %>% length() )),
  width = 18, height = 8)


### Count of cells per cluster by orig.ident
MSCneurons.reclust@meta.data %>%
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- table_orig.idents_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  # mutate(orig.ident = factor(orig.ident, levels = levels(MSCneurons.reclust@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable),
           position = 'stack', 
           stat = 'identity') +
  geom_text(data = temp_labels,
            aes(x = orig.ident, 
                y = Inf, 
                label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), 
                vjust = -1),
            color = 'black', 
            size = 2.8) +
  # scale_fill_manual(name = 'Cluster', values = custom_colors$discrete) +
  scale_y_continuous(name = 'Number of cells', 
                     labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(legend.position = 'left',
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

p2 <- table_clusters_by_orig.idents %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  # mutate(cluster = factor(cluster, levels = levels(MSCneurons.reclust@meta.data$integrated_snn_res.0.4))) %>%
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

ggsave('neurons/composition_orig.idents_ClustersByNum.png',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      MSCneurons.reclust@meta.data$orig.ident %>% 
        unique() %>% length(),
      MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
        unique() %>% length() )),
  width = 18, height = 8)

#### Build neuron hierarchical tree ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

### create tree of cell types
BuildClusterTree(MSCneurons.reclust,
                 dims = 1:15,
                 reorder = FALSE,
                 reorder.numeric = FALSE,
                 slot = 'data',
                 assay = "SCT") -> MSCneurons.reclust

## create tree
neuron.tree <- MSCneurons.reclust@tools$BuildClusterTree
# tree$tip.label <- paste0("Cluster ", tree$tip.label)

# # convert tree to tibble
# neuron.tree = as_tibble(neuron.tree)
# 
# ## add cell type data
# neuron.tree = full_join(neuron.tree,
#                         MSCneurons.reclust@meta.data %>% 
#                    dplyr::select(integrated_snn_res.0.4,
#                                  parent_id.exp.prob) %>% 
#                    distinct(),
#                  by = c("label" = "integrated_snn_res.0.4"))
# 
# # convert back to tree
# neuron.tree = treeio::as.treedata(neuron.tree)

### graph tree of celltypes
ggtree::ggtree(neuron.tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(shape = 16, 
                        size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))

ggsave('neurons/neuron.cluster.tree.png',
       height = 6, width = 6)


#### Neuron marker specificity score ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

### Prep seurat object for running find markers 
## need to run before FindAllMarkers
# needs to recorrect again after filtering down to neurons? 
PrepSCTFindMarkers(MSCneurons.reclust,
                   assay = "SCT") -> MSCneurons.reclust

### calculate specificity score with DEGs for each cluster
FindAllMarkers(MSCneurons.reclust, 
               logfc.threshold = 0.1, # increase to increase speed
               test.use = 'wilcox',
               verbose = TRUE,
               min.pct = 0.10, 
               assay = 'SCT') -> neuron.markers.df 

# create gene column
neuron.markers.df %>% 
  rownames_to_column('gene') -> neuron.markers.df

neuron.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2)) -> neuron.markers.df

## save data
write_csv(neuron.markers.df, 'neurons/stats/neuron_clusterMarkers.csv')

### create heatmap of 6 top markers per cluster
## reduce to top 6 genes per cluster
neuron.markers.df %>% 
  group_by(cluster) %>% 
  slice_max(specificity, n = 6) %>% 
  ungroup() -> neuron.markers.df.reduce

neuron.markers.df %>% 
  filter(gene %in% unique(neuron.markers.df.reduce$gene)) %>% 
  dplyr::select(cluster, gene, specificity) %>% 
  left_join(neuron.markers.df.reduce %>% 
              dplyr::select(cluster, gene, specificity) %>% 
              mutate(keep = cluster)) %>% 
  arrange(keep) %>%  
  pivot_wider(names_from = 'cluster',
              values_from = 'specificity',
              id_cols = 'gene') %>% 
  column_to_rownames('gene') %>% 
  as.matrix() -> neuron.markers.df.reduce.heatmap

# replace NA with 0
neuron.markers.df.reduce.heatmap[is.na(neuron.markers.df.reduce.heatmap)] = 0

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
                     filename = 'neurons/neuron_markers.heatmap.png',
                     width = 10,
                     height = 10)

#### Neuron cluster limma-trend (DGE) ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

### create a loop for each neuron cluster 
# limmatrend
# set idents
Idents(object = MSCneurons.reclust) <- "seurat_clusters"

## run across all neuron clusters
cluster_nums = levels(unique(MSCneurons.reclust@meta.data$seurat_clusters))

# remove cluster 13 - not enough cells for contrasts 
for (j in cluster_nums[-length(cluster_nums)]) {

  # create new folder
  dir.create(paste0('./neurons/stats/limma_trend/cluster_', j))
  
  #subset to neuron clusters
  tmp = subset(MSCneurons.reclust, idents = j)
  
  # subset with SCT data
  DefaultAssay(tmp) = 'SCT'
  
  ### extract gene expression and orig.idents from neuron_cluster
  ## create expression dataframe of neuropeptides
  ### calculate variable genes
  
  # use integrated assay for variable features
FindVariableFeatures(tmp, 
                     assay = 'integrated',
                     selection.method = "vst", 
                     verbose = F) -> tmp.var

  # identify top variable genes
VariableFeatures(tmp.var, 
                 assay = 'integrated') -> neuron_cluster.reduce.group.topgenes.prep
  
  # create dummy
full_join(full_join(tmp@reductions$umap@cell.embeddings %>% 
                      as.data.frame() %>% 
                      rownames_to_column("Cell.id"), 
                    tmp@meta.data %>%
                      rownames_to_column("Cell.id")),
          tmp@assays$SCT@data %>% 
            as.data.frame() %>% t() %>%
            as.data.frame() %>% 
            rownames_to_column('Cell.id')) -> tmp.expression
  
  
  ## create vector of factor
tmp.expression %>% 
    mutate(neuron_cluster.orig.ident = orig.ident %>% 
             as.factor()) %>% 
    pull(neuron_cluster.orig.ident) %>% 
    droplevels() -> neuron_cluster.reduce.DomvsSub.vector.list.sct.prep
  
  ### counts matrix
  ## raw read count matrix
  ## rows = genes, columns = cells
  # only keep 500 variable genes
  # set negative values to 0

GetAssayData(tmp, assay = 'SCT') %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column('gene') %>% 
    dplyr::select(c(gene, tmp.expression %>% 
                      pull(Cell.id))) %>% 
    filter(gene %in% neuron_cluster.reduce.group.topgenes.prep) %>% 
    column_to_rownames('gene') %>% 
    as.matrix() %>% 
    pmax(0) -> neuron_cluster.reduce.vector.count.sct.prep
  
  #create list
neuron_cluster.reduce.vector.limma.sct.prep = 
  list(count = neuron_cluster.reduce.vector.count.sct.prep,
       condt = neuron_cluster.reduce.DomvsSub.vector.list.sct.prep)
  
  ### run function  
neuron_cluster.limma.results.sct = run_limmatrend_neuron_cluster(neuron_cluster.reduce.vector.limma.sct.prep)
  
  # save results to dataframe
full_join(neuron_cluster.limma.results.sct$ttDvsS %>% 
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
                rownames_to_column("Gene")) -> neuron_cluster.limma.results.sct.df
  
  # add color for significance  
neuron_cluster.limma.results.sct.df %>% 
    mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05, 'Sig', 'Not Sig'),
           Direction.type_DvsS = ifelse(logFC_DvsS > 0, 'up', 'down'),
           Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig', 
                                       Direction.type_DvsS, 
                                       Sig_DvsS)) %>% 
    mutate(Sig_DvsS_M = ifelse(adj.P.Val_DvsS_M <= 0.05, 'Sig', 'Not Sig'),
           Direction.type_DvsS_M = ifelse(logFC_DvsS_M > 0, 'up', 'down'),
           Sig.direction_DvsS_M = ifelse(Sig_DvsS_M == 'Sig',
                                         Direction.type_DvsS_M,
                                         Sig_DvsS_M)) %>% 
    mutate(Sig_DvsS_F = ifelse(adj.P.Val_DvsS_F <= 0.05, 'Sig','Not Sig'),
           Direction.type_DvsS_F = ifelse(logFC_DvsS_F > 0,'up','down'),
           Sig.direction_DvsS_F = ifelse(Sig_DvsS_F == 'Sig',
                                         Direction.type_DvsS_F,
                                         Sig_DvsS_F)) %>% 
    mutate(Sig_MvsF = ifelse(adj.P.Val_MvsF <= 0.05,'Sig','Not Sig'),
           Direction.type_MvsF = ifelse(logFC_MvsF > 0, 'up', 'down'),
           Sig.direction_MvsF = ifelse(Sig_MvsF == 'Sig',
                                       Direction.type_MvsF,
                                       Sig_MvsF)) %>% 
    mutate(Sig_MvsFD = ifelse(adj.P.Val_MvsFD <= 0.05, 'Sig', 'Not Sig'),
           Direction.type_MvsFD = ifelse(logFC_MvsFD > 0, 'up', 'down'),
           Sig.direction_MvsFD = ifelse(Sig_MvsFD == 'Sig',
                                        Direction.type_MvsFD,
                                        Sig_MvsFD)) %>% 
    mutate(Sig_MvsFS = ifelse(adj.P.Val_MvsFS <= 0.05, 'Sig', 'Not Sig'),
           Direction.type_MvsFS = ifelse(logFC_MvsFS > 0, 'up', 'down'),
           Sig.direction_MvsFS = ifelse(Sig_MvsFS == 'Sig',
                                        Direction.type_MvsFS,
                                        Sig_MvsFS)) -> neuron_cluster.limma.results.sct.df
  
  # save data
write_csv(neuron_cluster.limma.results.sct.df,
          file = paste0('./neurons/stats/limma_trend/cluster_',
                        j, '/', j,
                        '_neuropeptides.neuron_cluster.limmaSCTdf.csv'))
  
  
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
      mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig', Gene, '')) %>% 
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
      ylab(paste0("-log10(adj.P.Val_", i,")")) +
      ggtitle(paste0(i, " volcano plot"))
    
    ggsave(paste0('./neurons/stats/limma_trend/cluster_',
                  j, '/', j,
                  '_limma.neuron_cluster', i, '.reduce.volcano.png'),
           width = 5, height = 5, create.dir = TRUE)
  }
}

#### Red ribbon neuron clusters ####
setwd(paste0(root.dir, "/DGE_CellTypes/neurons"))

# exclude cluster 13
for (j in c(0:12)) {
  
  ## load data 
  tmp = read.csv(paste0('./stats/limma_trend/cluster_', j, '/', j,
                        '_neuropeptides.neuron_cluster.limmaSCTdf.csv'))
  
  ### compare dom vs sub across sexes
  
  # males only
  neuron_cluster.DvsS.M.rrho2 = data.frame(gene = tmp$Gene,
                                           value = -log10(tmp$P.Value_DvsS_M),
                                           direction = tmp$Direction.type_DvsS_M, 
                                           stringsAsFactors = FALSE)
  # set positive negative
  neuron_cluster.DvsS.M.rrho2 = neuron_cluster.DvsS.M.rrho2 %>% 
    mutate(value = ifelse(direction == 'down', -1*value, value)) %>% 
    dplyr::select(-c(direction))
  
  # females only
  neuron_cluster.DvsS.F.rrho2 = data.frame(gene = tmp$Gene,
                                           value = -log10(tmp$P.Value_DvsS_F),
                                           direction = tmp$Direction.type_DvsS_F,
                                           stringsAsFactors = FALSE)
  # set positive negative
  neuron_cluster.DvsS.F.rrho2 = neuron_cluster.DvsS.F.rrho2 %>% 
    mutate(value = ifelse(direction == 'down', -1*value, value)) %>% 
    dplyr::select(-c(direction))

  RedRibbon.all(celltype = paste0("neuron_cluster_", j),
                dataset.a = neuron_cluster.DvsS.M.rrho2,
                dataset.b = neuron_cluster.DvsS.F.rrho2,
                dataset.a.type = "male",
                dataset.b.type = "female",
                a.variable = "value",
                b.variable = "value",
                new.max.log = NULL,
                file.name = './stats/RedRibbon/')
}

# count genes in each quadrant
for (j in c(0:12)) {
  
  ## load data 
  tmp = read.csv(paste0('./stats/RedRibbon/neuron_cluster_', j,
                        '.quadrant.genes.csv'))
  
  # graph count of genes per quadrant
  tmp %>% 
    group_by(RRquadrant) %>% 
    summarise(Total = n()) %>% 
    as.data.frame() %>% 
    ggplot(aes(x = reorder(RRquadrant, -Total),
               y = Total,
               label = Total)) +
    geom_label() +
    theme_classic() +
    ylab('Number of genes per quadrant') +
    xlab('') +
    ggtitle(paste0('Neuron cluster ', j,' RR quadrants'))
  
  ggsave(paste0('./stats/RedRibbon/neuron_cluster_', j,
                'RR_quad_genes_count.png'),
         width = 7, height = 4)
  
}

