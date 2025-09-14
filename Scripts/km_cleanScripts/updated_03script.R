# R-4.3.1, Seurat v.4.4.0
# Reclustering with neuronal nuclei only (exclude all glia and blood cells)

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
set.seed(12345)
compression = "xz" # slower, but usually smallest compression
setwd(root.dir) 

# functions
source("./Scripts/km_cleanScripts/functions/DGE_fun.R")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "Seurat", "limma", "edgeR", "RedRibbon", "RRHO2", "marginaleffects",
                 "modglm", "emmeans", "multcomp", "gridExtra", "gghalves", "patchwork"), 
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
  DimPlot() +
  ggtitle('Neuron Clusters') +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()

ggsave('./neurons/UMAP/neurons.recluster.clusters.png',
       height = 6, width = 6)

# broad cell types
MSCneurons.reclust %>%
  DimPlot(group.by = 'parent_id.broad.prob') +
  ggtitle('Excitatory vs. Inhibitory Neurons') +
  scale_color_manual(values = c('#1b9e77',
                                '#d95f02')) +
  theme_classic()

ggsave('./neurons/UMAP/neurons.recluster.parent_id.broad.png',
       height = 5, width = 6)

# specific sub-types
MSCneurons.reclust %>%
  DimPlot(group.by = 'predicted.prob') +
  ggtitle('Neuron Subtypes') +
  theme_classic()

ggsave('./neurons/UMAP/neurons.recluster_predicted.prob.png',
       height = 5, width = 11)

# sex/social status
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

ggsave('./neurons/UMAP/neurons.recluster.orig.ident.png',
       height = 5.5, width = 7)

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
       height = 5, width = 7)

# indiv_genotype (predicted by souporcell) for each group
# set idents
Idents(MSCneurons.reclust) = "orig.ident"
subset_idents <- unique(Idents(MSCneurons.reclust))

# subset Seurat obj
reclust <- subset_by_ident(MSCneurons.reclust,
                           subset_idents = subset_idents,
                           cluster = FALSE)

lapply(reclust, \(x) {
  x %>%
    DimPlot(group.by = 'indiv_genotype') +
    theme_classic() +
    ggtitle(paste0(x@project.name, ' - indiv_genotype'))
  ggsave(paste0('./neurons/UMAP/',
                sub("\\.data$", "", x@project.name),
                '_neurons.recluster_genotypes.png'),
         height = 5, width = 6)
  }
)

# dom candidate genes?
DotPlot(MSCneurons.reclust,
        assay = 'SCT', features = c("Lhx8", "Crym", "Ache", "Myb", "Mog")) +
  theme_classic()

ggsave('./neurons/dom.genes.dotplot.png',
       height = 5, width = 6)

# scaled counts of neurons in each cluster (as percent)
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


ggsave('./neurons/clustersVorig.ident.neurons.recluster.genotype.scaled.png',
       height = 7, width = 20)

#### Prep for stats ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# get cluster names?
select_neur %>%
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

# save neuron cluster names
write_csv(neuron.cluster.names, 'neurons/cluster_stats/neuron_clusterNames.csv')

# format data
sum_neur_indiv %>%
  pivot_longer(
    cols = matches("^(Percent|Freq)_"),
    names_to = c(".value", "seurat_clusters"),
    names_pattern = "(Percent|Freq)_(.*)"
  ) %>%
  dplyr::select(orig.ident,
                indiv_genotype,
                seurat_clusters,
                Freq,
                total.neuron.count,
                Percent) %>% 
  mutate(seurat_clusters = str_remove(seurat_clusters, 
                                      "cluster_") %>% 
           as.numeric()) %>% 
  mutate(Sex =
           ifelse(grepl("female", orig.ident),
                  "Female",
                  "Male")) %>%
  mutate(Status =
           ifelse(grepl("dom", orig.ident),
                  "Dom",
                  "Sub")) -> dat.for.stats

#### ANOVA ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

### run anova on every cluster
# create empty data frame
neuron.cluster.aov = data.frame(matrix(ncol = 6, nrow = 0))
# provide column names
colnames(neuron.cluster.aov) <- c("Df",
                                  "Sum Sq",
                                  "Mean Sq",
                                  "F value",
                                  "Pr(>F)",
                                  "seurat_clusters")
## loop through each cluster
for (i in unique(dat.for.stats$seurat_clusters)) {
  
  ## run two way anova on sex and status
  tmp <- aov(Percent ~ Sex*Status,
             data = dat.for.stats %>%
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

# adjust pvalue (FDR)
neuron.cluster.aov %>%
  mutate(p_adj = p.adjust(`Pr(>F)`,
                          method = 'fdr',
                          n = nrow(neuron.cluster.aov))) %>%
  rownames_to_column('term') -> neuron.cluster.aov

# clean up
neuron.cluster.aov %>%
  mutate(Predictor = case_when(grepl("Sex", term) ~ "Sex",
                               grepl("Status", term) ~ "Status")) %>%
  mutate(Predictor = ifelse(grepl("Sex:Status", term),
                            "Sex:Status", Predictor)) %>%
  dplyr::select(!(c(term, Df))) -> neuron.cluster.aov

neuron.cluster.aov <- neuron.cluster.aov[order(neuron.cluster.aov$seurat_clusters),]


# save anova table
write_csv(neuron.cluster.aov, 'neurons/cluster_stats/neuron_clusterANOVA.csv')


#### GLM binomial ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# GLM on neuron proportions
neuron.clusters.glm <- data.frame()
neuron.clusters.pairwise <- data.frame()
neuron.clusters.modglm <- data.frame()

for (i in unique(dat.for.stats$seurat_clusters)) {
  
  # binomial GLM on sex and status w/ interaction term
  tmp = glm(cbind(Freq,
                  total.neuron.count - Freq) ~ Sex*Status,
            data = subset(dat.for.stats,
                          seurat_clusters == i),
            family = binomial)
  
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  data.frame(Estimate = tmp.res$test$coefficients,
             Std.error = tmp.res$test$sigma,
             z.value = tmp.res$test$tstat,
             p.value = tmp.res$test$pvalues) %>%
    rownames_to_column('Predictor') %>%
    mutate(neuron.clusters = i) -> tmp.res.df
  
  # FDR correction for all comparisons
  tmp.res.df %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(p.adjust.fdr =
             round(p.adjust(p.value,
                            method = 'fdr'),
                   digits = 4) ) -> neuron.clusters.glm.fdr
  
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint) %>%
    mutate(Predictor = row.names(.)) %>%
    filter(!Predictor == "(Intercept)")
  neuron.clusters.glm.fdr %>%
    full_join(ci_df, by = c("Predictor",
                            "Estimate")) -> neuron.clusters.glm.fdr
  # graph results
  ci_df %>%
    ggplot(aes(x = Estimate,
               y = Predictor)) +
    geom_point() +
    geom_errorbar(aes(xmin = lwr,
                      xmax = upr),
                  width = 0.2) +
    geom_vline(xintercept = 0,
               linetype = 2) +
    labs(x = "Difference in Means",
         y = "Comparison") +
    theme_classic() +
    ggtitle(paste0("Cluster", i," neurons: binomial GLM"))
  
  ggsave(paste0('neurons/cluster_stats/neuron_clusterGLM/graphs/cluster_', i, '_binomGLM_cell.count.png'),
         width = 6, height = 5)
  
  # Average marginal effects (AMEs)
  ame <- slopes(tmp) %>%
    as.data.frame() %>%
    mutate(neuron.clusters = i,
           EffectType = "AME")
  
  # Marginal effects at representative values (MERs)
  mer <- slopes(tmp, type = "response") %>%
    as.data.frame() %>%
    mutate(neuron.clusters = i,
           EffectType = "MER") %>%
    dplyr::select(!c(rowid,
                     predicted_lo,
                     predicted_hi)) %>%
    distinct()
  
  merplots <- plot_mer_facet(mer, group = paste0("Cluster", i))
  
  ggsave(paste0('neurons/cluster_stats/neuron_clusterGLM/graphs/cluster_', i, '_MER_cell.count.png'),
         plot = merplots, width = 10, height = 5)
  
  write_csv(mer,
            file = paste0('neurons/cluster_stats/neuron_clusterGLM/MER/cluster_', i, '_MERtable.csv'))
  
  dat.for.mod <- subset(dat.for.stats,
                        seurat_clusters == i) %>%
    mutate(Sex = ifelse(Sex=="Female", 1, 0),
           Status = ifelse(Status=="Dom", 1, 0)) %>%
    dplyr::select(Freq, total.neuron.count, Sex, Status)
  
  tmp.for.mod <- glm(cbind(Freq,
                           total.neuron.count - Freq) ~ Sex*Status,
                     data = dat.for.mod,
                     family = binomial)
  
  # modglm interaction effects
  ints <- modglm(model = tmp.for.mod,
                 vars = c("Sex", "Status"),
                 data = dat.for.mod,
                 hyps = "means",
                 type = "dd")
  
  # Average interaction effect (AIE)
  aie_df <- data.frame(neuron.clusters = i,
                       EffectType = "AIE",
                       term = "Sex:Status",
                       Estimate = ints$aie$aie.est,
                       Std.error = ints$aie$aie.se.delta,
                       lwr = ints$aie$aie.ll,
                       upr = ints$aie$aie.ll,
                       t.val = NA)
  
  # Interaction at hypothetical values
  inthyp_df <- data.frame(
    neuron.clusters = i,
    EffectType = "HypotheticalInt",
    term = "Sex:Status",
    Estimate = ints$inthyp$hat,
    Std.error = ints$inthyp$se.int.est,
    lwr = ints$inthyp$hat - ints$inthyp$se.int.est,
    upr = ints$inthyp$hat + ints$inthyp$se.int.est,
    t.val = ints$inthyp$t.val)
  
  neuron.clusters.mod <- rbind(aie_df, inthyp_df)
  
  # Tukey for all pairwise comparisons
  tmp2 = glm(cbind(Freq,
                   total.neuron.count - Freq) ~ orig.ident,
             data = subset(dat.for.stats,
                           seurat_clusters == i),
             family = binomial)
  
  # use emmeans to get pairwise differences
  emm = emmeans(tmp2, ~ "orig.ident")
  emm.res = pairs(emm)
  data.frame(contrast = summary(emm.res)$contrast,
             Estimate = summary(emm.res)$estimate,
             Std.error = summary(emm.res)$SE,
             z.ratio = summary(emm.res)$z.ratio,
             adj.pval = round((summary(emm.res)$p.value),
                              digits = 4)) %>%
    mutate(neuron.clusters = i) -> emm.res.df
  
  # combine results and save
  neuron.clusters.glm = rbind(neuron.clusters.glm,
                              neuron.clusters.glm.fdr)
  neuron.clusters.pairwise = rbind(neuron.clusters.pairwise,
                                   emm.res.df)
  neuron.clusters.modglm = rbind(neuron.clusters.modglm,
                                 neuron.clusters.mod)
  
  rm(neuron.clusters.mod, emm.res.df, neuron.clusters.glm.fdr)
  cat("Processed cluster", i, "\n")
}


write_csv(neuron.clusters.glm,
          file = 'neurons/cluster_stats/neuron_clusterGLM/neuron_clusterGLM.csv')
write_csv(neuron.clusters.pairwise,
          file = 'neurons/cluster_stats/neuron_clusterGLM/neuron_clusters_pairwise.csv')
write_csv(neuron.clusters.modglm,
          file = 'neurons/cluster_stats/neuron_clusterGLM/neuron_clusters_modglm.csv')


#### Neuron cluster distributions ####
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
  

# Percent of cells per cluster by sample
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
  scale_y_continuous(name = 'Percentage [%]', 
                     labels = scales::percent_format(),
                     expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Group") +
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


# cells per cluster by orig.ident
MSCneurons.reclust@meta.data %>%
  group_by(orig.ident) %>%
  tally() %>% 
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) -> temp_labels

p3 <- table_orig.idents_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
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
  scale_y_continuous(name = 'Number of cells', 
                     labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Cluster") +
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

ggsave('neurons/composition_orig.idents_ClustersByNum.png',
  p3 + p4 +
    plot_layout(ncol = 2, widths = c(
      MSCneurons.reclust@meta.data$orig.ident %>% 
        unique() %>% length(),
      MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
        unique() %>% length() )),
  width = 18, height = 8)

p5 <- p4 + 
  # scale_fill_viridis_d(option = "magma",
  #                      # direction = -1,
  #                      alpha = 0.85) +
  scale_fill_manual(values = adjustcolor(c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"),
                                         alpha.f = 0.75)) +
  labs(x = "Neuron Cluster",
       y = "Total Nuclei") +
  theme(axis.title = element_text(size = 16))
p5

p6 <- p1 +
  scale_fill_viridis_d(option = "plasma",
                       alpha = 0.75)

p6


ggsave('neurons/composition_neuron.clusters_byPct_Num.png',
       p6 + p5 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 18, height = 8)


#### Build neuron hierarchical tree ####

# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

setwd(paste0(root.dir, "/DGE_CellTypes"))

# hierarchical cluster tree
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

ggsave('neurons/neuron.cluster.tree.png',
       height = 6, width = 6)


#### Neuron marker specificity score ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

Idents(MSCneurons.reclust) = "seurat_clusters"

# Prep seurat object for running find markers - need to run before FindAllMarkers
PrepSCTFindMarkers(MSCneurons.reclust,
                   assay = "SCT") -> MSCneurons.reclust

# calculate specificity score with DEGs for each cluster
FindAllMarkers(MSCneurons.reclust, 
               logfc.threshold = 0.1, # increase to increase speed
               test.use = 'wilcox',
               verbose = TRUE,
               min.pct = 0.10, 
               assay = 'SCT') -> neuron.markers.df 

# create gene column
neuron.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2)) -> neuron.markers.df

write_csv(neuron.markers.df, 'neurons/cluster_stats/neuron_clusterMarkers.csv')

# heatmap of 6 top markers per cluster - reduce to top 6 genes per cluster
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
                     breaks = c(-0.1, 0.1, 0.5, 2, 5, 15),
                     legend_breaks = c(1, 5, 15, 30),
                     angle_col = 45,
                     fontsize = 15,
                     fontsize_row = 15,
                     fontsize_col = 12,
                     fontsize_legend = 7,
                     col = c('white',
                             'lightblue',
                             'yellow',
                             'orange',
                             'red'),
                     border_color = "white",
                     cluster_rows = F,
                     cluster_cols = F,
                     main = 'Cluster Markers',
                     filename = 'neurons/neuron_markers.heatmap.png',
                     width = 5,
                     height = 15)


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

write_csv(markers.weighted.spec, 'neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')

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

ggsave(file = "./neurons/cluster_marker.test.png",
       width = 7.5, height = 15)


#### Neuron cluster limma-trend (DGE) ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# set idents
Idents(MSCneurons.reclust) = "seurat_clusters"
DefaultAssay(MSCneurons.reclust) = "SCT"

# get clusters
subset_idents <- sort(unique(Idents(MSCneurons.reclust)))

# exclude cluster 13 for now - not enough nuclei
subset_idents = subset_idents[-length(subset_idents)]

# subset Seurat obj
l.dfs <- subset_by_ident(MSCneurons.reclust, 
                         subset_idents = subset_idents,
                         cluster = TRUE)

# raw counts and norm.counts of variable genes
dge_data <- prep.for.DGE(l.dfs, 
                         selection.method = "vst", 
                         SCTransformed = TRUE, 
                         assay = 'integrated')

# get results and graph p-values
limma_results <- run_limmatrend(dge_data$results, "./neurons/cluster_stats/limma_trend")

save(limma_results, file = "./neurons/cluster_stats/limma_trend/neuronCluster_limma_results.rda")

#### RRHO/RedRibbon ####
setwd(paste0(root.dir, "/DGE_CellTypes/neurons"))

# how to determine max log scale before graphing? 
# max.log.scale = 100;

# make sure these match upper/lower case
sex = c("Female", "Male");
status = c("Dom", "Sub");

rrho_results <- get.RRHO(limma_results,
                         group.by = sex,
                         compare.across = status,
                         # new.max.log = max.log.scale,
                         outdir = "./cluster_stats/RRHO")

save(rrho_results, file = "./cluster_stats/RRHO/rrho_results.rda")


# graph gene counts in each RR quadrant
plot.RRHO.counts(rrho_results, "./cluster_stats/RRHO")
