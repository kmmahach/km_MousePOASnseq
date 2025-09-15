# R-4.3.1, Seurat v.4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# libraries 
load_packages(c("Seurat", "limma", "annotables", "patchwork", "tidyverse", "edgeR", "tidytext", "ggpol"),
              out_prefix = "4.75")

#### load data ####
load("data/integrated_seurat_onlyNeurons.rda")

# subset to Oxt/Avp expressing neurons
MSCneurons.reclust@active.assay = "RNA"

np.sub <- subset_by_gene(MSCneurons.reclust,
                          subset_genes = c("Oxt", "Avp"),
                          min_count = 2)

all.meta <- as.data.frame(MSCneurons.reclust@meta.data)

oxt.exp <- as.data.frame(np.sub$Oxt@assays$RNA@data) %>% 
  filter(rownames(.) == "Oxt") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Oxt_expr") 

oxt.meta <- all.meta %>% 
  filter(Cell_ID %in% oxt.exp$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, seurat_clusters, parent_id.broad, Cell_ID) 

oxt.exp = left_join(oxt.exp, oxt.meta,
                    by = "Cell_ID")

oxt.summary <- oxt.exp %>% 
  group_by(orig.ident, indiv_genotype) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))

oxt.summary %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  # facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/meanNormCounts_OXT_RNAdata.png"),
       width = 6, height = 4)

oxt.by.clust <- oxt.exp %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))

oxt.by.clust %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/cluster_NormCounts_OXT_RNAdata.png"),
       width = 10, height = 4)

avp.exp <- as.data.frame(np.sub$Avp@assays$RNA@data) %>% 
  filter(rownames(.) == "Avp") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Avp_expr")

avp.meta <- all.meta %>% 
  filter(Cell_ID %in% avp.exp$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, seurat_clusters, parent_id.broad, Cell_ID)

avp.exp = left_join(avp.exp, avp.meta,
                    by = "Cell_ID")


avp.summary <- avp.exp %>% 
  group_by(orig.ident, indiv_genotype) %>% 
  summarise(avg_AvpExpr = mean(Avp_expr))

avp.summary %>% 
  ggplot(aes(orig.ident, avg_AvpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  # facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/meanNormCounts_AVP_RNAdata.png"),
       width = 6, height = 4)

avp.by.clust <- avp.exp %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(avg_AvpExpr = mean(Avp_expr))

avp.by.clust %>% 
  ggplot(aes(orig.ident, avg_AvpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/cluster_NormCounts_AVP_RNAdata.png"),
       width = 10, height = 4)


#### Diffs in cell type? ####
cluster0_avp <- avp.exp %>% 
  filter(seurat_clusters == 0)

avp.clust0.sum <- cluster0_avp %>% 
  group_by(orig.ident, indiv_genotype, parent_id.broad) %>% 
  summarise(avg_AvpExpr = mean(Avp_expr))

avp.clust0.sum %>% 
  ggplot(aes(orig.ident, avg_AvpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ parent_id.broad, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in cluster0 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust0subtypes_NormCounts_AVP_RNAdata.png"),
       width = 7, height = 4)

cluster0_oxt <- oxt.exp %>% 
  filter(seurat_clusters == 0)

oxt.clust0.sum <- cluster0_oxt %>% 
  group_by(orig.ident, indiv_genotype, parent_id.broad) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))

oxt.clust0.sum %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ parent_id.broad, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in cluster0 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust0subtypes_NormCounts_OXT_RNAdata.png"),
       width = 7, height = 4)


avpclust0.meta <- all.meta %>% 
  filter(Cell_ID %in% avp.exp$Cell_ID &
           seurat_clusters == 0) %>% 
  select(predicted.prob, Cell_ID)

avp.cluster0.exp = left_join(cluster0_avp, avpclust0.meta,
                    by = "Cell_ID")

avp.clust0.subtypes <- avp.cluster0.exp %>% 
  group_by(orig.ident, indiv_genotype, predicted.prob) %>% 
  summarise(avg_AvpExpr = mean(Avp_expr))

avp.clust0.subtypes %>% 
  ggplot(aes(orig.ident, avg_AvpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ predicted.prob, nrow = 6) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in cluster0 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust0allSubtypes_NormCounts_AVP_RNAdata.png"),
       width = 15, height = 13)

oxtclust0.meta <- all.meta %>% 
  filter(Cell_ID %in% oxt.exp$Cell_ID &
           seurat_clusters == 0) %>% 
  select(predicted.prob, Cell_ID)

oxt.cluster0.exp = left_join(cluster0_oxt, oxtclust0.meta,
                             by = "Cell_ID")

oxt.clust0.subtypes <- oxt.cluster0.exp %>% 
  group_by(orig.ident, indiv_genotype, predicted.prob) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))

oxt.clust0.subtypes %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ predicted.prob, nrow = 5) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in cluster0 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust0allSubtypes_NormCounts_OXT_RNAdata.png"),
       width = 17, height = 13)


# GLU/GABA percent of clusters
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
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- neuron.props.ident %>%
  ggplot(aes(orig.ident, n_cells)) +
  geom_bar(aes(fill = parent_id.broad),
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

p2 <- neuron.props.ident %>% 
  ggplot(aes(cluster, n_cells)) +
  geom_bar(aes(fill = parent_id.broad), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Subtype",
       x = "Neuron Cluster",
       y = "") +
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

ggsave(paste0(root.dir, '/DGE_CellTypes/neurons/cluster_stats/composition_GABAvsGLU_ClustersByNum.png'),
       p1 + p2 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 18, height = 8)

### Just the Caprin2 GLU neurons ####
caprin2 <- subset(MSCneurons.reclust, 
                  MSCneurons.reclust$predicted.prob == "C66-22: Caprin2.GLU-6")

caprin2.exp.Avp <- as.data.frame(caprin2@assays$RNA@data) %>% 
  filter(rownames(.) == "Avp") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Avp_expr") 

caprin2.Avp.meta <- all.meta %>% 
  filter(Cell_ID %in% caprin2.exp.Avp$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, seurat_clusters, parent_id.broad, Cell_ID) 

caprin2.exp.Avp = left_join(caprin2.exp.Avp, caprin2.Avp.meta,
                    by = "Cell_ID")

caprin2.Avp.summary <- caprin2.exp.Avp %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(avg_avpExpr = mean(Avp_expr))


caprin2.Avp.summary %>% 
  ggplot(aes(orig.ident, avg_avpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in caprin2+ GLU neurons",
          subtitle = "(by cluster)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/Avp_normCounts_inCaprin2Neurons.png"),
       width = 6, height = 3.5)


caprin2_freq_Avp <- caprin2.exp.Avp %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(n_cells = n())

caprin2_freq_Avp %>% 
  ggplot(aes(orig.ident, n_cells, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Cell count", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("Number of caprin2+ GLU neurons",
          subtitle = "(that expr Avp)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/AVPcells_inCaprin2Neurons.png"),
       width = 6, height = 3.5)

caprin2.exp.Avp %>% 
ggplot(aes(orig.ident, Avp_expr, 
           color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("Avp expr in caprin2+ GLU neurons by cluster",
          subtitle = "(individual neurons)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/indiv_neuronAVPexpr_inCaprin2Neurons.png"),
       width = 6, height = 3.5)

####

caprin2.exp.Oxt <- as.data.frame(caprin2@assays$RNA@data) %>% 
  filter(rownames(.) == "Oxt") %>% 
  pivot_longer(cols = everything(),
               names_to = "Cell_ID",
               values_to = "Oxt_expr") 

caprin2.Oxt.meta <- all.meta %>% 
  filter(Cell_ID %in% caprin2.exp.Oxt$Cell_ID) %>% 
  select(orig.ident, indiv_genotype, seurat_clusters, parent_id.broad, Cell_ID) 

caprin2.exp.Oxt = left_join(caprin2.exp.Oxt, caprin2.Oxt.meta,
                            by = "Cell_ID")

caprin2.Oxt.summary <- caprin2.exp.Oxt %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))


caprin2.Oxt.summary %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in caprin2+ GLU neurons",
          subtitle = "(by cluster)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/Oxt_normCounts_inCaprin2Neurons.png"),
       width = 6, height = 3.5)


caprin2_freq_oxt <- caprin2.exp.Oxt %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(n_cells = n())

caprin2_freq_oxt %>% 
  ggplot(aes(orig.ident, n_cells, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Cell count", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("Number of caprin2+ GLU neurons",
          subtitle = "(that expr Oxt)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/OXTcells_inCaprin2Neurons.png"),
       width = 6, height = 3.5)

caprin2.exp.Oxt %>% 
  ggplot(aes(orig.ident, Oxt_expr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("Oxt expr in caprin2+ GLU neurons by cluster",
          subtitle = "(individual neurons)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/indiv_neuronOXTexpr_inCaprin2Neurons.png"),
       width = 6, height = 3.5)


####

caprin2_freq_all <- as.data.frame(caprin2@meta.data)
caprin2_freq_sum <- caprin2_freq_all %>% 
  group_by(orig.ident, indiv_genotype, seurat_clusters) %>% 
  summarise(n_cells = n())


caprin2_freq_sum %>% 
  ggplot(aes(orig.ident, n_cells, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Cell count", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("Total number of caprin2+ GLU neurons",
          subtitle = "(by cluster)")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/Caprin2Neurons_count.png"),
       width = 6, height = 3.5)


### Coexpression?

caprin2.coexpr <- as.data.frame(caprin2@assays$RNA@data) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp")) %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = !gene,
               names_to = "Cell_ID",
               values_to = "norm_counts") %>% 
  pivot_wider(names_from = "gene",
              values_from = norm_counts) %>% 
  mutate(coexpr = ifelse(Oxt > 0 & Avp > 0, "yes", "no"))

caprin2.coexpr.meta <- as.data.frame(caprin2@meta.data) %>% 
  filter(Cell_ID %in% caprin2.coexpr$Cell_ID)

caprin2.coexpr <- left_join(caprin2.coexpr, caprin2.coexpr.meta,
                            by = "Cell_ID")

caprin2.coexpr %>% 
  ggplot(aes(orig.ident, Avp, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ coexpr, nrow = 1) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 13)) + 
  ggtitle("Avp expr in caprin2+ GLU neurons",
          subtitle = "(indiv neurons - either coexpress Oxt or not)") -> coex1


caprin2.coexpr %>% 
  ggplot(aes(orig.ident, Oxt, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ coexpr, nrow = 1) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 13)) + 
  ggtitle("Oxt expr in caprin2+ GLU neurons",
          subtitle = "(indiv neurons - either coexpress Avp or not)") -> coex2


ggsave(plot = grid.arrange(coex1, coex2), 
       file = paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/coexpr_OxtAvp_inCaprin2Neurons.png"),
       width = 7, height = 6.5)


count_sum <- caprin2.coexpr %>% 
  group_by(orig.ident, indiv_genotype, coexpr) %>% 
  summarise(n_cells = n()) 


count_sum %>% 
  ggplot(aes(orig.ident, n_cells, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ coexpr, nrow = 1) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Cell count", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 13)) + 
  ggtitle("caprin2+ GLU neurons that coexpr Avp/Oxt")


ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/caprin2GLU/Caprin2Neurons_coexprOxtAvp_cellCount.png"),
       width = 6, height = 3.5)

caprin2.coexpr %>% 
  filter(coexpr == "yes") %>% 
  ggplot(aes(x = Avp,
             y = Oxt,
             color = orig.ident)) +
  geom_point()

# all neurons?
oxt.avp.coexpr_all <- as.data.frame(MSCneurons.reclust@assays$RNA@counts) %>% 
  filter(rownames(.) %in% c("Oxt", "Avp", "Caprin2")) %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(cols = !gene,
               names_to = "Cell_ID",
               values_to = "norm_counts") %>% 
  pivot_wider(names_from = "gene",
              values_from = norm_counts) %>% 
  mutate(coexpr = ifelse(Oxt > 0 & Avp > 0, "yes", "no"))

oxt.avp.coexpr_all.meta <- as.data.frame(MSCneurons.reclust@meta.data) 

oxt.avp.coexpr_all <- left_join(oxt.avp.coexpr_all, oxt.avp.coexpr_all.meta,
                            by = "Cell_ID")

oxt.avp.coexpr_all %>% 
filter(coexpr == "yes" & !(predicted.prob == "C66-22: Caprin2.GLU-6") & Oxt < 20) %>% 
  ggplot(aes(x = Avp,
             y = Oxt,
             color = orig.ident)) +
  geom_point(position = position_jitter(),
             alpha = 0.5)




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


# cells per cluster by orig.ident
predicted.prob.df %>%
  filter(!(predicted.prob == "Unknown")) %>% 
  group_by(orig.ident, seurat_clusters, Region_summarized) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         orig.ident = sub("\\.[^.]*$", "", orig.ident)) -> region.props.ident

MSCneurons.reclust@meta.data %>%
  mutate(orig.ident = sub("\\.[^.]*$", "", orig.ident)) %>% 
  group_by(orig.ident) %>%
  tally() -> temp_labels

p1 <- region.props.ident %>%
  ggplot(aes(orig.ident, n_cells)) +
  geom_bar(aes(fill = Region_summarized),
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

p2 <- region.props.ident %>% 
  ggplot(aes(cluster, n_cells)) +
  geom_bar(aes(fill = Region_summarized), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  labs(fill = "Predicted brain region",
       x = "Neuron Cluster",
       y = "") +
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

ggsave(paste0(root.dir, '/DGE_CellTypes/neurons/cluster_stats/composition_predictedRegion_ClustersByNum.png'),
       p1 + p2 +
         plot_layout(ncol = 2, widths = c(
           MSCneurons.reclust@meta.data$orig.ident %>% 
             unique() %>% length(),
           MSCneurons.reclust@meta.data$integrated_snn_res.0.4 %>%
             unique() %>% length() )),
       width = 18, height = 8)


# Region as percent of Caprin2 neurons
predicted.prob.df %>% 
  filter(!(predicted.prob == "Unknown") & 
           Cell_ID %in% caprin2@meta.data$Cell_ID) %>% 
  group_by(seurat_clusters, Region_summarized) %>% 
  summarise(n_cells = n()) %>% 
  rename(cluster = seurat_clusters) %>% 
  group_by(cluster) %>% 
  mutate(prop = n_cells / sum(n_cells),
         cluster = as.numeric(cluster) - 1) -> caprin2.region.props


### cluster 7 oxt/avp subtypes ####
cluster7_avp <- avp.exp %>% 
  filter(seurat_clusters == 7)

avpclust7.meta <- all.meta %>% 
  filter(Cell_ID %in% avp.exp$Cell_ID &
           seurat_clusters == 7) %>% 
  select(predicted.prob, Cell_ID)

avp.cluster7.exp = left_join(cluster7_avp, avpclust7.meta,
                             by = "Cell_ID")

avp.clust7.subtypes <- avp.cluster7.exp %>% 
  group_by(orig.ident, indiv_genotype, predicted.prob) %>% 
  summarise(avg_AvpExpr = mean(Avp_expr))

avp.clust7.subtypes %>% 
  ggplot(aes(orig.ident, avg_AvpExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ predicted.prob, nrow = 6) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("AVP expr in cluster7 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust7allSubtypes_NormCounts_AVP_RNAdata.png"),
       width = 15, height = 13)


cluster7_oxt <- oxt.exp %>% 
  filter(seurat_clusters == 7)

oxtclust7.meta <- all.meta %>% 
  filter(Cell_ID %in% oxt.exp$Cell_ID &
           seurat_clusters == 7) %>% 
  select(predicted.prob, Cell_ID)

oxt.cluster7.exp = left_join(cluster7_oxt, oxtclust7.meta,
                             by = "Cell_ID")

oxt.clust7.subtypes <- oxt.cluster7.exp %>% 
  group_by(orig.ident, indiv_genotype, predicted.prob) %>% 
  summarise(avg_oxtExpr = mean(Oxt_expr))

oxt.clust7.subtypes %>% 
  ggplot(aes(orig.ident, avg_oxtExpr, 
             color = orig.ident, fill = orig.ident))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ predicted.prob, nrow = 5) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 7.5)) + 
  ggtitle("OXT expr in cluster7 neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/clust7allSubtypes_NormCounts_OXT_RNAdata.png"),
       width = 17, height = 13)
