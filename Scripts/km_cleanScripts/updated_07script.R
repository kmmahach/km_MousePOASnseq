# R-4.3.1, Seurat v4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# load functions 
source("./functions/network_graph_fun.R")

load_packages(c("tidyverse", "Seurat", "patchwork", "ggrepel", "multcomp", "emmeans", "igraph"),
              out_prefix = "07")


# load data
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

#### get presence/absence (1/0) for neuropeptide.genes ####
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


#### GLM binomial ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# GLM on neuron proportions
neuron.np.glm <- data.frame()
neuron.np.pairwise <- data.frame()
neuron.np.modglm <- data.frame()

for (i in unique(dat.for.stats$neuropeptide)) {
  
  # binomial GLM on sex and status w/ interaction term
  tmp = glm(cbind(Freq,
                  total.neuron.count - Freq) ~ Sex*Status,
            data = subset(dat.for.stats,
                          neuropeptide == i),
            family = binomial)
  
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  data.frame(Estimate = tmp.res$test$coefficients,
             Std.error = tmp.res$test$sigma,
             z.value = tmp.res$test$tstat,
             p.value = tmp.res$test$pvalues) %>%
    rownames_to_column('Predictor') %>%
    mutate(neuron.np = i) -> tmp.res.df
  
  # FDR correction for all comparisons
  tmp.res.df %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(p.adjust.fdr =
             round(p.adjust(p.value,
                            method = 'fdr'),
                   digits = 4) ) -> neuron.np.glm.fdr
  
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint) %>%
    mutate(Predictor = row.names(.)) %>%
    filter(!Predictor == "(Intercept)")
  neuron.np.glm.fdr %>%
    full_join(ci_df, by = c("Predictor",
                            "Estimate")) -> neuron.np.glm.fdr
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
    ggtitle(paste0(i," neurons: binomial GLM"))
  
  ggsave(paste0('neurons/neuropeptides/stats/graphs/', i, '_binomGLM_cell.count.png'),
         width = 6, height = 5)
  
  # Average marginal effects (AMEs)
  ame <- slopes(tmp) %>%
    as.data.frame() %>%
    mutate(neuron.np = i,
           EffectType = "AME")
  
  # Marginal effects at representative values (MERs)
  mer <- slopes(tmp, type = "response") %>%
    as.data.frame() %>%
    mutate(neuron.np = i,
           EffectType = "MER") %>%
    dplyr::select(!c(rowid,
                     predicted_lo,
                     predicted_hi)) %>%
    distinct()
  
  merplots <- plot_mer_facet(mer, group = i)
  
  ggsave(paste0('neurons/neuropeptides/stats/graphs/', i, '_MER_cell.count.png'),
         plot = merplots, width = 10, height = 5)
  
  write_csv(mer,
            file = paste0('neurons/neuropeptides/stats/', i, '_MERtable.csv'))
  
  dat.for.mod <- subset(dat.for.stats,
                        neuropeptide == i) %>%
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
  aie_df <- data.frame(neuron.np = i,
                       EffectType = "AIE",
                       term = "Sex:Status",
                       Estimate = ints$aie$aie.est,
                       Std.error = ints$aie$aie.se.delta,
                       lwr = ints$aie$aie.ll,
                       upr = ints$aie$aie.ll,
                       t.val = NA)
  
  # Interaction at hypothetical values
  inthyp_df <- data.frame(
    neuron.np = i,
    EffectType = "HypotheticalInt",
    term = "Sex:Status",
    Estimate = ints$inthyp$hat,
    Std.error = ints$inthyp$se.int.est,
    lwr = ints$inthyp$hat - ints$inthyp$se.int.est,
    upr = ints$inthyp$hat + ints$inthyp$se.int.est,
    t.val = ints$inthyp$t.val)
  
  neuron.np.mod <- rbind(aie_df, inthyp_df)
  
  # Tukey for all pairwise comparisons
  tmp2 = glm(cbind(Freq,
                   total.neuron.count - Freq) ~ orig.ident,
             data = subset(dat.for.stats,
                           neuropeptide == i),
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
    mutate(neuron.np = i) -> emm.res.df
  
  # combine results and save
  neuron.np.glm = rbind(neuron.np.glm,
                              neuron.np.glm.fdr)
  neuron.np.pairwise = rbind(neuron.np.pairwise,
                                   emm.res.df)
  neuron.np.modglm = rbind(neuron.np.modglm,
                                 neuron.np.mod)
  
  rm(neuron.np.mod, emm.res.df, neuron.np.glm.fdr)
  cat("Processed GLM for", i, "\n")
}


write_csv(neuron.np.glm,
          file = 'neurons/neuropeptides/stats/neuron_npGLM.csv')
write_csv(neuron.np.pairwise,
          file = 'neurons/neuropeptides/stats/neuron_np_pairwise.csv')
write_csv(neuron.np.modglm,
          file = 'neurons/neuropeptides/stats/neuron_np_modglm.csv')

dat.for.stats %>% 
group_by(Sex,
         Status,
         neuropeptide) %>% 
  summarise(percent.avg = mean(Percent)) -> avg.NPs_per_neuron


#### Neuropeptide proportion bias ####
# get list of significant genes for males
neuron.np.pairwise %>% 
  filter(contrast == 'male.dom.data - male.sub.data' & adj.pval <= 0.05) %>% 
  mutate(Male.direction = ifelse(Estimate > 0,
                                 "Male.Dom.Bias",
                                 "Male.Sub.Bias")) %>% 
  dplyr::select(Male.direction,
                neuron.np) %>% 
  dplyr::rename(neuropeptide = neuron.np) -> male.np.results

# plot direction of bias
avg.NPs_per_neuron %>% 
  filter(Sex == 'Male') %>% 
  pivot_wider(id_cols = c('neuropeptide'),
              names_from = 'Status',
              values_from = 'percent.avg') %>% 
  full_join(male.np.results, by = "neuropeptide") %>% 
  mutate(label = ifelse(is.na(Male.direction) & neuropeptide != "Scg2", NA, neuropeptide),
         Male.direction = ifelse(is.na(Male.direction), 'no bias', Male.direction)) %>% 
  ggplot(aes(x = Sub,
             y = Dom,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(color = Male.direction),
             size = 3) +
  geom_text_repel(aes(
    # alpha = Male.color,
    size = 15),
    size = 3, max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = c("#ff7d00", 
                                "#7e38b7",
                                "grey"),
                     name = "Status Bias") +
  # scale_alpha_discrete(guide = "none",
  #                      range = c(1, 0.50)) +
  scale_size(guide = 'none') +
  theme(legend.position = c(0.8, 0.2)) +
  theme(text = element_text(size = 14)) +
  xlab('Subordinate Male neuropeptide proportion (%)') +
  ylab('Dominant Male neuropeptide proportion (%)') +
  ggtitle('Male Neuropeptide Cell Proportion') +
  xlim(0,100) +
  ylim(0,100)

ggsave('./neurons/neuropeptides/stats/graphs/male_countsGLM_statusBias.png',
       height = 5.25,
       width = 5.5)



# get list of significant genes for females
neuron.np.pairwise %>% 
  filter(contrast == 'female.dom.data - female.sub.data' & adj.pval <= 0.05) %>% 
  mutate(Female.direction = ifelse(Estimate > 0,
                                   "Female.Dom.Bias",
                                   "Female.Sub.Bias")) %>% 
  dplyr::select(Female.direction,
                neuron.np) %>% 
  dplyr::rename(neuropeptide = neuron.np) -> female.np.results

# label with significace 
avg.NPs_per_neuron %>% 
  filter(Sex == 'Female') %>% 
  pivot_wider(id_cols = c('neuropeptide'),
              names_from = 'Status',
              values_from = 'percent.avg') %>% 
  full_join((female.np.results %>% 
               add_row(Female.direction = "Female.Dom.Bias", neuropeptide = NA)), 
             by = "neuropeptide") %>% 
  mutate(label = ifelse(is.na(Female.direction), NA, neuropeptide),
         Female.direction = ifelse(is.na(Female.direction), 'no bias', Female.direction)) %>% 
  ggplot(aes(x = Sub,
             y = Dom,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(color = Female.direction),
             size = 3) +
  geom_text_repel(aes(
    # alpha = Female.color,
    size = 15),
    size = 3, max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "grey"),
                     name = "Status Bias") +
  # scale_alpha_discrete(guide = "none",
  #                      range = c(1, 0.50)) +
  scale_size(guide = 'none') +
  theme(legend.position = c(0.8, 0.2)) +
  theme(text = element_text(size = 14)) +
  xlab('Subordinate Female neuropeptide proportion (%)') +
  ylab('Dominant Female neuropeptide proportion (%)') +
  ggtitle('Female Neuropeptide Cell Proportion') +
  xlim(0,100) +
  ylim(0,100)

ggsave('./neurons/neuropeptides/stats/graphs/female_countsGLM_statusBias.png',
       height = 5.25,
       width = 5.5)


#### Distribution of co-expression levels ####
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
                 y = avg.Percent,
                 group = orig.ident,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3,
             shape = 3) +
  geom_smooth(aes(x = count.expr,
                  y = avg.Percent,
                  group = orig.ident,
                  color = orig.ident,
                  fill = orig.ident),
              se = FALSE,
              span = 0.4) + 
  # geom_line(aes(x = count.expr,
  #               y = avg.Percent,
  #               group = orig.ident,
  #               color = orig.ident)) +
  geom_point(data = np.coexpr.summary,
             aes(x = count.expr,
                 y = Percent,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3, alpha = 0.5) +
  geom_line(data = np.coexpr.summary,
            aes(x = count.expr,
                y = Percent,
                color = orig.ident)) +
  labs(y = "Percent of Neurons (exclusive)",
       x = "Number of Co-expressed NPs") +
  theme_classic() + 
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "#ff7d00", 
                                "#7e38b7"))

ggsave('./neurons/neuropeptides/coexprNPs_distribution.png',
       height = 10, width = 10)


# create cumulative freq plot
np.coexpr.summary %>% 
  group_by(orig.ident, count.expr) %>% 
  mutate(avg.Cum.Percent = mean(Cum.Percent)) %>% 
  ggplot() +
  geom_point(aes(x = count.expr,
                 y = avg.Cum.Percent,
                 group = orig.ident,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3,
             shape = 3) +
  geom_smooth(aes(x = count.expr,
                  y = avg.Cum.Percent,
                  group = orig.ident,
                  color = orig.ident,
                  fill = orig.ident),
              se = FALSE,
              span = 0.5) + 
  # geom_line(aes(x = count.expr,
  #               y = avg.Cum.Percent,
  #               group = orig.ident,
  #               color = orig.ident)) +
  geom_point(data = np.coexpr.summary,
             aes(x = count.expr,
                 y = Cum.Percent,
                 color = orig.ident,
                 fill = orig.ident),
             size = 3, alpha = 0.5) +
  geom_line(data = np.coexpr.summary,
            aes(x = count.expr,
                y = Cum.Percent,
                color = orig.ident)) +
  labs(y = "Percent of Neurons (cumulative)",
       x = "Number of Co-expressed NPs") +
  theme_classic() + 
  scale_color_manual(values = c("#f94449", 
                                "#408D8E",
                                "#ff7d00", 
                                "#7e38b7"))

ggsave('./neurons/neuropeptides/coexprNPs_distribution_cumulative.png',
       height = 10, width = 10)  



#### Presence/absence network of co-expression ####
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


sapply(unique(MSCneurons.reclust$orig.ident), 
       \(ident, df, genes, seurat_obj) {
  
  cells = df %>% 
    filter(orig.ident == ident) %>% 
    pull("Cell_ID")
  
  adj.mat = graph.adjacency(crossprod(t(as.matrix(seurat_obj@assays$RNA@data[genes, cells]))),
                            mode = "undirected",
                            weighted = TRUE,
                            diag = FALSE) 
  
  simp.adj = igraph::simplify(adj.mat,
                              remove.multiple = TRUE,
                              remove.loops = TRUE)
  
}, 
df = filtered_neur, genes = topNP, seurat_obj = MSCneurons.reclust) -> adj.matList



sapply(adj.matList, \(adj) {
  
  max.weight = max(E(adj)$weight)
  min.weight = min(E(adj)$weight)
  
  df = data.frame(min.weight, max.weight)
  sapply(df, as.numeric)
  
}) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(scale.max = 10*max.weight/max(max.weight),
         scale.min = 1*min.weight/max(min.weight)) -> net.weights


plotlist = list()

for(adj in names(adj.matList)) {
  
  group = sub("\\.data", "",  adj)
  title = sub("\\.", " ", group) %>% str_to_title()
  adj.mat = adj.matList[adj]
  scale.min = net.weights[net.weights$ID == adj, "scale.min"]
  scale.max = net.weights[net.weights$ID == adj, "scale.max"]
  
  lapply(adj.mat, \(adj) {
    adj %>% 
    ggraph(layout = 'linear',
           circular = TRUE,
           sort.by = vSizes.all) +
      geom_edge_parallel(aes(width = weight,
                             alpha = weight,
                             start_cap = label_rect(node1.name),
                             end_cap = label_rect(node2.name))) +
      geom_node_label(aes(label = toupper(name),
                          size = vSizes.all)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),
            title = element_text(size = 13,
                                 face = 'bold'),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks  = element_blank(),
            axis.text = element_blank(),
            legend.position = 'none') +
      ggtitle(title) +
      scale_edge_width(range = c(scale.min, scale.max)) +
      scale_edge_alpha(range = c(0, scale.max/10)) 
    
  }) -> plotlist[group]

}
  

ggsave(file = "./neurons/neuropeptides/coexpression_networks.png",
       wrap_plots(plotlist, ncol = 2),
       width = 8.5, height = 8)

save(plotlist,
     file = paste0(root.dir, "/manuscriptFigures/plotlist.rda"))
