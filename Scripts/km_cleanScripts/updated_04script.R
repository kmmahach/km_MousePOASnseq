# R-4.3.1, Seurat v.4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# libraries 
load_packages(c("Seurat", "RRHO2", "DEsingle", "dendextend", "RedRibbon", "emmeans", "limma",
              "pheatmap", "ggalluvial", "multcomp", "clustree", "patchwork", "tidyverse", "edgeR"),
              out_prefix = "04")

#### load data ####
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# subset with SCT data
DefaultAssay(int.ldfs) = "SCT"

#### Neuropeptide candidates ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# load neuropeptide list - where did this come from? 
neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))

# genes
neuropeptides.genes = neuropeptides %>% 
  pull(Gene.name)

select_genes <- c("Avp", 
                  "Avpr1a", 
                  "Avpr2",
                  "Oxt",
                  "Oxtr",
                  "Esr1",
                  "Esr2",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2",
                  "Crh")

# subset genes
subset_genes <- c("Esr1",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2")

# graph overall neurons
DotPlot(int.ldfs,
        assay = 'SCT',
        features = select_genes,
        group.by = 'orig.ident',
        col.min = 0) +
  theme_classic()

ggsave('neurons/neuroendocrine_genes/selectGenes_neurons.png',
       height = 4, width = 9)

#### Candidate neuroendocrine genes ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# subset by select_genes
sub.MSCneurons <- subset_by_gene(int.ldfs,
                                 select_genes,
                                 slot = "data",
                                 min_count = 0.5)

# get presence/absence (1/0) for select_genes
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

# counts & percentages
sum_neur_indiv <- select_neur %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(cell_count = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

# percent of neurons expr select genes by individual
sum_neur_indiv %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "gene",
    names_prefix = "pct_",
    values_to = "percentage"
  ) %>%
  ggplot(aes(x = reorder(gene, -percentage),
             y = percentage,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF")) +
  xlab('') + ylab('% of neurons') +
  ggtitle('Percent nuclei expressing neur-endo genes')

ggsave('neurons/neuroendocrine_genes/pctNeurons_exprSelectGenes_by.indiv_genotype.png',
       height = 5, width = 10)

# percent of neurons expr select genes by group
sum_neur_group <- select_neur %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

sum_neur_group %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "gene",
    names_prefix = "pct_",
    values_to = "percentage"
  ) %>%
  ggplot(aes(x = reorder(gene, -percentage),
             y = percentage,
             color = orig.ident)) +
  geom_point(size = 3) +
  theme_classic() +
  ggtitle('Percent nuclei expressing neur-endo genes')

ggsave('neurons/neuroendocrine_genes/pctNeurons_exprSelectGenes_by_orig.ident.png',
       height = 5, width = 10)

# total - all groups combined
sum_neur_group %>%
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "sum"
  ) %>%
  dplyr::select(orig.ident,
                gene,
                cell_count,
                sum) %>%
  group_by(gene) %>%
  summarise(percentage = sum(sum)/sum(cell_count) * 100) %>%
  filter(percentage < 99) %>%
  mutate(percentage = signif(percentage, 2)) %>%
  ggplot(aes(x = reorder(gene, -percentage),
             y = percentage,
             label = percentage)) +
  geom_label() +
  theme_classic() +
  ggtitle('Select genes - neurons only')

ggsave('neurons/neuroendocrine_genes/pctNeurons_exprSelectGenes_all.png',
       height = 5, width = 7)

# subset genes
subset_genes <- c("Esr1",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2")
# format data
sum_neur_indiv %>%
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "Freq"
  ) %>%
  dplyr::select(orig.ident,
                indiv_genotype,
                gene,
                Freq,
                cell_count) %>%
  rename(total.neuron.count = cell_count) %>% 
  mutate(Sex =
           ifelse(grepl("female", orig.ident),
                  "Female",
                  "Male")) %>%
  mutate(Status =
           ifelse(grepl("dom", orig.ident),
                  "Dom",
                  "Sub")) -> dat.for.stats


#### Neuroendocrine genes; neuron stats ####

# GLM on neuron proportions
neuron.genes.glm <- data.frame()
neuron.genes.pairwise <- data.frame()
neuron.genes.modglm <- data.frame()

for (gene in subset_genes) {
  
  # binomial GLM on sex and status w/ interaction term
  tmp = glm(cbind(Freq,
                  total.neuron.count - Freq) ~ Sex*Status,
            data = subset(dat.for.stats,
                          gene == gene),
            family = binomial)
  
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  data.frame(Estimate = tmp.res$test$coefficients,
             Std.error = tmp.res$test$sigma,
             z.value = tmp.res$test$tstat,
             p.value = tmp.res$test$pvalues) %>%
    rownames_to_column('Predictor') %>%
    mutate(gene = gene) -> tmp.res.df
  
  # FDR correction for all comparisons
  tmp.res.df %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(p.adjust.fdr =
             round(p.adjust(p.value,
                            method = 'fdr'),
                   digits = 4) ) -> neuron.genes.glm.fdr
  
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint) %>%
    mutate(Predictor = row.names(.)) %>%
    filter(!Predictor == "(Intercept)")
  neuron.genes.glm.fdr %>%
    full_join(ci_df, by = c("Predictor",
                            "Estimate")) -> neuron.genes.glm.fdr
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
    ggtitle(paste0(gene," neurons: binomial GLM"))
  
  ggsave(paste0('neurons/neuroendocrine_genes/selectGenesGLM/graphs', gene, '_binomGLM_cell.count.png'),
         width = 6, height = 5)
  
  # Average marginal effects (AMEs)
  ame <- slopes(tmp) %>%
    as.data.frame() %>%
    mutate(neuron.genes = gene,
           EffectType = "AME")
  
  # Marginal effects at representative values (MERs)
  mer <- slopes(tmp, type = "response") %>%
    as.data.frame() %>%
    mutate(neuron.genes = gene,
           EffectType = "MER") %>%
    dplyr::select(!c(rowid,
                     predicted_lo,
                     predicted_hi)) %>%
    distinct()
  
  merplots <- plot_mer_facet(mer, group = gene)
  
  ggsave(paste0('neurons/neuroendocrine_genes/selectGenesGLM/graphs/', gene, '_MER_cell.count.png'),
         plot = merplots, width = 12, height = 5)
  
  write_csv(mer,
            file = paste0('neurons/neuroendocrine_genes/selectGenesGLM/MER/', gene, '_MERtable.csv'))
  
  dat.for.mod <- subset(dat.for.stats,
                        gene == gene) %>%
    mutate(Sex = ifelse(Sex=="Female", 1, 0),
           Status = ifelse(Status=="Dom", 1, 0)) %>%
    dplyr::select(Freq, total.neuron.count, Sex, Status)
  
  tmp.for.mod <- glm(cbind(Freq,
                           total.neuron.count - Freq) ~ Sex*Status,
                     data = dat.for.mod,
                     family = binomial)
  
  # modglm interaction effects
  ints <- suppressWarnings(modglm(model = tmp.for.mod,
                                  vars = c("Sex", "Status"),
                                  data = dat.for.mod,
                                  hyps = "means",
                                  type = "dd"))
  
  # Average interaction effect (AIE)
  aie_df <- data.frame(neuron.genes = gene,
                       EffectType = "AIE",
                       term = "Sex:Status",
                       Estimate = ints$aie$aie.est,
                       Std.error = ints$aie$aie.se.delta,
                       lwr = ints$aie$aie.ll,
                       upr = ints$aie$aie.ll,
                       t.val = NA)
  
  # Interaction at hypothetical values
  inthyp_df <- data.frame(
    neuron.genes = gene,
    EffectType = "HypotheticalInt",
    term = "Sex:Status",
    Estimate = ints$inthyp$hat,
    Std.error = ints$inthyp$se.int.est,
    lwr = ints$inthyp$hat - ints$inthyp$se.int.est,
    upr = ints$inthyp$hat + ints$inthyp$se.int.est,
    t.val = ints$inthyp$t.val)
  
  neuron.genes.mod <- rbind(aie_df, inthyp_df)
  
  # Tukey for all pairwise comparisons
  tmp2 = glm(cbind(Freq,
                   total.neuron.count - Freq) ~ orig.ident,
             data = subset(dat.for.stats,
                           gene == gene),
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
    mutate(neuron.genes = gene) -> emm.res.df
  
  # combine results and save
  neuron.genes.glm = rbind(neuron.genes.glm,
                              neuron.genes.glm.fdr)
  neuron.genes.pairwise = rbind(neuron.genes.pairwise,
                                   emm.res.df)
  neuron.genes.modglm = rbind(neuron.genes.modglm,
                                 neuron.genes.mod)
  
  rm(neuron.genes.mod, emm.res.df, neuron.genes.glm.fdr)
  cat("Processed gene", gene, "\n")
}


write_csv(neuron.genes.glm,
          file = 'neurons/neuroendocrine_genes/selectGenesGLM/selectGenes_neuronGLM.csv')
write_csv(neuron.genes.pairwise,
          file = 'neurons/neuroendocrine_genes/selectGenesGLM/selectGenes_neuron_pairwise.csv')
write_csv(neuron.genes.modglm,
          file = 'neurons/neuroendocrine_genes/selectGenesGLM/selectGenes_neuron_modglm.csv')


#### DGE with limma ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# # load if starting here: 
# load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))
# 
# # set idents
# Idents(object = int.ldfs) <- "parent_id.broad.prob"
# 
# # subset to neurons
# int.ldfs = subset(int.ldfs,
#                   idents = c("C7-2: GABA", "C7-1: GLU"))
# 
# # subset genes
# subset_genes <- c("Esr1",
#                   "Ar",
#                   "Pgr",
#                   "Nr3c1",
#                   "Nr3c2")

l.dfs <- subset_by_gene(int.ldfs, 
                        subset_genes, 
                        slot = "counts",
                        min_count = 2)


# raw counts and norm.counts of variable genes
dge_data <- prep.for.DGE(l.dfs, 
                         selection.method = "vst", 
                         SCTransformed = TRUE, 
                         assay = 'integrated')

# # pseudo-bulk expr by sample
# dge_data_bulk <- prep.for.DGE(l.dfs,
#                               selection.method = "vst",
#                               SCTransformed = TRUE,
#                               pseudo_bulk = TRUE,
#                               group.by = c("indiv_genotype", "orig.ident")


# get results and graph p-values
limma_results <- run_limmatrend(dge_data$results, "./neurons/neuroendocrine_genes/limma_trend")

  save(limma_results, file = "./neurons/neuroendocrine_genes/limma_trend/limma_results.rda")

  
#### RRHO/RedRibbon ####
# how to determine max log scale before graphing? 
  max.log.scale = 235;
# make sure these match upper/lower case
  sex = c("Female", "Male");
  status = c("Dom", "Sub");

rrho_results <- get.RRHO(limma_results,
                         group.by = sex,
                         compare.across = status,
                         new.max.log = max.log.scale,
                         outdir = "./neurons/neuroendocrine_genes/RRHO")

  save(rrho_results, file = "./neurons/neuroendocrine_genes/RRHO/rrho_results.rda")
