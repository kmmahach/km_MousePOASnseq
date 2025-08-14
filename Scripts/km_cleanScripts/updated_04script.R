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

ggsave('neurons/neuropeptides/selectGenes_neurons.png',
       height = 4, width = 9)


# counts per sample
select_neur = full_join(full_join(int.ldfs@reductions$umap@cell.embeddings %>% 
                                    as.data.frame() %>% 
                                    rownames_to_column("Cell.id"),
                                  int.ldfs@meta.data %>%
                                    rownames_to_column("Cell.id")),
                        int.ldfs@assays$SCT@data %>% 
                          as.data.frame() %>% 
                          filter(rownames(int.ldfs@assays$SCT@data) %in% select_genes) %>% 
                          t() %>% as.data.frame() %>% 
                          rownames_to_column('Cell.id'))


# graph counts
# only want to count cells with read above 2
select_neur.table = select_neur %>% 
  dplyr::select(c(select_genes)) %>% 
  as.matrix() %>% 
  pmin(log(2)) 


# set value to 1 if read above 2
select_neur.table[select_neur.table > 0.5] <- 1

# create count table
select_neur.table = select_neur.table %>% 
  as.data.frame() %>% 
  mutate(All = 1) %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename(Counts = '.') %>% 
  mutate(Total.neurons = max(Counts),
         Percentage = 100*Counts/Total.neurons)

# compare counts across samples
# set expression to present/absent
tmp = int.ldfs@assays$SCT@data %>% 
  as.data.frame() %>% 
  filter(rownames(int.ldfs@assays$SCT@data) %in% select_genes) %>% 
  t() %>% 
  pmin(log(2)) 

tmp[tmp > 0.5] <- 1

tmp = tmp %>% 
  as.data.frame() %>% 
  rownames_to_column('Cell.id')

full_join(full_join(int.ldfs@reductions$umap@cell.embeddings %>% 
                      as.data.frame() %>% 
                      rownames_to_column("Cell.id"),
                    int.ldfs@meta.data %>%
                      rownames_to_column("Cell.id")), tmp) -> select_neur.table_orig.ident
rm(tmp)

#### Neuroendocrine genes; neuron stats ####

# subset genes by sample
select_neur.table_orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident,
           indiv_genotype) %>% 
  summarise_at(c(subset_genes, 'All'), 
               sum) %>% 
  pivot_longer(cols = c(subset_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident,
           indiv_genotype) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count, 2)) %>% 
  filter(Genes != 'All') -> select_neur.counts

# add variable for status and sex
select_neur.counts %>% 
  mutate(Sex = ifelse(grepl("female", orig.ident),
                      "female", "male")) %>% 
  mutate(Status = ifelse(grepl("dom", orig.ident),
                         "dom", "sub")) -> select_neur.counts

# graph percentage neuroendocrine per sample
select_neur.counts %>% 
  ggplot(aes(x = reorder(Genes, -Percentage),
             y = Percentage,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 3) +
  theme_classic() + 
  scale_color_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF")) +
  xlab('') + ylab('% of neurons') +
  ggtitle('Percent nuclei expressing neur-endo genes')

ggsave('neurons/neuropeptides/pctNeurons_exprSelectGenes_by.indiv_genotype.png',
       height = 5, width = 5)

# GLM (binom) on cluster proportions
neuron.genes.glm = data.frame(matrix(ncol = 7, nrow = 0))

colnames(neuron.genes.glm) <- c("contrast",
                                "Estimate",
                                "Std.error",
                                "z.value",
                                "adj.pvalue",
                                "z.ratio",
                                "Genes")

# loop through each cluster to generate stats
for (i in unique(select_neur.counts$Genes)) {
  
  # binomial GLM on sex and status w/ interaction term
  tmp = glm(cbind(Counts, Others) ~ Sex*Status, 
            data = select_neur.counts %>% 
              filter(Genes == i) %>% 
              mutate(Others = Total.count - Counts,
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
    mutate(Genes = i)
  
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
    ggtitle(paste0("Genes ", i,": binomial GLM cell count"))
  
  ggsave(paste0('neurons/neuropeptides/', i, '_binomialGLM_cell.count.png'),
         width = 6, height = 5)
  
  # Tukey for all pairwise comparisons
  tmp2 = glm(cbind(Counts, Others) ~ orig.ident, 
             data = select_neur.counts %>% 
               filter(Genes == i) %>% 
               mutate(Others = Total.count - Counts,
                      Sex = as.factor(Sex),
                      Status = as.factor(Status)), 
             family = binomial)
  
  # use emmeans to get pairwise differences
  emm = emmeans(tmp2, ~ "orig.ident")
  emm.res = pairs(emm)
  
  
  # save results
  emm.res.df = data.frame(Estimate = summary(emm.res)$estimate,
                          Std.error = summary(emm.res)$SE,
                          z.ratio = summary(emm.res)$z.ratio,
                          adj.pvalue = summary(emm.res)$p.value,
                          contrast = summary(emm.res)$contrast,
                          z.value = NA) %>%
    mutate(Genes = i)
  
  # combine results
  neuron.genes.glm = neuron.genes.glm %>% 
    rbind(tmp.res.df) %>% 
    rbind(emm.res.df)
}


# FDR correction for all comparisons
neuron.genes.glm.fdr = neuron.genes.glm %>%
  filter(is.na(z.ratio)) %>%
  filter(contrast != "(Intercept)") %>%
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                                 method = 'fdr'))

# combine with data frame
neuron.genes.glm = neuron.genes.glm %>%
  left_join(neuron.genes.glm.fdr)


# round p-vals to .0001
neuron.genes.glm = neuron.genes.glm %>%
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

write_csv(neuron.genes.glm,
          file = 'neurons/stats/neuron_selectGenesGLM.csv')

#### graph select neur-endo genes ####

select_neur.table %>% 
  rownames_to_column('Genes') %>% 
  filter(Percentage < 99) %>% 
  mutate(Percentage = signif(Percentage, 2)) %>% 
  ggplot(aes(x = reorder(Genes, -Percentage),
             y = Percentage,
             label = Percentage)) +
  geom_label() +
  theme_classic() +
  ggtitle('Select genes - neurons only')

ggsave('neurons/neuropeptides/pctNeurons_exprSelectGenes_all.png',
       height = 5, width = 7)


# graph table 
select_neur.table_orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident) %>% 
  summarise_at(c(select_genes, 'All'),
               sum) %>% 
  pivot_longer(cols = c(select_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count, 2)) %>% 
  filter(Genes != 'All') %>% 
  ggplot(aes(x = reorder(Genes, -Percentage),
             y= Percentage,
             color = orig.ident)) +
  geom_point() + 
  theme_classic()

ggsave('neurons/neuropeptides/pctNeurons_exprSelectGenes_by_orig.ident.png',
       height = 5, width = 9)

# reduced gene subset
select_neur.table_orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident) %>% 
  summarise_at(c(select_genes, 'All'),
               sum) %>% 
  pivot_longer(cols = c(select_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count, 2)) %>% 
  filter(Genes != 'All') %>% 
  filter(Percentage >= 5) %>% 
  ggplot(aes(x = reorder(Genes, -Percentage),
             y = Percentage,
             color = orig.ident)) +
  geom_point() +
  theme_classic()

ggsave('neurons/neuropeptides/pctNeurons_exprSelectGenes_by_orig.ident_filtered.png',
       height = 5, width = 9)

# reduced subset again
select_neur.table_orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident) %>% 
  summarise_at(c(subset_genes,
                 'All'),
               sum) %>% 
  pivot_longer(cols = c(subset_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count,
                            2)) %>% 
  filter(Genes != 'All') %>% 
  filter(Percentage >= 5) %>% 
  ggplot(aes(x = reorder(Genes,
                         -Percentage),
             y= Percentage,
             color = orig.ident)) +
  geom_point() +
  theme_classic()

ggsave('neurons/neuropeptides/pctNeurons_exprSubset_SelectGenes_by_orig.ident_filtered.png',
       height = 5, width = 9)


#### DGE with limma ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# load if starting here: 
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# subset genes
subset_genes <- c("Esr1",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2")

l.dfs <- subset_by_gene(int.ldfs, subset_genes, min_count = 2)


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
# may want to filter to top 500 most variable genes? 
limma_results <- run_limmatrend(dge_data$results, "./neurons/neuropeptides/limma_trend")

  save(limma_results, file = "./neurons/neuropeptides/limma_trend/limma_results.rda")

  
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
                         outdir = "./neurons/neuropeptides/RRHO")

  save(rrho_results, file = "./neurons/neuropeptides/RRHO/rrho_results.rda")
