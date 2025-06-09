#### Mouse snseq seurat analysis
### DEG analysis
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3

### set working directory
setwd("/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/")

#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")
# install.packages("patchwork",
#                  repos = "https://cloud.r-project.org")
#install.packages('BiocManager')
#BiocManager::install('multtest')
# install.packages('metap',
#                  repos = "https://cloud.r-project.org")
# install.packages('clustree')
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DEsingle")
# install.packages("ggalluvial")
# library(devtools)
# install_github("RRHO2/RRHO2")

#load libraries
library(Seurat)
library(patchwork)
library(clustree)
library(pheatmap)
library(DEsingle)
library(dendextend)
library(tidyverse)
library(RRHO2)
library(ggalluvial)
library(multcomp)
library(emmeans)
library(RedRibbon)

#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#### load data ####
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

## select neurons
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.broad.prob"

#subset to neurons
mouse.snseq.combined.sct.neurons = subset(mouse.snseq.combined.sct,
                                          idents = c("C7-2: GABA",
                                                     "C7-1: GLU"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.neurons) = 'SCT'

### load neuropeptide list 
neuropeptides = read.csv('./gene.lists/neuropeptides.list.csv')
# genes
neuropeptides.genes = neuropeptides %>% 
  pull(Gene.name)


##### select neuropeptide analysis ####
### create dotplot of select genes
## select genes
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

## subset genes
subset_genes <- c("Esr1",
                  "Ar",
                  "Pgr",
                  "Nr3c1",
                  "Nr3c2")

# ## graph
# DotPlot(mouse.snseq.combined.sct.neurons,
#         assay = 'SCT',
#         features = select_genes,
#         group.by = 'orig.ident',
#         col.min = 0)
# ggsave('neurons/neuropeptides/Select genes neurons dotplot.png',
#        height = 10,
#        width = 10)
# 
# VlnPlot(mouse.snseq.combined.sct.neurons,
#         assay = 'SCT',
#         features = select_genes,
#         group.by = 'orig.ident',
#         slot = 'counts') +
#   ylim(2,20)

## counts per sample
mouse.snseq.combined.sct.neurons.expression.select = full_join(full_join(mouse.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
                                                                                     as.data.frame() %>% 
                                                                                     rownames_to_column("Cell.id"),
                                                                                   mouse.snseq.combined.sct.neurons@meta.data %>%
                                                                                     rownames_to_column("Cell.id")),
                                                                         mouse.snseq.combined.sct.neurons@assays$SCT@data %>% 
                                                                           as.data.frame() %>% 
                                                                           filter(rownames(mouse.snseq.combined.sct.neurons@assays$SCT@data) %in% select_genes) %>% 
                                                                           t() %>% as.data.frame() %>% 
                                                                           rownames_to_column('Cell.id'))


## graph counts
# only want to count cells with read above 2
mouse.snseq.combined.sct.neurons.expression.select.table = mouse.snseq.combined.sct.neurons.expression.select %>% 
  dplyr::select(c(select_genes)) %>% 
  as.matrix() %>% 
  pmin(log(2)) 
# set value to 1 if read above 2
mouse.snseq.combined.sct.neurons.expression.select.table[mouse.snseq.combined.sct.neurons.expression.select.table > 0.5] <-1
# create count table
mouse.snseq.combined.sct.neurons.expression.select.table = mouse.snseq.combined.sct.neurons.expression.select.table %>% 
  as.data.frame() %>% 
  mutate(All = 1) %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename(Counts = '.') %>% 
  mutate(Total.neurons = max(Counts),
         Percentage = 100*Counts/Total.neurons)


### compare counts across samples
## set expression to present/absent
tmp = mouse.snseq.combined.sct.neurons@assays$SCT@data %>% 
  as.data.frame() %>% 
  filter(rownames(mouse.snseq.combined.sct.neurons@assays$SCT@data) %in% select_genes) %>% 
  t() %>% 
  pmin(log(2)) 
tmp[tmp > 0.5] <-1
tmp = tmp %>% 
  as.data.frame() %>% 
  rownames_to_column('Cell.id')

mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident = full_join(full_join(mouse.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
                                                                                            as.data.frame() %>% 
                                                                                            rownames_to_column("Cell.id"),
                                                                                          mouse.snseq.combined.sct.neurons@meta.data %>%
                                                                                            rownames_to_column("Cell.id")),
                                                                                tmp)
rm(tmp)


#### neuroendocrine per sample count stats ####
## subset genes by sample
mouse.snseq.neuroendocrine.counts = mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident,
           Genotype) %>% 
  summarise_at(c(subset_genes,
                 'All'),
               sum) %>% 
  pivot_longer(cols = c(subset_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident,
           Genotype) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count,
                            2)) %>% 
  filter(Genes != 'All') 

# add variable for status and sex
mouse.snseq.neuroendocrine.counts = mouse.snseq.neuroendocrine.counts %>% 
  mutate(Sex = ifelse(grepl("Female",
                            orig.ident),
                      "Female",
                      "Male")) %>% 
  mutate(Status = ifelse(grepl("Dom",
                               orig.ident),
                         "Dom",
                         "Sub"))

# graph percentage neuroendocrine per sample
mouse.snseq.neuroendocrine.counts %>% 
  ggplot(aes(x = reorder(Genes,
                         -Percentage),
             y= Percentage,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = Genotype,
                 group = orig.ident),
             size = 3)+
  # facet_grid(~Genes,
  #            scales = "free_x") +
  theme_classic()+ 
  scale_color_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF")) +
  xlab('') +
  ggtitle('Percentage neurons expressing neuroendocrine genes')
ggsave('neurons/neuropeptides/Subset genes neurons percentage orig.ident genotype.png',
       height = 5,
       width = 5)


#### GLM binomial
### run glm on every cluster
## use with proportion data
# create empty data frame
neuron.genes.glm = data.frame(matrix(ncol = 7, nrow = 0))

#provide column names
colnames(neuron.genes.glm) <- c("contrast",
                                  "Estimate",
                                  "Std.error",
                                  "z.value",
                                  "adj.pvalue",
                                  "z.ratio",
                                  "Genes")
## loop through each cluster
for (i in unique(mouse.snseq.neuroendocrine.counts$Genes)) {
  ## run interaction binomial GLM on sex and status
  tmp = glm(cbind(Counts, Others) ~ Sex*Status, 
            data = mouse.snseq.neuroendocrine.counts %>% 
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
    ggtitle(paste0("Genes ",
                   i,
                   ": binomial GLM cell count"))
  ggsave(paste0('neurons/neuron_Genes/',
                i,
                ' binomial GLM cell count.png'))
  
  ## run all pairwise comparison
  # Tukey needed for all pairwise comparisons
  tmp2 = glm(cbind(Counts, Others) ~ orig.ident, 
             data = mouse.snseq.neuroendocrine.counts %>% 
               filter(Genes == i) %>% 
               mutate(Others = Total.count - Counts,
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
    mutate(Genes = i)
  
  ## combine in data frame
  neuron.genes.glm = neuron.genes.glm %>% 
    rbind(tmp.res.df) %>% 
    rbind(emm.res.df)
}


## FDR correction for all comparisons
neuron.genes.glm.fdr = neuron.genes.glm %>%
  filter(is.na(z.ratio)) %>%
  filter(contrast != "(Intercept)") %>%
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                                 method = 'fdr'))

# combine with data frame
neuron.genes.glm = neuron.genes.glm %>%
  left_join(neuron.genes.glm.fdr)



## create rounded p.value to make it easier to read
neuron.genes.glm = neuron.genes.glm %>%
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

## save glm table
write_csv(neuron.genes.glm,
          file = 'neurons/neuron_seurat_genes_glm.csv')

## create p value dataframe
neuron.genes.glm.p = neuron.genes.glm %>%
  dplyr::select(Genes,
                adj.pvalue.round,
                contrast) %>%
  pivot_wider(names_from = contrast,
              values_from = adj.pvalue.round)











#### graph select ####
# ggplot
mouse.snseq.combined.sct.neurons.expression.select.table %>% 
  rownames_to_column('Genes') %>% 
  filter(Percentage < 99) %>% 
  mutate(Percentage = signif(Percentage, 2)) %>% 
  ggplot(aes(x = reorder(Genes,
                         -Percentage),
             y = Percentage,
             label = Percentage)) +
  geom_label() +
  theme_classic() +
  ggtitle('Select genes neurons')
ggsave('neurons/neuropeptides/Select genes neurons percentage.png',
       height = 10,
       width = 10)


## graph table 
mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident) %>% 
  summarise_at(c(select_genes,
                 'All'),
               sum) %>% 
  pivot_longer(cols = c(select_genes, 'All'),
               names_to = 'Genes',
               values_to = 'Counts') %>% 
  group_by(orig.ident) %>% 
  mutate(Total.count = max(Counts)) %>% 
  ungroup() %>% 
  mutate(Percentage = round(100*Counts/Total.count,
                            2)) %>% 
  filter(Genes != 'All') %>% 
  ggplot(aes(x = reorder(Genes,
                         -Percentage),
             y= Percentage,
             color = orig.ident)) +
  geom_point()
ggsave('neurons/neuropeptides/Select genes neurons percentage orig.ident.png',
       height = 10,
       width = 10)

## graph table 
mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident %>% 
  mutate(All = 1) %>% 
  group_by(orig.ident) %>% 
  summarise_at(c(select_genes,
                 'All'),
               sum) %>% 
  pivot_longer(cols = c(select_genes, 'All'),
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
  geom_point()
ggsave('neurons/neuropeptides/Select genes neurons percentage orig.ident filter.png',
       height = 10,
       width = 10)

# subset genes
mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident %>% 
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
ggsave('neurons/neuropeptides/Subset genes neurons percentage orig.ident filter.png',
       height = 10,
       width = 10)


# 
# ## graph table 
# mouse.snseq.combined.sct.neurons.expression.select.table.orig.ident %>% 
#   mutate(All = 1) %>% 
#   group_by(orig.ident) %>% 
#   summarise_at(c(select_genes,
#                  'All'),
#                sum) %>% 
#   pivot_longer(cols = c(select_genes, 'All'),
#                names_to = 'Genes',
#                values_to = 'Counts') %>% 
#   group_by(orig.ident) %>% 
#   mutate(Total.count = max(Counts)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = round(100*Counts/Total.count,
#                             2)) %>% 
#   filter(Genes == 'All') %>% 
#   ggplot(aes(x = reorder(orig.ident,
#                          -Counts),
#              y= Counts,
#              color = orig.ident)) +
#   geom_point()
# ggsave('neurons/neuropeptides/Neurons counts orig.ident.png',
#        height = 10,
#        width = 10)



# ### neuropeptides
# ## counts per sample
# mouse.snseq.combined.sct.neurons.expression.np = full_join(full_join(mouse.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
#                                                                            as.data.frame() %>% 
#                                                                            rownames_to_column("Cell.id"),
#                                                                          mouse.snseq.combined.sct.neurons@meta.data %>%
#                                                                            rownames_to_column("Cell.id")),
#                                                                mouse.snseq.combined.sct.neurons@assays$SCT@data %>% 
#                                                                  as.data.frame() %>% 
#                                                                  filter(toupper(rownames(mouse.snseq.combined.sct.neurons@assays$SCT@data)) %in% neuropeptides.genes) %>% 
#                                                                  t() %>% as.data.frame() %>% 
#                                                              rename_all(., .funs = toupper) %>% 
#                                                                  rownames_to_column('Cell.id'))
# 
# ## graph counts
# # only want to count cells with read above 2
# mouse.snseq.combined.sct.neurons.expression.np.table = mouse.snseq.combined.sct.neurons.expression.np %>% 
#   dplyr::select(c(intersect(neuropeptides.genes, colnames(mouse.snseq.combined.sct.neurons.expression.np)))) %>% 
#   as.matrix() %>% 
#   pmin(log(2)) 
# # set value to 1 if read above 2
# mouse.snseq.combined.sct.neurons.expression.np.table[mouse.snseq.combined.sct.neurons.expression.np.table > 0.5] <-1
# # create count table
# mouse.snseq.combined.sct.neurons.expression.np.table = mouse.snseq.combined.sct.neurons.expression.np.table %>% 
#   as.data.frame() %>% 
#   mutate(All = 1) %>% 
#   colSums() %>% 
#   as.data.frame() %>% 
#   dplyr::rename(Counts = '.') %>% 
#   mutate(Total.neurons = max(Counts),
#          Percentage = 100*Counts/Total.neurons)


#### Ar limmatrend ####
### subset into Ar

# subset with count above 2 reads 
mouse.snseq.combined.sct.neurons.Ar = subset(mouse.snseq.combined.sct.neurons,
                                                    subset = Ar >= 2,
                                             slot = 'counts')


### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.Ar <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.Ar, 
                                                                                              assay = 'integrated',
                                                                                              selection.method = "vst", 
                                                                                              verbose = F)
# identify top variable genes
Ar.neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.Ar,
                                                               assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.Ar.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.Ar@reductions$umap@cell.embeddings %>% 
                                                                                           as.data.frame() %>% 
                                                                                           rownames_to_column("Cell.id"),
                                                                     mouse.snseq.combined.sct.neurons.Ar@meta.data %>%
                                                                                           rownames_to_column("Cell.id")),
                                                           mouse.snseq.combined.sct.neurons.Ar@assays$SCT@data %>% 
                                                                                 as.data.frame() %>% 
                                                                                 filter(rownames(mouse.snseq.combined.sct.neurons.Ar@assays$SCT@data) %in% c('Ar')) %>% 
                                                                                 t() %>% as.data.frame() %>% 
                                                                                 rownames_to_column('Cell.id'))


## create vector of factor
Ar.neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.Ar.expression %>% 
  mutate(Ar.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Ar.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Ar.neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons.Ar,
                                                                assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.Ar.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Ar.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Ar.neuron.reduce.vector.limma.sct.prep = list(count = Ar.neuron.reduce.vector.count.sct.prep,
                                               condt = Ar.neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Ar <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
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
  })
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Ar/limmatrend/limmatrend.histograms.pdf" )
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
  dev.off()
  
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Ar/limmatrend/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(1,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF)
}


###run function  
Ar.neuron.limma.results.sct = run_limmatrend_Ar(Ar.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
Ar.neuron.limma.results.sct.df = full_join(Ar.neuron.limma.results.sct$ttDvsS %>% 
                                                                                rename_with(~paste0(.,"_DvsS")) %>% 
                                                                                rownames_to_column("Gene"),
                                                                              Ar.neuron.limma.results.sct$ttDvsSM %>% 
                                                                                rename_with(~paste0(.,"_DvsS_M")) %>% 
                                                                                rownames_to_column("Gene")) %>% 
  full_join(Ar.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Ar.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Ar.neuron.limma.results.sct.df = Ar.neuron.limma.results.sct.df %>% 
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
                                       Sig_MvsF)) 

##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  Ar.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/Ar/limmatrend/limma.Ar.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# save results
save(Ar.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Ar/limmatrend/Ar.neuron.limma.results.sct.df.RData')

load('./neurons/neuropeptides/Ar/limmatrend/Ar.neuron.limma.results.sct.df.RData')

write.csv(Ar.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Ar/limmatrend/Ar.neuron.limma.results.sct.df.csv')




#### Pgr limmatrend ####
### subset into Pgr

# subset with count above 2 reads 
mouse.snseq.combined.sct.neurons.Pgr = subset(mouse.snseq.combined.sct.neurons,
                                             subset = Pgr >= 2,
                                             slot = 'counts')


### calculate variable genes
## identify top variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.Pgr <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.Pgr, 
                                                            assay = 'integrated',
                                                            selection.method = "vst", 
                                                            verbose = F)
# identify top  variable genes
Pgr.neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.Pgr,
                                                              assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.Pgr.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.Pgr@reductions$umap@cell.embeddings %>% 
                                                                       as.data.frame() %>% 
                                                                       rownames_to_column("Cell.id"),
                                                                     mouse.snseq.combined.sct.neurons.Pgr@meta.data %>%
                                                                       rownames_to_column("Cell.id")),
                                                           mouse.snseq.combined.sct.neurons.Pgr@assays$SCT@data %>% 
                                                             as.data.frame() %>% 
                                                             filter(rownames(mouse.snseq.combined.sct.neurons.Pgr@assays$SCT@data) %in% c('Pgr')) %>% 
                                                             t() %>% as.data.frame() %>% 
                                                             rownames_to_column('Cell.id'))


## create vector of factor
Pgr.neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.Pgr.expression %>% 
  mutate(Pgr.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Pgr.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Pgr.neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons.Pgr,
                                                      assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.Pgr.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Pgr.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Pgr.neuron.reduce.vector.limma.sct.prep = list(count = Pgr.neuron.reduce.vector.count.sct.prep,
                                              condt = Pgr.neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Pgr <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
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
  })
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Pgr/limmatrend/limmatrend.histograms.pdf" )
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
  dev.off()
  
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Pgr/limmatrend/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(1,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF)
}


###run function  
Pgr.neuron.limma.results.sct = run_limmatrend_Pgr(Pgr.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
Pgr.neuron.limma.results.sct.df = full_join(Pgr.neuron.limma.results.sct$ttDvsS %>% 
                                             rename_with(~paste0(.,"_DvsS")) %>% 
                                             rownames_to_column("Gene"),
                                           Pgr.neuron.limma.results.sct$ttDvsSM %>% 
                                             rename_with(~paste0(.,"_DvsS_M")) %>% 
                                             rownames_to_column("Gene")) %>% 
  full_join(Pgr.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Pgr.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Pgr.neuron.limma.results.sct.df = Pgr.neuron.limma.results.sct.df %>% 
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
                                     Sig_MvsF)) 

##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  Pgr.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/Pgr/limmatrend/limma.Pgr.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# save results
save(Pgr.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Pgr/limmatrend/Pgr.neuron.limma.results.sct.df.RData')

load('./neurons/neuropeptides/Pgr/limmatrend/Pgr.neuron.limma.results.sct.df.RData')

write.csv(Pgr.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Pgr/limmatrend/Pgr.neuron.limma.results.sct.df.csv')





#### Esr1 limmatrend ####
### subset into Esr1

# subset with count above 2 reads 
mouse.snseq.combined.sct.neurons.Esr1 = subset(mouse.snseq.combined.sct.neurons,
                                              subset = Esr1 >= 2,
                                              slot = 'counts')


### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.Esr1 <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.Esr1, 
                                                             assay = 'integrated',
                                                             selection.method = "vst", 
                                                             verbose = F)
# identify top variable genes
Esr1.neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.Esr1,
                                                               assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.Esr1.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.Esr1@reductions$umap@cell.embeddings %>% 
                                                                        as.data.frame() %>% 
                                                                        rownames_to_column("Cell.id"),
                                                                      mouse.snseq.combined.sct.neurons.Esr1@meta.data %>%
                                                                        rownames_to_column("Cell.id")),
                                                            mouse.snseq.combined.sct.neurons.Esr1@assays$SCT@data %>% 
                                                              as.data.frame() %>% 
                                                              filter(rownames(mouse.snseq.combined.sct.neurons.Esr1@assays$SCT@data) %in% c('Esr1')) %>% 
                                                              t() %>% as.data.frame() %>% 
                                                              rownames_to_column('Cell.id'))


## create vector of factor
Esr1.neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.Esr1.expression %>% 
  mutate(Esr1.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Esr1.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Esr1.neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons.Esr1,
                                                       assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.Esr1.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Esr1.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Esr1.neuron.reduce.vector.limma.sct.prep = list(count = Esr1.neuron.reduce.vector.count.sct.prep,
                                               condt = Esr1.neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Esr1 <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
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
  })
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Esr1/limmatrend/limmatrend.histograms.pdf" )
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
  dev.off()
  
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Esr1/limmatrend/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(1,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF)
}


###run function  
Esr1.neuron.limma.results.sct = run_limmatrend_Esr1(Esr1.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
Esr1.neuron.limma.results.sct.df = full_join(Esr1.neuron.limma.results.sct$ttDvsS %>% 
                                              rename_with(~paste0(.,"_DvsS")) %>% 
                                              rownames_to_column("Gene"),
                                            Esr1.neuron.limma.results.sct$ttDvsSM %>% 
                                              rename_with(~paste0(.,"_DvsS_M")) %>% 
                                              rownames_to_column("Gene")) %>% 
  full_join(Esr1.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Esr1.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Esr1.neuron.limma.results.sct.df = Esr1.neuron.limma.results.sct.df %>% 
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
                                     Sig_MvsF)) 

##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  Esr1.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/Esr1/limmatrend/limma.Esr1.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# save results
save(Esr1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Esr1/limmatrend/Esr1.neuron.limma.results.sct.df.RData')

load('./neurons/neuropeptides/Esr1/limmatrend/Esr1.neuron.limma.results.sct.df.RData')

write.csv(Esr1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Esr1/limmatrend/Esr1.neuron.limma.results.sct.df.csv')





#### Nr3c1 limmatrend ####
### subset into Nr3c1

# subset with count above 2 reads 
mouse.snseq.combined.sct.neurons.Nr3c1 = subset(mouse.snseq.combined.sct.neurons,
                                               subset = Nr3c1 >= 2,
                                               slot = 'counts')


### calculate variable genes
## identify top variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.Nr3c1 <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.Nr3c1, 
                                                              assay = 'integrated',
                                                              selection.method = "vst", 
                                                              verbose = F)
# identify top variable genes
Nr3c1.neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.Nr3c1,
                                                                assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.Nr3c1.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.Nr3c1@reductions$umap@cell.embeddings %>% 
                                                                         as.data.frame() %>% 
                                                                         rownames_to_column("Cell.id"),
                                                                       mouse.snseq.combined.sct.neurons.Nr3c1@meta.data %>%
                                                                         rownames_to_column("Cell.id")),
                                                             mouse.snseq.combined.sct.neurons.Nr3c1@assays$SCT@data %>% 
                                                               as.data.frame() %>% 
                                                               filter(rownames(mouse.snseq.combined.sct.neurons.Nr3c1@assays$SCT@data) %in% c('Nr3c1')) %>% 
                                                               t() %>% as.data.frame() %>% 
                                                               rownames_to_column('Cell.id'))


## create vector of factor
Nr3c1.neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.Nr3c1.expression %>% 
  mutate(Nr3c1.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Nr3c1.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Nr3c1.neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons.Nr3c1,
                                                        assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.Nr3c1.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Nr3c1.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Nr3c1.neuron.reduce.vector.limma.sct.prep = list(count = Nr3c1.neuron.reduce.vector.count.sct.prep,
                                                condt = Nr3c1.neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Nr3c1 <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
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
  })
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Nr3c1/limmatrend/limmatrend.histograms.pdf" )
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
  dev.off()
  
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Nr3c1/limmatrend/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(1,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF)
}


###run function  
Nr3c1.neuron.limma.results.sct = run_limmatrend_Nr3c1(Nr3c1.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
Nr3c1.neuron.limma.results.sct.df = full_join(Nr3c1.neuron.limma.results.sct$ttDvsS %>% 
                                               rename_with(~paste0(.,"_DvsS")) %>% 
                                               rownames_to_column("Gene"),
                                             Nr3c1.neuron.limma.results.sct$ttDvsSM %>% 
                                               rename_with(~paste0(.,"_DvsS_M")) %>% 
                                               rownames_to_column("Gene")) %>% 
  full_join(Nr3c1.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Nr3c1.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Nr3c1.neuron.limma.results.sct.df = Nr3c1.neuron.limma.results.sct.df %>% 
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
                                     Sig_MvsF)) 

##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  Nr3c1.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/Nr3c1/limmatrend/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# save results
save(Nr3c1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c1/limmatrend/Nr3c1.neuron.limma.results.sct.df.RData')

load('./neurons/neuropeptides/Nr3c1/limmatrend/Nr3c1.neuron.limma.results.sct.df.RData')

write.csv(Nr3c1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c1/limmatrend/Nr3c1.neuron.limma.results.sct.df.csv')





#### Nr3c2 limmatrend ####
### subset into Nr3c2

# subset with count above 2 reads 
mouse.snseq.combined.sct.neurons.Nr3c2 = subset(mouse.snseq.combined.sct.neurons,
                                               subset = Nr3c2 >= 2,
                                               slot = 'counts')


### calculate variable genes
## identify top  variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.Nr3c2 <- FindVariableFeatures(mouse.snseq.combined.sct.neurons.Nr3c2, 
                                                              assay = 'integrated',
                                                              selection.method = "vst", 
                                                              verbose = F)
# identify top variable genes
Nr3c2.neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.Nr3c2,
                                                                assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.Nr3c2.expression = full_join(full_join(mouse.snseq.combined.sct.neurons.Nr3c2@reductions$umap@cell.embeddings %>% 
                                                                         as.data.frame() %>% 
                                                                         rownames_to_column("Cell.id"),
                                                                       mouse.snseq.combined.sct.neurons.Nr3c2@meta.data %>%
                                                                         rownames_to_column("Cell.id")),
                                                             mouse.snseq.combined.sct.neurons.Nr3c2@assays$SCT@data %>% 
                                                               as.data.frame() %>% 
                                                               filter(rownames(mouse.snseq.combined.sct.neurons.Nr3c2@assays$SCT@data) %in% c('Nr3c2')) %>% 
                                                               t() %>% as.data.frame() %>% 
                                                               rownames_to_column('Cell.id'))


## create vector of factor
Nr3c2.neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.Nr3c2.expression %>% 
  mutate(Nr3c2.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Nr3c2.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Nr3c2.neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons.Nr3c2,
                                                        assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.Nr3c2.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Nr3c2.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Nr3c2.neuron.reduce.vector.limma.sct.prep = list(count = Nr3c2.neuron.reduce.vector.count.sct.prep,
                                                condt = Nr3c2.neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Nr3c2 <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS = (treatMale.Dom + treatFemale.Dom)/2 - (treatMale.Sub + treatFemale.Sub)/2, 
                               DvsSM = treatMale.Dom - treatMale.Sub, 
                               DvsSF = treatFemale.Dom - treatFemale.Sub,
                               MvsF = (treatMale.Dom + treatMale.Sub)/2 - (treatFemale.Dom + treatFemale.Sub)/2,
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
  })
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Nr3c2/limmatrend/limmatrend.histograms.pdf" )
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
  dev.off()
  
  # Open pdf file
  pdf(file= "./neurons/neuropeptides/Nr3c2/limmatrend/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(1,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDvsS = ttDvsS,
       ttDvsSM = ttDvsSM,
       ttDvsSF = ttDvsSF,
       ttMvsF = ttMvsF)
}


###run function  
Nr3c2.neuron.limma.results.sct = run_limmatrend_Nr3c2(Nr3c2.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
Nr3c2.neuron.limma.results.sct.df = full_join(Nr3c2.neuron.limma.results.sct$ttDvsS %>% 
                                               rename_with(~paste0(.,"_DvsS")) %>% 
                                               rownames_to_column("Gene"),
                                             Nr3c2.neuron.limma.results.sct$ttDvsSM %>% 
                                               rename_with(~paste0(.,"_DvsS_M")) %>% 
                                               rownames_to_column("Gene")) %>% 
  full_join(Nr3c2.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Nr3c2.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Nr3c2.neuron.limma.results.sct.df = Nr3c2.neuron.limma.results.sct.df %>% 
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
                                     Sig_MvsF)) 

##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS",
                          "DvsS_M",
                          "DvsS_F",
                          "MvsF")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  Nr3c2.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/Nr3c2/limmatrend/limma.Nr3c2.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# save results
save(Nr3c2.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c2/limmatrend/Nr3c2.neuron.limma.results.sct.df.RData')

load('./neurons/neuropeptides/Nr3c2/limmatrend/Nr3c2.neuron.limma.results.sct.df.RData')

write.csv(Nr3c2.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c2/limmatrend/Nr3c2.neuron.limma.results.sct.df.csv')





#### save and load receptor limmatrend data ####
save(Ar.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Ar/Ar.neuron.limma.results.sct.df.Rdata')
save(Pgr.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Pgr/Pgr.neuron.limma.results.sct.df.Rdata')
save(Esr1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Esr1/Esr1.neuron.limma.results.sct.df.Rdata')
save(Nr3c1.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c1/Nr3c1.neuron.limma.results.sct.df.Rdata')
save(Nr3c2.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/Nr3c2/Nr3c2.neuron.limma.results.sct.df.Rdata')

# load('./neurons/neuropeptides/Ar/Ar.neuron.limma.results.sct.df.Rdata')
# load('./neurons/neuropeptides/Pgr/Pgr.neuron.limma.results.sct.df.Rdata')
# load('./neurons/neuropeptides/Esr1/Esr1.neuron.limma.results.sct.df.Rdata')
# load('./neurons/neuropeptides/Nr3c1/Nr3c1.neuron.limma.results.sct.df.Rdata')
# load('./neurons/neuropeptides/Nr3c2/Nr3c2.neuron.limma.results.sct.df.Rdata')

#### RRHO2 receptors ####
#### compare dom vs sub females across sexes
### Ar
##male data
Ar.neuron.DvsS.M.rrho2 = data.frame(gene = Ar.neuron.limma.results.sct.df$Gene,
                                    value = -log10(Ar.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                    direction = Ar.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                    stringsAsFactors = FALSE)
#set positive negative
Ar.neuron.DvsS.M.rrho2 = Ar.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Ar.neuron.DvsS.F.rrho2 = data.frame(gene = Ar.neuron.limma.results.sct.df$Gene,
                                    value = -log10(Ar.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                    direction = Ar.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                    stringsAsFactors = FALSE)
#set positive negative
Ar.neuron.DvsS.F.rrho2 = Ar.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Ar.neuron.RRHO_obj <-  RRHO2_initialize(Ar.neuron.DvsS.M.rrho2, 
                              Ar.neuron.DvsS.F.rrho2, 
                              boundary = 0.05,
                              labels = c('males',
                                         'females'))
## graph
png('./neurons/neuropeptides/Ar/limmatrend/RRHO2.Ar.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Ar.neuron.RRHO_obj,
              main = 'Ar DvsS')
dev.off()

#poster
Ar.neuron.RRHO_obj.poster <-  RRHO2_initialize(Ar.neuron.DvsS.M.rrho2, 
                                        Ar.neuron.DvsS.F.rrho2, 
                                        boundary = 0.05)
pdf('./neurons/neuropeptides/Ar/limmatrend/RRHO2.Ar.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(Ar.neuron.RRHO_obj.poster)
dev.off()

### Pgr
##male data
Pgr.neuron.DvsS.M.rrho2 = data.frame(gene = Pgr.neuron.limma.results.sct.df$Gene,
                                    value = -log10(Pgr.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                    direction = Pgr.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                    stringsAsFactors = FALSE)
#set positive negative
Pgr.neuron.DvsS.M.rrho2 = Pgr.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Pgr.neuron.DvsS.F.rrho2 = data.frame(gene = Pgr.neuron.limma.results.sct.df$Gene,
                                    value = -log10(Pgr.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                    direction = Pgr.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                    stringsAsFactors = FALSE)
#set positive negative
Pgr.neuron.DvsS.F.rrho2 = Pgr.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Pgr.neuron.RRHO_obj <-  RRHO2_initialize(Pgr.neuron.DvsS.M.rrho2, 
                                        Pgr.neuron.DvsS.F.rrho2, 
                                        boundary = 0.05,
                                        labels = c('males',
                                                   'females'))
## graph
png('./neurons/neuropeptides/Pgr/limmatrend/RRHO2.Pgr.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Pgr.neuron.RRHO_obj,
              main = 'Pgr DvsS')
dev.off()

#poster
Pgr.neuron.RRHO_obj.poster <-  RRHO2_initialize(Pgr.neuron.DvsS.M.rrho2, 
                                               Pgr.neuron.DvsS.F.rrho2, 
                                               boundary = 0.05)
pdf('./neurons/neuropeptides/Pgr/limmatrend/RRHO2.Pgr.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(Pgr.neuron.RRHO_obj.poster)
dev.off()

### Esr1
##male data
Esr1.neuron.DvsS.M.rrho2 = data.frame(gene = Esr1.neuron.limma.results.sct.df$Gene,
                                     value = -log10(Esr1.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = Esr1.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
Esr1.neuron.DvsS.M.rrho2 = Esr1.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Esr1.neuron.DvsS.F.rrho2 = data.frame(gene = Esr1.neuron.limma.results.sct.df$Gene,
                                     value = -log10(Esr1.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = Esr1.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
Esr1.neuron.DvsS.F.rrho2 = Esr1.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Esr1.neuron.RRHO_obj <-  RRHO2_initialize(Esr1.neuron.DvsS.M.rrho2, 
                                         Esr1.neuron.DvsS.F.rrho2, 
                                         boundary = 0.05,
                                         labels = c('males',
                                                    'females'))
## graph
png('./neurons/neuropeptides/Esr1/limmatrend/RRHO2.Esr1.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Esr1.neuron.RRHO_obj,
              main = 'Esr1 DvsS')
dev.off()

#poster
Esr1.neuron.RRHO_obj.poster <-  RRHO2_initialize(Esr1.neuron.DvsS.M.rrho2, 
                                                Esr1.neuron.DvsS.F.rrho2, 
                                                boundary = 0.05)
pdf('./neurons/neuropeptides/Esr1/limmatrend/RRHO2.Esr1.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(Esr1.neuron.RRHO_obj.poster)
dev.off()

### Nr3c1
##male data
Nr3c1.neuron.DvsS.M.rrho2 = data.frame(gene = Nr3c1.neuron.limma.results.sct.df$Gene,
                                      value = -log10(Nr3c1.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                      direction = Nr3c1.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                      stringsAsFactors = FALSE)
#set positive negative
Nr3c1.neuron.DvsS.M.rrho2 = Nr3c1.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Nr3c1.neuron.DvsS.F.rrho2 = data.frame(gene = Nr3c1.neuron.limma.results.sct.df$Gene,
                                      value = -log10(Nr3c1.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                      direction = Nr3c1.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                      stringsAsFactors = FALSE)
#set positive negative
Nr3c1.neuron.DvsS.F.rrho2 = Nr3c1.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Nr3c1.neuron.RRHO_obj <-  RRHO2_initialize(Nr3c1.neuron.DvsS.M.rrho2, 
                                          Nr3c1.neuron.DvsS.F.rrho2, 
                                          boundary = 0.05,
                                          labels = c('males',
                                                     'females'))
## graph
png('./neurons/neuropeptides/Nr3c1/limmatrend/RRHO2.Nr3c1.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Nr3c1.neuron.RRHO_obj,
              main = 'Nr3c1 DvsS')
dev.off()

#poster
Nr3c1.neuron.RRHO_obj.poster <-  RRHO2_initialize(Nr3c1.neuron.DvsS.M.rrho2, 
                                           Nr3c1.neuron.DvsS.F.rrho2, 
                                           boundary = 0.05)
pdf('./neurons/neuropeptides/Nr3c1/limmatrend/RRHO2.Nr3c1.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(Nr3c1.neuron.RRHO_obj.poster)
dev.off()

### Nr3c2
##male data
Nr3c2.neuron.DvsS.M.rrho2 = data.frame(gene = Nr3c2.neuron.limma.results.sct.df$Gene,
                                       value = -log10(Nr3c2.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                       direction = Nr3c2.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                       stringsAsFactors = FALSE)
#set positive negative
Nr3c2.neuron.DvsS.M.rrho2 = Nr3c2.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Nr3c2.neuron.DvsS.F.rrho2 = data.frame(gene = Nr3c2.neuron.limma.results.sct.df$Gene,
                                       value = -log10(Nr3c2.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                       direction = Nr3c2.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                       stringsAsFactors = FALSE)
#set positive negative
Nr3c2.neuron.DvsS.F.rrho2 = Nr3c2.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Nr3c2.neuron.RRHO_obj <-  RRHO2_initialize(Nr3c2.neuron.DvsS.M.rrho2, 
                                           Nr3c2.neuron.DvsS.F.rrho2, 
                                           boundary = 0.05,
                                           labels = c('males',
                                                      'females'))
## graph
png('./neurons/neuropeptides/Nr3c2/limmatrend/RRHO2.Nr3c2.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Nr3c2.neuron.RRHO_obj,
              main = 'Nr3c2 DvsS')
dev.off()


#poster
Nr3c2.neuron.RRHO_obj.poster <-  RRHO2_initialize(Nr3c2.neuron.DvsS.M.rrho2, 
                                                 Nr3c2.neuron.DvsS.F.rrho2, 
                                                 boundary = 0.05)
pdf('./neurons/neuropeptides/Nr3c2/limmatrend/RRHO2.Nr3c2.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(Nr3c2.neuron.RRHO_obj.poster)
dev.off()



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
  RR.list <- data.frame(Gene = df[quad$upup$positions,]$gene,
                        RRquadrant = 'upup') %>% 
    rbind(data.frame(Gene = df[quad$downdown$positions,]$gene,
                     RRquadrant = 'downdown')) %>% 
    rbind(data.frame(Gene = df[quad$updown$positions,]$gene,
                     RRquadrant = 'updown')) %>% 
    rbind(data.frame(Gene = df[quad$downup$positions,]$gene,
                     RRquadrant = 'downup')) %>% 
    mutate(Sample = celltype)
  
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


#### RedRibbon all cell types ####
## set scale value 
scale.value.RR = 217

## Ar
RedRibbon.all(celltype = "Ar",
              dataset.a = Ar.neuron.DvsS.M.rrho2,
              dataset.b = Ar.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'neurons/neuropeptides/RedRibbon/')

# scaled
RedRibbon.all(celltype = "Ar",
              dataset.a = Ar.neuron.DvsS.M.rrho2,
              dataset.b = Ar.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = scale.value.RR,
              file.name = 'neurons/neuropeptides/RedRibbon/')

## Pgr
RedRibbon.all(celltype = "Pgr",
              dataset.a = Pgr.neuron.DvsS.M.rrho2,
              dataset.b = Pgr.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'neurons/neuropeptides/RedRibbon/')

# scaled
RedRibbon.all(celltype = "Pgr",
              dataset.a = Pgr.neuron.DvsS.M.rrho2,
              dataset.b = Pgr.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = scale.value.RR,
              file.name = 'neurons/neuropeptides/RedRibbon/')

## Esr1
RedRibbon.all(celltype = "Esr1",
              dataset.a = Esr1.neuron.DvsS.M.rrho2,
              dataset.b = Esr1.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'neurons/neuropeptides/RedRibbon/')

# scaled
RedRibbon.all(celltype = "Esr1",
              dataset.a = Esr1.neuron.DvsS.M.rrho2,
              dataset.b = Esr1.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = scale.value.RR,
              file.name = 'neurons/neuropeptides/RedRibbon/')

## Nr3c1
RedRibbon.all(celltype = "Nr3c1",
              dataset.a = Nr3c1.neuron.DvsS.M.rrho2,
              dataset.b = Nr3c1.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'neurons/neuropeptides/RedRibbon/')

# scaled
RedRibbon.all(celltype = "Nr3c1",
              dataset.a = Nr3c1.neuron.DvsS.M.rrho2,
              dataset.b = Nr3c1.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = scale.value.RR,
              file.name = 'neurons/neuropeptides/RedRibbon/')

## Nr3c2
RedRibbon.all(celltype = "Nr3c2",
              dataset.a = Nr3c2.neuron.DvsS.M.rrho2,
              dataset.b = Nr3c2.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'neurons/neuropeptides/RedRibbon/')

# scaled
RedRibbon.all(celltype = "Nr3c2",
              dataset.a = Nr3c2.neuron.DvsS.M.rrho2,
              dataset.b = Nr3c2.neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = scale.value.RR,
              file.name = 'neurons/neuropeptides/RedRibbon/')


# # create data frame for red ribbon
# # needs to have an id (gene) col and one called 'a' and one called 'b'
# df = Nr3c2.neuron.DvsS.M.rrho2 %>% 
#   rename('a' =  'value') %>% 
#   full_join(Nr3c2.neuron.DvsS.F.rrho2 %>% 
#               rename('b' = 'value'))
# 
# ## Create RedRibbon object
# rr <- RedRibbon(df, 
#                 enrichment_mode="hyper-two-tailed")
# 
# ## Run the overlap using evolutionnary algorithm,
# ## computing permutation adjusted p-value for the four quadrants
# quad <- quadrants(rr, 
#                   algorithm="ea",
#                   permutation=TRUE, 
#                   whole=FALSE)
# 
# ## Plots the results
# ggRedRibbon(rr,
#             quadrants=quad) + 
#   coord_fixed(ratio = 1, 
#               clip = "off") +
#   xlab('Males') +
#   ylab('Females') +
#   ggtitle('neurons DvsS')
# 
# # scaled
# ggRedRibbon.rrho.scale(rr, 
#                        quadrants=quad,
#                        new.max.log = 155) + 
#   coord_fixed(ratio = 1, 
#               clip = "off") +
#   xlab('Males') +
#   ylab('Females') +
#   ggtitle('neurons DvsS')
# ggsave(paste0('./neurons/neuropeptides/RedRibbon/',
#               i,
#               ' RedRibbon.neurons.DvsS.sexes.png'),
#        height = 10,
#        width = 10)
# 
# ### compare RRHO2 to Redribbon
# # create list of RRHO outcomes
# RR.list <- list(
#   RRuu = df[quad$upup$positions,]$gene, 
#   RRdd = df[quad$downdown$positions,]$gene, 
#   RRud = df[quad$updown$positions,]$gene, 
#   RRdu = df[quad$downup$positions,]$gene)


#### create neuropeptide expression dataframe ####
### extract neuropeptide expression and orig.idents from neurons
## create expression dataframe of neuropeptides
mouse.snseq.combined.sct.neurons.np.expression = mouse.snseq.combined.sct.neurons@assays$SCT@data %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  mutate(gene = toupper(gene)) %>%
  filter(gene %in% neuropeptides.genes) %>% 
  column_to_rownames('gene') %>% 
  t() %>% 
  as.data.frame()

## create counts of presence absence
# use reads above 2 (in this case log(2) = 0.69)
mouse.snseq.combined.sct.neurons.np.expression.counts = mouse.snseq.combined.sct.neurons.np.expression %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < .68, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= .68, 1, x))) 

# create counts total matrix 
mouse.snseq.combined.sct.neurons.np.expression.counts.total = mouse.snseq.combined.sct.neurons.np.expression.counts %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename('Counts' = ".") %>% 
  rownames_to_column('gene')

## graph histogram of neuropeptides
mouse.snseq.combined.sct.neurons.np.expression.counts.total %>%
  ggplot(aes(Counts)) +
  geom_histogram(binwidth = 10) +
  theme_bw()
ggsave('./neurons/neuropeptides/comparison/histogram all neuropeptides.png',
       height = 10,
       width = 10)

### add orig.ident information
mouse.snseq.combined.sct.neurons.np.expression.sample = mouse.snseq.combined.sct.neurons.np.expression %>% 
  rownames_to_column('Cell.id') %>% 
  full_join(mouse.snseq.combined.sct.neurons@meta.data %>% 
              rownames_to_column('Cell.id') %>% 
              dplyr::select(orig.ident,
                            Cell.id))

### create matrix for each orig.ident
##male
#dom
mouse.snseq.combined.sct.neurons.np.expression.male.dom = mouse.snseq.combined.sct.neurons.np.expression.sample %>% 
  filter(orig.ident == 'Male.Dom') %>% 
  dplyr::select(-c(orig.ident,
                   Cell.id))
#sub
mouse.snseq.combined.sct.neurons.np.expression.male.sub = mouse.snseq.combined.sct.neurons.np.expression.sample %>% 
  filter(orig.ident == 'Male.Sub') %>% 
  dplyr::select(-c(orig.ident,
                   Cell.id))

##female
#dom
mouse.snseq.combined.sct.neurons.np.expression.Female.Dom = mouse.snseq.combined.sct.neurons.np.expression.sample %>% 
  filter(orig.ident == 'Female.Dom') %>% 
  dplyr::select(-c(orig.ident,
                   Cell.id))
#sub
mouse.snseq.combined.sct.neurons.np.expression.female.sub = mouse.snseq.combined.sct.neurons.np.expression.sample %>% 
  filter(orig.ident == 'Female.Dub') %>% 
  dplyr::select(-c(orig.ident,
                   Cell.id))

### create counts of presence absence for each sample
# use reads above 2 (in this case log(2) = 0.69)
## male
#dom
mouse.snseq.combined.sct.neurons.np.expression.male.dom.counts = mouse.snseq.combined.sct.neurons.np.expression.male.dom %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < .68, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= .68, 1, x))) 
#sub
mouse.snseq.combined.sct.neurons.np.expression.male.sub.counts = mouse.snseq.combined.sct.neurons.np.expression.male.sub %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < .68, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= .68, 1, x))) 
## female
#dom
mouse.snseq.combined.sct.neurons.np.expression.Female.Dom.counts = mouse.snseq.combined.sct.neurons.np.expression.Female.Dom %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < .68, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= .68, 1, x))) 
## male
#sub
mouse.snseq.combined.sct.neurons.np.expression.female.sub.counts = mouse.snseq.combined.sct.neurons.np.expression.female.sub %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < .68, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= .68, 1, x))) 

### create counts total matrix for each sample
##male
#dom
mouse.snseq.combined.sct.neurons.np.expression.male.dom.counts.total = mouse.snseq.combined.sct.neurons.np.expression.male.dom.counts %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename('Counts' = ".") %>% 
  rownames_to_column('gene') %>% 
  mutate(orig.ident = 'Male.Dom')
#sub
mouse.snseq.combined.sct.neurons.np.expression.male.sub.counts.total = mouse.snseq.combined.sct.neurons.np.expression.male.sub.counts %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename('Counts' = ".") %>% 
  rownames_to_column('gene') %>% 
  mutate(orig.ident = 'Male.Sub')
##female
#dom
mouse.snseq.combined.sct.neurons.np.expression.Female.Dom.counts.total = mouse.snseq.combined.sct.neurons.np.expression.Female.Dom.counts %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename('Counts' = ".") %>% 
  rownames_to_column('gene') %>% 
  mutate(orig.ident = 'Female.Dom')
#sub
mouse.snseq.combined.sct.neurons.np.expression.female.sub.counts.total = mouse.snseq.combined.sct.neurons.np.expression.female.sub.counts %>% 
  colSums() %>% 
  as.data.frame() %>% 
  dplyr::rename('Counts' = ".") %>% 
  rownames_to_column('gene') %>% 
  mutate(orig.ident = 'Female.Sub')

### combine dataframes together
mouse.snseq.combined.sct.neurons.np.expression.samples.counts.total = rbind(mouse.snseq.combined.sct.neurons.np.expression.male.dom.counts.total,
                                                                            mouse.snseq.combined.sct.neurons.np.expression.male.sub.counts.total) %>% 
  rbind(mouse.snseq.combined.sct.neurons.np.expression.Female.Dom.counts.total) %>% 
  rbind(mouse.snseq.combined.sct.neurons.np.expression.female.sub.counts.total)

#### neuropeptide genes limmatrend ####
### get list of neuropeptide genes to filter
filter.50.genes.subject = mouse.snseq.combined.sct.neurons.np.expression.samples.counts.total %>% 
  filter(Counts >= 50) %>% 
  pull(gene) %>% 
  unique() 

### extract neuropeptide expression and orig.idents from neurons
## create expression dataframe of neuropeptides
# remove neurons with no neuropeptides
mouse.snseq.combined.sct.neurons.np.expression = mouse.snseq.combined.sct.neurons@assays$SCT@counts %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  mutate(gene = toupper(gene)) %>%
  filter(gene %in% neuropeptides.genes) %>% 
  column_to_rownames('gene') %>% 
  select_if(colSums(.) != 0) %>% 
  t() %>% 
  as.data.frame()

## add a single count to PRL in sub male
mouse.snseq.combined.sct.neurons.np.expression[rownames(mouse.snseq.combined.sct.neurons.np.expression) %in% c('AAACCCAGTCGACTTA-1_3'),'PRL'] = 1


### add orig.ident information
mouse.snseq.combined.sct.neurons.np.expression.sample = mouse.snseq.combined.sct.neurons.np.expression %>% 
  rownames_to_column('Cell.id') %>% 
  left_join(mouse.snseq.combined.sct.neurons@meta.data %>% 
              rownames_to_column('Cell.id') %>% 
              dplyr::select(orig.ident,
                            Cell.id))

## create vector of factor
neuropeptides.neuron.reduce.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.np.expression.sample %>% 
  mutate(neuropeptides.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(neuropeptides.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# set negative values to 0
neuropeptides.neuron.reduce.vector.count.sct.prep = mouse.snseq.combined.sct.neurons.np.expression %>% 
  t() 

#create list
neuropeptides.neuron.reduce.vector.limma.sct.prep = list(count = neuropeptides.neuron.reduce.vector.count.sct.prep,
                                                 condt = neuropeptides.neuron.reduce.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub


# MDS not working

#create function
run_limmatrend_neuropeptides <- function(L) {
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
  pdf(file= "./neurons/neuropeptides/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./neurons/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
neuropeptides.neuron.limma.results.sct = run_limmatrend_neuropeptides(neuropeptides.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
neuropeptides.neuron.limma.results.sct.df = full_join(neuropeptides.neuron.limma.results.sct$ttDvsS %>% 
                                                rename_with(~paste0(.,"_DvsS")) %>% 
                                                rownames_to_column("Gene"),
                                              neuropeptides.neuron.limma.results.sct$ttDvsSM %>% 
                                                rename_with(~paste0(.,"_DvsS_M")) %>% 
                                                rownames_to_column("Gene")) %>% 
  full_join(neuropeptides.neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuropeptides.neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuropeptides.neuron.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuropeptides.neuron.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
neuropeptides.neuron.limma.results.sct.df = neuropeptides.neuron.limma.results.sct.df %>% 
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
  neuropeptides.neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/neuropeptides/limmatrend/limma.neuropeptides.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
#### save results neuropeptide limmatrend ####
save(neuropeptides.neuron.limma.results.sct.df,
     file = './neurons/neuropeptides/limmatrend/neuropeptides.neuron.limma.results.sct.df.RData')

# load('./neurons/neuropeptides/limmatrend/neuropeptides.neuron.limma.results.sct.df.RData')
#### RRHO2 neuropeptides ####
#### compare dom vs sub across sexes
##male data
neuropeptides.neuron.DvsS.M.rrho2 = data.frame(gene = neuropeptides.neuron.limma.results.sct.df$Gene,
                                       value = -log10(neuropeptides.neuron.limma.results.sct.df$P.Value_DvsS_M),
                                       direction = neuropeptides.neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                       stringsAsFactors = FALSE)
#set positive negative
neuropeptides.neuron.DvsS.M.rrho2 = neuropeptides.neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
neuropeptides.neuron.DvsS.F.rrho2 = data.frame(gene = neuropeptides.neuron.limma.results.sct.df$Gene,
                                       value = -log10(neuropeptides.neuron.limma.results.sct.df$P.Value_DvsS_F),
                                       direction = neuropeptides.neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                       stringsAsFactors = FALSE)
#set positive negative
neuropeptides.neuron.DvsS.F.rrho2 = neuropeptides.neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
neuropeptides.neuron.RRHO_obj.DvsS <-  RRHO2_initialize(neuropeptides.neuron.DvsS.M.rrho2, 
                                           neuropeptides.neuron.DvsS.F.rrho2, 
                                           boundary = 0.1,
                                           labels = c('males',
                                                      'females'),
                                           )
## graph
png('./neurons/neuropeptides/limmatrend/RRHO2.neuropeptides.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(neuropeptides.neuron.RRHO_obj.DvsS,
              main = 'neuropeptides DvsS')
dev.off()

#poster
neuropeptides.neuron.RRHO_obj.DvsS.poster <-  RRHO2_initialize(neuropeptides.neuron.DvsS.M.rrho2, 
                                                        neuropeptides.neuron.DvsS.F.rrho2, 
                                                        boundary = 0.1)

pdf('./neurons/neuropeptides/limmatrend/RRHO2.neuropeptides.neurons.DvsS.sexes.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(neuropeptides.neuron.RRHO_obj.DvsS.poster)
dev.off()



#### compare sexes across status
##dom data
neuropeptides.neuron.MvsF.D.rrho2 = data.frame(gene = neuropeptides.neuron.limma.results.sct.df$Gene,
                                               value = -log10(neuropeptides.neuron.limma.results.sct.df$P.Value_MvsFD),
                                               direction = neuropeptides.neuron.limma.results.sct.df$Direction.type_MvsFD, 
                                               stringsAsFactors = FALSE)
#set positive negative
neuropeptides.neuron.MvsF.D.rrho2 = neuropeptides.neuron.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
neuropeptides.neuron.MvsF.S.rrho2 = data.frame(gene = neuropeptides.neuron.limma.results.sct.df$Gene,
                                               value = -log10(neuropeptides.neuron.limma.results.sct.df$P.Value_MvsFS),
                                               direction = neuropeptides.neuron.limma.results.sct.df$Direction.type_MvsFS,
                                               stringsAsFactors = FALSE)
#set positive negative
neuropeptides.neuron.MvsF.S.rrho2 = neuropeptides.neuron.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
neuropeptides.neuron.RRHO_obj.MvsF <-  RRHO2_initialize(neuropeptides.neuron.MvsF.D.rrho2, 
                                                   neuropeptides.neuron.MvsF.S.rrho2, 
                                                   boundary = 0.1,
                                                   labels = c('Doms',
                                                              'Subs')
                                                   )
## graph
png('./neurons/neuropeptides/limmatrend/RRHO2.neuropeptides.neurons.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(neuropeptides.neuron.RRHO_obj.MvsF,
              main = 'neuropeptides MvsF')
dev.off()

#poster
neuropeptides.neuron.RRHO_obj.MvsF.poster <-  RRHO2_initialize(neuropeptides.neuron.MvsF.D.rrho2, 
                                                        neuropeptides.neuron.MvsF.S.rrho2, 
                                                        boundary = 0.1)

pdf('./neurons/neuropeptides/limmatrend/RRHO2.neuropeptides.neurons.MvsF.status.pdf',
    height = 10,
    width = 11)
RRHO2_heatmap(neuropeptides.neuron.RRHO_obj.MvsF.poster)
dev.off()


