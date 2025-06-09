#### Mouse snseq seurat analysis
### DEG/RRHO2 analysis
###Note: Seurat requires R version > 4
## use ccbbcomp1 to run R
# > R-4.3.1

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
# library(clustree)
library(pheatmap)
library(DEsingle)
library(dendextend)
library(tidyverse)
library(RRHO2)
# library(ggalluvial)
library(limma)
library(edgeR)

#### load data ####
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

#check cell types
#broad
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::pull(parent_id.broad.prob) %>% 
  unique()

# expanded
mouse.snseq.combined.sct@meta.data %>% 
  dplyr::pull(parent_id.exp.prob) %>% 
  unique()

### load neuropeptide list 
neuropeptides = read.csv('./gene.lists/neuropeptides.list.csv')
# genes
neuropeptides.genes = neuropeptides %>% 
  pull(Gene.name)


#### neurons limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.broad.prob"

#subset to neurons
mouse.snseq.combined.sct.neurons = subset(mouse.snseq.combined.sct,
                                          idents = c("C7-2: GABA",
                                                     "C7-1: GLU"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.neurons) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.neurons@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 3085       1240       1667       2190 


### extract gene expression and orig.idents from neurons
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.neurons.var <- FindVariableFeatures(mouse.snseq.combined.sct.neurons, 
                                                             assay = 'integrated',
                                                             selection.method = "vst", 
                                                             verbose = F)
# identify top variable genes
neuron.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.neurons.var,
                                                           assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.neurons.expression = full_join(full_join(mouse.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
                                                                    as.data.frame() %>% 
                                                                    rownames_to_column("Cell.id"),
                                                                  mouse.snseq.combined.sct.neurons@meta.data %>%
                                                                    rownames_to_column("Cell.id")),
                                                        mouse.snseq.combined.sct.neurons@assays$SCT@data %>% 
                                                          as.data.frame()  %>% 
                                                          t() %>% as.data.frame() %>% 
                                                          rownames_to_column('Cell.id'))


## create vector of factor
neuron.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.neurons.expression %>% 
  mutate(neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
neuron.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.neurons,
                                                   assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.neurons.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
neuron.reduce.vector.limma.sct.prep = list(count = neuron.reduce.vector.count.sct.prep,
                                           condt = neuron.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_neurons <- function(L) {
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
  pdf(file= "./neurons/limmatrend/limmatrend.histograms.pdf" )
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
neuron.limma.results.sct = run_limmatrend_neurons(neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
neuron.limma.results.sct.df = full_join(neuron.limma.results.sct$ttDvsS %>% 
                                          rename_with(~paste0(.,"_DvsS")) %>% 
                                          rownames_to_column("Gene"),
                                        neuron.limma.results.sct$ttDvsSM %>% 
                                          rename_with(~paste0(.,"_DvsS_M")) %>% 
                                          rownames_to_column("Gene")) %>% 
  full_join(neuron.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(neuron.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
neuron.limma.results.sct.df = neuron.limma.results.sct.df %>% 
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
  neuron.limma.results.sct.df %>% 
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
  ggsave(paste0('./neurons/limmatrend/limma.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

# highlight neuropeptides
neuron.limma.results.sct.df.add.np = neuron.limma.results.sct.df %>% 
  mutate(GENE = toupper(Gene)) %>% 
  left_join(neuropeptides.genes %>% 
              as.data.frame() %>% 
              rename('GENE' = '.') %>% 
              mutate(Neuropeptide = 1)) %>% 
  mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
                               0,
                               Neuropeptide))


for (i in limma.genotpye.vector) {
  # graph volcano plot
  neuron.limma.results.sct.df.add.np %>% 
    mutate(sig.label = ifelse(Neuropeptide == 1,
                              GENE,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = factor(Neuropeptide)))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("grey","red")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('./neurons/limmatrend/limma.neurons.reduce.volcano.', i, '.neuropeptides.png'),
         width = 5,
         height = 5)
}

#### save results neuron limmatrend ####
save(neuron.limma.results.sct.df,
     file = './neurons/limmatrend/neuropeptides.neuron.limma.results.sct.df.RData')

write.csv(neuron.limma.results.sct.df,
     file = './neurons/limmatrend/neuron.limma.results.sct.df.csv')

# load('./neurons/limmatrend/neuropeptides.neuron.limma.results.sct.df.RData')
#### RRHO2 neurons ####
#### compare dom vs sub across sexes
##male data
neuron.DvsS.M.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                 value = -log10(neuron.limma.results.sct.df$P.Value_DvsS_M),
                                 direction = neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                 stringsAsFactors = FALSE)
#set positive negative
neuron.DvsS.M.rrho2 = neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
neuron.DvsS.F.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                 value = -log10(neuron.limma.results.sct.df$P.Value_DvsS_F),
                                 direction = neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                 stringsAsFactors = FALSE)
#set positive negative
neuron.DvsS.F.rrho2 = neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

#### create rrho2 object ####
neuron.RRHO_obj.DvsS <-  RRHO2_initialize(neuron.DvsS.M.rrho2, 
                                          neuron.DvsS.F.rrho2, 
                                          boundary = 0.1,
                                          labels = c('males',
                                                     'females'),
)
## graph
png('./neurons/limmatrend/RRHO2.neurons.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(neuron.RRHO_obj.DvsS,
              main = 'neurons DvsS')
dev.off()

### save out list 
# need to add NA for updown quadrant 
# because no genes present
neuron.RRHO_obj.DvsS.df = data.frame(type = 'upup.overlap',
                                     gene = neuron.RRHO_obj.DvsS$genelist_uu$gene_list_overlap_uu) %>%
  rbind(data.frame(type = 'downdown.overlap',
                   gene = neuron.RRHO_obj.DvsS$genelist_dd$gene_list_overlap_dd)) %>% 
  rbind(data.frame(type = 'updown.overlap',
                   gene = c(neuron.RRHO_obj.DvsS$genelist_ud$gene_list_overlap_ud, 
                            NA))) %>% 
  rbind(data.frame(type = 'downup.overlap',
                   gene = neuron.RRHO_obj.DvsS$genelist_du$gene_list_overlap_du)) %>% 
  separate_wider_delim(type,
                        delim = '.',
                        names = c('direction',
                                  'RRHO2_type'))
# save file to neuron folder  
write_csv(neuron.RRHO_obj.DvsS.df,
          file = 'neurons/neuron_RRHO2_results.csv')  
# # load file
# neuron.RRHO_obj.DvsS.df = read_csv('neurons/neuron_RRHO2_results.csv')

## compare results with cluster marker specificity list
# load cluster marker specificity list
# from: neurons_snseq_mouse_IMC.R
neuron.markers.df = read.csv('neurons/neuron_seurat_cluster_markers.csv')

# check length of lists
# take top ~1000 markers
# calculate percentage of RRHO genes
# 43.56%
100*(neuron.markers.df %>% 
  filter(p_val_adj <= 0.05,
         abs(specificity) >= 0.5) %>% 
  pull(gene) %>% 
  intersect(neuron.RRHO_obj.DvsS.df %>% 
              pull(gene)) %>% 
  length())/length(neuron.RRHO_obj.DvsS.df %>% 
                     pull(gene))

# calculate for just up/up
# 80.29%
100*(neuron.markers.df %>% 
       filter(p_val_adj <= 0.05,
              specificity >= 0.5) %>% 
       pull(gene) %>% 
       intersect(neuron.RRHO_obj.DvsS.df %>% 
                   filter(direction == 'upup') %>% 
                   pull(gene)) %>% 
       length())/length(neuron.RRHO_obj.DvsS.df %>% 
                          filter(direction == 'upup') %>% 
                          pull(gene))

# create graph to show number of RRHO up up genes per cluster
neuron.markers.df %>% 
  filter(p_val_adj <= 0.05,
         specificity >= 0.5) %>% 
  filter(gene %in% c(neuron.RRHO_obj.DvsS.df %>%
           filter(direction == 'upup') %>% 
           pull(gene))) %>% 
  group_by(cluster) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = reorder(as.character(cluster),
                         -count),
             y = count)) +
  geom_point() +
  theme_classic() +
  xlab('neuron cluster') +
  ylab('number of genes') +
  ggtitle('Overlap of RRHO2 upup genes and neuron cluster markers',
          subtitle = 'adj. p-value < 0.05; specificity > 0.58')
ggsave('neurons/neurons_RRHO2_upup_cluster_markers.png',
       height = 10,
       width = 10)

# graph count of genes per quadrant
neuron.RRHO_obj.DvsS.df %>% 
  group_by(direction) %>% 
  summarise(Total = n()) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = reorder(direction,
                         -Total),
             y = Total)) +
  geom_point() +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RRHO2 neuron quadrant results')
ggsave('neurons/neurons_RRHO2_number_gene_quad.png')

### run GO analysis on enriched genes

## load GO libraries 
# load libraries 
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

## upup create list of genes to test with GO
neuron.RRHO.DvsS.upup.genes = neuron.RRHO_obj.DvsS.df %>%
  filter(direction == 'upup') %>% 
  pull(gene)

# ## convert to ENSEMBL ID
# # load biomart to conovert to gene ID
# library(biomaRt)
# # get mouse data mart
# mouse.mart <- useEnsembl("ensembl",
#                          "mmusculus_gene_ensembl",
#                          mirror = 'useast')
# # check attributes list
# # listAttributes(mouse.mart) %>% view()
# # convert from names to ensembl ID and GO ID
# neuron.RRHO.DvsS.upup.genes.ensemble <- getBM(c("ensembl_gene_id",
#                                                 "go_id",
#                                                 'external_gene_name'),
#                                               filters = "external_gene_name", 
#                                               neuron.RRHO.DvsS.upup.genes,
#                                               mouse.mart)
# 
# 
# ## run GO enrichment
# neuron.RRHO.DvsS.upup.genes.go <- enrichGO(gene = neuron.RRHO.DvsS.upup.genes.ensemble %>% 
#                                              pull(ensembl_gene_id) %>% 
#                                              unique(), 
#                        OrgDb = "org.Mm.eg.db",
#                        keyType = "ENSEMBL", 
#                        ont = "BP")

# can use gene symbol, don't need ensembl ids
neuron.RRHO.DvsS.upup.genes.go <- enrichGO(gene = neuron.RRHO.DvsS.upup.genes %>% 
                                             unique(), 
                                           OrgDb = "org.Mm.eg.db",
                                           keyType = "SYMBOL", 
                                           ont = "BP")

## graph GO results
# create plot
neuron.RRHO.DvsS.upup.genes.go.fit <- plot(barplot(neuron.RRHO.DvsS.upup.genes.go,
                                                   showCategory = 15))

neuron.RRHO.DvsS.upup.genes.go.fit

# save plot
png("neurons/neuron_RRHO2_results_upup.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(neuron.RRHO.DvsS.upup.genes.go.fit)
dev.off()

## create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

# graph tree
png("neurons/neuron_RRHO2_results_upup_tree.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
neuron.RRHO.DvsS.upup.genes.go %>% 
  pairwise_termsim() %>% 
  treeplot()
dev.off()

## down down create list of genes to test with GO
neuron.RRHO.DvsS.downdown.genes = neuron.RRHO_obj.DvsS.df %>%
  filter(direction == 'downdown') %>% 
  pull(gene)


# can use gene symbol, don't need ensembl ids
neuron.RRHO.DvsS.downdown.genes.go <- enrichGO(gene = neuron.RRHO.DvsS.downdown.genes %>% 
                                             unique(), 
                                           OrgDb = "org.Mm.eg.db",
                                           keyType = "SYMBOL", 
                                           ont = "BP")

## graph GO results
# create plot
neuron.RRHO.DvsS.downdown.genes.go.fit <- plot(barplot(neuron.RRHO.DvsS.downdown.genes.go %>% 
                                                         filter(p.adjust <= 0.001) %>% 
                                                         filter(Count >= 5),
                                                   showCategory = 10))

neuron.RRHO.DvsS.downdown.genes.go.fit

# save plot
png("neurons/neuron_RRHO2_results_downdown.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(neuron.RRHO.DvsS.downdown.genes.go.fit)
dev.off()

## create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot

# graph tree
png("neurons/neuron_RRHO2_results_downdown_tree.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
neuron.RRHO.DvsS.downdown.genes.go %>% 
  enrichplot::pairwise_termsim() %>% 
  enrichplot::treeplot()
dev.off()



#### compare sexes across status
##dom data
neuron.MvsF.D.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                 value = -log10(neuron.limma.results.sct.df$P.Value_MvsFD),
                                 direction = neuron.limma.results.sct.df$Direction.type_MvsFD, 
                                 stringsAsFactors = FALSE)
#set positive negative
neuron.MvsF.D.rrho2 = neuron.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
neuron.MvsF.S.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                 value = -log10(neuron.limma.results.sct.df$P.Value_MvsFS),
                                 direction = neuron.limma.results.sct.df$Direction.type_MvsFS,
                                 stringsAsFactors = FALSE)
#set positive negative
neuron.MvsF.S.rrho2 = neuron.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
neuron.RRHO_obj.MvsF <-  RRHO2_initialize(neuron.MvsF.D.rrho2, 
                                          neuron.MvsF.S.rrho2, 
                                          boundary = 0.1,
                                          labels = c('Doms',
                                                     'Subs')
)
## graph
png('./neurons/limmatrend/RRHO2.neurons.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(neuron.RRHO_obj.MvsF,
              main = 'neurons MvsF')
dev.off()




#### RedRibbon neurons ####
# https://github.com/antpiron/RedRibbon
# install
# devtools::install_github("antpiron/RedRibbon")
# install library
library(RedRibbon)

# create data frame for red ribbon
# needs to have an id (gene) col and one called 'a' and one called 'b'
df = neuron.DvsS.M.rrho2 %>% 
  rename('a' =  'value') %>% 
  full_join(neuron.DvsS.F.rrho2 %>% 
              rename('b' = 'value'))

## Create RedRibbon object
rr <- RedRibbon(df, 
                enrichment_mode="hyper-two-tailed")

## Run the overlap using evolutionnary algorithm,
## computing permutation adjusted p-value for the four quadrants
quad <- quadrants(rr, 
                  algorithm="ea",
                  permutation=TRUE, 
                  whole=FALSE)

## Plots the results
ggRedRibbon(rr, 
            quadrants=quad) + 
  coord_fixed(ratio = 1, 
              clip = "off") +
  xlab('Males') +
  ylab('Females') +
  ggtitle('neurons DvsS')
ggsave('./neurons/limmatrend/RedRibbon.neurons.DvsS.sexes.png',
       height = 10,
       width = 10)

## Get the down-down quandrant list of genes
# df[quad$upup$positions,] %>% View()
### save out list 
# need to add NA for updown quadrant 
# because no genes present
neuron.rr_obj.DvsS.df = data.frame(type = 'upup.overlap',
                                     gene = df[quad$upup$positions,]$gene) %>%
  rbind(data.frame(type = 'downdown.overlap',
                   gene = df[quad$downdown$positions,]$gene)) %>% 
  rbind(data.frame(type = 'updown.overlap',
                   gene = df[quad$updown$positions,]$gene)) %>% 
  rbind(data.frame(type = 'downup.overlap',
                   gene = df[quad$downup$positions,]$gene)) %>% 
  separate_wider_delim(type,
                       delim = '.',
                       names = c('direction',
                                 'RRHO2_type'))
# save file to neuron folder  
write_csv(neuron.rr_obj.DvsS.df,
          file = 'neurons/neuron_redribbon_results.csv')  
# # load file
# neuron.rr_obj.DvsS.df = read_csv('neurons/neuron_redribbon_results.csv')

## compare results with cluster marker specificity list
# load cluster marker specificity list
# from: neurons_snseq_mouse_IMC.R
neuron.markers.df = read.csv('neurons/neuron_seurat_cluster_markers.csv')

# check length of lists
# take top ~1185 markers
# calculate percentage of RRHO genes
# 26.53%
100*(neuron.markers.df %>% 
       filter(p_val_adj <= 0.05,
              abs(specificity) >= 0.5) %>% 
       pull(gene) %>% 
       intersect(neuron.rr_obj.DvsS.df %>% 
                   pull(gene)) %>% 
       length())/length(neuron.rr_obj.DvsS.df %>% 
                          pull(gene))

# calculate for just up/up
# 79.02%
100*(neuron.markers.df %>% 
       filter(p_val_adj <= 0.05,
              specificity >= 0.5) %>% 
       pull(gene) %>% 
       intersect(neuron.rr_obj.DvsS.df %>% 
                   filter(direction == 'upup') %>% 
                   pull(gene)) %>% 
       length())/length(neuron.rr_obj.DvsS.df %>% 
                          filter(direction == 'upup') %>% 
                          pull(gene))

# create list of overlapping genes between RRHO upup and marker genes
neuron.markers.df.rr.upup = neuron.markers.df %>% 
  filter(p_val_adj <= 0.05,
         specificity >= 0.5) %>% 
  filter(gene %in% c(neuron.rr_obj.DvsS.df %>%
                       filter(direction == 'upup') %>% 
                       pull(gene)))

# save
write_csv(neuron.markers.df.rr.upup,
          file = 'neurons/neuron_seurat_cluster_markers_rr_upup.csv')

# create graph to show number of RRHO up up genes per cluster
neuron.markers.df %>% 
  filter(p_val_adj <= 0.05,
         specificity >= 0.5) %>% 
  filter(gene %in% c(neuron.rr_obj.DvsS.df %>%
                       filter(direction == 'upup') %>% 
                       pull(gene))) %>% 
  group_by(cluster) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x = reorder(as.character(cluster),
                         -count),
             y = count)) +
  geom_point() +
  theme_classic() +
  xlab('neuron cluster') +
  ylab('number of genes') +
  ggtitle('Overlap of RR upup genes and neuron cluster markers',
          subtitle = 'adj. p-value < 0.05; specificity > 0.58')
ggsave('neurons/neurons_rr_upup_cluster_markers.png',
       height = 10,
       width = 10)

# graph count of genes per quadrant
neuron.rr_obj.DvsS.df %>% 
  group_by(direction) %>% 
  summarise(Total = n()) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = reorder(direction,
                         -Total),
             y = Total)) +
  geom_point() +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RR neuron quadrant results')
ggsave('neurons/neurons_rr_number_gene_quad.png')

### run GO analysis on enriched genes

## load GO libraries 
# load libraries 
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

## upup create list of genes to test with GO
neuron.rr.DvsS.upup.genes = neuron.rr_obj.DvsS.df %>%
  filter(direction == 'upup') %>% 
  pull(gene)

# can use gene symbol, don't need ensembl ids
neuron.rr.DvsS.upup.genes.go <- enrichGO(gene = neuron.rr.DvsS.upup.genes %>% 
                                             unique(), 
                                           OrgDb = "org.Mm.eg.db",
                                           keyType = "SYMBOL", 
                                           ont = "BP")

## graph GO results
# create plot
neuron.rr.DvsS.upup.genes.go.fit <- plot(barplot(neuron.rr.DvsS.upup.genes.go,
                                                   showCategory = 15))

neuron.rr.DvsS.upup.genes.go.fit

# save plot
png("neurons/neuron_rr_results_upup.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(neuron.rr.DvsS.upup.genes.go.fit)
dev.off()

## create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

library(enrichplot)
# graph tree
png("neurons/neuron_rr_results_upup_tree.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
neuron.rr.DvsS.upup.genes.go %>% 
  enrichplot::pairwise_termsim() %>% 
  enrichplot::treeplot()
dev.off()

## down down create list of genes to test with GO
neuron.rr.DvsS.downdown.genes = neuron.rr_obj.DvsS.df %>%
  filter(direction == 'downdown') %>% 
  pull(gene)


# can use gene symbol, don't need ensembl ids
neuron.rr.DvsS.downdown.genes.go <- enrichGO(gene = neuron.rr.DvsS.downdown.genes %>% 
                                                 unique(), 
                                               OrgDb = "org.Mm.eg.db",
                                               keyType = "SYMBOL", 
                                               ont = "BP")

## graph GO results
# create plot
neuron.rr.DvsS.downdown.genes.go.fit <- plot(barplot(neuron.rr.DvsS.downdown.genes.go %>% 
                                                         filter(p.adjust <= 0.001) %>% 
                                                         filter(Count >= 5),
                                                       showCategory = 10))

neuron.rr.DvsS.downdown.genes.go.fit

# save plot
png("neurons/neuron_rr_results_downdown.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(neuron.rr.DvsS.downdown.genes.go.fit)
dev.off()

## create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot

# graph tree
png("neurons/neuron_rr_results_downdown_tree.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
neuron.rr.DvsS.downdown.genes.go %>% 
  enrichplot::pairwise_termsim() %>% 
  enrichplot::treeplot()
dev.off()

# # compare with RRHO2
# neuron.RRHO_obj.DvsS[["genelist_ii"]][["gene_list_overlap_uu"]] %>% 
#   intersect(df[quad$upup$positions,] %>% 
#               pull(gene)) 
# 
# ### compare RRHO2 to Redribbon
# # create list of RRHO outcomes
# RRHO.list <- list(
#   RRHO2uu = neuron.RRHO_obj.DvsS[["genelist_uu"]][["gene_list_overlap_uu"]], 
#   RRHO2dd = neuron.RRHO_obj.DvsS[["genelist_dd"]][["gene_list_overlap_dd"]], 
#   RRHO2ud = neuron.RRHO_obj.DvsS[["genelist_"]][["gene_list_overlap_ud"]], 
#   RRHO2du = neuron.RRHO_obj.DvsS[["genelist_du"]][["gene_list_overlap_du"]], 
#   RRuu = df[quad$upup$positions,]$gene, 
#   RRdd = df[quad$downdown$positions,]$gene, 
#   RRud = df[quad$updown$positions,]$gene, 
#   RRdu = df[quad$downup$positions,]$gene 
# )
# 
# # load UpsetR library
# library(UpSetR)
# 
# ## Graph upset plot to compare results of RedRibbon with RRHO2
# png('./neurons/limmatrend/RedRibbon vs RRHO2 DvsS.png',
#     height = 10,
#     width = 10,
#     units = 'in',
#     res = 300)
# upset(fromList(RRHO.list), 
#       order.by = "freq")
# dev.off()



#### GLU limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.broad.prob"

#subset to GLU
mouse.snseq.combined.sct.GLU = subset(mouse.snseq.combined.sct,
                                             idents = c("C7-1: GLU"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.GLU) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.GLU@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 748        298        582        740 

### extract gene expression and orig.idents from GLU
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.GLU.var <- FindVariableFeatures(mouse.snseq.combined.sct.GLU, 
                                                                assay = 'integrated',
                                                                selection.method = "vst", 
                                                                verbose = F)
# identify top variable genes
GLU.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.GLU.var,
                                                          assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.GLU.expression = full_join(full_join(mouse.snseq.combined.sct.GLU@reductions$umap@cell.embeddings %>% 
                                                                       as.data.frame() %>% 
                                                                       rownames_to_column("Cell.id"),
                                                                     mouse.snseq.combined.sct.GLU@meta.data %>%
                                                                       rownames_to_column("Cell.id")),
                                                           mouse.snseq.combined.sct.GLU@assays$SCT@data %>% 
                                                             as.data.frame()  %>% 
                                                             t() %>% as.data.frame() %>% 
                                                             rownames_to_column('Cell.id'))


## create vector of factor
GLU.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.GLU.expression %>% 
  mutate(GLU.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(GLU.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
GLU.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.GLU,
                                                       assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.GLU.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% GLU.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
GLU.reduce.vector.limma.sct.prep = list(count = GLU.reduce.vector.count.sct.prep,
                                               condt = GLU.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_GLU <- function(L) {
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
  pdf(file= "./GLU/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./GLU/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
GLU.limma.results.sct = run_limmatrend_GLU(GLU.reduce.vector.limma.sct.prep)

# save results to dataframe
GLU.limma.results.sct.df = full_join(GLU.limma.results.sct$ttDvsS %>% 
                                              rename_with(~paste0(.,"_DvsS")) %>% 
                                              rownames_to_column("Gene"),
                                            GLU.limma.results.sct$ttDvsSM %>% 
                                              rename_with(~paste0(.,"_DvsS_M")) %>% 
                                              rownames_to_column("Gene")) %>% 
  full_join(GLU.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GLU.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GLU.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GLU.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
GLU.limma.results.sct.df = GLU.limma.results.sct.df %>% 
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
  GLU.limma.results.sct.df %>% 
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
  ggsave(paste0('./GLU/limmatrend/limma.GLU.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# GLU.limma.results.sct.df.add.np = GLU.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   GLU.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./GLU/limmatrend/limma.GLU.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results GLU limmatrend ####
save(GLU.limma.results.sct.df,
     file = './GLU/limmatrend/neuropeptides.GLU.limma.results.sct.df.RData')

# load('./GLU/limmatrend/neuropeptides.GLU.limma.results.sct.df.RData')
#### RRHO2 GLU ####
#### compare dom vs sub across sexes
##male data
GLU.DvsS.M.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                                     value = -log10(GLU.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = GLU.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
GLU.DvsS.M.rrho2 = GLU.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
GLU.DvsS.F.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                                     value = -log10(GLU.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = GLU.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
GLU.DvsS.F.rrho2 = GLU.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
GLU.RRHO_obj.DvsS <-  RRHO2_initialize(GLU.DvsS.M.rrho2, 
                                              GLU.DvsS.F.rrho2, 
                                              boundary = 0.1,
                                              labels = c('males',
                                                         'females'),
)
## graph
png('./GLU/limmatrend/RRHO2.GLU.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(GLU.RRHO_obj.DvsS,
              main = 'GLU DvsS')
dev.off()

#### compare sexes across status
##dom data
GLU.MvsF.D.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                                     value = -log10(GLU.limma.results.sct.df$P.Value_MvsFD),
                                     direction = GLU.limma.results.sct.df$Direction.type_MvsFD, 
                                     stringsAsFactors = FALSE)
#set positive negative
GLU.MvsF.D.rrho2 = GLU.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
GLU.MvsF.S.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                                     value = -log10(GLU.limma.results.sct.df$P.Value_MvsFS),
                                     direction = GLU.limma.results.sct.df$Direction.type_MvsFS,
                                     stringsAsFactors = FALSE)
#set positive negative
GLU.MvsF.S.rrho2 = GLU.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
GLU.RRHO_obj.MvsF <-  RRHO2_initialize(GLU.MvsF.D.rrho2, 
                                              GLU.MvsF.S.rrho2, 
                                              boundary = 0.1,
                                              labels = c('Doms',
                                                         'Subs')
)
## graph
png('./GLU/limmatrend/RRHO2.GLU.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(GLU.RRHO_obj.MvsF,
              main = 'GLU MvsF')
dev.off()





#### GABA limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.broad.prob"

#subset to GABA
mouse.snseq.combined.sct.GABA = subset(mouse.snseq.combined.sct,
                                      idents = c("C7-2: GABA"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.GABA) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.GABA@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 748        298        582        740 

### extract gene expression and orig.idents from GABA
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.GABA.var <- FindVariableFeatures(mouse.snseq.combined.sct.GABA, 
                                                         assay = 'integrated',
                                                         selection.method = "vst", 
                                                         verbose = F)
# identify top variable genes
GABA.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.GABA.var,
                                                   assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.GABA.expression = full_join(full_join(mouse.snseq.combined.sct.GABA@reductions$umap@cell.embeddings %>% 
                                                                as.data.frame() %>% 
                                                                rownames_to_column("Cell.id"),
                                                              mouse.snseq.combined.sct.GABA@meta.data %>%
                                                                rownames_to_column("Cell.id")),
                                                    mouse.snseq.combined.sct.GABA@assays$SCT@data %>% 
                                                      as.data.frame()  %>% 
                                                      t() %>% as.data.frame() %>% 
                                                      rownames_to_column('Cell.id'))


## create vector of factor
GABA.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.GABA.expression %>% 
  mutate(GABA.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(GABA.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
GABA.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.GABA,
                                                assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.GABA.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% GABA.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
GABA.reduce.vector.limma.sct.prep = list(count = GABA.reduce.vector.count.sct.prep,
                                        condt = GABA.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_GABA <- function(L) {
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
  pdf(file= "./GABA/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./GABA/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
GABA.limma.results.sct = run_limmatrend_GABA(GABA.reduce.vector.limma.sct.prep)

# save results to dataframe
GABA.limma.results.sct.df = full_join(GABA.limma.results.sct$ttDvsS %>% 
                                       rename_with(~paste0(.,"_DvsS")) %>% 
                                       rownames_to_column("Gene"),
                                     GABA.limma.results.sct$ttDvsSM %>% 
                                       rename_with(~paste0(.,"_DvsS_M")) %>% 
                                       rownames_to_column("Gene")) %>% 
  full_join(GABA.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GABA.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GABA.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(GABA.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
GABA.limma.results.sct.df = GABA.limma.results.sct.df %>% 
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
  GABA.limma.results.sct.df %>% 
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
  ggsave(paste0('./GABA/limmatrend/limma.GABA.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# GABA.limma.results.sct.df.add.np = GABA.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   GABA.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./GABA/limmatrend/limma.GABA.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results GABA limmatrend ####
save(GABA.limma.results.sct.df,
     file = './GABA/limmatrend/neuropeptides.GABA.limma.results.sct.df.RData')

# load('./GABA/limmatrend/neuropeptides.GABA.limma.results.sct.df.RData')
#### RRHO2 GABA ####
#### compare dom vs sub across sexes
##male data
GABA.DvsS.M.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                              value = -log10(GABA.limma.results.sct.df$P.Value_DvsS_M),
                              direction = GABA.limma.results.sct.df$Direction.type_DvsS_M, 
                              stringsAsFactors = FALSE)
#set positive negative
GABA.DvsS.M.rrho2 = GABA.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
GABA.DvsS.F.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                              value = -log10(GABA.limma.results.sct.df$P.Value_DvsS_F),
                              direction = GABA.limma.results.sct.df$Direction.type_DvsS_F,
                              stringsAsFactors = FALSE)
#set positive negative
GABA.DvsS.F.rrho2 = GABA.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
GABA.RRHO_obj.DvsS <-  RRHO2_initialize(GABA.DvsS.M.rrho2, 
                                       GABA.DvsS.F.rrho2, 
                                       boundary = 0.1,
                                       labels = c('males',
                                                  'females'),
)
## graph
png('./GABA/limmatrend/RRHO2.GABA.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(GABA.RRHO_obj.DvsS,
              main = 'GABA DvsS')
dev.off()

#### compare sexes across status
##dom data
GABA.MvsF.D.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                              value = -log10(GABA.limma.results.sct.df$P.Value_MvsFD),
                              direction = GABA.limma.results.sct.df$Direction.type_MvsFD, 
                              stringsAsFactors = FALSE)
#set positive negative
GABA.MvsF.D.rrho2 = GABA.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
GABA.MvsF.S.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                              value = -log10(GABA.limma.results.sct.df$P.Value_MvsFS),
                              direction = GABA.limma.results.sct.df$Direction.type_MvsFS,
                              stringsAsFactors = FALSE)
#set positive negative
GABA.MvsF.S.rrho2 = GABA.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
GABA.RRHO_obj.MvsF <-  RRHO2_initialize(GABA.MvsF.D.rrho2, 
                                       GABA.MvsF.S.rrho2, 
                                       boundary = 0.1,
                                       labels = c('Doms',
                                                  'Subs')
)
## graph
png('./GABA/limmatrend/RRHO2.GABA.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(GABA.RRHO_obj.MvsF,
              main = 'GABA MvsF')
dev.off()





#### astrocytes limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

#subset to astrocytes
mouse.snseq.combined.sct.astrocytes = subset(mouse.snseq.combined.sct,
                                          idents = c("C25-18: Astrocytes"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.astrocytes) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.astrocytes@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 1175       1443       1479       1256 

### extract gene expression and orig.idents from astrocytes
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.astrocytes.var <- FindVariableFeatures(mouse.snseq.combined.sct.astrocytes, 
                                                             assay = 'integrated',
                                                             selection.method = "vst", 
                                                             verbose = F)
# identify top variable genes
astrocytes.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.astrocytes.var,
                                                           assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.astrocytes.expression = full_join(full_join(mouse.snseq.combined.sct.astrocytes@reductions$umap@cell.embeddings %>% 
                                                                    as.data.frame() %>% 
                                                                    rownames_to_column("Cell.id"),
                                                                  mouse.snseq.combined.sct.astrocytes@meta.data %>%
                                                                    rownames_to_column("Cell.id")),
                                                        mouse.snseq.combined.sct.astrocytes@assays$SCT@data %>% 
                                                          as.data.frame()  %>% 
                                                          t() %>% as.data.frame() %>% 
                                                          rownames_to_column('Cell.id'))


## create vector of factor
astrocytes.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.astrocytes.expression %>% 
  mutate(astrocytes.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(astrocytes.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
astrocytes.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.astrocytes,
                                                   assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.astrocytes.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% astrocytes.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
astrocytes.reduce.vector.limma.sct.prep = list(count = astrocytes.reduce.vector.count.sct.prep,
                                           condt = astrocytes.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_astrocytes <- function(L) {
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
  pdf(file= "./astrocytes/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./astrocytes/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
astrocytes.limma.results.sct = run_limmatrend_astrocytes(astrocytes.reduce.vector.limma.sct.prep)

# save results to dataframe
astrocytes.limma.results.sct.df = full_join(astrocytes.limma.results.sct$ttDvsS %>% 
                                          rename_with(~paste0(.,"_DvsS")) %>% 
                                          rownames_to_column("Gene"),
                                        astrocytes.limma.results.sct$ttDvsSM %>% 
                                          rename_with(~paste0(.,"_DvsS_M")) %>% 
                                          rownames_to_column("Gene")) %>% 
  full_join(astrocytes.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(astrocytes.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(astrocytes.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(astrocytes.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
astrocytes.limma.results.sct.df = astrocytes.limma.results.sct.df %>% 
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
  astrocytes.limma.results.sct.df %>% 
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
  ggsave(paste0('./astrocytes/limmatrend/limma.astrocytes.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# astrocytes.limma.results.sct.df.add.np = astrocytes.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   astrocytes.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./astrocytes/limmatrend/limma.astrocytes.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results astrocytes limmatrend ####
save(astrocytes.limma.results.sct.df,
     file = './astrocytes/limmatrend/neuropeptides.astrocytes.limma.results.sct.df.RData')

# load('./astrocytes/limmatrend/neuropeptides.astrocytes.limma.results.sct.df.RData')
#### RRHO2 astrocytes ####
#### compare dom vs sub across sexes
##male data
astrocytes.DvsS.M.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                 value = -log10(astrocytes.limma.results.sct.df$P.Value_DvsS_M),
                                 direction = astrocytes.limma.results.sct.df$Direction.type_DvsS_M, 
                                 stringsAsFactors = FALSE)
#set positive negative
astrocytes.DvsS.M.rrho2 = astrocytes.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
astrocytes.DvsS.F.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                 value = -log10(astrocytes.limma.results.sct.df$P.Value_DvsS_F),
                                 direction = astrocytes.limma.results.sct.df$Direction.type_DvsS_F,
                                 stringsAsFactors = FALSE)
#set positive negative
astrocytes.DvsS.F.rrho2 = astrocytes.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
astrocytes.RRHO_obj.DvsS <-  RRHO2_initialize(astrocytes.DvsS.M.rrho2, 
                                          astrocytes.DvsS.F.rrho2, 
                                          boundary = 0.1,
                                          labels = c('males',
                                                     'females'),
)
## graph
png('./astrocytes/limmatrend/RRHO2.astrocytes.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(astrocytes.RRHO_obj.DvsS,
              main = 'astrocytes DvsS')
dev.off()

#### compare sexes across status
##dom data
astrocytes.MvsF.D.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                 value = -log10(astrocytes.limma.results.sct.df$P.Value_MvsFD),
                                 direction = astrocytes.limma.results.sct.df$Direction.type_MvsFD, 
                                 stringsAsFactors = FALSE)
#set positive negative
astrocytes.MvsF.D.rrho2 = astrocytes.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
astrocytes.MvsF.S.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                 value = -log10(astrocytes.limma.results.sct.df$P.Value_MvsFS),
                                 direction = astrocytes.limma.results.sct.df$Direction.type_MvsFS,
                                 stringsAsFactors = FALSE)
#set positive negative
astrocytes.MvsF.S.rrho2 = astrocytes.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
astrocytes.RRHO_obj.MvsF <-  RRHO2_initialize(astrocytes.MvsF.D.rrho2, 
                                          astrocytes.MvsF.S.rrho2, 
                                          boundary = 0.1,
                                          labels = c('Doms',
                                                     'Subs')
)
## graph
png('./astrocytes/limmatrend/RRHO2.astrocytes.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(astrocytes.RRHO_obj.MvsF,
              main = 'astrocytes MvsF')
dev.off()




#### oligodendrocytes limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

#subset to oligodendrocytes
mouse.snseq.combined.sct.oligodendrocytes = subset(mouse.snseq.combined.sct,
                                             idents = c("C25-19: Oligodendrocytes"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.oligodendrocytes) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.oligodendrocytes@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 845       1300       1126        990 

### extract gene expression and orig.idents from oligodendrocytes
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.oligodendrocytes.var <- FindVariableFeatures(mouse.snseq.combined.sct.oligodendrocytes, 
                                                                assay = 'integrated',
                                                                selection.method = "vst", 
                                                                verbose = F)
# identify top  variable genes
oligodendrocytes.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.oligodendrocytes.var,
                                                               assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.oligodendrocytes.expression = full_join(full_join(mouse.snseq.combined.sct.oligodendrocytes@reductions$umap@cell.embeddings %>% 
                                                                       as.data.frame() %>% 
                                                                       rownames_to_column("Cell.id"),
                                                                     mouse.snseq.combined.sct.oligodendrocytes@meta.data %>%
                                                                       rownames_to_column("Cell.id")),
                                                           mouse.snseq.combined.sct.oligodendrocytes@assays$SCT@data %>% 
                                                             as.data.frame()  %>% 
                                                             t() %>% as.data.frame() %>% 
                                                             rownames_to_column('Cell.id'))


## create vector of factor
oligodendrocytes.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.oligodendrocytes.expression %>% 
  mutate(oligodendrocytes.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(oligodendrocytes.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
oligodendrocytes.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.oligodendrocytes,
                                                       assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.oligodendrocytes.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% oligodendrocytes.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
oligodendrocytes.reduce.vector.limma.sct.prep = list(count = oligodendrocytes.reduce.vector.count.sct.prep,
                                               condt = oligodendrocytes.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_oligodendrocytes <- function(L) {
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
  pdf(file= "./oligodendrocytes/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./oligodendrocytes/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
oligodendrocytes.limma.results.sct = run_limmatrend_oligodendrocytes(oligodendrocytes.reduce.vector.limma.sct.prep)

# save results to dataframe
oligodendrocytes.limma.results.sct.df = full_join(oligodendrocytes.limma.results.sct$ttDvsS %>% 
                                              rename_with(~paste0(.,"_DvsS")) %>% 
                                              rownames_to_column("Gene"),
                                            oligodendrocytes.limma.results.sct$ttDvsSM %>% 
                                              rename_with(~paste0(.,"_DvsS_M")) %>% 
                                              rownames_to_column("Gene")) %>% 
  full_join(oligodendrocytes.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(oligodendrocytes.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(oligodendrocytes.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(oligodendrocytes.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
oligodendrocytes.limma.results.sct.df = oligodendrocytes.limma.results.sct.df %>% 
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
  oligodendrocytes.limma.results.sct.df %>% 
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
  ggsave(paste0('./oligodendrocytes/limmatrend/limma.oligodendrocytes.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# oligodendrocytes.limma.results.sct.df.add.np = oligodendrocytes.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   oligodendrocytes.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./oligodendrocytes/limmatrend/limma.oligodendrocytes.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results oligodendrocytes limmatrend ####
save(oligodendrocytes.limma.results.sct.df,
     file = './oligodendrocytes/limmatrend/neuropeptides.oligodendrocytes.limma.results.sct.df.RData')

# load('./oligodendrocytes/limmatrend/neuropeptides.oligodendrocytes.limma.results.sct.df.RData')
#### RRHO2 oligodendrocytes ####
#### compare dom vs sub across sexes
##male data
oligodendrocytes.DvsS.M.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.DvsS.M.rrho2 = oligodendrocytes.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
oligodendrocytes.DvsS.F.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.DvsS.F.rrho2 = oligodendrocytes.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
oligodendrocytes.RRHO_obj.DvsS <-  RRHO2_initialize(oligodendrocytes.DvsS.M.rrho2, 
                                              oligodendrocytes.DvsS.F.rrho2, 
                                              boundary = 0.1,
                                              labels = c('males',
                                                         'females'),
)
## graph
png('./oligodendrocytes/limmatrend/RRHO2.oligodendrocytes.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(oligodendrocytes.RRHO_obj.DvsS,
              main = 'oligodendrocytes DvsS')
dev.off()

#### compare sexes across status
##dom data
oligodendrocytes.MvsF.D.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_MvsFD),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_MvsFD, 
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.MvsF.D.rrho2 = oligodendrocytes.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
oligodendrocytes.MvsF.S.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_MvsFS),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_MvsFS,
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.MvsF.S.rrho2 = oligodendrocytes.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
oligodendrocytes.RRHO_obj.MvsF <-  RRHO2_initialize(oligodendrocytes.MvsF.D.rrho2, 
                                              oligodendrocytes.MvsF.S.rrho2, 
                                              boundary = 0.1,
                                              labels = c('Doms',
                                                         'Subs')
)
## graph
png('./oligodendrocytes/limmatrend/RRHO2.oligodendrocytes.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(oligodendrocytes.RRHO_obj.MvsF,
              main = 'oligodendrocytes MvsF')
dev.off()




#### OPCs limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

#subset to OPCs
mouse.snseq.combined.sct.OPCs = subset(mouse.snseq.combined.sct,
                                                   idents = c("C25-20: OPC"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.OPCs) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.OPCs@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 265        249        223        255 

### extract gene expression and orig.idents from OPCs
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.OPCs.var <- FindVariableFeatures(mouse.snseq.combined.sct.OPCs, 
                                                                      assay = 'integrated',
                                                                      selection.method = "vst", 
                                                                      verbose = F)
# identify top variable genes
OPCs.reduce.group.topgenes.prep <-  VariableFeatures(mouse.snseq.combined.sct.OPCs.var,
                                                     assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.OPCs.expression = full_join(full_join(mouse.snseq.combined.sct.OPCs@reductions$umap@cell.embeddings %>% 
                                                                             as.data.frame() %>% 
                                                                             rownames_to_column("Cell.id"),
                                                                           mouse.snseq.combined.sct.OPCs@meta.data %>%
                                                                             rownames_to_column("Cell.id")),
                                                                 mouse.snseq.combined.sct.OPCs@assays$SCT@data %>% 
                                                                   as.data.frame()  %>% 
                                                                   t() %>% as.data.frame() %>% 
                                                                   rownames_to_column('Cell.id'))


## create vector of factor
OPCs.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.OPCs.expression %>% 
  mutate(OPCs.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(OPCs.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
OPCs.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.OPCs,
                                                             assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.OPCs.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% OPCs.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
OPCs.reduce.vector.limma.sct.prep = list(count = OPCs.reduce.vector.count.sct.prep,
                                                     condt = OPCs.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_OPCs <- function(L) {
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
  pdf(file= "./OPCs/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./OPCs/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
OPCs.limma.results.sct = run_limmatrend_OPCs(OPCs.reduce.vector.limma.sct.prep)

# save results to dataframe
OPCs.limma.results.sct.df = full_join(OPCs.limma.results.sct$ttDvsS %>% 
                                                    rename_with(~paste0(.,"_DvsS")) %>% 
                                                    rownames_to_column("Gene"),
                                                  OPCs.limma.results.sct$ttDvsSM %>% 
                                                    rename_with(~paste0(.,"_DvsS_M")) %>% 
                                                    rownames_to_column("Gene")) %>% 
  full_join(OPCs.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(OPCs.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(OPCs.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(OPCs.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
OPCs.limma.results.sct.df = OPCs.limma.results.sct.df %>% 
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
  OPCs.limma.results.sct.df %>% 
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
  ggsave(paste0('./OPCs/limmatrend/limma.OPCs.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# OPCs.limma.results.sct.df.add.np = OPCs.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   OPCs.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./OPCs/limmatrend/limma.OPCs.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results OPCs limmatrend ####
save(OPCs.limma.results.sct.df,
     file = './OPCs/limmatrend/neuropeptides.OPCs.limma.results.sct.df.RData')

# load('./OPCs/limmatrend/neuropeptides.OPCs.limma.results.sct.df.RData')
#### RRHO2 OPCs ####
#### compare dom vs sub across sexes
##male data
OPCs.DvsS.M.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                           value = -log10(OPCs.limma.results.sct.df$P.Value_DvsS_M),
                                           direction = OPCs.limma.results.sct.df$Direction.type_DvsS_M, 
                                           stringsAsFactors = FALSE)
#set positive negative
OPCs.DvsS.M.rrho2 = OPCs.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
OPCs.DvsS.F.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                           value = -log10(OPCs.limma.results.sct.df$P.Value_DvsS_F),
                                           direction = OPCs.limma.results.sct.df$Direction.type_DvsS_F,
                                           stringsAsFactors = FALSE)
#set positive negative
OPCs.DvsS.F.rrho2 = OPCs.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
OPCs.RRHO_obj.DvsS <-  RRHO2_initialize(OPCs.DvsS.M.rrho2, 
                                                    OPCs.DvsS.F.rrho2, 
                                                    boundary = 0.1,
                                                    labels = c('males',
                                                               'females'),
)
## graph
png('./OPCs/limmatrend/RRHO2.OPCs.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(OPCs.RRHO_obj.DvsS,
              main = 'OPCs DvsS')
dev.off()

#### compare sexes across status
##dom data
OPCs.MvsF.D.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                           value = -log10(OPCs.limma.results.sct.df$P.Value_MvsFD),
                                           direction = OPCs.limma.results.sct.df$Direction.type_MvsFD, 
                                           stringsAsFactors = FALSE)
#set positive negative
OPCs.MvsF.D.rrho2 = OPCs.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
OPCs.MvsF.S.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                           value = -log10(OPCs.limma.results.sct.df$P.Value_MvsFS),
                                           direction = OPCs.limma.results.sct.df$Direction.type_MvsFS,
                                           stringsAsFactors = FALSE)
#set positive negative
OPCs.MvsF.S.rrho2 = OPCs.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
OPCs.RRHO_obj.MvsF <-  RRHO2_initialize(OPCs.MvsF.D.rrho2, 
                                                    OPCs.MvsF.S.rrho2, 
                                                    boundary = 0.1,
                                                    labels = c('Doms',
                                                               'Subs')
)
## graph
png('./OPCs/limmatrend/RRHO2.OPCs.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(OPCs.RRHO_obj.MvsF,
              main = 'OPCs MvsF')
dev.off()




#### Immune limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

#subset to Immune
mouse.snseq.combined.sct.Immune = subset(mouse.snseq.combined.sct,
                                       idents = c("C25-21: Immune"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.Immune) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.Immune@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 300        726        482        499 

# use variable genes approx equal to the lowest cell count

### extract gene expression and orig.idents from Immune
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.Immune.var <- FindVariableFeatures(mouse.snseq.combined.sct.Immune, 
                                                          assay = 'integrated',
                                                          selection.method = "vst", 
                                                          verbose = F)
# identify top variable genes
Immune.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.Immune.var,
                                                         assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.Immune.expression = full_join(full_join(mouse.snseq.combined.sct.Immune@reductions$umap@cell.embeddings %>% 
                                                                 as.data.frame() %>% 
                                                                 rownames_to_column("Cell.id"),
                                                               mouse.snseq.combined.sct.Immune@meta.data %>%
                                                                 rownames_to_column("Cell.id")),
                                                     mouse.snseq.combined.sct.Immune@assays$SCT@data %>% 
                                                       as.data.frame()  %>% 
                                                       t() %>% as.data.frame() %>% 
                                                       rownames_to_column('Cell.id'))


## create vector of factor
Immune.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.Immune.expression %>% 
  mutate(Immune.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Immune.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Immune.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.Immune,
                                                 assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.Immune.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Immune.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Immune.reduce.vector.limma.sct.prep = list(count = Immune.reduce.vector.count.sct.prep,
                                         condt = Immune.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Immune <- function(L) {
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
  pdf(file= "./Immune/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./Immune/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
Immune.limma.results.sct = run_limmatrend_Immune(Immune.reduce.vector.limma.sct.prep)

# save results to dataframe
Immune.limma.results.sct.df = full_join(Immune.limma.results.sct$ttDvsS %>% 
                                        rename_with(~paste0(.,"_DvsS")) %>% 
                                        rownames_to_column("Gene"),
                                      Immune.limma.results.sct$ttDvsSM %>% 
                                        rename_with(~paste0(.,"_DvsS_M")) %>% 
                                        rownames_to_column("Gene")) %>% 
  full_join(Immune.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Immune.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Immune.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Immune.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Immune.limma.results.sct.df = Immune.limma.results.sct.df %>% 
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
  Immune.limma.results.sct.df %>% 
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
  ggsave(paste0('./Immune/limmatrend/limma.Immune.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# Immune.limma.results.sct.df.add.np = Immune.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   Immune.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./Immune/limmatrend/limma.Immune.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results Immune limmatrend ####
save(Immune.limma.results.sct.df,
     file = './Immune/limmatrend/neuropeptides.Immune.limma.results.sct.df.RData')

# load('./Immune/limmatrend/neuropeptides.Immune.limma.results.sct.df.RData')
#### RRHO2 Immune ####
#### compare dom vs sub across sexes
##male data
Immune.DvsS.M.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                               value = -log10(Immune.limma.results.sct.df$P.Value_DvsS_M),
                               direction = Immune.limma.results.sct.df$Direction.type_DvsS_M, 
                               stringsAsFactors = FALSE)
#set positive negative
Immune.DvsS.M.rrho2 = Immune.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Immune.DvsS.F.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                               value = -log10(Immune.limma.results.sct.df$P.Value_DvsS_F),
                               direction = Immune.limma.results.sct.df$Direction.type_DvsS_F,
                               stringsAsFactors = FALSE)
#set positive negative
Immune.DvsS.F.rrho2 = Immune.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Immune.RRHO_obj.DvsS <-  RRHO2_initialize(Immune.DvsS.M.rrho2, 
                                        Immune.DvsS.F.rrho2, 
                                        boundary = 0.1,
                                        labels = c('males',
                                                   'females'),
)
## graph
png('./Immune/limmatrend/RRHO2.Immune.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Immune.RRHO_obj.DvsS,
              main = 'Immune DvsS')
dev.off()

#### compare sexes across status
##dom data
Immune.MvsF.D.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                               value = -log10(Immune.limma.results.sct.df$P.Value_MvsFD),
                               direction = Immune.limma.results.sct.df$Direction.type_MvsFD, 
                               stringsAsFactors = FALSE)
#set positive negative
Immune.MvsF.D.rrho2 = Immune.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
Immune.MvsF.S.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                               value = -log10(Immune.limma.results.sct.df$P.Value_MvsFS),
                               direction = Immune.limma.results.sct.df$Direction.type_MvsFS,
                               stringsAsFactors = FALSE)
#set positive negative
Immune.MvsF.S.rrho2 = Immune.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Immune.RRHO_obj.MvsF <-  RRHO2_initialize(Immune.MvsF.D.rrho2, 
                                        Immune.MvsF.S.rrho2, 
                                        boundary = 0.1,
                                        labels = c('Doms',
                                                   'Subs')
)
## graph
png('./Immune/limmatrend/RRHO2.Immune.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Immune.RRHO_obj.MvsF,
              main = 'Immune MvsF')
dev.off()




#### Endothelial limmatrend ####
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

#subset to Endothelial
mouse.snseq.combined.sct.Endothelial = subset(mouse.snseq.combined.sct,
                                         idents = c("C25-24: Mural+Endothelial"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.Endothelial) = 'SCT'

#check number of cells 
mouse.snseq.combined.sct.Endothelial@meta.data %>% 
  dplyr::select(orig.ident) %>% 
  table()

# orig.ident
# Female.Dom Female.Sub   Male.Dom   Male.Sub 
# 393        521        521        432 

# use variable genes approx equal to the lowest cell count

### extract gene expression and orig.idents from Endothelial
## create expression dataframe of neuropeptides
### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.Endothelial.var <- FindVariableFeatures(mouse.snseq.combined.sct.Endothelial, 
                                                            assay = 'integrated',
                                                            selection.method = "vst", 
                                                            verbose = F)
# identify top variable genes
Endothelial.reduce.group.topgenes.prep <- VariableFeatures(mouse.snseq.combined.sct.Endothelial.var,
                                                           assay = 'integrated')

# create dummy
mouse.snseq.combined.sct.Endothelial.expression = full_join(full_join(mouse.snseq.combined.sct.Endothelial@reductions$umap@cell.embeddings %>% 
                                                                   as.data.frame() %>% 
                                                                   rownames_to_column("Cell.id"),
                                                                 mouse.snseq.combined.sct.Endothelial@meta.data %>%
                                                                   rownames_to_column("Cell.id")),
                                                       mouse.snseq.combined.sct.Endothelial@assays$SCT@data %>% 
                                                         as.data.frame()  %>% 
                                                         t() %>% as.data.frame() %>% 
                                                         rownames_to_column('Cell.id'))


## create vector of factor
Endothelial.reduce.DomvsSub.vector.list.sct.prep = mouse.snseq.combined.sct.Endothelial.expression %>% 
  mutate(Endothelial.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(Endothelial.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
Endothelial.reduce.vector.count.sct.prep = GetAssayData(mouse.snseq.combined.sct.Endothelial,
                                                   assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  mouse.snseq.combined.sct.Endothelial.expression %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% Endothelial.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

#create list
Endothelial.reduce.vector.limma.sct.prep = list(count = Endothelial.reduce.vector.count.sct.prep,
                                           condt = Endothelial.reduce.DomvsSub.vector.list.sct.prep)

# [1] Male.Dom   Female.Dom
# [3] Male.Sub   Female.Sub

#create function
run_limmatrend_Endothelial <- function(L) {
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
  pdf(file= "./Endothelial/limmatrend/limmatrend.histograms.pdf" )
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
  # pdf(file= "./Endothelial/neuropeptides/limmatrend/limmatrend.MDS.pdf" )
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
Endothelial.limma.results.sct = run_limmatrend_Endothelial(Endothelial.reduce.vector.limma.sct.prep)

# save results to dataframe
Endothelial.limma.results.sct.df = full_join(Endothelial.limma.results.sct$ttDvsS %>% 
                                          rename_with(~paste0(.,"_DvsS")) %>% 
                                          rownames_to_column("Gene"),
                                        Endothelial.limma.results.sct$ttDvsSM %>% 
                                          rename_with(~paste0(.,"_DvsS_M")) %>% 
                                          rownames_to_column("Gene")) %>% 
  full_join(Endothelial.limma.results.sct$ttDvsSF %>% 
              rename_with(~paste0(.,"_DvsS_F")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Endothelial.limma.results.sct$ttMvsF %>% 
              rename_with(~paste0(.,"_MvsF")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Endothelial.limma.results.sct$ttMvsFD %>% 
              rename_with(~paste0(.,"_MvsFD")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(Endothelial.limma.results.sct$ttMvsFS %>% 
              rename_with(~paste0(.,"_MvsFS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
Endothelial.limma.results.sct.df = Endothelial.limma.results.sct.df %>% 
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
  Endothelial.limma.results.sct.df %>% 
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
  ggsave(paste0('./Endothelial/limmatrend/limma.Endothelial.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}
# 
# # highlight neuropeptides
# Endothelial.limma.results.sct.df.add.np = Endothelial.limma.results.sct.df %>% 
#   mutate(GENE = toupper(Gene)) %>% 
#   left_join(neuropeptides.genes %>% 
#               as.data.frame() %>% 
#               rename('GENE' = '.') %>% 
#               mutate(Neuropeptide = 1)) %>% 
#   mutate(Neuropeptide = ifelse(is.na(Neuropeptide),
#                                0,
#                                Neuropeptide))
# 
# 
# for (i in limma.genotpye.vector) {
#   # graph volcano plot
#   Endothelial.limma.results.sct.df.add.np %>% 
#     mutate(sig.label = ifelse(Neuropeptide == 1,
#                               GENE,
#                               '')) %>% 
#     ggplot(aes(x = get(paste0("logFC_", i)),
#                y = -log10(get(paste0("adj.P.Val_", i))),
#                color = factor(Neuropeptide)))+
#     geom_hline(yintercept = -log10(0.05),
#                linetype = 'dotted') +
#     geom_point(size = 5) +
#     theme_classic() + 
#     scale_color_manual(values=c("grey","red")) +
#     geom_text(aes(label = sig.label),
#               vjust = 0, 
#               nudge_y = 0.10,
#               size = 5) +
#     theme(text = element_text(size = 20),
#           legend.position = 'none') +
#     xlab(paste0("logFC_", i)) +
#     ylab( paste0("-log10(adj.P.Val_", i,")")) +
#     ggtitle(paste0(i, " volcano plot"))
#   ggsave(paste0('./Endothelial/limmatrend/limma.Endothelial.reduce.volcano.', i, '.neuropeptides.png'),
#          width = 5,
#          height = 5)
# }

#### save results Endothelial limmatrend ####
save(Endothelial.limma.results.sct.df,
     file = './Endothelial/limmatrend/neuropeptides.Endothelial.limma.results.sct.df.RData')

# load('./Endothelial/limmatrend/neuropeptides.Endothelial.limma.results.sct.df.RData')
#### RRHO2 Endothelial ####
#### compare dom vs sub across sexes
##male data
Endothelial.DvsS.M.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                 value = -log10(Endothelial.limma.results.sct.df$P.Value_DvsS_M),
                                 direction = Endothelial.limma.results.sct.df$Direction.type_DvsS_M, 
                                 stringsAsFactors = FALSE)
#set positive negative
Endothelial.DvsS.M.rrho2 = Endothelial.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Endothelial.DvsS.F.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                 value = -log10(Endothelial.limma.results.sct.df$P.Value_DvsS_F),
                                 direction = Endothelial.limma.results.sct.df$Direction.type_DvsS_F,
                                 stringsAsFactors = FALSE)
#set positive negative
Endothelial.DvsS.F.rrho2 = Endothelial.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Endothelial.RRHO_obj.DvsS <-  RRHO2_initialize(Endothelial.DvsS.M.rrho2, 
                                          Endothelial.DvsS.F.rrho2, 
                                          boundary = 0.1,
                                          labels = c('males',
                                                     'females'),
)
## graph
png('./Endothelial/limmatrend/RRHO2.Endothelial.DvsS.sexes.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Endothelial.RRHO_obj.DvsS,
              main = 'Endothelial DvsS')
dev.off()

#### compare sexes across status
##dom data
Endothelial.MvsF.D.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                 value = -log10(Endothelial.limma.results.sct.df$P.Value_MvsFD),
                                 direction = Endothelial.limma.results.sct.df$Direction.type_MvsFD, 
                                 stringsAsFactors = FALSE)
#set positive negative
Endothelial.MvsF.D.rrho2 = Endothelial.MvsF.D.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##sub data
Endothelial.MvsF.S.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                 value = -log10(Endothelial.limma.results.sct.df$P.Value_MvsFS),
                                 direction = Endothelial.limma.results.sct.df$Direction.type_MvsFS,
                                 stringsAsFactors = FALSE)
#set positive negative
Endothelial.MvsF.S.rrho2 = Endothelial.MvsF.S.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))
## create rrho2 object
Endothelial.RRHO_obj.MvsF <-  RRHO2_initialize(Endothelial.MvsF.D.rrho2, 
                                          Endothelial.MvsF.S.rrho2, 
                                          boundary = 0.1,
                                          labels = c('Doms',
                                                     'Subs')
)
## graph
png('./Endothelial/limmatrend/RRHO2.Endothelial.MvsF.status.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 300)
RRHO2_heatmap(Endothelial.RRHO_obj.MvsF,
              main = 'Endothelial MvsF')
dev.off()




#### Get GO terms ####
## load ensembl
library(biomaRt)
###selecting biomart database
ensembl = useMart('ensembl')
# #list datasets
# dataset = listDatasets(ensembl)
##select dataset
#mouse
ensembl.mouse = useDataset('mmusculus_gene_ensembl',
                             mart=ensembl)

# # get list of all attributes
# listAttributes(ensembl.mouse) %>% View


#create attributes lists
mouse.attributes = c('external_gene_name',
                       'ensembl_gene_id',
                       'go_id')

# create gene list
mouse.gene.list = c(neuron.limma.results.sct.df$Gene,
                    astrocytes.limma.results.sct.df$Gene,
                    oligodendrocytes.limma.results.sct.df$Gene,
                    OPCs.limma.results.sct.df$Gene,
                    Immune.limma.results.sct.df$Gene,
                    Endothelial.limma.results.sct.df$Gene) %>% 
  unique()

##identify GO terms for WGCNA mouse genes
# use gene names
mouse.wgcna.go.terms.gene = getBM(attributes = mouse.attributes,
                                    mart = ensembl.mouse,
                                    values = mouse.gene.list,
                                    filter = 'external_gene_name',
                                    useCache = FALSE) # useCache has to do with version of R not being up to date?
# use ensembl gene IDs

# single gene column 
mouse.wgcna.go.terms = mouse.wgcna.go.terms.gene %>% 
  mutate(gene_name = ifelse(is.na(external_gene_name),
                            ensembl_gene_id,
                            external_gene_name)) %>% 
  dplyr::select(c(gene_name,
                  go_id))

#check length
#25030
mouse.wgcna.go.terms %>%
  nrow()

#check number of mouse genes?
#1295
mouse.wgcna.go.terms %>%
  pull(gene_name) %>%
  unique() %>%
  length()

#check number of GO IDs?
#5913
mouse.wgcna.go.terms %>%
  pull(go_id) %>%
  unique() %>%
  length()

## check duplicates
# 25030 - 25012 = 18 duplicates
mouse.wgcna.go.terms %>%
  distinct() %>% 
  nrow()

# remove duplicates
mouse.wgcna.go.terms = mouse.wgcna.go.terms %>% 
  distinct()

## collapse go_id into gene list
mouse.wgcna.go.terms = mouse.wgcna.go.terms %>%
  group_by(gene_name) %>%
  summarize(go_id = str_c(go_id, collapse = ";"))

# add 'unknown' 
mouse.wgcna.go.terms = mouse.wgcna.go.terms %>% 
  mutate(go_id = ifelse(go_id == '',
                        'unknown',
                        go_id))

# save to csv 
write.csv(mouse.wgcna.go.terms,
          file = './GO_terms/mouse.wgcna.go.terms.csv')
# save to tab delimited with no column names
write_tsv(mouse.wgcna.go.terms,
          file = './GO_terms/mouse.wgcna.go.terms.tsv',
          col_names = FALSE)

#### prepare for GO_MWU ####
### create data input list for each cell type 
# save files
#neuron
write.csv(neuron.limma.results.sct.df,
          file ='./GO_terms/neuron.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)

#astrocytes
write.csv(astrocytes.limma.results.sct.df,
          file ='./GO_terms/astrocytes.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)

#oligodendrocytes
write.csv(oligodendrocytes.limma.results.sct.df,
          file ='./GO_terms/oligodendrocytes.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)

#OPCs
write.csv(OPCs.limma.results.sct.df,
          file ='./GO_terms/OPCs.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)

#Immune
write.csv(Immune.limma.results.sct.df,
          file ='./GO_terms/Immune.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)

#Endothelial
write.csv(Endothelial.limma.results.sct.df,
          file ='./GO_terms/Endothelial.limma.results.sct.df.csv',
          row.names = FALSE,
          quote = FALSE)




## need wgcna.DMEs.neurons.log2FC.csv and mouse.wgcna.go.terms.tsv for MWU_GO
# need to run on desktop 
# https://github.com/z0on/GO_MWU
# https://github.com/schmidte10/GO-MWU-automation-

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


#### load data for RedRibbon ####
#### RedRibbon neurons 
load('./neurons/limmatrend/neuropeptides.neuron.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
neuron.DvsS.M.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                     value = -log10(neuron.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = neuron.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
neuron.DvsS.M.rrho2 = neuron.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
neuron.DvsS.F.rrho2 = data.frame(gene = neuron.limma.results.sct.df$Gene,
                                     value = -log10(neuron.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = neuron.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
neuron.DvsS.F.rrho2 = neuron.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

### RedRibbon astrocytes 
load('./astrocytes/limmatrend/neuropeptides.astrocytes.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
astrocytes.DvsS.M.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                     value = -log10(astrocytes.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = astrocytes.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
astrocytes.DvsS.M.rrho2 = astrocytes.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
astrocytes.DvsS.F.rrho2 = data.frame(gene = astrocytes.limma.results.sct.df$Gene,
                                     value = -log10(astrocytes.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = astrocytes.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
astrocytes.DvsS.F.rrho2 = astrocytes.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

### RedRibbon oligodendrocytes 
load('./oligodendrocytes/limmatrend/neuropeptides.oligodendrocytes.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
oligodendrocytes.DvsS.M.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.DvsS.M.rrho2 = oligodendrocytes.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
oligodendrocytes.DvsS.F.rrho2 = data.frame(gene = oligodendrocytes.limma.results.sct.df$Gene,
                                     value = -log10(oligodendrocytes.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = oligodendrocytes.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
oligodendrocytes.DvsS.F.rrho2 = oligodendrocytes.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

### RedRibbon OPCs 
load('./OPCs/limmatrend/neuropeptides.OPCs.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
OPCs.DvsS.M.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                     value = -log10(OPCs.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = OPCs.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
OPCs.DvsS.M.rrho2 = OPCs.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
OPCs.DvsS.F.rrho2 = data.frame(gene = OPCs.limma.results.sct.df$Gene,
                                     value = -log10(OPCs.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = OPCs.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
OPCs.DvsS.F.rrho2 = OPCs.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))


### RedRibbon Immune 
load('./Immune/limmatrend/neuropeptides.Immune.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
Immune.DvsS.M.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                                     value = -log10(Immune.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = Immune.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
Immune.DvsS.M.rrho2 = Immune.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Immune.DvsS.F.rrho2 = data.frame(gene = Immune.limma.results.sct.df$Gene,
                                     value = -log10(Immune.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = Immune.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
Immune.DvsS.F.rrho2 = Immune.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))


### RedRibbon Endothelial 
load('./Endothelial/limmatrend/neuropeptides.Endothelial.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
Endothelial.DvsS.M.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                     value = -log10(Endothelial.limma.results.sct.df$P.Value_DvsS_M),
                                     direction = Endothelial.limma.results.sct.df$Direction.type_DvsS_M, 
                                     stringsAsFactors = FALSE)
#set positive negative
Endothelial.DvsS.M.rrho2 = Endothelial.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
Endothelial.DvsS.F.rrho2 = data.frame(gene = Endothelial.limma.results.sct.df$Gene,
                                     value = -log10(Endothelial.limma.results.sct.df$P.Value_DvsS_F),
                                     direction = Endothelial.limma.results.sct.df$Direction.type_DvsS_F,
                                     stringsAsFactors = FALSE)
#set positive negative
Endothelial.DvsS.F.rrho2 = Endothelial.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

### RedRibbon GABA 
load('./GABA/limmatrend/neuropeptides.GABA.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
GABA.DvsS.M.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                                      value = -log10(GABA.limma.results.sct.df$P.Value_DvsS_M),
                                      direction = GABA.limma.results.sct.df$Direction.type_DvsS_M, 
                                      stringsAsFactors = FALSE)
#set positive negative
GABA.DvsS.M.rrho2 = GABA.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
GABA.DvsS.F.rrho2 = data.frame(gene = GABA.limma.results.sct.df$Gene,
                                      value = -log10(GABA.limma.results.sct.df$P.Value_DvsS_F),
                                      direction = GABA.limma.results.sct.df$Direction.type_DvsS_F,
                                      stringsAsFactors = FALSE)
#set positive negative
GABA.DvsS.F.rrho2 = GABA.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

### RedRibbon GLU 
load('./GLU/limmatrend/neuropeptides.GLU.limma.results.sct.df.RData')
### compare dom vs sub across sexes
##male data
GLU.DvsS.M.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                               value = -log10(GLU.limma.results.sct.df$P.Value_DvsS_M),
                               direction = GLU.limma.results.sct.df$Direction.type_DvsS_M, 
                               stringsAsFactors = FALSE)
#set positive negative
GLU.DvsS.M.rrho2 = GLU.DvsS.M.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))

##female data
GLU.DvsS.F.rrho2 = data.frame(gene = GLU.limma.results.sct.df$Gene,
                               value = -log10(GLU.limma.results.sct.df$P.Value_DvsS_F),
                               direction = GLU.limma.results.sct.df$Direction.type_DvsS_F,
                               stringsAsFactors = FALSE)
#set positive negative
GLU.DvsS.F.rrho2 = GLU.DvsS.F.rrho2 %>% 
  mutate(value = ifelse(direction == 'down',
                        -1*value,
                        value)) %>% 
  dplyr::select(-c(direction))


#### RedRibbon all cell types ####
## neuron
RedRibbon.all(celltype = "neuron",
              dataset.a = neuron.DvsS.M.rrho2,
              dataset.b = neuron.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## astrocytes
RedRibbon.all(celltype = "astrocytes",
              dataset.a = astrocytes.DvsS.M.rrho2,
              dataset.b = astrocytes.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## oligodendrocytes
RedRibbon.all(celltype = "oligodendrocytes",
              dataset.a = oligodendrocytes.DvsS.M.rrho2,
              dataset.b = oligodendrocytes.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## OPCs
RedRibbon.all(celltype = "OPCs",
              dataset.a = OPCs.DvsS.M.rrho2,
              dataset.b = OPCs.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## Immune
RedRibbon.all(celltype = "Immune",
              dataset.a = Immune.DvsS.M.rrho2,
              dataset.b = Immune.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## Endothelial
RedRibbon.all(celltype = "Endothelial",
              dataset.a = Endothelial.DvsS.M.rrho2,
              dataset.b = Endothelial.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## GABA
RedRibbon.all(celltype = "GABA",
              dataset.a = GABA.DvsS.M.rrho2,
              dataset.b = GABA.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')

## GLU
RedRibbon.all(celltype = "GLU",
              dataset.a = GLU.DvsS.M.rrho2,
              dataset.b = GLU.DvsS.F.rrho2,
              dataset.a.type = "male",
              dataset.b.type = "female",
              a.variable = "value",
              b.variable = "value",
              new.max.log = NULL,
              file.name = 'RedRibbon/')





