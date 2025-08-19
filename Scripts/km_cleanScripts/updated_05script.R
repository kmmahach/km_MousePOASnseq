# R-4.3.1, Seurat v4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# libraries 
load_packages(c("Seurat", "RRHO2", "DEsingle", "dendextend", "RedRibbon", "pheatmap", "limma",
                "ggalluvial", "multcomp", "clustree", "patchwork", "tidyverse", "edgeR", 
                "enrichplot", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi"),
              out_prefix = "05")

#### load data ####
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# subset with SCT data
DefaultAssay(int.ldfs) = "SCT"


#### DGE with limma - all neurons #### 
setwd(paste0(root.dir, "/DGE_CellTypes"))

# raw counts and norm.counts of variable genes
dge_data <- prep.for.DGE(int.ldfs, 
                         selection.method = "vst", 
                         SCTransformed = TRUE, 
                         assay = 'integrated')

names(dge_data[[1]]) <- "integrated_neuronSeurat"

# get results and graph p-values
limma_results <- run_limmatrend(dge_data$results, "./neurons/all_neurons/limma_trend")

  save(limma_results, file = "./neurons/all_neurons/limma_trend/limma_results.rda")

  
#### RRHO/RedRibbon - all neurons ####
# make sure these match upper/lower case
  sex = c("Female", "Male");
  status = c("Dom", "Sub");
  
rrho_results <- get.RRHO(limma_results,
                         group.by = sex,
                         compare.across = status,
                         outdir = "./neurons/all_neurons/RRHO")

  save(rrho_results, file = "./neurons/all_neurons/RRHO/rrho_results.rda")
  
# graph count of genes per quadrant
  plot.RRHO.counts(rrho_results, "./neurons/all_neurons/RRHO")

  
#### Compare concordant genes to neuron cluster markers ####
  # cluster marker genes from updated_03script.R
setwd(paste0(root.dir, "/DGE_CellTypes"))
  
neuron.markers <- read.csv("./neurons/cluster_stats/neuron_clusterMarkers.csv")

# check % overlaps with uu and dd gene lists for RRHO2 & RedRibbon
neuron.markers %>% 
  filter(p_val_adj <= 0.05, 
         abs(specificity) >= 0.5) -> neuron_clusterMarkers

check.quadrants(rrho_results, 
                neuron_clusterMarkers, 
                quadrants_to_check = NULL, 
                outdir = "./neurons/all_neurons/RRHO")

  plot.overlaps(rrho_results,
                neuron_clusterMarkers,
                quadrants_to_check = NULL,
                subtitle = "adj. p-value < 0.05; specificity > 0.5",
                group.by = "cluster",
                outdir = "./neurons/all_neurons/RRHO")


#### GO Enrichment Analysis for upup genes (all neurons) ####
# note: top 15 GO terms are the same for RedRibbon & RRHO2 gene lists. 
# plots are identical

# RRHO2 list 
# get list of upup genes to run GO analysis
neuron.RRHO.DvsS.upup.genes = rrho_results$integrated_neuronSeurat$genelist_uu$gene_list_overlap_uu

upupGO <- enrichGO(gene = unique(neuron.RRHO.DvsS.upup.genes), 
                   OrgDb = "org.Mm.eg.db",
                   keyType = "SYMBOL", 
                   ont = "BP")


GO15.plot <- plot(barplot(upupGO,
                          showCategory = 15))

png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RRHO2_results_upup.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(GO15.plot)
dev.off()

# graph tree
png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RRHO2_results_upup_tree.png", 
    res = 300, 
    width = 13.5, 
    height = 10,
    units = 'in')

upupGO %>% 
  pairwise_termsim() %>%  # note - warning does not affect results:
  treeplot()              # ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned. 
dev.off()

# RedRibbon list 
# get list of upup genes to run GO analysis
neuron.RR.DvsS.upup.genes = rrho_results$integrated_neuronSeurat$df[gene_lists$upup$positions,1]

upupGO <- enrichGO(gene = unique(neuron.RR.DvsS.upup.genes), 
                   OrgDb = "org.Mm.eg.db",
                   keyType = "SYMBOL", 
                   ont = "BP")


GO15.plot <- plot(barplot(upupGO,
                          showCategory = 15))

png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RedRibbon_results_upup.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(GO15.plot)
dev.off()

# graph tree
png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RedRibbon_results_upup_tree.png", 
    res = 300, 
    width = 13.5, 
    height = 10,
    units = 'in')

upupGO %>% 
  pairwise_termsim() %>%  # note - warning does not affect results:
  treeplot()              # ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned. 
dev.off()

#### GO Enrichment Analysis for downdown genes (all neurons) ####
# note: top 15 GO terms are the same for RedRibbon & RRHO2 gene lists. 
# plots are identical

# RRHO2 list 
# get list of downdown genes to run GO analysis
neuron.RRHO.DvsS.downdown.genes = rrho_results$integrated_neuronSeurat$genelist_dd$gene_list_overlap_dd

downdownGO <- enrichGO(gene = unique(neuron.RRHO.DvsS.downdown.genes), 
                   OrgDb = "org.Mm.eg.db",
                   keyType = "SYMBOL", 
                   ont = "BP")


GO15.plot <- plot(barplot(downdownGO,
                          showCategory = 15))

png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RRHO2_results_downdown.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(GO15.plot)
dev.off()

# graph tree
png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RRHO2_results_downdown_tree.png", 
    res = 300, 
    width = 13.5, 
    height = 10,
    units = 'in')

downdownGO %>% 
  pairwise_termsim() %>%  # note - warning does not affect results:
  treeplot()              # ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned. 
dev.off()

# RedRibbon list 
# get list of downdown genes to run GO analysis
neuron.RR.DvsS.downdown.genes = rrho_results$integrated_neuronSeurat$df[gene_lists$downdown$positions,1]

downdownGO <- enrichGO(gene = unique(neuron.RR.DvsS.downdown.genes), 
                   OrgDb = "org.Mm.eg.db",
                   keyType = "SYMBOL", 
                   ont = "BP")


GO15.plot <- plot(barplot(downdownGO,
                          showCategory = 15))

png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RedRibbon_results_downdown.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(GO15.plot)
dev.off()

# graph tree
png("neurons/all_neurons/RRHO/GO_enrich/GOneuron_RedRibbon_results_downdown_tree.png", 
    res = 300, 
    width = 13.5, 
    height = 10,
    units = 'in')

downdownGO %>% 
  pairwise_termsim() %>%  # note - warning does not affect results:
  treeplot()              # ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned. 
dev.off()

#### analyses for glut & gaba neurons separately in updated_5.5script.R