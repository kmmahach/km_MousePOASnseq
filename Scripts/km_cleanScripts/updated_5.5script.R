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
              out_prefix = "5.5")

#### load data ####
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# subset with SCT data
DefaultAssay(int.ldfs) = "SCT"

# subset Seurat obj
subset_idents = c("C7-2: GABA", "C7-1: GLU")

l.dfs <- subset_by_ident(int.ldfs, 
                         subset_idents = subset_idents,
                         cluster = FALSE)


#### DGE with limma - glut & gaba neurons #### 
setwd(paste0(root.dir, "/DGE_CellTypes"))

# raw counts and norm.counts of variable genes
dge_data <- prep.for.DGE(l.dfs, 
                         selection.method = "vst", 
                         SCTransformed = TRUE, 
                         assay = 'integrated')

names(dge_data$results) <- lapply(names(dge_data$results), \(nm) { 
  nm = gsub("^.{6}", "", nm) })

# get results and graph p-values
limma_results <- run_limmatrend(dge_data$results, "./neurons/all_neurons/limma_trend")

save(limma_results, file = "./neurons/all_neurons/limma_trend/GLU_GABA.limma_results.rda")

#### RRHO/RedRibbon - glut & gaba neurons ####
# make sure these match upper/lower case
sex = c("Female", "Male");
status = c("Dom", "Sub");

rrho_results <- get.RRHO(limma_results,
                         group.by = sex,
                         compare.across = status,
                         outdir = "./neurons/all_neurons/RRHO")

save(rrho_results, file = "./neurons/all_neurons/RRHO/GLU_GABA.rrho_results.rda")

#### Compare concordant genes to neuron cluster markers - GLU ####
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

# graph count of genes per quadrant (RRHO2)
gene_lists <- rrho_results$GLU[grep("^genelist", names(rrho_results$GLU))]

overlap_df <- map_dfr(names(gene_lists), function(nm) {
  overlap_name <- paste0("gene_list_overlap_", str_remove(nm, "genelist_"))
  tibble(
    list_id = str_remove(nm, "genelist_"),
    overlap_length = length(gene_lists[[nm]][[overlap_name]])
  )
})

ggplot(overlap_df, aes(x = reorder(list_id, -overlap_length), 
                       y = overlap_length)) +
  geom_point(size = 3) +
  geom_text(aes(label = overlap_length), 
            vjust = -1, 
            size = 3) +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RRHO2 GLU quadrant results')


ggsave('neurons/all_neurons/RRHO/GLU_RRHO2_number_gene_quad.png',
       width = 5, height = 5)


# graph count of genes per quadrant (RedRibbon)
gene_lists <- rrho_results$GLU$RedRibbon.quads[names(rrho_results$GLU$RedRibbon.quads)]

overlap_df <- map_dfr(names(gene_lists), function(nm) {
  overlap_name <- paste(nm)
  tibble(
    list_id = overlap_name,
    overlap_length = length(rrho_results$GLU$df[gene_lists[[overlap_name]]$positions,1])
  )
})

ggplot(overlap_df, aes(x = reorder(list_id, -overlap_length), 
                       y = overlap_length)) +
  geom_point(size = 3) +
  geom_text(aes(label = overlap_length), 
            vjust = -1, 
            size = 3) +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RedRibbon GLU quadrant results')

ggsave('neurons/all_neurons/RRHO/GLU_RedRibbon_number_gene_quad.png',
       width = 5, height = 5)


#### Compare concordant genes to neuron cluster markers - GABA ####
# cluster marker genes from updated_03script.R
setwd(paste0(root.dir, "/DGE_CellTypes"))

neuron.markers <- read.csv("./neurons/cluster_stats/neuron_clusterMarkers.csv")

# check % overlaps with uu and dd gene lists for RRHO2
100*(neuron.markers %>% 
       filter(p_val_adj <= 0.05,
              abs(specificity) >= 0.5) %>% 
       pull(gene) %>% 
       intersect(rrho_results$GABA$genelist_uu$gene_list_overlap_uu) %>% 
       length())/length(rrho_results$GABA$genelist_uu$gene_list_overlap_uu)

# [1] 

100*(neuron.markers %>% 
       filter(p_val_adj <= 0.05,
              abs(specificity) >= 0.5) %>% 
       pull(gene) %>% 
       intersect(rrho_results$GABA$genelist_dd$gene_list_overlap_dd) %>% 
       length())/length(rrho_results$GABA$genelist_dd$gene_list_overlap_dd)

# [1] 

# overlaps: % of marker genes in concordance (uu + dd) 
# 

# check % overlaps with uu and dd gene lists for RedRibbon 
100*(neuron.markers %>% 
       filter(p_val_adj <= 0.05,
              abs(specificity) >= 0.5) %>% 
       pull(gene) %>% 
       intersect(rrho_results$GABA$df[rrho_results$GABA$RedRibbon.quads$upup$positions,1]) %>% 
       length())/length(rrho_results$GABA$df[rrho_results$GABA$RedRibbon.quads$upup$positions,1])

# [1] 

100*(neuron.markers %>% 
       filter(p_val_adj <= 0.05,
              abs(specificity) >= 0.5) %>% 
       pull(gene) %>% 
       intersect(rrho_results$GABA$df[rrho_results$GABA$RedRibbon.quads$downdown$positions,1]) %>% 
       length())/length(rrho_results$GABA$df[rrho_results$GABA$RedRibbon.quads$downdown$positions,1])

# [1] 

# overlaps: % of marker genes in concordance (uu + dd) 
# 

# graph count of genes per quadrant (RRHO2)
gene_lists <- rrho_results$GABA[grep("^genelist", names(rrho_results$GABA))]

overlap_df <- map_dfr(names(gene_lists), function(nm) {
  overlap_name <- paste0("gene_list_overlap_", str_remove(nm, "genelist_"))
  tibble(
    list_id = str_remove(nm, "genelist_"),
    overlap_length = length(gene_lists[[nm]][[overlap_name]])
  )
})

ggplot(overlap_df, aes(x = reorder(list_id, -overlap_length), 
                       y = overlap_length)) +
  geom_point(size = 3) +
  geom_text(aes(label = overlap_length), 
            vjust = -1, 
            size = 3) +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RRHO2 GABA quadrant results')


ggsave('neurons/all_neurons/RRHO/GABA_RRHO2_number_gene_quad.png',
       width = 5, height = 5)


# graph count of genes per quadrant (RedRibbon)
gene_lists <- rrho_results$GABA$RedRibbon.quads[names(rrho_results$GABA$RedRibbon.quads)]

overlap_df <- map_dfr(names(gene_lists), function(nm) {
  overlap_name <- paste(nm)
  tibble(
    list_id = overlap_name,
    overlap_length = length(rrho_results$GABA$df[gene_lists[[overlap_name]]$positions,1])
  )
})

ggplot(overlap_df, aes(x = reorder(list_id, -overlap_length), 
                       y = overlap_length)) +
  geom_point(size = 3) +
  geom_text(aes(label = overlap_length), 
            vjust = -1, 
            size = 3) +
  theme_classic() +
  ylab('Number of genes per quadrant') +
  xlab('') +
  ggtitle('RedRibbon GABA quadrant results')

ggsave('neurons/all_neurons/RRHO/GABA_RedRibbon_number_gene_quad.png',
       width = 5, height = 5)


