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
  
# graph count of genes per quadrant
  plot.RRHO.counts(rrho_results, "./neurons/all_neurons/RRHO")

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

  plot.overlaps(rrho_results,
                neuron_clusterMarkers,
                quadrants_to_check = NULL,
                subtitle = "adj. p-value < 0.05; specificity > 0.5",
                group.by = "cluster",
                outdir = "./neurons/all_neurons/RRHO")
  


