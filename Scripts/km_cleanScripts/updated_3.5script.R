# R-4.3.1, Seurat v.4.4.0

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
set.seed(12345)
compression = "xz" # slower, but usually smallest compression
setwd(root.dir) 

# functions
source("./Scripts/km_cleanScripts/functions/network_graph_fun.R")

# set up python dependencies for hdWGCNA install (with miniconda3, will install automatically or update if already installed)
system2("bash", "./hdWGCNA/activate_hdWGCNA_env.sh")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "Seurat", "WGCNA", "hdWGCNA", "limma", "edgeR", "reticulate",
                "igraph", "GeneOverlap", "RedRibbon", "corrplot", "ggpattern",
                "emmeans", "multcomp", "cowplot", "ggrepel", "patchwork", "gghalves"), 
              out_prefix = "3.5", folder = "./Scripts/km_cleanScripts")

use_condaenv("hdwgcna_env", required = TRUE)
py_config()

# load data
load(paste0(root.dir, "/Scripts/km_cleanScripts/data/integrated_seurat_onlyNeurons.rda"))

#### Neurons only WGCNA ####
# note - this analysis was not very informative; only produced 3 modules
# possibly change network construction params for smaller modules? 
setwd(paste0(root.dir, "/hdWGCNA"))

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# create new variable features list for neuron data
MSCneurons.reclustNET <- FindVariableFeatures(MSCneurons.reclust,
                                              selection.method = "vst", 
                                              nfeatures = 3000, 
                                              verbose = FALSE,
                                              assay = 'integrated')

# gets only 2240 genes
# set up object for analysis
MSCneurons.reclustNET <- SetupForWGCNA(MSCneurons.reclustNET,
                                       features = VariableFeatures(MSCneurons.reclustNET),
                                       wgcna_name = "neurons")


# construct metacells  in each group
MSCneurons.reclustNET <- MetacellsByGroups(seurat_obj = MSCneurons.reclustNET, 
                                           group.by = c("indiv_genotype", "orig.ident"), 
                                           k = 25, # nearest-neighbors parameter
                                           max_shared = 10, # max number of shared nuclei between 2 metacells
                                           ident.group = 'orig.ident',
                                           slot = "data",
                                           assay = "SCT", 
                                           mode = "average", # metacell expr profile determined by avg of constituent nuclei
                                           min_cells = 50) # set min number of cells per metacell

# set up expr matrix
MSCneurons.reclustNET <- SetDatExpr(MSCneurons.reclustNET,
                                    assay = 'SCT', # using SCT assay
                                    slot = 'data') # using normalized data

# select soft-power threshold
MSCneurons.reclustNET <- TestSoftPowers(MSCneurons.reclustNET,
                                        networkType = 'signed')

# plot the results:
plot_list.neuron <- PlotSoftPowers(MSCneurons.reclustNET)

# assemble with patchwork
png('./neurons/softPowerThresholds.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)

wrap_plots(plot_list.neuron, ncol = 2)
dev.off()

# use soft power threshold 7
# construct co-expression network:
MSCneurons.reclustNET <- ConstructNetwork(MSCneurons.reclustNET, 
                                          soft_power = 7,
                                          tom_outdir = paste0(root.dir, '/hdWGCNA/neurons/TOM/'),
                                          tom_name = 'neuron',
                                          overwrite_tom = TRUE)

# graph dendrogram
png('./neurons/HypoMap_dendrogram.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)

PlotDendrogram(MSCneurons.reclustNET,
               main = 'Neurons hdWGCNA Dendro')
dev.off()

# compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
MSCneurons.reclustNET <- ScaleData(MSCneurons.reclustNET,
                                   features = VariableFeatures(MSCneurons.reclustNET))

# compute all MEs in the full single-cell dataset
MSCneurons.reclustNET <- ModuleEigengenes(MSCneurons.reclustNET,
                                          exclude_grey = TRUE, 
                                          group.by.vars = 'indiv_genotype')

# module eigengenes:
MEs.neurons <- GetMEs(MSCneurons.reclustNET,
                      harmonized = FALSE)

# compute eigengene-based connectivity (kME):
MSCneurons.reclustNET <- ModuleConnectivity(MSCneurons.reclustNET)

#### WGCNA graphs ####
# plot genes ranked by kME for each module
png('./neurons/neuron_kMEs.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)

PlotKMEs(MSCneurons.reclustNET,
         ncol = 2)
dev.off()

# correlation between modules
# again, not very informative because of only small module number

png('./neurons/neuron_modulesCorrelogram.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)

ModuleCorrelogram(MSCneurons.reclustNET,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001),
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.neurons.me <- ModuleFeaturePlot(MSCneurons.reclustNET,
                                          features = 'MEs', # plot the MEs
                                          order = TRUE) # order so the points with highest MEs are on top

# stitch together with patchwork
png('./neurons/UMAP_neuronMEs.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)

wrap_plots(plot_list.neurons.me, ncol=2)
dev.off()

# get mods from object
mods.neurons <- colnames(MEs.neurons)
mods.neurons <- mods.neurons[mods.neurons != 'grey']

# add MEs to Seurat meta-data:
MSCneurons.reclustNET@meta.data <- cbind(MSCneurons.reclustNET@meta.data, MEs.neurons)

# neuropeptides
modules.neurons <- GetModules(MSCneurons.reclustNET)

#### Stop here ####
# line 2017 - 3130 in neurons_snseq_mouse_IMC.R removed, not included in manuscript
# doesn't run - also seems like cluster resolution = 0.8 was used for initial analysis
# but cluster res = 0.4 is final; can potentially revisit? 

# heatmap of kME and neuropeptides
neuropeptides.df <- data.frame(gene_name = neuropeptides.list %>% 
                                filter(!is.na(Gene.name.nile.tilapia)) %>% 
                                pull(Gene.name.nile.tilapia) %>% 
                                unique(), neuropeptide = TRUE)

