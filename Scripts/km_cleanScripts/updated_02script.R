# R-4.3.1, Seurat v.5.3.0
# HypoMap reference data from https://github.com/lsteuernagel/mapscvi

# net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq" 
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts"))
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/QC_filtering_fun.R")

# libraries (save dependencies and package versions)
load_packages(c("tidyverse", "Seurat", "ggalluvial", "mapscvi", "SeuratDisk"),
              out_prefix = "02")


#### load data ####
# load integrated Seurat object
# load('./data/mouse.snseq.combined.sct.RData') # IMC final copy
load('./data/integrated_seurat_withScType.rda')

# load HypoMap reference data
reference_hypoMap_downsample <- mapscvi::reference_hypoMap_downsample

# check intersection of data with reference
# 25982 overlap out of 26007 
# reference_hypoMap_downsample@assays[["RNA"]]@data@Dimnames[[1]] %>% 
#   int.ldfs@assays[["RNA"]]@data@Dimnames[[1]] %>% 
#   # intersect(dimnames(int.ldfs[["RNA"]])[[1]]) %>% # if v5
#   length()

## trying to add back in genes/cell names (works but not necessary?)
# int.ldfs[["RNA"]]@layers$counts@Dimnames = dimnames(int.ldfs[["RNA"]]) # if v5


# someone else had similar prob w/ RNA assay: 
# https://github.com/mojaveazure/seurat-disk/issues/27
# here was their solution
DefaultAssay(int.ldfs) <-  'RNA' # temporarily making 'RNA' active assay
int.ldfs <- FindVariableFeatures(int.ldfs)
DefaultAssay(int.ldfs) <-  'SCT' # returning 'SCT' as the default assay

save(int.ldfs, file = "./data/integrated_seurat_withScType.rda",
     compress = compression)

#### Map data to HypoMap ####

## Run bash script ./data/km_rscvi_docker/run_scvi_docker.sh
## Requires docker install! If running on lambcomp02 or other rental POD, comes pre-installed

## Do not try to run docker desktop on Windows if using VPN to access data files via ssh;
## can't bind mount to remote server without SFTP enabled by sys admin - not worth the headache, also:
## If running on local machine (not recommended w/o large RAM), check all filepaths & setwd() or clone this repo

system2("bash", "./data/km_rscvi_docker/run_scvi_docker.sh")

# This pulls the docker image from lsteuernagel/mapscvi and modifies it to pre-load all dependencies
# (see ./data/km_rscvi_docker/Dockerfile),
# then uses a modified map_new_seurat_hypoMap wrapper function 
# (./data/km_rscvi_docker/mapscvi_wrapper.R) to perform appropriate data prep and mapping,
# and finally saves the new Seurat object: ./data/integrated_seurat_withHypoMap.rda

#### Plot UMAP ####
load('./data/integrated_seurat_withHypoMap.rda') # is compressed as usual, may take a while

# reset working directory to save plots
setwd(paste0(root.dir, "/HypoMap"))

# graph umap
plot_query_labels(query_seura_object = int.ldfs,
                  reference_seurat = mapscvi::reference_hypoMap_downsample,
                  label_col = "C66_named",
                  overlay = FALSE,
                  labelonplot = FALSE)  

ggsave('projection.umap.predicted.png',
       height = 5,
       width = 8)

# projection onto total
plot_query_labels(query_seura_object = int.ldfs,
                  reference_seurat = reference_hypoMap_downsample,
                  label_col = "C66_named",
                  overlay = TRUE,
                  query_pt_size = 0.4,
                  labelonplot = TRUE,
                  label.size = 1.5)

ggsave('predicted.umap.projection.png',
       height = 5, width = 8)

# prediction probability
FeaturePlot(int.ldfs,
            features = "prediction_probability") + 
  NoAxes()

ggsave('predicted.probability.umap.png',
       height = 5, width = 8)

# prediction probability
DimPlot(int.ldfs,
        group.by = "predicted") + 
  NoAxes() 

ggsave('predicted.umap.png',
       height = 10, width = 20)  

#### Compare cell type results ####

## data from Steuernagel et al., 2022;
## download excel file: https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-022-00657-y/MediaObjects/42255_2022_657_MOESM3_ESM.xlsx
cell.type.hypomap.data = readxl::read_excel('./data/42255_2022_657_MOESM3_ESM.xlsx', sheet = 6)

## subset into the cell types used for hypomap
cell.type.hypomap.data %>% 
  dplyr::select(c(cluster_name,
                  parent_id)) %>% 
  droplevels() %>% 
  distinct() -> cell.type.hypomap.data

## pull out cell types of interest
# create list of cells
int.ldfs@meta.data %>% 
  pull(predicted) %>% 
  unique() %>% 
  c() -> mouse.snseq.cell.hypomap.list

# filter down to those cells
cell.type.hypomap.data %>% 
  filter(cluster_name %in% mouse.snseq.cell.hypomap.list) -> cell.type.hypomap.data

#### STOP and find where HypoMap_cluster_parent_ids.csv came from #### 
#load parent id information
cell.type.hypomap.data.parent.id = read.csv('./data/HypoMap_cluster_parent_ids.csv')

## add parent id
cell.type.hypomap.data %>% 
  left_join(cell.type.hypomap.data.parent.id,
            by = "parent_id") %>% 
  dplyr::select(-c(parent_id)) -> cell.type.hypomap.data

## create metadata of broad cell type
metadata.mouse.snseq = int.ldfs@meta.data

# add broad cell data
metadata.mouse.snseq %>% 
  left_join(cell.type.hypomap.data %>% 
              dplyr::rename(predicted = cluster_name),
            by = "predicted") -> metadata.mouse.snseq

## select relevant variables
metadata.mouse.snseq %>% 
  column_to_rownames('Cell_ID') %>% 
  dplyr::select(c(predicted,
                  parent_id.exp,
                  parent_id.broad)) -> metadata.mouse.snseq

## add metadata to seurat obj
int.ldfs <- AddMetaData(int.ldfs, 
                        metadata.mouse.snseq)

#### add unknown cells ####
#### project_query prob cell type score 
### specify labels and reduction  
cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C66_named
reference_reduction = "scvi"

### already have 'Batch_ID' and 'Cell_ID' in metadata 
# make sure to use 'Cell_ID' for cell identifiers instead of 'cell.id' -
# 'cell.id' is from pre-integration and may contain duplicates
query_seurat_object = int.ldfs

### project data into UMAP space from reference
## predict cell types

query_seurat_object = project_query(query_seurat_object = query_seurat_object,
                                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                    query_reduction = "scvi",
                                    label_vec = cluster_labels,
                                    result_type = "all") 

## add unknown 
### use prediction probability to assign cell type
predicted.prob.df = query_seurat_object@meta.data

# predicted tmp
tmp = project_query(query_seurat_object = query_seurat_object,
                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                    query_reduction = "scvi",
                    label_vec = cluster_labels)


tmp@meta.data %>% 
  mutate(predicted.prob.tmp = ifelse(prediction_probability >= 0.75,
                                     predicted,
                                     'Unknown')) -> tmp.df

## combine dataframes
predicted.prob.df %>% 
  full_join(tmp.df %>% select(c(Cell_ID, 
                                predicted.prob.tmp,
                                prediction_probability, 
                                predicted)), 
            join_by(Cell_ID)) %>% 
  column_to_rownames("Cell_ID") -> predicted.prob.df

# create cell type score
predicted.prob.df %>%
  mutate(neuron = rowSums(select(., contains(c('GABA', 'GLU'))))) %>%
  mutate(astrocytes = rowSums(select(., contains(c('Astrocytes', 'Tanycytes','Ependymal'))))) %>%
  mutate(oligodendrocytes = rowSums(select(., contains(c('OPC', 'Oligodendrocytes'))))) %>%
  mutate(immune = rowSums(select(., contains(c('Immune'))))) %>%
  mutate(vascular = rowSums(select(., contains(c('ParsTuber'))))) %>%
  mutate(parstuber = rowSums(select(., contains(c('Fibroblasts', 'Endothelial', 'Mural'))))) %>% 
# create column for comparison
  mutate(neuron.keep = ifelse(neuron >= 0.75, 1, 0)) %>% 
  mutate(astrocytes.keep = ifelse(astrocytes >= 0.75, 1, 0)) %>% 
  mutate(oligodendrocytes.keep = ifelse(oligodendrocytes >= 0.75, 1, 0)) %>% 
  mutate(immune.keep = ifelse(immune >= 0.75, 1, 0)) %>% 
  mutate(vascular.keep = ifelse(vascular >= 0.75, 1, 0)) %>% 
  mutate(parstuber.keep = ifelse(parstuber >= 0.75, 1, 0)) %>% 
# check counts
  mutate(unknown.count = rowSums(select(., contains(c('.keep'))))) %>% 
  mutate(unknown.count = as.numeric(unknown.count)) -> test.df

# plot probability distributions
test.df %>%
  select(c(predicted,
           prediction_probability)) %>% 
  droplevels() %>% 
  ggplot(aes(x = prediction_probability,
             fill = prediction_probability > 0.75)) +
  geom_vline(xintercept = 0.75) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(predicted~.) +
  theme_classic() +
  ggtitle('Hypomap predicted probability')

ggsave("cell.type.predicted.prob.png",
       width = 25, height = 20)

# just neurons
test.df %>% 
  select(c(neuron,
           predicted.prob.tmp,
           prediction_probability,
           predicted)) %>% 
  droplevels() %>% 
  ggplot(aes(x = neuron,
             color = neuron > 0.75,
             y = prediction_probability)) +
  geom_vline(xintercept = 0.75) +
  geom_hline(yintercept = 0.75) +
  geom_point() +
  facet_wrap(predicted~.) +
  theme_classic() +
  ylim(0,1)+
  ggtitle('Hypomap neuron score')

ggsave("neuron.score.predicted.prob.png",
       width = 25, height = 20)

# unknown
test.df %>% 
  select(c(neuron,
           predicted.prob.tmp,
           prediction_probability,
           predicted)) %>% 
  filter(prediction_probability < 0.75) %>% 
  droplevels() %>% 
  ggplot(aes(x = neuron,
             fill = neuron > 0.75)) +
  geom_vline(xintercept = 0.75) +
  geom_hline(yintercept = 0.75) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(predicted~.) +
  theme_classic()+
  ggtitle('Hypomap Unknown Neuron score')

ggsave("neuron.score.predicted.prob.histo.png",
       width = 25, height = 20)

### use prediction probability to assign cell type
# celltype call
test.df %>% 
  mutate(predicted.prob = ifelse(unknown.count > 0,
                                 predicted,
                                 'Unknown'),
         parent_id.exp.prob = ifelse(unknown.count > 0,
                                     parent_id.exp,
                                     'Unknown'),
         parent_id.broad.prob = ifelse(unknown.count > 0,
                                       parent_id.broad,
                                       'Unknown')) %>% 
  select(predicted.prob, parent_id.exp.prob, parent_id.broad.prob) -> test.df.meta

## add metadata to object
# check
  identical(rownames(test.df.meta), colnames(int.ldfs))
# [1] TRUE
int.ldfs = AddMetaData(int.ldfs, 
                       test.df.meta)

  save(int.ldfs, file = "./data/integrated_seurat_withHypoMap_predictions.rda",
       compress = compression)

#### graph comparison ####
### graph parent umaps
## Prediction broad
DimPlot(int.ldfs,
        group.by = "parent_id.broad.prob",
        cols = c('#1b9e77',
                 '#d95f02',
                 '#7570b3',
                 '#e7298a',
                 '#66a61e',
                 '#e6ab02',
                 '#a6761d',
                 'grey')) +
  NoAxes() 

ggsave('predicted.broad.UMAP.png',
       height = 10, width = 12)  


## predicted exp
DimPlot(int.ldfs, 
        group.by = "parent_id.exp.prob") +
  NoAxes() 

ggsave('predicted.expandedUMAP.png',
       height = 10, width = 15)  

### graph cell counts across samples
## broad
int.ldfs@meta.data %>% 
  select(c(orig.ident,
           parent_id.broad.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = parent_id.broad.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('predicted.cell.counts.broad.png',
       height = 5, width = 8)  

## expanded
int.ldfs@meta.data %>% 
  select(c(orig.ident,
           parent_id.exp.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = parent_id.exp.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('predicted.cell.counts.expanded.png',
       height = 5, width = 8)  

#reduce 
int.ldfs@meta.data %>% 
  select(c(orig.ident,
           parent_id.exp.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  ggplot(aes(x = reorder(parent_id.exp.prob,
                         -Freq),
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('predicted.cell.counts.expanded.reduced.png',
       height = 5,
       width = 8)  

## predicted
int.ldfs@meta.data %>% 
  select(c(orig.ident,
           predicted.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = predicted.prob,
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('predicted.cell.counts.predicted.png',
       height = 5, width = 12)

#reduce 
int.ldfs@meta.data %>% 
  select(c(orig.ident,
           predicted.prob)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  ggplot(aes(x = reorder(predicted.prob, -Freq),
             y = Freq,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('predicted.cell.counts.predicted.reduced.png',
       height = 5, width = 8)

#reduce 
#neurons
#genotype
##get list of cell types with at least 50 in one sample
int.ldfs@meta.data %>% 
  filter(parent_id.broad.prob %in% 
           c("C7-2: GABA", "C7-1: GLU")) %>% 
  select(c(orig.ident,
           predicted.prob,
           indiv_genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  pull(predicted.prob) %>% 
  unique() -> neuron.cell.types.predicted.50

#graph
int.ldfs@meta.data %>% 
  filter(parent_id.broad.prob %in%
           c("C7-2: GABA", "C7-1: GLU")) %>% 
  select(c(orig.ident,
           predicted.prob,
           indiv_genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(predicted.prob %in% 
           neuron.cell.types.predicted.50) %>%
  ggplot(aes(x = reorder(predicted.prob, -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_hline(yintercept = 50) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width = 0.75),
             aes(shape = indiv_genotype, group = orig.ident), 
             size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('neurons.predicted.cell.counts.predicted.reduced.png',
       height = 5, width = 8)

#scale
#reduce 
#neurons
#genotype
##get list of cell types with at least 50 in one sample
int.ldfs@meta.data %>% 
  filter(parent_id.broad.prob %in% 
           c("C7-2: GABA", "C7-1: GLU")) %>% 
  select(c(orig.ident,
           predicted.prob,
           indiv_genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq >= 50) %>% 
  pull(predicted.prob) %>% 
  unique() -> neuron.cell.types.predicted.50

#dataframe of scale
int.ldfs@meta.data %>% 
  filter(parent_id.broad.prob %in% 
           c("C7-2: GABA", "C7-1: GLU")) %>%
  select(c(orig.ident,
           indiv_genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  rename(total.neuron.count = Freq) -> neuron.genotype.scale

#dataframe of cluster counts
int.ldfs@meta.data %>% 
  filter(parent_id.broad.prob %in% 
           c("C7-2: GABA", "C7-1: GLU")) %>%
  select(c(orig.ident,
           predicted.prob,
           indiv_genotype)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(predicted.prob %in% neuron.cell.types.predicted.50) %>% 
  as.data.frame() -> neuron.cluster.genotype.count

#combine counts and scale
neuron.cluster.genotype.count %>% 
  full_join(neuron.genotype.scale) %>% 
  mutate(Percent = 100*Freq/
           total.neuron.count) -> neuron.cluster.genotype.count

#graph
neuron.cluster.genotype.count %>% 
  ggplot(aes(x = reorder(predicted.prob, -Percent),
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position=position_dodge(width = 0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident), size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1))

ggsave('neurons.predicted.cell.counts.predicted.reduced.scale.png',
       height = 5, width = 8)

###  sctype expanded vs hypomap broad
int.ldfs@meta.data %>%
  select(c(sctype.integrated,
           parent_id.broad.prob)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis2 = reorder(sctype.integrated, -Freq.scale),
             axis3 = reorder(parent_id.broad.prob, -Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = parent_id.broad.prob)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.integrated",
                              "parent_id.broad.prob"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()

ggsave('alluvial.sctype.integratedVparent_id.broad.png', 
       height = 10, width = 15)

# projection with broad cell types
plot_query_labels(query_seura_object = int.ldfs,
                  reference_seurat = reference_hypoMap_downsample,
                  label_col_query = "parent_id.broad",
                  label_col = "C7_named",
                  overlay = TRUE,
                  query_pt_size = 0.4,
                  labelonplot = TRUE,
                  label.size = 5) 

ggsave('predicted.UMAPprojection.broad.label.png',
       height = 10, width = 12)

#without label
plot_query_labels_test(query_seura_object = int.ldfs,
                       reference_seurat = reference_hypoMap_downsample,
                       label_col_query = "parent_id.broad.prob",
                       label_col = "C7_named",
                       overlay = TRUE,
                       query_pt_size = 0.4,
                       labelonplot = FALSE,
                       label.size = 5) + 
  theme(legend.position = "none")

ggsave('predicted.UMAPprojection.broad.unlabeled.png',
       height = 10, width = 10)

#without label
# add unknown 
plot_query_labels_test(query_seura_object = int.ldfs,
                       reference_seurat = reference_hypoMap_downsample,
                       label_col_query = "parent_id.broad.prob",
                       label_col = "C7_named",
                       overlay = TRUE,
                       query_pt_size = 0.4,
                       labelonplot = FALSE,
                       label.size = 5,
                       cols = c('#1b9e77',
                                '#d95f02',
                                '#7570b3',
                                '#e7298a',
                                '#66a61e',
                                '#e6ab02',
                                '#a6761d',
                                'grey'))

ggsave('predictedUMAPprojection.broad.unknown.png',
       height = 10, width = 12)

#### check unknown cell types ####
### compare known to unknown cell counts
## broad cell types
int.ldfs@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  select(c(parent_id.broad, orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(int.ldfs@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              select(c(parent_id.broad, orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(parent_id.broad, -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  xlab('')

ggsave('knownVunknown.cell.counts.broad.png',
       height = 5, width = 8)  

## expanded cell types
int.ldfs@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  select(c(parent_id.exp, orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(int.ldfs@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              select(c(parent_id.exp, orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(parent_id.exp, -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width = 0.75)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  xlab('')

ggsave('knownVunkown.cell.counts.exp.png',
       height = 5, width = 10)  

## expanded cell types reduce
int.ldfs@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  select(c(parent_id.exp, orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(int.ldfs@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              select(c(parent_id.exp, orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  filter(Freq > 50) %>% 
  ggplot(aes(x = reorder(parent_id.exp, -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.5,
                                   hjust = 0.5)) +
  xlab('')

ggsave('knownVunkown.cell.counts.exp.reduced.png',
       height = 5, width = 8)  

## predicted cell types
int.ldfs@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  select(c(predicted, orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(int.ldfs@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              select(c(predicted, orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  ggplot(aes(x = reorder(predicted, -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 1,
                                   hjust = 0)) +
  xlab('')

ggsave('knownVunknown.cell.counts.png',
       height = 5, width = 10)  

## predicted cell types reduce
int.ldfs@meta.data %>% 
  filter(predicted.prob != 'Unknown') %>% 
  select(c(predicted, orig.ident)) %>%
  table() %>% 
  data.frame() %>% 
  mutate(Known = 'Known') %>% 
  full_join(int.ldfs@meta.data %>% 
              filter(predicted.prob == 'Unknown') %>% 
              select(c(predicted, orig.ident)) %>%
              table() %>% 
              data.frame() %>% 
              mutate(Known = 'Unknown')) %>% 
  filter(Freq > 50) %>% 
  ggplot(aes(x = reorder(predicted, -Freq),
             y = Freq,
             color = Known)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 1,
                                   hjust = 0)) +
  xlab('')

ggsave('knownVunknown.cell.counts.reduced.png',
       height = 5, width = 8)  


### check distribution of probability for cell types
int.ldfs@meta.data %>% 
  select(c(predicted,
           orig.ident,
           prediction_probability,
           parent_id.broad)) %>% 
  droplevels() %>% 
  ggplot(aes(x = prediction_probability,
             fill = prediction_probability > 0.25)) +
  geom_vline(xintercept = 0.25) +
  geom_histogram(binwidth = 0.02) +
  facet_wrap(parent_id.broad~., nrow = 2) +
  theme_classic()

ggsave('predicted.prob.histo.png',
       height = 8, width = 10)  


#### project_query spatial ####
### specify labels and reduction  
cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C286_named
reference_reduction = "scvi"
query_seurat_object = int.ldfs

### project data into UMAP space from reference
## predict cell types
query_seurat_object = project_query(query_seurat_object = query_seurat_object,
                                    reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                    reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                    query_reduction = "scvi",
                                    label_vec = cluster_labels)  


## graph
plot_query_labels(query_seura_object = query_seurat_object,
                  reference_seurat = mapscvi::reference_hypoMap_downsample,
                  label_col = "C286_named",
                  overlay = FALSE,
                  labelonplot = FALSE)

Seurat::DimPlot(query_seurat_object,
                group.by = "predicted") +
  Seurat::NoLegend()

ggsave('predicted.umap.unlabeled.png',
       width = 9, height = 7)

## add unknown 
### use prediction probability to assign cell type
predicted.prob.df = query_seurat_object@meta.data

# celltype call
predicted.prob.df <- predicted.prob.df %>% 
  mutate(predicted.prob = ifelse(prediction_probability >= 0.25, 
                                 predicted, 'Unknown'))

# combine with spatial data
data.spatial = read.csv('./data/clusterC286_brain_annotation_HypoMap.csv')

predicted.prob.df.spatial = predicted.prob.df %>% 
  select(c(predicted)) %>% 
  left_join(data.spatial %>% 
              select(c(cluster_name, Region_summarized)) %>% 
              rename(predicted = cluster_name))

## add metadata to object
query_seurat_object <- AddMetaData(query_seurat_object, 
                                   predicted.prob.df.spatial$Region_summarized,
                                   col.name = 'Region_summarized')

# save metadata
meta <- as.data.frame(query_seurat_object@meta.data)
write.csv(meta, "./data/metadata_with_predictedRegions.csv")

# Seurat::DimPlot(query_seurat_object,
#                 group.by = "predicted.prob")

Seurat::DimPlot(query_seurat_object,
                group.by = "Region_summarized")

ggsave('predicted.brain.region.umap.png',
       width = 12, height = 8)

## graph counts
query_seurat_object@meta.data %>% 
  select(c(orig.ident, Region_summarized)) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(y = reorder(Region_summarized, -Freq),
             x = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1)) +
  ylab('')

ggsave('brain.region.count.all.png',
       height = 8, width = 6.5)  

## graph counts
# known
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  select(c(orig.ident, Region_summarized)) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(y = reorder(Region_summarized, -Freq),
             x = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1)) +
  ylab('')

ggsave('brain.region.count.known.png',
       height = 8, width = 6.5)  

# known
# reduce
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  select(c(orig.ident,
           Region_summarized,
           indiv_genotype)) %>% 
  table() %>% 
  data.frame() %>% 
  filter(Region_summarized %in% c('Medial preoptic area',
                                  'Lateral preoptic area',
                                  'Lateral hypotalamic area',
                                  'Zona incerta',
                                  'Paraventricular hypothalamic nucleus')) %>% 
  ggplot(aes(x = reorder(Region_summarized, -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position = position_dodge(width = 0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5, 
                                   hjust = 0.5)) +
  xlab('') 

ggsave('brain.region.count.known.reducedByGenotype.png',
       height = 7, width = 7) 

# known
# reduce
# scale
query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  select(c(orig.ident,
           Region_summarized,
           indiv_genotype)) %>% 
  table() %>% 
  data.frame() %>% 
  filter(Region_summarized %in% 
           c('Medial preoptic area',
             'Lateral preoptic area',
             'Lateral hypotalamic area',
             'Zona incerta',
             'Paraventricular hypothalamic nucleus')) -> nuclei.genotype.count

#dataframe of scale
nuclei.genotype.scale = query_seurat_object@meta.data %>% 
  filter(parent_id.broad.prob != 'Unknown') %>% 
  select(c(orig.ident, indiv_genotype)) %>% 
  table() %>%
  as.data.frame() %>% 
  rename(total.count = Freq)

#combine counts and scale
nuclei.genotype.count = nuclei.genotype.count %>% 
  full_join(nuclei.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.count)

#graph
nuclei.genotype.count %>% 
  ggplot(aes(x = reorder(Region_summarized, -Percent),
             y = Percent,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75)) +
  geom_point(position = position_dodge(width = 0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5, 
                                   hjust = 0.5)) +
  xlab('') 

ggsave('brain.region.count.known.reducedByGenotype.scale.png',
       height = 7, width = 7) 

# ~fin~