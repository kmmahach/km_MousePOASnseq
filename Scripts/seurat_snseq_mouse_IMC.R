#### Mouse snseq seurat analysis
### Integrating dom and sub, male and female snseq data
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

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
#load libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(HGNChelper)
library(clustree)
library(ggalluvial)
library(ggrepel)
library(scCustomize)

# load libraries for bubble plot
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

#### load data ####
### Load data
## male
# male dom data
male.dom.data <- Read10X(data.dir = "../count/MaleDom/outs/filtered_feature_bc_matrix/")

# male sub data
male.sub.data <- Read10X(data.dir = "../count/MaleSub/outs/filtered_feature_bc_matrix/")

## female
# female dom data
female.dom.data <- Read10X(data.dir = "../count/FemaleDom/outs/filtered_feature_bc_matrix/")

# male sub data
female.sub.data <- Read10X(data.dir = "../count/FemaleSub/outs/filtered_feature_bc_matrix/")

### load souporcell data
##dom
# male
dom.male.clusters = read_tsv('../souporcell/dommale_soupercell/clusters.tsv')
#rename barcode to Cell.id
dom.male.clusters = dom.male.clusters %>%
  dplyr::rename(Cell.id = barcode) %>% 
  dplyr::rename(Doublet = status) %>% 
  dplyr::rename(Genotype = assignment) %>% 
  dplyr::select(Cell.id,
                Doublet,
                Genotype) 
  
# female
dom.female.clusters = read_tsv('../souporcell/domfemale_soupercell/clusters.tsv')
#rename barcode to Cell.id
dom.female.clusters = dom.female.clusters %>%
  dplyr::rename(Cell.id = barcode) %>% 
  dplyr::rename(Doublet = status) %>% 
  dplyr::rename(Genotype = assignment) %>% 
  dplyr::select(Cell.id,
                Doublet,
                Genotype)

##sub
# male
sub.male.clusters = read_tsv('../souporcell/submale_soupercell/clusters.tsv')
#rename barcode to Cell.id
sub.male.clusters = sub.male.clusters %>%
  dplyr::rename(Cell.id = barcode) %>% 
  dplyr::rename(Doublet = status) %>% 
  dplyr::rename(Genotype = assignment) %>% 
  dplyr::select(Cell.id,
                Doublet,
                Genotype)

# female
sub.female.clusters = read_tsv('../souporcell/subfemale_soupercell/clusters.tsv')
#rename barcode to Cell.id
sub.female.clusters = sub.female.clusters %>%
  dplyr::rename(Cell.id = barcode) %>% 
  dplyr::rename(Doublet = status) %>% 
  dplyr::rename(Genotype = assignment) %>% 
  dplyr::select(Cell.id,
                Doublet,
                Genotype)


# load combined data
# load('mouse.snseq.combined.sct.RData')

#### functions ####
### modify gene_sets_prepare function to integrate with user added data
## do not require xlsx for gene_sets_prepare
## remove checkGeneSymbols function 
# remove forcing genes to uppercase
gene_sets_prepare.df = function(db_file, cell_type){
  
  cell_markers = db_file #change to dataframe
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}
#### prepare sctype ####
#### cell type assignment
### load functions
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

### prepare gene set list
## download gene lists
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# set tissue type
tissue = "Brain" 
# prepare gene sets
gs_list = gene_sets_prepare(db_,
                            tissue)

#### finding min.cells threshold ####
### compare the number of genes for different threshold of min.cells
#create range of min.cells to check
gene.min.cells.range = c(3,
                         seq(from = 10, 
                             to = 100, 
                             by = 10),
                         seq(from = 200, 
                             to = 1000, 
                             by = 100))

##DomMale
# crate empty dataframe
gene.min.cells.dom.male = data.frame()

## count genes for range of min.cells
for (i in gene.min.cells.range) {
  temp <- CreateSeuratObject(counts = male.dom.data,
                             project = "Male.Dom", 
                             min.cells = i,
                             min.features = 700)
  
  temp.data.frame = data.frame('min.cells' = i,
                               'gene.count' = temp@assays[["RNA"]]@data %>% 
                                 nrow())
  
  gene.min.cells.dom.male = rbind(gene.min.cells.dom.male,
                                  temp.data.frame)
}

# delete temporary dataframes
rm(temp)
rm(temp.data.frame)

# graph 
gene.min.cells.dom.male %>% 
  ggplot(aes(x = min.cells,
             y = gene.count,
             label = paste(min.cells,
                           gene.count,
                           sep = ','))) +
  geom_line()+
  geom_label()+
  ylim(0,
       max(gene.min.cells.dom.male$gene.count)) +
  theme_classic()+
  ggtitle('DomMale') 
ggsave('./QC_filtering/MaleDom/min.cells.male.dom.png',
       height = 10,
       width = 10)

##DomFemale
# crate empty dataframe
gene.min.cells.dom.female = data.frame()

## count genes for range of min.cells
for (i in gene.min.cells.range) {
  temp <- CreateSeuratObject(counts = female.dom.data,
                             project = "Male.Dom", 
                             min.cells = i,
                             min.features = 700)
  
  temp.data.frame = data.frame('min.cells' = i,
                               'gene.count' = temp@assays[["RNA"]]@data %>% 
                                 nrow())
  
  gene.min.cells.dom.female = rbind(gene.min.cells.dom.female,
                                  temp.data.frame)
}

# delete temporary dataframes
rm(temp)
rm(temp.data.frame)

# graph 
gene.min.cells.dom.female %>% 
  ggplot(aes(x = min.cells,
             y = gene.count,
             label = paste(min.cells,
                           gene.count,
                           sep = ','))) +
  geom_line()+
  geom_label()+
  ylim(0,
       max(gene.min.cells.dom.female$gene.count)) +
  theme_classic()+
  ggtitle('DomFemale') 
ggsave('./QC_filtering/FemaleDom/min.cells.female.dom.png',
       height = 10,
       width = 10)

##SubMale
# crate empty dataframe
gene.min.cells.sub.male = data.frame()

## count genes for range of min.cells
for (i in gene.min.cells.range) {
  temp <- CreateSeuratObject(counts = male.sub.data,
                             project = "Male.Sub", 
                             min.cells = i,
                             min.features = 700)
  
  temp.data.frame = data.frame('min.cells' = i,
                               'gene.count' = temp@assays[["RNA"]]@data %>% 
                                 nrow())
  
  gene.min.cells.sub.male = rbind(gene.min.cells.sub.male,
                                  temp.data.frame)
}

# delete temporary dataframes
rm(temp)
rm(temp.data.frame)

# graph 
gene.min.cells.sub.male %>% 
  ggplot(aes(x = min.cells,
             y = gene.count,
             label = paste(min.cells,
                           gene.count,
                           sep = ','))) +
  geom_line()+
  geom_label()+
  ylim(0,
       max(gene.min.cells.sub.male$gene.count)) +
  theme_classic()+
  ggtitle('SubMale') 
ggsave('./QC_filtering/MaleSub/min.cells.male.sub.png',
       height = 10,
       width = 10)

##SubFemale
# crate empty dataframe
gene.min.cells.sub.female = data.frame()

## count genes for range of min.cells
for (i in gene.min.cells.range) {
  temp <- CreateSeuratObject(counts = female.sub.data,
                             project = "Male.Sub", 
                             min.cells = i,
                             min.features = 700)
  
  temp.data.frame = data.frame('min.cells' = i,
                               'gene.count' = temp@assays[["RNA"]]@data %>% 
                                 nrow())
  
  gene.min.cells.sub.female = rbind(gene.min.cells.sub.female,
                                  temp.data.frame)
}

# delete temporary dataframes
rm(temp)
rm(temp.data.frame)

# graph 
gene.min.cells.sub.female %>% 
  ggplot(aes(x = min.cells,
             y = gene.count,
             label = paste(min.cells,
                           gene.count,
                           sep = ','))) +
  geom_line()+
  geom_label()+
  ylim(0,
       max(gene.min.cells.sub.female$gene.count)) +
  theme_classic()+
  ggtitle('SubFemale') 
ggsave('./QC_filtering/FemaleSub/min.cells.female.sub.png',
       height = 10,
       width = 10)

#### MaleDom clustering ####
###filter cells first
##create object
male.dom.object <- CreateSeuratObject(counts = male.dom.data,
                                      project = "Male.Dom", 
                                      min.cells = 3,
                                      min.features = 200)

## add souporcell data
male.dom.object = AddMetaData(object = male.dom.object,
                              metadata = dom.male.clusters %>% 
                                column_to_rownames('Cell.id')) 

## label mitochondria genes
male.dom.object <- PercentageFeatureSet(male.dom.object, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")

## QC graph
VlnPlot(male.dom.object,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleDom/QCmetrics.male.dom.png',
       height = 10,
       width = 10)

#doublet
VlnPlot(subset(male.dom.object, 
               subset = Doublet == 'doublet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleDom/QCmetrics.male.dom.doublet.png',
       height = 10,
       width = 10)

#singlet
VlnPlot(subset(male.dom.object, 
               subset = Doublet == 'singlet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleDom/QCmetrics.male.dom.singlet.png',
       height = 10,
       width = 10)


#unassigned
VlnPlot(subset(male.dom.object, 
               subset = Doublet == 'unassigned'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleDom/QCmetrics.male.dom.unassigned.png',
       height = 10,
       width = 10)

### graph
##features and counts
FeatureScatter(male.dom.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/MaleDom/featurescatter.male.dom.png',
       height = 10,
       width = 10)

#filter
FeatureScatter(subset(male.dom.object,
                      subset = nFeature_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 700)+
  geom_hline(yintercept = 3000) +
  annotate("text", x=0, y=700, label="700") +
  annotate("text", x=0, y=3000, label="3000")
ggsave('./QC_filtering/MaleDom/featurescatter.male.dom.limits.png',
       height = 10,
       width = 10)

##histogram
male.dom.object@meta.data %>%
  subset(nFeature_RNA < 6000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 700 & nFeature_RNA < 3000,
                       'Kept',
                       'Remove')) %>%
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('DomMale UMI counts')  +
  scale_fill_manual(values = c('red',
                               'black'))
ggsave('./QC_filtering/MaleDom/UMI.counts.filtered.male.dom.png',
       height = 10,
       width = 10)

## check filters
# want ~10,000 cells per sample
#6559
male.dom.object@meta.data %>%
  subset(nFeature_RNA <= 3000) %>%
  subset(nFeature_RNA >= 700) %>%
  subset(percent.mt < 5) %>% 
  subset(Doublet == 'singlet') %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##create new object with filters
male.dom.object = subset(male.dom.object, 
       subset = Doublet == 'singlet' & nFeature_RNA <= 3000 & nFeature_RNA >= 700 & percent.mt < 5)

#check filters worked
# should have ~6559 cells
male.dom.object@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##features and counts
FeatureScatter(male.dom.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/MaleDom/featurescatter.male.dom.filter.png',
       height = 10,
       width = 10)

##histogram
male.dom.object@meta.data %>%
  ggplot(aes(nFeature_RNA))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('DomMale UMI counts')  
ggsave('./QC_filtering/MaleDom/UMI.counts.filtered.male.dom.filter.png',
       height = 10,
       width = 10)

## run SCTransform, PCA, UMAP, and cluster 
#use resolution from clustree below: 1
male.dom.object = male.dom.object %>% 
  SCTransform() %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

## check PC's in elbow plot
ElbowPlot(male.dom.object)
ggsave('./QC_filtering/MaleDom/Elbowplot.dom.male.png',
       height = 10,
       width = 10)

##graph umap
DimPlot(male.dom.object,
        reduction = "umap",
        label = TRUE,
        repel = TRUE)+
  ggtitle('DomMale')
ggsave('./QC_filtering/MaleDom/Clusters.dimplot.male.dom.png',
       height = 10,
       width = 10)

##graph UMI counts on umap
FeaturePlot(male.dom.object,
        reduction = "umap",
        features = "nFeature_RNA",
        label = TRUE,
        repel = TRUE)+
  ggtitle('DomMale') +
  labs(color = "nFeature_RNA")
ggsave('./QC_filtering/MaleDom/nFeature_RNA.dimplot.male.dom.png',
       height = 10,
       width = 10)

## variable features
VariableFeaturePlot_scCustom(male.dom.object)
ggsave('./QC_filtering/MaleDom/Variable.features.male.dom.png',
       height = 10,
       width = 10)

### Clustree
### Find clusters using a range of resolutions
# cluster across range
male.dom.object.clustree <- Seurat::FindClusters(object = male.dom.object,
                                                               resolution = seq(from = 0, to = 2, by = 0.2))
#check data
# head(male.dom.object.clustree[[]])
#clustree
clustree(male.dom.object.clustree,
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions')+
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/MaleDom/clustree.dom.male.png',
       height = 10,
       width = 10)
#reduced range
male.dom.object.clustree.reduce <- Seurat::FindClusters(object = male.dom.object,
                                                               resolution = seq(from = 0, to = 1.4, by = 0.2))
#clustree
clustree(male.dom.object.clustree.reduce,
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions reduced') +
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/MaleDom/clustree.reduced.dom.male.png',
       height = 10,
       width = 10)

#remove clustree
rm(male.dom.object.clustree)
rm(male.dom.object.clustree.reduce)

#### FemaleDom clustering ####
###filter cells first
##create object
female.dom.object <- CreateSeuratObject(counts = female.dom.data,
                                      project = "Female.Dom", 
                                      min.cells = 3,
                                      min.features = 200)

## add souporcell data
female.dom.object = AddMetaData(object = female.dom.object,
                              metadata = dom.female.clusters %>% 
                                column_to_rownames('Cell.id')) 

## label mitochondria genes
female.dom.object <- PercentageFeatureSet(female.dom.object, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")

## QC graph
VlnPlot(female.dom.object,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleDom/QCmetrics.female.dom.png',
       height = 10,
       width = 10)

#doublet
VlnPlot(subset(female.dom.object, 
               subset = Doublet == 'doublet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleDom/QCmetrics.female.dom.doublet.png',
       height = 10,
       width = 10)

#singlet
VlnPlot(subset(female.dom.object, 
               subset = Doublet == 'singlet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleDom/QCmetrics.female.dom.singlet.png',
       height = 10,
       width = 10)


#unassigned
VlnPlot(subset(female.dom.object, 
               subset = Doublet == 'unassigned'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleDom/QCmetrics.female.dom.unassigned.png',
       height = 10,
       width = 10)

### graph
##features and counts
FeatureScatter(female.dom.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/FemaleDom/featurescatter.female.dom.png',
       height = 10,
       width = 10)

#filter
FeatureScatter(subset(female.dom.object,
                      subset = nFeature_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 700)+
  geom_hline(yintercept = 1400) +
  annotate("text", x=0, y=700, label="700") +
  annotate("text", x=0, y=1400, label="1400")
ggsave('./QC_filtering/FemaleDom/featurescatter.female.dom.limits.png',
       height = 10,
       width = 10)

##histogram
female.dom.object@meta.data %>%
  subset(nFeature_RNA < 6000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 700 & nFeature_RNA < 1400,
                       'Kept',
                       'Remove')) %>%
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('DomFemale UMI counts')  +
  scale_fill_manual(values = c('red',
                               'black'))
ggsave('./QC_filtering/FemaleDom/UMI.counts.filtered.female.dom.png',
       height = 10,
       width = 10)

## check filters
# want ~10,000 cells per sample
#6869
female.dom.object@meta.data %>%
  subset(nFeature_RNA <= 1400) %>%
  subset(nFeature_RNA >= 700) %>%
  subset(percent.mt < 5) %>% 
  subset(Doublet == 'singlet') %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##create new object with filters
female.dom.object = subset(female.dom.object, 
                         subset = Doublet == 'singlet' & nFeature_RNA <= 1400 & nFeature_RNA >= 700 & percent.mt < 5)

#check filters worked
# should have ~6869 cells
female.dom.object@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##features and counts
FeatureScatter(female.dom.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/FemaleDom/featurescatter.female.dom.filter.png',
       height = 10,
       width = 10)

##histogram
female.dom.object@meta.data %>%
  ggplot(aes(nFeature_RNA))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('DomFemale UMI counts')  
ggsave('./QC_filtering/FemaleDom/UMI.counts.filtered.female.dom.filter.png',
       height = 10,
       width = 10)

## run SCTransform, PCA, UMAP, and cluster 
#use resolution from clustree below: 1
female.dom.object = female.dom.object %>% 
  SCTransform() %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

## check PC's in elbow plot
ElbowPlot(female.dom.object)
ggsave('./QC_filtering/FemaleDom/Elbowplot.dom.female.png',
       height = 10,
       width = 10)

##graph umap
DimPlot(female.dom.object,
        reduction = "umap",
        label = TRUE,
        repel = TRUE)+
  ggtitle('DomFemale')
ggsave('./QC_filtering/FemaleDom/Clusters.dimplot.female.dom.png',
       height = 10,
       width = 10)

##graph UMI counts on umap
FeaturePlot(female.dom.object,
            reduction = "umap",
            features = "nFeature_RNA",
            label = TRUE,
            repel = TRUE)+
  ggtitle('DomFemale') +
  labs(color = "nFeature_RNA")
ggsave('./QC_filtering/FemaleDom/nFeature_RNA.dimplot.female.dom.png',
       height = 10,
       width = 10)

## variable features
VariableFeaturePlot_scCustom(female.dom.object)
ggsave('./QC_filtering/FemaleDom/Variable.features.female.dom.png',
       height = 10,
       width = 10)

### Clustree
### Find clusters using a range of resolutions
# cluster across range
female.dom.object.clustree <- Seurat::FindClusters(object = female.dom.object,
                                                 resolution = seq(from = 0, to = 2, by = 0.2))
#check data
# head(female.dom.object.clustree[[]])
#clustree
clustree(female.dom.object.clustree,
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions')+
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/FemaleDom/clustree.dom.female.png',
       height = 10,
       width = 10)
#reduced range
female.dom.object.clustree.reduce <- Seurat::FindClusters(object = female.dom.object,
                                                        resolution = seq(from = 0, to = 1.4, by = 0.2))
#clustree
clustree(female.dom.object.clustree.reduce,
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions reduced') +
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/FemaleDom/clustree.reduced.dom.female.png',
       height = 10,
       width = 10)

#remove clustree
rm(female.dom.object.clustree)
rm(female.dom.object.clustree.reduce)

#### MaleSub clustering ####
###filter cells first
##create object
male.sub.object <- CreateSeuratObject(counts = male.sub.data,
                                      project = "Male.Sub", 
                                      min.cells = 3,
                                      min.features = 200)

## add souporcell data
male.sub.object = AddMetaData(object = male.sub.object,
                              metadata = sub.male.clusters %>% 
                                column_to_rownames('Cell.id')) 

## label mitochondria genes
male.sub.object <- PercentageFeatureSet(male.sub.object, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")

## QC graph
VlnPlot(male.sub.object,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleSub/QCmetrics.male.sub.png',
       height = 10,
       width = 10)

#doublet
VlnPlot(subset(male.sub.object, 
               subset = Doublet == 'doublet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleSub/QCmetrics.male.sub.doublet.png',
       height = 10,
       width = 10)

#singlet
VlnPlot(subset(male.sub.object, 
               subset = Doublet == 'singlet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleSub/QCmetrics.male.sub.singlet.png',
       height = 10,
       width = 10)


#unassigned
VlnPlot(subset(male.sub.object, 
               subset = Doublet == 'unassigned'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/MaleSub/QCmetrics.male.sub.unassigned.png',
       height = 10,
       width = 10)

### graph
##features and counts
FeatureScatter(male.sub.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/MaleSub/featurescatter.male.sub.png',
       height = 10,
       width = 10)

#filter
FeatureScatter(subset(male.sub.object,
                      subset = nFeature_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 700)+
  geom_hline(yintercept = 2500) +
  annotate("text", x=0, y=700, label="700") +
  annotate("text", x=0, y=2500, label="2500")
ggsave('./QC_filtering/MaleSub/featurescatter.male.sub.limits.png',
       height = 10,
       width = 10)

##histogram
male.sub.object@meta.data %>%
  subset(nFeature_RNA < 6000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 700 & nFeature_RNA < 2500,
                       'Kept',
                       'Remove')) %>%
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('SubMale UMI counts')  +
  scale_fill_manual(values = c('red',
                               'black'))
ggsave('./QC_filtering/MaleSub/UMI.counts.filtered.male.sub.png',
       height = 10,
       width = 10)

## check filters
# want ~10,000 cells per sample
#6543
male.sub.object@meta.data %>%
  subset(nFeature_RNA <= 2500) %>%
  subset(nFeature_RNA >= 700) %>%
  subset(percent.mt < 5) %>% 
  subset(Doublet == 'singlet') %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##create new object with filters
male.sub.object = subset(male.sub.object, 
                         subset = Doublet == 'singlet' & nFeature_RNA <= 2500 & nFeature_RNA >= 700 & percent.mt < 5)

#check filters worked
# should have ~6543 cells
male.sub.object@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##features and counts
FeatureScatter(male.sub.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/MaleSub/featurescatter.male.sub.filter.png',
       height = 10,
       width = 10)

##histogram
male.sub.object@meta.data %>%
  ggplot(aes(nFeature_RNA))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('SubMale UMI counts')  
ggsave('./QC_filtering/MaleSub/UMI.counts.filtered.male.sub.filter.png',
       height = 10,
       width = 10)

## run SCTransform, PCA, UMAP, and cluster 
#use resolution from clustree below: 1
male.sub.object = male.sub.object %>% 
  SCTransform() %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

## check PC's in elbow plot
ElbowPlot(male.sub.object)
ggsave('./QC_filtering/MaleSub/Elbowplot.sub.male.png',
       height = 10,
       width = 10)

##graph umap
DimPlot(male.sub.object,
        reduction = "umap",
        label = TRUE,
        repel = TRUE)+
  ggtitle('SubMale')
ggsave('./QC_filtering/MaleSub/Clusters.dimplot.male.sub.png',
       height = 10,
       width = 10)

##graph UMI counts on umap
FeaturePlot(male.sub.object,
            reduction = "umap",
            features = "nFeature_RNA",
            label = TRUE,
            repel = TRUE)+
  ggtitle('SubMale') +
  labs(color = "nFeature_RNA")
ggsave('./QC_filtering/MaleSub/nFeature_RNA.dimplot.male.sub.png',
       height = 10,
       width = 10)

## variable features
VariableFeaturePlot_scCustom(male.sub.object)
ggsave('./QC_filtering/MaleSub/Variable.features.male.sub.png',
       height = 10,
       width = 10)

### Clustree
### Find clusters using a range of resolutions
# cluster across range
male.sub.object.clustree <- Seurat::FindClusters(object = male.sub.object,
                                                 resolution = seq(from = 0, to = 2, by = 0.2))
#check data
# head(male.sub.object.clustree[[]])
#clustree
clustree(male.sub.object.clustree,
         prefix = "SCT_snn_res.")+
  ggtitle('Sub clusters across resolutions')+
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/MaleSub/clustree.sub.male.png',
       height = 10,
       width = 10)
#reduced range
male.sub.object.clustree.reduce <- Seurat::FindClusters(object = male.sub.object,
                                                        resolution = seq(from = 0, to = 1.4, by = 0.2))
#clustree
clustree(male.sub.object.clustree.reduce,
         prefix = "SCT_snn_res.")+
  ggtitle('Sub clusters across resolutions reduced') +
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/MaleSub/clustree.reduced.sub.male.png',
       height = 10,
       width = 10)

#remove clustree
rm(male.sub.object.clustree)
rm(male.sub.object.clustree.reduce)

#### FemaleSub clustering ####
###filter cells first
##create object
female.sub.object <- CreateSeuratObject(counts = female.sub.data,
                                        project = "Female.Sub", 
                                        min.cells = 3,
                                        min.features = 200)

## add souporcell data
female.sub.object = AddMetaData(object = female.sub.object,
                                metadata = sub.female.clusters %>% 
                                  column_to_rownames('Cell.id')) 

## label mitochondria genes
female.sub.object <- PercentageFeatureSet(female.sub.object, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

## QC graph
VlnPlot(female.sub.object,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleSub/QCmetrics.female.sub.png',
       height = 10,
       width = 10)

#doublet
VlnPlot(subset(female.sub.object, 
               subset = Doublet == 'doublet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleSub/QCmetrics.female.sub.doublet.png',
       height = 10,
       width = 10)

#singlet
VlnPlot(subset(female.sub.object, 
               subset = Doublet == 'singlet'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleSub/QCmetrics.female.sub.singlet.png',
       height = 10,
       width = 10)


#unassigned
VlnPlot(subset(female.sub.object, 
               subset = Doublet == 'unassigned'),
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3)
ggsave('./QC_filtering/FemaleSub/QCmetrics.female.sub.unassigned.png',
       height = 10,
       width = 10)

### graph
##features and counts
FeatureScatter(female.sub.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/FemaleSub/featurescatter.female.sub.png',
       height = 10,
       width = 10)

#filter
FeatureScatter(subset(female.sub.object,
                      subset = nFeature_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 700)+
  geom_hline(yintercept = 1500) +
  annotate("text", x=0, y=700, label="700") +
  annotate("text", x=0, y=1500, label="1500")
ggsave('./QC_filtering/FemaleSub/featurescatter.female.sub.limits.png',
       height = 10,
       width = 10)

##histogram
female.sub.object@meta.data %>%
  subset(nFeature_RNA < 6000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 700 & nFeature_RNA < 1500,
                       'Kept',
                       'Remove')) %>%
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('SubFemale UMI counts')  +
  scale_fill_manual(values = c('red',
                               'black'))
ggsave('./QC_filtering/FemaleSub/UMI.counts.filtered.female.sub.png',
       height = 10,
       width = 10)

## check filters
# want ~10,000 cells per sample
#6848
female.sub.object@meta.data %>%
  subset(nFeature_RNA <= 1500) %>%
  subset(nFeature_RNA >= 700) %>%
  subset(percent.mt < 5) %>% 
  subset(Doublet == 'singlet') %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##create new object with filters
female.sub.object = subset(female.sub.object, 
                           subset = Doublet == 'singlet' & nFeature_RNA <= 1500 & nFeature_RNA >= 700 & percent.mt < 5)

#check filters worked
# should have ~6848 cells
female.sub.object@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            median = median(nFeature_RNA),
            count = n())

##features and counts
FeatureScatter(female.sub.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('./QC_filtering/FemaleSub/featurescatter.female.sub.filter.png',
       height = 10,
       width = 10)

##histogram
female.sub.object@meta.data %>%
  ggplot(aes(nFeature_RNA))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('SubFemale UMI counts')  
ggsave('./QC_filtering/FemaleSub/UMI.counts.filtered.female.sub.filter.png',
       height = 10,
       width = 10)

## run SCTransform, PCA, UMAP, and cluster 
#use resolution from clustree below: 1
female.sub.object = female.sub.object %>% 
  SCTransform() %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

## check PC's in elbow plot
ElbowPlot(female.sub.object)
ggsave('./QC_filtering/FemaleSub/Elbowplot.sub.female.png',
       height = 10,
       width = 10)

##graph umap
DimPlot(female.sub.object,
        reduction = "umap",
        label = TRUE,
        repel = TRUE)+
  ggtitle('SubFemale')
ggsave('./QC_filtering/FemaleSub/Clusters.dimplot.female.sub.png',
       height = 10,
       width = 10)

##graph UMI counts on umap
FeaturePlot(female.sub.object,
            reduction = "umap",
            features = "nFeature_RNA",
            label = TRUE,
            repel = TRUE)+
  ggtitle('SubFemale') +
  labs(color = "nFeature_RNA")
ggsave('./QC_filtering/FemaleSub/nFeature_RNA.dimplot.female.sub.png',
       height = 10,
       width = 10)

## variable features
VariableFeaturePlot_scCustom(female.sub.object)
ggsave('./QC_filtering/FemaleSub/Variable.features.female.sub.png',
       height = 10,
       width = 10)

### Clustree
### Find clusters using a range of resolutions
# cluster across range
female.sub.object.clustree <- Seurat::FindClusters(object = female.sub.object,
                                                   resolution = seq(from = 0, to = 2, by = 0.2))
#check data
# head(female.sub.object.clustree[[]])
#clustree
clustree(female.sub.object.clustree,
         prefix = "SCT_snn_res.")+
  ggtitle('Sub clusters across resolutions')+
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/FemaleSub/clustree.sub.female.png',
       height = 10,
       width = 10)
#reduced range
female.sub.object.clustree.reduce <- Seurat::FindClusters(object = female.sub.object,
                                                          resolution = seq(from = 0, to = 1.4, by = 0.2))
#clustree
clustree(female.sub.object.clustree.reduce,
         prefix = "SCT_snn_res.")+
  ggtitle('Sub clusters across resolutions reduced') +
  guides(edge_colour = FALSE,
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black",
                              high = "black")
ggsave('./QC_filtering/FemaleSub/clustree.reduced.sub.female.png',
       height = 10,
       width = 10)

#remove clustree
rm(female.sub.object.clustree)
rm(female.sub.object.clustree.reduce)

#### sctype DomMale ####

### assign cell types to clusters
# get cell-type by cell matrix
es.max.dom.male = sctype_score(scRNAseqData = male.dom.object[["SCT"]]@scale.data, 
                      scaled = TRUE, 
                      gs = gs_list$gs_positive, 
                      gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.dom.male = do.call("rbind", lapply(unique(male.dom.object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max.dom.male[ ,rownames(male.dom.object@meta.data[male.dom.object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(male.dom.object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.dom.male = cL_results.dom.male %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.dom.male$type[as.numeric(as.character(sctype_scores.dom.male$scores)) < sctype_scores.dom.male$ncells/4] = "Unknown"
print(sctype_scores.dom.male[,1:3])

#graph umap
male.dom.object@meta.data$sctype.ind = ""
for(j in unique(sctype_scores.dom.male$cluster)){
  cl_type.dom.male = sctype_scores.dom.male[sctype_scores.dom.male$cluster==j,]; 
  male.dom.object@meta.data$sctype.ind[male.dom.object@meta.data$seurat_clusters == j] = as.character(cl_type.dom.male$type[1])
}

#graph umap
DimPlot(male.dom.object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.ind') +
  ggtitle('DomMale')
ggsave('./QC_filtering/MaleDom/sctype.dimplot.male.dom.png',
       height = 10,
       width = 10)

#' ### bubble network
#' # prepare edges
#' cL_results.dom.male=cL_results.dom.male[order(cL_results.dom.male$cluster),]; edges.dom.male = cL_results.dom.male; edges.dom.male$type = paste0(edges.dom.male$type,"_",edges.dom.male$cluster); edges.dom.male$cluster = paste0("cluster ", edges.dom.male$cluster); edges.dom.male = edges.dom.male[,c("cluster", "type")]; colnames(edges.dom.male) = c("from", "to"); rownames(edges.dom.male) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1.dom.male = sctype_scores.dom.male[,c("cluster", "ncells")]; nodes_lvl1.dom.male$cluster = paste0("cluster ", nodes_lvl1.dom.male$cluster); nodes_lvl1.dom.male$Colour = "#f1f1ef"; nodes_lvl1.dom.male$ord = 1; nodes_lvl1.dom.male$realname = nodes_lvl1.dom.male$cluster; nodes_lvl1.dom.male = as.data.frame(nodes_lvl1.dom.male); nodes_lvl2.dom.male = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results.dom.male$cluster))){
#'   dt_tmp.dom.male = cL_results.dom.male[cL_results.dom.male$cluster == unique(cL_results.dom.male$cluster)[i], ]; nodes_lvl2.dom.male = rbind(nodes_lvl2.dom.male, data.frame(cluster = paste0(dt_tmp.dom.male$type,"_",dt_tmp.dom.male$cluster), ncells = dt_tmp.dom.male$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp.dom.male$type))
#' }
#' nodes.dom.male = rbind(nodes_lvl1.dom.male, nodes_lvl2.dom.male); nodes.dom.male$ncells[nodes.dom.male$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes.dom.male = merge(nodes.dom.male, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes.dom.male$shortName[is.na(nodes.dom.male$shortName)] = nodes.dom.male$realname[is.na(nodes.dom.male$shortName)]; nodes.dom.male = nodes.dom.male[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph.dom.male <- graph_from_data_frame(edges.dom.male, 
#'                                  vertices=nodes.dom.male %>% 
#'                                    filter(shortName != 'Immune system',
#'                                           shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr.dom.male<- ggraph(mygraph.dom.male, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP.dom.male = DimPlot(male.dom.object, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP.dom.male,
#'                    gggr.dom.male) +
#'   ggtitle('Dom.Male')
#' ggsave('./QC_filtering/MaleDom/sctype.bubbleplot.male.dom.png',
#'        height = 10,
#'        width = 10)


#### sctype DomFemale ####

### assign cell types to clusters
# get cell-type by cell matrix
es.max.dom.female = sctype_score(scRNAseqData = female.dom.object[["SCT"]]@scale.data, 
                               scaled = TRUE, 
                               gs = gs_list$gs_positive, 
                               gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.dom.female = do.call("rbind", lapply(unique(female.dom.object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max.dom.female[ ,rownames(female.dom.object@meta.data[female.dom.object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(female.dom.object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.dom.female = cL_results.dom.female %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.dom.female$type[as.numeric(as.character(sctype_scores.dom.female$scores)) < sctype_scores.dom.female$ncells/4] = "Unknown"
print(sctype_scores.dom.female[,1:3])

#graph umap
female.dom.object@meta.data$sctype.ind = ""
for(j in unique(sctype_scores.dom.female$cluster)){
  cl_type.dom.female = sctype_scores.dom.female[sctype_scores.dom.female$cluster==j,]; 
  female.dom.object@meta.data$sctype.ind[female.dom.object@meta.data$seurat_clusters == j] = as.character(cl_type.dom.female$type[1])
}

#graph umap
DimPlot(female.dom.object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.ind') +
  ggtitle('DomFemale')
ggsave('./QC_filtering/FemaleDom/sctype.dimplot.female.dom.png',
       height = 10,
       width = 10)

#' ### bubble network
#' # prepare edges
#' cL_results.dom.female=cL_results.dom.female[order(cL_results.dom.female$cluster),]; edges.dom.female = cL_results.dom.female; edges.dom.female$type = paste0(edges.dom.female$type,"_",edges.dom.female$cluster); edges.dom.female$cluster = paste0("cluster ", edges.dom.female$cluster); edges.dom.female = edges.dom.female[,c("cluster", "type")]; colnames(edges.dom.female) = c("from", "to"); rownames(edges.dom.female) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1.dom.female = sctype_scores.dom.female[,c("cluster", "ncells")]; nodes_lvl1.dom.female$cluster = paste0("cluster ", nodes_lvl1.dom.female$cluster); nodes_lvl1.dom.female$Colour = "#f1f1ef"; nodes_lvl1.dom.female$ord = 1; nodes_lvl1.dom.female$realname = nodes_lvl1.dom.female$cluster; nodes_lvl1.dom.female = as.data.frame(nodes_lvl1.dom.female); nodes_lvl2.dom.female = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results.dom.female$cluster))){
#'   dt_tmp.dom.female = cL_results.dom.female[cL_results.dom.female$cluster == unique(cL_results.dom.female$cluster)[i], ]; nodes_lvl2.dom.female = rbind(nodes_lvl2.dom.female, data.frame(cluster = paste0(dt_tmp.dom.female$type,"_",dt_tmp.dom.female$cluster), ncells = dt_tmp.dom.female$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp.dom.female$type))
#' }
#' nodes.dom.female = rbind(nodes_lvl1.dom.female, nodes_lvl2.dom.female); nodes.dom.female$ncells[nodes.dom.female$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes.dom.female = merge(nodes.dom.female, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes.dom.female$shortName[is.na(nodes.dom.female$shortName)] = nodes.dom.female$realname[is.na(nodes.dom.female$shortName)]; nodes.dom.female = nodes.dom.female[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph.dom.female <- graph_from_data_frame(edges.dom.female, 
#'                                           vertices=nodes.dom.female %>% 
#'                                             filter(shortName != 'Immune system',
#'                                                    shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr.dom.female<- ggraph(mygraph.dom.female, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP.dom.female = DimPlot(female.dom.object, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP.dom.female,
#'                    gggr.dom.female) +
#'   ggtitle('Dom.Female')
#' ggsave('./QC_filtering/FemaleDom/sctype.bubbleplot.female.dom.png',
#'        height = 10,
#'        width = 10)

#### sctype SubMale ####

### assign cell types to clusters
# get cell-type by cell matrix
es.max.sub.male = sctype_score(scRNAseqData = male.sub.object[["SCT"]]@scale.data, 
                               scaled = TRUE, 
                               gs = gs_list$gs_positive, 
                               gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.sub.male = do.call("rbind", lapply(unique(male.sub.object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max.sub.male[ ,rownames(male.sub.object@meta.data[male.sub.object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(male.sub.object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.sub.male = cL_results.sub.male %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.sub.male$type[as.numeric(as.character(sctype_scores.sub.male$scores)) < sctype_scores.sub.male$ncells/4] = "Unknown"
print(sctype_scores.sub.male[,1:3])

#graph umap
male.sub.object@meta.data$sctype.ind = ""
for(j in unique(sctype_scores.sub.male$cluster)){
  cl_type.sub.male = sctype_scores.sub.male[sctype_scores.sub.male$cluster==j,]; 
  male.sub.object@meta.data$sctype.ind[male.sub.object@meta.data$seurat_clusters == j] = as.character(cl_type.sub.male$type[1])
}

#graph umap
DimPlot(male.sub.object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.ind')  +
  ggtitle('SubMale')
ggsave('./QC_filtering/MaleSub/sctype.dimplot.male.sub.png',
       height = 10,
       width = 10)

#' ### bubble network
#' # prepare edges
#' cL_results.sub.male=cL_results.sub.male[order(cL_results.sub.male$cluster),]; edges.sub.male = cL_results.sub.male; edges.sub.male$type = paste0(edges.sub.male$type,"_",edges.sub.male$cluster); edges.sub.male$cluster = paste0("cluster ", edges.sub.male$cluster); edges.sub.male = edges.sub.male[,c("cluster", "type")]; colnames(edges.sub.male) = c("from", "to"); rownames(edges.sub.male) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1.sub.male = sctype_scores.sub.male[,c("cluster", "ncells")]; nodes_lvl1.sub.male$cluster = paste0("cluster ", nodes_lvl1.sub.male$cluster); nodes_lvl1.sub.male$Colour = "#f1f1ef"; nodes_lvl1.sub.male$ord = 1; nodes_lvl1.sub.male$realname = nodes_lvl1.sub.male$cluster; nodes_lvl1.sub.male = as.data.frame(nodes_lvl1.sub.male); nodes_lvl2.sub.male = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results.sub.male$cluster))){
#'   dt_tmp.sub.male = cL_results.sub.male[cL_results.sub.male$cluster == unique(cL_results.sub.male$cluster)[i], ]; nodes_lvl2.sub.male = rbind(nodes_lvl2.sub.male, data.frame(cluster = paste0(dt_tmp.sub.male$type,"_",dt_tmp.sub.male$cluster), ncells = dt_tmp.sub.male$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp.sub.male$type))
#' }
#' nodes.sub.male = rbind(nodes_lvl1.sub.male, nodes_lvl2.sub.male); nodes.sub.male$ncells[nodes.sub.male$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes.sub.male = merge(nodes.sub.male, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes.sub.male$shortName[is.na(nodes.sub.male$shortName)] = nodes.sub.male$realname[is.na(nodes.sub.male$shortName)]; nodes.sub.male = nodes.sub.male[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph.sub.male <- graph_from_data_frame(edges.sub.male, 
#'                                           vertices=nodes.sub.male %>% 
#'                                             filter(shortName != 'Immune system',
#'                                                    shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr.sub.male<- ggraph(mygraph.sub.male, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP.sub.male = DimPlot(male.sub.object, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP.sub.male,
#'                    gggr.sub.male) +
#'   ggtitle('Sub.Male')
#' ggsave('./QC_filtering/MaleSub/sctype.bubbleplot.male.sub.png',
#'        height = 10,
#'        width = 10)

#### sctype SubFemale ####

### assign cell types to clusters
# get cell-type by cell matrix
es.max.sub.female = sctype_score(scRNAseqData = female.sub.object[["SCT"]]@scale.data, 
                                 scaled = TRUE, 
                                 gs = gs_list$gs_positive, 
                                 gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.sub.female = do.call("rbind", lapply(unique(female.sub.object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max.sub.female[ ,rownames(female.sub.object@meta.data[female.sub.object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(female.sub.object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.sub.female = cL_results.sub.female %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.sub.female$type[as.numeric(as.character(sctype_scores.sub.female$scores)) < sctype_scores.sub.female$ncells/4] = "Unknown"
print(sctype_scores.sub.female[,1:3])

#graph umap
female.sub.object@meta.data$sctype.ind = ""
for(j in unique(sctype_scores.sub.female$cluster)){
  cl_type.sub.female = sctype_scores.sub.female[sctype_scores.sub.female$cluster==j,]; 
  female.sub.object@meta.data$sctype.ind[female.sub.object@meta.data$seurat_clusters == j] = as.character(cl_type.sub.female$type[1])
}

#graph umap
DimPlot(female.sub.object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.ind') +
  ggtitle('SubFemale')
ggsave('./QC_filtering/FemaleSub/sctype.dimplot.female.sub.png',
       height = 10,
       width = 10)
#' 
#' ### bubble network
#' # prepare edges
#' cL_results.sub.female=cL_results.sub.female[order(cL_results.sub.female$cluster),]; edges.sub.female = cL_results.sub.female; edges.sub.female$type = paste0(edges.sub.female$type,"_",edges.sub.female$cluster); edges.sub.female$cluster = paste0("cluster ", edges.sub.female$cluster); edges.sub.female = edges.sub.female[,c("cluster", "type")]; colnames(edges.sub.female) = c("from", "to"); rownames(edges.sub.female) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1.sub.female = sctype_scores.sub.female[,c("cluster", "ncells")]; nodes_lvl1.sub.female$cluster = paste0("cluster ", nodes_lvl1.sub.female$cluster); nodes_lvl1.sub.female$Colour = "#f1f1ef"; nodes_lvl1.sub.female$ord = 1; nodes_lvl1.sub.female$realname = nodes_lvl1.sub.female$cluster; nodes_lvl1.sub.female = as.data.frame(nodes_lvl1.sub.female); nodes_lvl2.sub.female = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results.sub.female$cluster))){
#'   dt_tmp.sub.female = cL_results.sub.female[cL_results.sub.female$cluster == unique(cL_results.sub.female$cluster)[i], ]; nodes_lvl2.sub.female = rbind(nodes_lvl2.sub.female, data.frame(cluster = paste0(dt_tmp.sub.female$type,"_",dt_tmp.sub.female$cluster), ncells = dt_tmp.sub.female$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp.sub.female$type))
#' }
#' nodes.sub.female = rbind(nodes_lvl1.sub.female, nodes_lvl2.sub.female); nodes.sub.female$ncells[nodes.sub.female$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes.sub.female = merge(nodes.sub.female, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes.sub.female$shortName[is.na(nodes.sub.female$shortName)] = nodes.sub.female$realname[is.na(nodes.sub.female$shortName)]; nodes.sub.female = nodes.sub.female[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph.sub.female <- graph_from_data_frame(edges.sub.female, 
#'                                             vertices=nodes.sub.female %>% 
#'                                               filter(shortName != 'Immune system',
#'                                                      shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr.sub.female<- ggraph(mygraph.sub.female, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP.sub.female = DimPlot(female.sub.object, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP.sub.female,
#'                    gggr.sub.female) +
#'   ggtitle('Sub.Female')
#' ggsave('./QC_filtering/FemaleSub/sctype.bubbleplot.female.sub.png',
#'        height = 10,
#'        width = 10)
#### integrate all data ####
### integrate data
## https://satijalab.org/seurat/articles/integration_introduction.html

## create list
mouse.snseq.list <- list(Male.Dom = male.dom.object, 
                         Female.Dom = female.dom.object,
                         Male.Sub = male.sub.object, 
                         Female.Sub = female.sub.object)

# ##create all gene list
# #19,518 genes
# mouse.snseq.genes <- lapply(mouse.snseq.list,
#                             row.names) %>%
#   Reduce(intersect, .)
# # length(mouse.snseq.genes)

##select features
mouse.snseq.features <- SelectIntegrationFeatures(object.list = mouse.snseq.list,
                                                  nfeatures = 3000)
## identify anchors
mouse.snseq.list <- PrepSCTIntegration(object.list = mouse.snseq.list,
                                       anchor.features = mouse.snseq.features)
## find anchors
# takes a long time to do every pairwise comparison
mouse.snseq.anchors <- FindIntegrationAnchors(object.list = mouse.snseq.list,
                                              normalization.method = "SCT",
                                              anchor.features = mouse.snseq.features,
                                              reduction = "rpca", #runs faster integration
                                              k.anchor = 20) # increase strength of alignment



## integrate data
#retain all genes
mouse.snseq.combined.sct <- IntegrateData(anchorset = mouse.snseq.anchors,
                                          normalization.method = "SCT",
                                          # features.to.integrate = mouse.snseq.genes,
                                          dims = 1:30)
# Integrating data
# Merging dataset 3 into 4 1
# Extracting anchors for merged samples
# Finding integration vectors
# Error in .T2C(newTMat(i = c(ij1[, 1], ij2[, 1]), j = c(ij1[, 2], ij2[,  : 
#                                                                        unable to coerce from TsparseMatrix to [CR]sparseMatrixwhen length of 'i' slot exceeds 2^31-1
#                                                                      Attaching SeuratObject

## Run PCA
mouse.snseq.combined.sct <- RunPCA(mouse.snseq.combined.sct, 
                                   verbose = FALSE)
#check elbowplot
ElbowPlot(mouse.snseq.combined.sct)
ggsave('./QC_filtering/Integrated/Elbowplot.integrated.png',
       height = 10,
       width = 10)

## run UMAP
mouse.snseq.combined.sct <- RunUMAP(mouse.snseq.combined.sct,
                                    reduction = "pca", 
                                    dims = 1:30)


### normalize SCT across samples
## rescale SCT assay for all samples
## subset genes
mouse.snseq.combined.sct = PrepSCTFindMarkers(mouse.snseq.combined.sct,
                                              assay = 'SCT',
                                              verbose = T)

# find neighbors
mouse.snseq.combined.sct <- FindNeighbors(mouse.snseq.combined.sct, 
                                          reduction = "pca", 
                                          dims = 1:30)
# find clusters
mouse.snseq.combined.sct <- FindClusters(mouse.snseq.combined.sct, 
                                         resolution = 0.4)

### graph
##dom vs sub
DimPlot(mouse.snseq.combined.sct,
        reduction = "umap",
        group.by = "orig.ident")
ggsave('./QC_filtering/Integrated/Dimplot.ident.integrated.png',
       height = 10,
       width = 10)

## clusters
DimPlot(mouse.snseq.combined.sct, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('./QC_filtering/Integrated/Clusters.dimplot.integrated.png',
       height = 10,
       width = 10)

#### clustree integrated ####
# Select a range of resolutions
#reduced range
resolution.range.reduced.2 <- seq(from = 0, to = 1, by = 0.2)
# cluster across resolutions
mouse.snseq.combined.sct.clustree <- Seurat::FindClusters(object = mouse.snseq.combined.sct,                                                  resolution = resolution.range.reduced.2)
#check data
head(mouse.snseq.combined.sct.clustree[[]])

#clustree
clustree(mouse.snseq.combined.sct.clustree, 
         prefix = "integrated_snn_res.",
         node_size_range = c(10,20),
         node_colour = 'cluster') +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('./QC_filtering/Integrated/Combined.clustree.png',
       width = 10,
       height = 10)

## check cells per group per cluster
#get table
Cells.per.cluster.table = table(mouse.snseq.combined.sct@active.ident, 
                                mouse.snseq.combined.sct@meta.data$orig.ident) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>%
  pivot_longer(cols = c("Male.Dom",
                        "Male.Sub",
                        "Female.Dom",
                        "Female.Sub"),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total)

#graph
Cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('./QC_filtering/Integrated/Cells.per.cluster.DomvsSub.png',
       width = 10,
       height = 10)

rm(mouse.snseq.combined.sct.clustree)

#### create umap plots for presentation  ####
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#empty cluster
DimPlot(mouse.snseq.combined.sct,
        reduction = "umap",
        cols = c(rep('grey',
                     43)),
        pt.size = 1) +
  theme(legend.position = 'none')
ggsave('./QC_filtering/Integrated/Empty.dimplot.comparison.all.png',
       width = 10,
       height = 10)
## presentation simple cluster
DimPlot(mouse.snseq.combined.sct,
        reduction = "umap",
        label = T,
        pt.size = 1,
        label.size = 10,
        repel = T) +
  theme(legend.position = 'none')
ggsave('./QC_filtering/Integrated/Clusters.dimplot.comparison.all.png',
       width = 10,
       height = 10)

## presentation across samples cluster
DimPlot(mouse.snseq.combined.sct,
        reduction = "umap",
        split.by = 'orig.ident',
        pt.size = 1) +
  theme(legend.position = 'none')
ggsave('./QC_filtering/Integrated/DomvsSub.clusters.dimplot.comparison.all.png',
       width = 20,
       height = 10)
## presentation across samples genotype
DimPlot(mouse.snseq.combined.sct,
        reduction = "umap",
        split.by = 'orig.ident',
        group.by = 'Genotype',
        pt.size = 1) 
ggsave('./QC_filtering/Integrated/DomvsSub.genotype.dimplot.comparison.all.png',
       width = 20,
       height = 10)


## presentation across samples nfeature
FeaturePlot(mouse.snseq.combined.sct,
        reduction = "umap",
        split.by = 'orig.ident',
        features = 'nFeature_SCT',
        pt.size = 1,
        max.cutoff = ) 
ggsave('./QC_filtering/Integrated/DomvsSub.nfeature.dimplot.comparison.all.png',
       width = 20,
       height = 10)

#### sctype Integrated ####

### assign cell types to clusters
# get cell-type by cell matrix
es.max.integrated = sctype_score(scRNAseqData = mouse.snseq.combined.sct[["integrated"]]@scale.data, 
                               scaled = TRUE, 
                               gs = gs_list$gs_positive, 
                               gs2 = gs_list$gs_negative,
                               gene_names_to_uppercase = TRUE) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.integrated = do.call("rbind", lapply(unique(mouse.snseq.combined.sct@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max.integrated[ ,rownames(mouse.snseq.combined.sct@meta.data[mouse.snseq.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(mouse.snseq.combined.sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.integrated = cL_results.integrated %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.integrated$type[as.numeric(as.character(sctype_scores.integrated$scores)) < sctype_scores.integrated$ncells/4] = "Unknown"
print(sctype_scores.integrated[,1:3])

#graph umap
mouse.snseq.combined.sct@meta.data$sctype.all = ""
for(j in unique(sctype_scores.integrated$cluster)){
  cl_type.integrated = sctype_scores.integrated[sctype_scores.integrated$cluster==j,]; 
  mouse.snseq.combined.sct@meta.data$sctype.all[mouse.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type.integrated$type[1])
}

#graph umap
DimPlot(mouse.snseq.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.all') +
  ggtitle('Integrated')
ggsave('./QC_filtering/Integrated/sctype.dimplot.integrated.png',
       height = 10,
       width = 10)

#' ### bubble network
#' # prepare edges
#' cL_results.integrated=cL_results.integrated[order(cL_results.integrated$cluster),]; edges.integrated = cL_results.integrated; edges.integrated$type = paste0(edges.integrated$type,"_",edges.integrated$cluster); edges.integrated$cluster = paste0("cluster ", edges.integrated$cluster); edges.integrated = edges.integrated[,c("cluster", "type")]; colnames(edges.integrated) = c("from", "to"); rownames(edges.integrated) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1.integrated = sctype_scores.integrated[,c("cluster", "ncells")]; nodes_lvl1.integrated$cluster = paste0("cluster ", nodes_lvl1.integrated$cluster); nodes_lvl1.integrated$Colour = "#f1f1ef"; nodes_lvl1.integrated$ord = 1; nodes_lvl1.integrated$realname = nodes_lvl1.integrated$cluster; nodes_lvl1.integrated = as.data.frame(nodes_lvl1.integrated); nodes_lvl2.integrated = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results.integrated$cluster))){
#'   dt_tmp.integrated = cL_results.integrated[cL_results.integrated$cluster == unique(cL_results.integrated$cluster)[i], ]; nodes_lvl2.integrated = rbind(nodes_lvl2.integrated, data.frame(cluster = paste0(dt_tmp.integrated$type,"_",dt_tmp.integrated$cluster), ncells = dt_tmp.integrated$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp.integrated$type))
#' }
#' nodes.integrated = rbind(nodes_lvl1.integrated, nodes_lvl2.integrated); nodes.integrated$ncells[nodes.integrated$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes.integrated = merge(nodes.integrated, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes.integrated$shortName[is.na(nodes.integrated$shortName)] = nodes.integrated$realname[is.na(nodes.integrated$shortName)]; nodes.integrated = nodes.integrated[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph.integrated <- graph_from_data_frame(edges.integrated, 
#'                                  vertices=nodes.integrated %>% 
#'                                    filter(shortName != 'Immune system',
#'                                           shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr.integrated<- ggraph(mygraph.integrated, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP.integrated = DimPlot(mouse.snseq.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP.integrated,
#'                    gggr.integrated) +
#'   ggtitle('integrated')
#' ggsave('./QC_filtering/Integrated/sctype.bubbleplot.male.dom.png',
#'        height = 10,
#'        width = 10)
#'        


### compare individual to all
## graph alluvial plot
mouse.snseq.combined.sct@meta.data %>%
  select(c(sctype.all,
           sctype.ind)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(sctype.all, -Freq.scale),
             axis2 = reorder(sctype.ind,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.all)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all",
                              "sctype.ind"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./QC_filtering/Integrated/Alluvial cells by sctype.ind vs sctype.all.png',
       height = 10,
       width = 10)

## graph by individual 
for (i in unique(mouse.snseq.combined.sct@meta.data$orig.ident)) {
  mouse.snseq.combined.sct@meta.data %>%
    filter(orig.ident == i) %>% 
    select(c(sctype.all,
             sctype.ind)) %>%
    table() %>%
    as.data.frame() %>%
    mutate(Count = sum(Freq)) %>%
    ungroup() %>%
    mutate(Freq.scale = 100*Freq/Count) %>%
    ggplot(aes(axis1 = reorder(sctype.all, -Freq.scale),
               axis2 = reorder(sctype.ind,-Freq.scale),
               y = Freq.scale)) +
    geom_alluvium(aes(fill = sctype.all)) +
    geom_stratum() +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("sctype.all",
                                "sctype.ind"),
                     expand = c(0.15, 0.05)) +
    scale_fill_viridis_d() +
    theme_classic() +
    ggtitle(paste(i))
  ggsave(paste('./QC_filtering/Integrated/Alluvial cells by sctype.ind vs sctype.all ',
               i,
               '.png',
               sep = ''),
         height = 10,
         width = 10)
}

## graph alluvial plot
mouse.snseq.combined.sct@meta.data %>% 
  select(c(integrated_snn_res.0.4,
           sctype.all,
           sctype.ind)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(sctype.all, -Freq.scale),
             axis2 = reorder(integrated_snn_res.0.4, -Freq.scale),
             axis3 = reorder(sctype.ind,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.all)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all",
                              "sctype.ind"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./QC_filtering/Integrated/Alluvial cells by sctype.ind vs sctype.all cluster.png',
       height = 10,
       width = 10)


## counts per sample
mouse.snseq.combined.sct@meta.data %>% 
  select(c(orig.ident,
           sctype.all)) %>%
  table() %>%
  as.data.frame() %>% 
  ggplot(aes(x = reorder(sctype.all,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point() +
  theme_classic() +
  xlab('')+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave('./QC_filtering/Integrated/Cell count sctype.all per sample.png',
       height = 10,
       width = 10)
#### save data point ####
##just single cell object
# save(mouse.snseq.combined.sct,
#      file = "mouse.snseq.combined.sct.RData")
# load('mouse.snseq.combined.sct.RData')

#### sctype HypoMap Integrated [stopped here] ####
### load HypoMap annotation
annotation.data.hypomap = read.csv("./Cell.type.gene.markers.HypoMap.csv")

# check genes in data
# 193/306
mouse.snseq.combined.sct@assays$SCT@scale.data %>% 
  rownames() %>% 
intersect(annotation.data.hypomap$gene) %>% 
  length()

## reduce to genes in data
annotation.data.hypomap.reduce= annotation.data.hypomap %>% 
  filter(gene %in% rownames(mouse.snseq.combined.sct@assays$SCT@scale.data))

# check reduced gene list
annotation.data.hypomap %>% 
  dplyr::select(cluster_name) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_col(aes(x = reorder(cluster_name, -Freq),
               y = Freq),
           fill = 'black') + 
  geom_col(data = annotation.data.hypomap.reduce %>% 
             dplyr::select(cluster_name) %>% 
             table() %>% 
             as.data.frame() ,
           aes(x = reorder(cluster_name, -Freq),
               y = Freq),
           fill = 'red') +
  theme_classic() +
  ggtitle("HypoMap genes per cell type") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('./QC_filtering/Integrated/HypoMap genes per cell type.png',
       height = 10,
       width = 10)


## remove NA
#get distinct 
annotation.data.hypomap = annotation.data.hypomap %>% 
  na.omit() %>% 
  dplyr::select(c(gene,
                  cluster_name)) %>% 
  dplyr::rename('Type' = 'cluster_name') %>% 
  dplyr::rename('Marker' = 'gene') %>% 
  distinct() 

## collapse dataframe to match format needed
annotation.data.hypomap = annotation.data.hypomap %>% 
  group_by(Type) %>% 
  summarise(geneSymbolmore1 = paste0(Marker, 
                                     collapse = ",")) %>% 
  rename(cellName = Type) %>% 
  mutate(tissueType = 'Brain',
         geneSymbolmore2 = 0) %>% 
  relocate(cellName, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore1, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore2, 
           .after = last_col()) %>% 
  as.data.frame()

# set tissue type
tissue = "Brain"
# prepare gene sets
gs_list.hypomap = gene_sets_prepare.df(annotation.data.hypomap,
                               tissue)

### assign cell types to clusters
# get cell-type by cell matrix
es.max.integrated.hypomap = sctype_score(scRNAseqData = mouse.snseq.combined.sct[["SCT"]]@scale.data, 
                                 scaled = TRUE, 
                                 gs = gs_list.hypomap$gs_positive) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.integrated.hypomap = do.call("rbind", lapply(unique(mouse.snseq.combined.sct@meta.data$seurat_clusters), function(cl){
  es.max.cl.hypomap = sort(rowSums(es.max.integrated.hypomap[ ,rownames(mouse.snseq.combined.sct@meta.data[mouse.snseq.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl.hypomap), scores = es.max.cl.hypomap, ncells = sum(mouse.snseq.combined.sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.integrated.hypomap = cL_results.integrated.hypomap %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.integrated.hypomap$type[as.numeric(as.character(sctype_scores.integrated.hypomap$scores)) < sctype_scores.integrated.hypomap$ncells/4] = "Unknown"
print(sctype_scores.integrated.hypomap[,1:3])

#graph umap
mouse.snseq.combined.sct@meta.data$sctype.all.hypomap = ""
for(j in unique(sctype_scores.integrated.hypomap$cluster)){
  cl_type.integrated.hypomap = sctype_scores.integrated.hypomap[sctype_scores.integrated.hypomap$cluster==j,]; 
  mouse.snseq.combined.sct@meta.data$sctype.all.hypomap[mouse.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type.integrated.hypomap$type[1])
}

#graph umap
DimPlot(mouse.snseq.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.all.hypomap') +
  ggtitle('Integrated')
ggsave('./QC_filtering/Integrated/sctype.dimplot.integrated hypomap.png',
       height = 10,
       width = 10)


#### sctype HypoMap Integrated expanded ####
### load HypoMap annotation
annotation.data.hypomap = read.csv("./Cell.type.gene.markers.HypoMap.expanded.csv")

# check genes in data
# 458/1044
mouse.snseq.combined.sct@assays$SCT@scale.data %>% 
  rownames() %>% 
  intersect(annotation.data.hypomap$gene) %>% 
  length()

## reduce to genes in data
annotation.data.hypomap.reduce= annotation.data.hypomap %>% 
  filter(gene %in% rownames(mouse.snseq.combined.sct@assays$SCT@scale.data))

# check reduced gene list
annotation.data.hypomap %>% 
  dplyr::select(cluster_name) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename("Cell.type" = ".") %>% 
  ggplot() +
  geom_col(aes(x = reorder(Cell.type, -Freq),
               y = Freq),
           fill = 'black') + 
  geom_col(data = annotation.data.hypomap.reduce %>% 
             dplyr::select(cluster_name) %>% 
             table() %>% 
             as.data.frame() %>% 
             dplyr::rename("Cell.type" = "."),
           aes(x = reorder(Cell.type, -Freq),
               y = Freq),
           fill = 'red') +
  theme_classic() +
  ggtitle("HypoMap genes per cell type expanded") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('./QC_filtering/Integrated/HypoMap genes per cell type expanded.png',
       height = 10,
       width = 10)

## remove NA
#get distinct 
annotation.data.hypomap = annotation.data.hypomap %>% 
  na.omit() %>% 
  dplyr::select(c(gene,
                  cluster_name)) %>% 
  dplyr::rename('Type' = 'cluster_name') %>% 
  dplyr::rename('Marker' = 'gene') %>% 
  distinct() 

## collapse dataframe to match format needed
annotation.data.hypomap = annotation.data.hypomap %>% 
  group_by(Type) %>% 
  summarise(geneSymbolmore1 = paste0(Marker, 
                                     collapse = ",")) %>% 
  rename(cellName = Type) %>% 
  mutate(tissueType = 'Brain',
         geneSymbolmore2 = 0) %>% 
  relocate(cellName, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore1, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore2, 
           .after = last_col()) %>% 
  as.data.frame()

# set tissue type
tissue = "Brain"
# prepare gene sets
gs_list.hypomap = gene_sets_prepare.df(annotation.data.hypomap,
                                       tissue)

### assign cell types to clusters
# get cell-type by cell matrix
es.max.integrated.hypomap = sctype_score(scRNAseqData = mouse.snseq.combined.sct[["SCT"]]@scale.data, 
                                         scaled = TRUE, 
                                         gs = gs_list.hypomap$gs_positive,
                                         gs2 = gs_list.hypomap$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results.integrated.hypomap = do.call("rbind", lapply(unique(mouse.snseq.combined.sct@meta.data$seurat_clusters), function(cl){
  es.max.cl.hypomap = sort(rowSums(es.max.integrated.hypomap[ ,rownames(mouse.snseq.combined.sct@meta.data[mouse.snseq.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl.hypomap), scores = es.max.cl.hypomap, ncells = sum(mouse.snseq.combined.sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores.integrated.hypomap = cL_results.integrated.hypomap %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores.integrated.hypomap$type[as.numeric(as.character(sctype_scores.integrated.hypomap$scores)) < sctype_scores.integrated.hypomap$ncells/4] = "Unknown"
print(sctype_scores.integrated.hypomap[,1:3])

#graph umap
mouse.snseq.combined.sct@meta.data$sctype.all.hypomap.exp = ""
for(j in unique(sctype_scores.integrated.hypomap$cluster)){
  cl_type.integrated.hypomap = sctype_scores.integrated.hypomap[sctype_scores.integrated.hypomap$cluster==j,]; 
  mouse.snseq.combined.sct@meta.data$sctype.all.hypomap.exp[mouse.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type.integrated.hypomap$type[1])
}

#graph umap
DimPlot(mouse.snseq.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.all.hypomap.exp') +
  ggtitle('Integrated expanded')
ggsave('./QC_filtering/Integrated/sctype.dimplot.integrated hypomap expanded.png',
       height = 10,
       width = 10)


#### graph alluvial plot ####
### hypomap vs hypomap expanded
mouse.snseq.combined.sct@meta.data %>%
  select(c(sctype.all.hypomap,
           sctype.all.hypomap.exp)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis2 = reorder(sctype.all.hypomap, -Freq.scale),
             axis3 = reorder(sctype.all.hypomap.exp,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.all.hypomap)) +
  geom_stratum() +
  geom_text_repel(direction = "x",
                  stat = "stratum",
                  aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all.hypomap",
                              "sctype.all.hypomap.exp"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./QC_filtering/Integrated/Alluvial cells by sctype.all.hypomap vs sctype.all.hypomap.exp.png',
       height = 10,
       width = 10)


### hypomap vs sctype.all
mouse.snseq.combined.sct@meta.data %>%
  select(c(sctype.all.hypomap,
           sctype.all)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis2 = reorder(sctype.all.hypomap, -Freq.scale),
             axis3 = reorder(sctype.all,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.all.hypomap)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all.hypomap",
                              "sctype.all"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./QC_filtering/Integrated/Alluvial cells by sctype.all.hypomap vs sctype.all.png',
       height = 10,
       width = 10)

### hypomap.exp vs sctype.all
mouse.snseq.combined.sct@meta.data %>%
  select(c(sctype.all.hypomap.exp,
           sctype.all)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis2 = reorder(sctype.all.hypomap.exp, -Freq.scale),
             axis3 = reorder(sctype.all,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.all.hypomap.exp)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctype.all.hypomap.exp",
                              "sctype.all"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./QC_filtering/Integrated/Alluvial cells by sctype.all.hypomap.exp vs sctype.all.png',
       height = 10,
       width = 10)

#### markers per clusters ####
### top markers per cluster
##find all markers
mouse.snseq.combined.sct.markers <- FindAllMarkers(mouse.snseq.combined.sct, 
                                                   only.pos = TRUE, 
                                                   min.pct = 0.25, 
                                                   logfc.threshold = 0.25)
# check top 2 for each cluster
# mouse.snseq.combined.sct.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, 
#             order_by = avg_log2FC)

##graph
#identify top 10 per cluster
mouse.snseq.combined.sct.markers.top10 = mouse.snseq.combined.sct.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, 
        wt = avg_log2FC) %>% 
  mutate(cluster.position = as.numeric(as.character(cluster))) %>% 
  arrange(cluster.position)
#graph heatmap
DoHeatmap(mouse.snseq.combined.sct,
          features = mouse.snseq.combined.sct.markers.top10$gene) +
  NoLegend()
ggsave('comparison/Top 10 marker genes per cluster.png')

#ordered?
DoHeatmap(mouse.snseq.combined.sct,
          features = mouse.snseq.combined.sct.markers.top10 %>% 
            pull(gene)) +
  NoLegend()
ggsave('comparison/Top 10 marker genes per cluster order.png')


#### identify clusters all genes ####

###View full gene list
# mouse.snseq.all.genes %>%
#   View()

### check overlap of cluster markers and marker genes
## merge cell markers with cluster markers
#   mouse.snseq.combined.sct.markers.celltype = full_join(mouse.snseq.combined.sct.markers,
#                                                         marker.genes.list %>%
#                                                           rename(gene = gene.id.tilapia))
# 
# ## remove all NA
# mouse.snseq.combined.sct.markers.celltype.subset = mouse.snseq.combined.sct.markers.celltype %>%
#   na.omit()
#   
# 
# #check number of cells expressing gene above threshold
# sum(GetAssayData(object = mouse.snseq.combined.sct.all, 
#                  slot = "data")["avp",]>5)
# 
# ###check top 10 markers per cluster
# #remove ensemble gene names
# mouse.snseq.combined.sct.markers.top10 %>% 
#   filter(!str_detect(gene,
#                      "ENSONIG")) %>% 
#   View()
# 
# full_join(mouse.snseq.combined.sct.markers.top10,
#           marker.genes.list %>% 
#             rename(gene = gene.id.tilapia)) %>% 
#   na.omit() %>% View()
# 
# ## heat map of markers per cluster
# #reorder levels
# levels(mouse.snseq.combined.sct.all) <- factor(0:23)
# #reorder markers
# top20 <- mouse.snseq.combined.sct.markers.celltype %>% 
#   group_by(cluster) %>% 
#   top_n(20, 
#         avg_logFC)
# 
# # across all clusters
# # DoHeatmap(mouse.snseq.combined.sct.all, 
# #           features = mouse.snseq.combined.sct.markers.celltype %>% 
# #             group_by(cluster) %>% 
# #             arrange(desc(avg_log2FC)) %>% 
# #             slice_head(n = 20) %>% 
# #             pull(gene), 
# #           size = 3)
# 
# #cluster specific
# DoHeatmap(mouse.snseq.combined.sct.all, 
#           features = mouse.snseq.combined.sct.markers.celltype.subset %>% 
#             pull(gene), 
#           size = 3)
# 
# #all markers
# DoHeatmap(mouse.snseq.combined.sct.all, 
#           features = marker.genes.list %>% 
#             pull(gene.id.tilapia), 
#           size = 3)


## make list of broad cell types
cell.types = marker.genes.list %>% 
  filter(cell.subtype.category == 'broad') %>% 
  pull(cell.subtype) %>% 
  unique()

## create marker score for each broad cell type
## Takes a while to fully run
for (x in cell.types) {
  #create marker score 
  mouse.snseq.combined.sct.all <- AddModuleScore(object = mouse.snseq.combined.sct.all, 
                                                 features = list(c(marker.genes.list %>% 
                                                                     filter(cell.subtype.category == 'broad') %>% 
                                                                     filter(cell.subtype == x) %>% 
                                                                     pull(gene.id.tilapia))), 
                                                 name = paste(x,
                                                              "score",
                                                              sep = '.'))
}



## graph cell type marker score
for (x in cell.types) {
  #graph cell type umap
  FeaturePlot(object = mouse.snseq.combined.sct.all, 
              features = paste(x,
                               ".score1",
                               sep =''),
              min.cutoff = 1,
              max.cutoff = 2,
              pt.size = 0.5,
              order = TRUE,
              cols = c("grey",
                       "red")) + 
    ggtitle(paste(x,
                  "marker score",
                  sep = ' '), 
            subtitle = paste("genes =", 
                             marker.genes.list %>% 
                               filter(cell.subtype.category == 'broad') %>% 
                               filter(cell.subtype == x) %>% 
                               pull(gene.id.tilapia) %>% 
                               length))
  ggsave(paste('cell.markers/marker.umap/',
               x,
               ' marker score umap.png',
               sep = ''),
         width = 10,
         height = 10)
  
  #create dotplot for markers of each cell type 
  DotPlot(object = mouse.snseq.combined.sct.all, 
          features = marker.genes.list %>% 
            filter(cell.subtype.category == 'broad') %>% 
            filter(cell.subtype == x) %>% 
            pull(gene.id.tilapia),
          col.min = 1,
          dot.min = 0.1) +
    ylab('Cluster ID') + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    ggtitle(paste(x,
                  ' marker genes',
                  sep = ''))
  ggsave(paste('cell.markers/marker.dotplot/cell type marker enrichment ',
               x,
               ' dotplot.png'),
         width = 10,
         height = 10)
}

# ###hoverplot
# FeaturePlot(object = mouse.snseq.combined.sct.all, 
#             features = "astrocytes.score1",
#             min.cutoff = 1,
#             max.cutoff = 2,
#             pt.size = 0.5,
#             order = TRUE,
#             cols = c("grey",
#                      "red")) %>% 
#   HoverLocator(information = FetchData(object = mouse.snseq.combined.sct.all, 
#                                        vars = 'ident'))

### dotplot
## graph
DotPlot(object = mouse.snseq.combined.sct.all, 
        features = mouse.snseq.combined.sct.all@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames(),
        col.min = 0,
        dot.min = .1) +
  ylab('Cluster ID') + 
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave('cell.markers/cell type marker enrichment per cluster dotplot.png',
       width = 10,
       height = 10)

##check data
##use marker score
cell.types.marker.per.cluster = DotPlot(object = mouse.snseq.combined.sct.all, 
                                        features = mouse.snseq.combined.sct.all@meta.data %>% 
                                          select(ends_with(".score1")) %>% 
                                          colnames())

#get counts of cell per cluster
cell.types.marker.per.cluster.count = table(mouse.snseq.combined.sct.all@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100)

##use marker gene score
cell.types.gene.marker.per.cluster = DotPlot(object = mouse.snseq.combined.sct.all, 
                                             features = marker.genes.list %>% 
                                               filter(cell.subtype.category == 'broad') %>% 
                                               pull(gene.id.tilapia))

#get counts of cell per cluster
cell.types.gene.marker.per.cluster.count = table(mouse.snseq.combined.sct.all@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.gene.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100) %>% 
  rename(gene.id.tilapia = features.plot) %>% 
  left_join(marker.genes.list)

## combine gene marker and marker score
cell.types.per.cluster.count = full_join(cell.types.gene.marker.per.cluster.count, 
                                         cell.types.marker.per.cluster.count) %>% 
  mutate(feature.id = ifelse(is.na(gene.id.tilapia),
                             features.plot %>% 
                               as.character(),
                             cell.subtype),
         weight.pct.avg = avg.exp.scaled*pct.exp/100)

#graph
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>% 
  filter(id == '5') %>% 
  ggplot(aes(x= avg.exp.scaled,
             y= pct.exp,
             label = feature.id,
             fill = weight.pct.avg)) +
  geom_label() +
  theme_bw()

#avg weight
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>%
  ggplot(aes(x= id,
             y= weight.pct.avg,
             label = features.plot)) +
  geom_text(angle = 90) +
  theme_bw()



# ### remap data to clusters
# #make dummy dataset
# mouse.snseq.combined.sct.all.label = mouse.snseq.combined.sct.all
# #create cluster ids
# new.cluster.ids <- c("0", 
#                      '1',
#                      '2',
#                      'neuron.3',
#                      '4',
#                      'astrocyte',
#                      '6',
#                      '7',
#                      'oligodendrocyte',
#                      'neuron.9',
#                      '10',
#                      'inhibitory.11',
#                      '12',
#                      '13',
#                      '14',
#                      'macrophage',
#                      '16',
#                      'endothelial.17',
#                      'endothelial.18',
#                      '19',
#                      'inhibitory.20',
#                      'neuron.21',
#                      '22',
#                      '23')
# #add names as levels
# names(new.cluster.ids) <- levels(mouse.snseq.combined.sct.all.label)
# #rename clusters
# mouse.snseq.combined.sct.all.label <- RenameIdents(mouse.snseq.combined.sct.all.label, 
#                      new.cluster.ids)
# #remap
# DimPlot(mouse.snseq.combined.sct.all.label,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) +
#   theme(legend.position = 'none')
# ggsave('comparison/Labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)
# #no label
# DimPlot(mouse.snseq.combined.sct.all.label, 
#         reduction = "umap",
#         pt.size = 1) +
#   theme(legend.position = 'none') 
# ggsave('comparison/No.labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)

# ## hover plot
# DimPlot(mouse.snseq.combined.sct.all,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) %>%
#   HoverLocator(information = FetchData(object = mouse.snseq.combined.sct.all,
#                                        vars = 'ident'))

#### save data ####
save.image('mouse.snseq.SCTransform.RData')
# load('mouse.snseq.SCTransform.RData')
##just single cell objact
# save(mouse.snseq.combined.sct.all,
# file = "mouse.snseq.combined.sct.all.RData")
# load('mouse.snseq.combined.sct.all.RData')

## save out barcodes
mouse.snseq.selected.barcodes = mouse.snseq.combined.sct.all@assays$RNA@data@Dimnames[[2]] %>% 
  as.data.frame()
#rename column
colnames(mouse.snseq.selected.barcodes) = 'Barcode.id'
#seperate barcode
mouse.snseq.selected.barcodes = mouse.snseq.selected.barcodes %>% 
  separate(Barcode.id, 
           into = c("Barcode","orig.ident"),
           sep = '_') %>% 
  mutate(orig.ident = ifelse(orig.ident == 1,
                             'dom',
                             'sub'))

#dom
dom.mouse.snseq.selected.barcodes = mouse.snseq.selected.barcodes %>% 
  filter(orig.ident == 'dom') %>% 
  select(Barcode)
#sub
sub.mouse.snseq.selected.barcodes = mouse.snseq.selected.barcodes %>% 
  filter(orig.ident == 'sub') %>% 
  select(Barcode)

##save tsv
#dom
write_tsv(dom.mouse.snseq.selected.barcodes, 
          'dom.mouse.snseq.selected.barcodes.tsv',
          col_names = FALSE)

#sub
write_tsv(dom.mouse.snseq.selected.barcodes, 
          'sub.mouse.snseq.selected.barcodes.tsv',
          col_names = FALSE)

#### QC data ####
### calculate percent mitochondria
tmp <- PercentageFeatureSet(mouse.snseq.combined.sct, pattern = "^mt-", col.name = 'percent.mt', assay = 'RNA')
## graph
# vlnplot
VlnPlot(tmp, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        group.by = 'orig.ident')
ggsave('QC_filtering/Integrated/mitochondira/Violinplot mito all.png',
       height = 5,
       width = 10)

# histogram
tmp@meta.data %>% 
  dplyr::select(percent.mt,
                orig.ident) %>% 
  ggplot(aes(x = percent.mt,
             fill = orig.ident)) +
  geom_histogram() +
  facet_grid(orig.ident~.) +
  theme_classic() +
  geom_vline(xintercept = 5)
ggsave('QC_filtering/Integrated/mitochondira/Histogram mito all.png',
       height = 10,
       width = 10)

# count
tmp@meta.data %>% 
  dplyr::select(percent.mt,
                orig.ident) %>% 
  mutate(mt.threshold = ifelse(percent.mt > 5,
                               'remove',
                               'keep')) %>% 
  dplyr::select(-c(percent.mt)) %>% 
  table() %>% 
  data.frame() %>%
  ggplot(aes(x = orig.ident,
             y = Freq,
             label = Freq,
             fill = mt.threshold)) +
  geom_label() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.5, 
                                   hjust=0.5))
ggsave('QC_filtering/Integrated/mitochondira/Count threshold mito all.png',
       height = 5,
       width = 5)

## compare percent mt to other features
plot1 <- FeatureScatter(tmp, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt",
                        group.by = 'orig.ident')
plot2 <- FeatureScatter(tmp, 
                        feature1 = "nFeature_RNA", 
                        feature2 = "percent.mt",
                        group.by = 'orig.ident')
plot1 + plot2 
ggsave('QC_filtering/Integrated/mitochondira/Compare features mito all.png',
       height = 5,
       width = 10)

### compare normalization methods
DefaultAssay(tmp) = 'SCT'

tmp_logtransform <- tmp
DefaultAssay(tmp_logtransform) = 'RNA'
tmp_logtransform <- NormalizeData(tmp_logtransform, verbose = FALSE)
tmp_logtransform <- FindVariableFeatures(tmp_logtransform, verbose = FALSE)
tmp_logtransform <- ScaleData(tmp_logtransform, verbose = FALSE)
tmp_logtransform <- RunPCA(tmp_logtransform, verbose = FALSE)
tmp_logtransform <- RunUMAP(tmp_logtransform, dims = 1:20, verbose = FALSE)

## graph
# clusters
Idents(tmp) = "seurat_clusters"
tmp_logtransform$clusterID <- Idents(tmp)
Idents(tmp_logtransform) <- "clusterID"
plot1 <- DimPlot(object = tmp, label = TRUE) + NoLegend() + ggtitle("sctransform")
#> Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
#> Please use `as_label()` or `as_name()` instead.
#> This warning is displayed once per session.
plot2 <- DimPlot(object = tmp_logtransform, label = TRUE)
plot2 <- plot2 + NoLegend() + ggtitle("Log-normalization") + scale_y_reverse()
plot1 + plot2
ggsave('QC_filtering/Integrated/normalization/Compare sctransform to log normalization cluters.png',
       height = 5,
       width = 10)

# cell types
Idents(tmp) = "parent_id.broad.prob"
tmp_logtransform$parent_id.broad.prob <- Idents(tmp)
Idents(tmp_logtransform) <- "parent_id.broad.prob"
plot1 <- DimPlot(object = tmp, label = TRUE) + NoLegend() + ggtitle("sctransform")
#> Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
#> Please use `as_label()` or `as_name()` instead.
#> This warning is displayed once per session.
plot2 <- DimPlot(object = tmp_logtransform, label = TRUE)
plot2 <- plot2  + ggtitle("Log-normalization") + scale_y_reverse()
plot1 + plot2
ggsave('QC_filtering/Integrated/normalization/Compare sctransform to log normalization cell types.png',
       height = 5,
       width = 10)

## compare variable features
library(scCustomize)

# get variable features for integrated assay
DefaultAssay(tmp) = 'integrated'
tmp2 = tmp
tmp2 <- FindVariableFeatures(tmp2, verbose = FALSE)

# graph
plot1 = VariableFeaturePlot_scCustom(tmp2,
                             assay = 'integrated',
                             num_features = 20, 
                             repel = TRUE)  +
  ggtitle('integrated')
plot2 =   VariableFeaturePlot_scCustom(tmp_logtransform,
                      assay = 'RNA',
                      num_features = 20, 
                      repel = TRUE) +
  ggtitle('RNA logtransform')
plot1 + plot2
ggsave('QC_filtering/Integrated/normalization/Compare sctransform to log normalization variance.png',
       height = 10,
       width = 15)

#### compare variable features across cell types
#set idents
Idents(object = mouse.snseq.combined.sct) <- "parent_id.exp.prob"

##subset to astrocytes
mouse.snseq.combined.sct.astrocytes = subset(mouse.snseq.combined.sct,
                                             idents = c("C25-18: Astrocytes"))

# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.astrocytes) = 'SCT'

### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.astrocytes <- FindVariableFeatures(mouse.snseq.combined.sct.astrocytes, 
                                                                assay = 'integrated',
                                                                selection.method = "vst", 
                                                                verbose = F)

# graph
VariableFeaturePlot_scCustom(mouse.snseq.combined.sct.astrocytes,
                             assay = 'integrated',
                             num_features = 20, 
                             repel = TRUE) +
  theme_classic() +
  ggtitle('astrocytes')
ggsave('QC_filtering/Integrated/normalization/astrocytes variance.png',
       height = 10,
       width = 10)

#subset to oligodendrocytes
mouse.snseq.combined.sct.oligodendrocytes = subset(mouse.snseq.combined.sct,
                                                   idents = c("C25-19: Oligodendrocytes"))
# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.oligodendrocytes) = 'SCT'

### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.oligodendrocytes <- FindVariableFeatures(mouse.snseq.combined.sct.oligodendrocytes, 
                                                            assay = 'integrated',
                                                            selection.method = "vst", 
                                                            verbose = F)

# graph
VariableFeaturePlot_scCustom(mouse.snseq.combined.sct.oligodendrocytes,
                             assay = 'integrated',
                             num_features = 20, 
                             repel = TRUE) +
  theme_classic() +
  ggtitle('oligodendrocytes')
ggsave('QC_filtering/Integrated/normalization/oligodendrocytes variance.png',
       height = 10,
       width = 10)

# subset to OPCs
mouse.snseq.combined.sct.OPCs = subset(mouse.snseq.combined.sct,
                                       idents = c("C25-20: OPC"))
# subset with SCT data
DefaultAssay(mouse.snseq.combined.sct.OPCs) = 'SCT'

### calculate variable genes
# use integrated assay for variable features
mouse.snseq.combined.sct.OPCs <- FindVariableFeatures(mouse.snseq.combined.sct.OPCs, 
                                                                  assay = 'integrated',
                                                                  selection.method = "vst", 
                                                                  verbose = F)

# graph
VariableFeaturePlot_scCustom(mouse.snseq.combined.sct.OPCs,
                             assay = 'integrated',
                             num_features = 20, 
                             repel = TRUE) +
  theme_classic() +
  ggtitle('OPCs')
ggsave('QC_filtering/Integrated/normalization/OPCs variance.png',
       height = 10,
       width = 10)


