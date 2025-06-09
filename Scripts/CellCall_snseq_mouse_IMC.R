#### Burtoni snseq seurat analysis
### Cell to cell communication
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###CellCall
#https://github.com/ShellyCoder/cellcall
## cellinker
#http://www.rna-society.org/cellinker/index.html

### set working directory
setwd("/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/")


#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")

#load libraries
# library(devtools)
## need to install vctrs and usethis
# clusterProfiler, ComplexHeatmap

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# 
# BiocManager::install("clusterProfiler")
# 
library(clusterProfiler)
# 
# test = warnings()

# install.packages("Matrix", INSTALL_opts = '--no-lock')
# install.packages("SeuratObject", INSTALL_opts = '--no-lock')
# install.packages("spatstat.utils", INSTALL_opts = '--no-lock')
# install.packages("tibble", INSTALL_opts = '--no-lock')
# install.packages("spatstat.data", INSTALL_opts = '--no-lock')
# install.packages("spatstat.sparse", INSTALL_opts = '--no-lock')
# install.packages("spatstat.geom", INSTALL_opts = '--no-lock')
# install.packages("Seurat", INSTALL_opts = '--no-lock')
# install.packages("spatstat.random", INSTALL_opts = '--no-lock')
# devtools::install_github("ShellyCoder/cellcall")
## need version 4.2 to install clusterProfiler

library(cellcall)
library(Seurat)
library(tidyverse)

#### load data ####
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

# #### get receptor, ligand, and TF genes #### 
# ### get human dataframe
# human.genes.ligand.receptor.TF.df <- read.table(system.file("extdata", "new_ligand_receptor_TFs.txt", 
#                                             package = "cellcall"), 
#                                 header = TRUE, 
#                                 quote = "",
#                                 sep = "\t", 
#                                 stringsAsFactors = FALSE)
# 
# # create list of distinct genes
# human.genes.ligand.receptor.TF.list = data.frame(Human.genes = c(human.genes.ligand.receptor.TF.df$Ligand_Symbol,
#                                                                  human.genes.ligand.receptor.TF.df$Receptor_Symbol,
#                                                                  human.genes.ligand.receptor.TF.df$TF_Symbol)) %>% 
#   distinct()
# 
# #### Map homologous genes from human ####
# library(biomaRt)
# 
# ###selecting biomart database
# ensembl = useMart('ensembl')
# #list datasets
# dataset = listDatasets(ensembl)
# ##select dataset
# #human
# ensembl.human = useDataset('hsapiens_gene_ensembl',
#                              mart=ensembl)
# 
# # get list of all human attributes in tilapia 
# listAttributes(ensembl.human) %>%
#   filter(grepl('onil',
#                name)) %>%
#   pull(name)
# 
# #create human attributes
# human.to.tilapia.attributes = c('external_gene_name',
#                                 'ensembl_gene_id',
#                                 'oniloticus_homolog_ensembl_gene',
#                                 'oniloticus_homolog_associated_gene_name',
#                                 'oniloticus_homolog_orthology_type',
#                                 'oniloticus_homolog_orthology_confidence')
# 
# #identify human to tilapia homology
# human.to.tilapia.homology = getBM(attributes = human.to.tilapia.attributes,
#                                   mart = ensembl.human,
#                                   values = human.genes.ligand.receptor.TF.list$Human.genes,
#                                   filter = 'external_gene_name',
#                                   useCache = FALSE) # useCache has to do with version of R not being up to date?
# 
# #count number of human genes?
# #1287
# human.to.tilapia.homology %>% 
#   na.omit() %>% 
#   nrow() 
# 
# #check number of tilapia genes?
# #497
# human.to.tilapia.homology %>% 
#   na.omit() %>% 
#   pull(ensembl_gene_id) %>% 
#   unique() %>% 
#   length()
# 
# # add tilapia gene names
# human.to.tilapia.homology = human.to.tilapia.homology %>% 
#   right_join(features.tsv %>%  
#                dplyr::select(-c(V3)) %>% 
#                dplyr::rename('oniloticus_homolog_ensembl_gene' = 'V1')%>% 
#                dplyr::rename('tilapia_gene_name' = 'V2')) %>% 
#   na.omit()
# 
# # keep only 1:1 orthologs
# # create list of single genes
# one2one.human <- human.to.tilapia.homology %>% 
#   group_by(external_gene_name) %>% 
#   summarise(n()) %>% 
#   filter(`n()` <= 1) %>%
#   pull(external_gene_name)
# 
# # filter to only orthologs 
# # remove NA 
# human.to.tilapia.orthologs <- human.to.tilapia.homology %>% 
#   filter(external_gene_name %in% one2one.human) %>% 
#   na.omit()
# 
# 
# # remove multiple human genes 
# one2many.human <- human.to.tilapia.orthologs %>% 
#   group_by(oniloticus_homolog_associated_gene_name) %>% 
#   summarise(n()) %>% 
#   filter(`n()` <= 1) %>%
#   pull(oniloticus_homolog_associated_gene_name)
# 
# # remove duplicate human genes
# human.to.tilapia.orthologs <- human.to.tilapia.orthologs %>% 
#   filter(oniloticus_homolog_associated_gene_name %in% one2many.human)
# 
# #### Map homologous genes ####
# library(biomaRt)
# 
# 
# 
# ###selecting biomart database
# ensembl = useMart('ensembl')
# #list datasets
# dataset = listDatasets(ensembl)
# ##select dataset
# ensembl.tilapia = useDataset('oniloticus_gene_ensembl',
#                              mart=ensembl)
# 
# # get list of all human attributes in tilapia 
# listAttributes(ensembl.tilapia) %>%
#   filter(grepl('hsapiens',
#                name)) %>%
#   pull(name)
# 
# #create tilapia attributes
# tilapia.to.human.attributes = c('ensembl_gene_id',
#                                    'hsapiens_homolog_ensembl_gene',
#                                    'hsapiens_homolog_associated_gene_name',
#                                    'hsapiens_homolog_orthology_type',
#                                    'hsapiens_homolog_orthology_confidence')
# 
# #identify human to tilapia homology
# tilapia.to.human.homology = getBM(attributes = tilapia.to.human.attributes,
#                                      mart = ensembl.tilapia,
#                                      values = features.tsv$V1,
#                                      filter = 'ensembl_gene_id',
#                                      useCache = FALSE) # useCache has to do with version of R not being up to date?
# 
# #count number of human genes?
# #24698
# tilapia.to.human.homology %>% 
#   nrow()
# 
# #check number of tilapia genes?
# #22089
# tilapia.to.human.homology %>% 
#   pull(ensembl_gene_id) %>% 
#   unique() %>% 
#   length()
# 
# # add tilapia gene names
# tilapia.to.human.homology = tilapia.to.human.homology %>% 
#   right_join(features.tsv %>% 
#                dplyr::select(-c(V3)) %>% 
#                dplyr::rename('ensembl_gene_id' = 'V1')%>% 
#                dplyr::rename('tilapia_gene_name' = 'V2'))
# 
# # keep only 1:1 orthologs
# # create list of single genes
# one2one <- tilapia.to.human.homology %>% 
#   group_by(tilapia_gene_name) %>% 
#   summarise(n()) %>% 
#   filter(`n()` <= 1) %>%
#   pull(tilapia_gene_name)
# 
# # filter to only orthologs 
# # remove NA 
# tilapia.to.human.orthologs <- tilapia.to.human.homology %>% 
#   filter(tilapia_gene_name %in% one2one) %>% 
#   na.omit()
# 
# 
# # remove multiple human genes 
# one2many <- tilapia.to.human.orthologs %>% 
#   group_by(hsapiens_homolog_associated_gene_name) %>% 
#   summarise(n()) %>% 
#   filter(`n()` <= 1) %>%
#   pull(hsapiens_homolog_associated_gene_name)
# 
# # remove duplicate human genes
# tilapia.to.human.orthologs <- tilapia.to.human.orthologs %>% 
#   filter(hsapiens_homolog_associated_gene_name %in% one2many)

#### neuron  ####
mouse.snseq.combined.sct.neurons = mouse.snseq.combined.sct

#set idents
Idents(object = mouse.snseq.combined.sct.neurons) <- "parent_id.broad"

#subset to neurons
mouse.snseq.combined.sct.neurons = subset(mouse.snseq.combined.sct.neurons,
                                          idents = c("C7-2: GABA",
                                                     "C7-1: GLU"))

# need to set to integrated for clustering
DefaultAssay(mouse.snseq.combined.sct.neurons) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
mouse.snseq.combined.sct.neurons.recluster = mouse.snseq.combined.sct.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.4)

## graph 
# idents to new clusters
Idents(object = mouse.snseq.combined.sct.neurons.recluster) <- "seurat_clusters"

# remove neuron only object
rm(mouse.snseq.combined.sct.neurons)

#### cell call dom males  ####
## subset to dom males
mouse.snseq.combined.sct.neurons.recluster.dom.males = mouse.snseq.combined.sct.neurons.recluster

#set idents
Idents(object = mouse.snseq.combined.sct.neurons.recluster.dom.males) <- "orig.ident"

mouse.snseq.combined.sct.neurons.recluster.dom.males = subset(mouse.snseq.combined.sct.neurons.recluster.dom.males,  
                                                    idents = c("Male.Dom"))

# get count data 
CellCall.mouse.males.dom <-  as.matrix(Seurat::GetAssayData(object = mouse.snseq.combined.sct.neurons.recluster.dom.males, 
                                                          slot = "counts", 
                                                          assay = "SCT"))

## set negative values to 0 
CellCall.mouse.males.dom = CellCall.mouse.males.dom %>%
  pmax(0)

# convert to dataframe
CellCall.mouse.males.dom <- as.data.frame(CellCall.mouse.males.dom)

# get grouping information
# use cluster information
myCelltype.neuron.males.dom <- as.character(mouse.snseq.combined.sct.neurons.recluster.dom.males@meta.data %>% 
                                              mutate(parent_id.broad.chr = case_when(parent_id.broad == "C7-2: GABA" ~ "GABA",
                                                                                     parent_id.broad == "C7-1: GLU" ~ "GLU",
                                                                                     TRUE ~ 'NA')) %>% 
                                              pull(parent_id.broad.chr))

# rename dataframe columns to match CellCall format
colnames(CellCall.mouse.males.dom) <- paste(1:length(myCelltype.neuron.males.dom), 
                                        myCelltype.neuron.males.dom, 
                                          sep = "_")



## filter out low expressing genes across all cells
CellCall.mouse.males.dom = CellCall.mouse.males.dom %>% 
  mutate(rowsum = rowSums(.)) %>%  
  filter(rowsum >= 1) %>% 
  dplyr::select(-c(rowsum))



# create object for cell communication object
CellCall.mouse.males.dom.obj <- CreateNichConObject(data=CellCall.mouse.males.dom,
                                                  min.feature = 3,
                                                  names.field = 2,
                                                  names.delim = "_",
                                                  source = "UMI",
                                                  scale.factor = 10^6,
                                                  Org = "Mus musculus",
                                                  project = "males.dom")

# Infer the cell-cell communication score
CellCall.mouse.males.dom.obj <- TransCommuProfile(object = CellCall.mouse.males.dom.obj,
                                                pValueCor = 0.05,
                                                CorValue = 0.1,
                                                topTargetCor=1,
                                                p.adjust = 0.05,
                                                use.type="median",
                                                probs = 0.9,
                                                method="weighted",
                                                IS_core = TRUE,
                                                Org = 'Mus musculus')


### fails to run


#Pathway activity analysis
n.males.dom <- CellCall.mouse.males.dom.obj@data$expr_l_r_log2_scale

pathway.hyper.list.males.dom <- lapply(colnames(n.males.dom), function(i){
  print(i)
  tmp <- getHyperPathway(data = n.males.dom, 
                         object = CellCall.mouse.males.dom.obj, 
                         cella_cellb = i, 
                         Org="Mus musculus")
  return(tmp)
})


#### graph dom  ####
### bubble plot
myPub.df.males.dom <- getForBubble(pathway.hyper.list.males.dom, 
                               cella_cellb=colnames(n.males.dom))
png('./cellcall/males.dom/Bubbleplot.males.dom.png')
plotBubble(myPub.df.males.dom)
dev.off()

## circle plot
# make cell color
cell_color <- data.frame(color=c('#e31a1c',
                                 '#1f78b4',
                                 '#e78ac3'), 
                         stringsAsFactors = FALSE)

rownames(cell_color) <- c("GABA",
                          "GLU")
#graph
png('./cellcall/males.dom/Circleplot.males.dom.png')
ViewInterCircos(object = CellCall.mouse.males.dom.obj,
                font = 2,  
                cellColor = cell_color,
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",
                arr.length = 0.04,
                trackhight1 = 0.05,
                slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, 
                gap.degree = 2,
                trackhight2 = 0.032,
                track.margin2 = c(0.01,0.12), 
                DIY = FALSE)
dev.off()

## pheatmap
png('./cellcall/males.dom/pHeatmap.males.dom.png')
viewPheatmap(object = CellCall.mouse.males.dom.obj, 
             slot="expr_l_r_log2_scale", 
             show_rownames = T,
             show_colnames = T,
             treeheight_row=0, 
             treeheight_col=10,
             cluster_rows = T,
             cluster_cols = F,
             fontsize = 12,
             angle_col = "45", 	
             main="score")
dev.off()
