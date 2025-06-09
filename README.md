### MousePOASnseq

### Project Description
These scripts relate to the analysis and figures for mouse POA snseq across dominant and subordinate male and female mice in manuscript snseq_mouse_poa_IMC_5.docx. 
All files are stored and R scripts were run on the Hofmann lab server: //lambstor01.ccbb.utexas.edu/HOFMANN/All_projects/Mouse_poa_snseq. 

## Data Details

-  Mouse snseq object Seurat – output of Seurat integration/analysis used for most downstream analyses and figures
   - /seurat/mouse.snseq.combined.sct.RData

-  Fastq files – raw sequence data files
   - /Data/basespace
-  Cellranger output files – single nuclei data after running 10x Cellranger and used as input for Seurat
   - /Data/count
-  Souporcell output files – genotype demultiplexing with souporcell
   - /souporcell/
-  Neuropeptide gene list – List of neuropeptide genes 
   - /gene.lists/neuropeptides.list.csv

## Script Details
All scripts can be found here: //lambstor01.ccbb.utexas.edu/HOFMANN/All_projects/Mouse_poa_snseq/Scripts

# Command line scripts - used for initial processing
-  basespace_poa_IMC.txt
   - Downloading data from basespace with command line
   - Dataset produced: fastq files
-  cellranger_poa_IMC.txt
   - Converting fastq files into snseq datasets with command line
   - Input data: /Data/basespace
   - Dataset produced: See Cellranger output files
-  souporcell_poa_IMC.txt
   - Genotype demultiplexing with souporcell
   - Input data: /Data/count
   - Dataset produced: See Souporcell output files

# R scripts – Seurat analysis and figures 
> Each R script uses collapsable commented headings to help with navigation which work best when viewed in Rstudio. 
> All figures and output from R scripts found in directory: /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/

-  seurat_snseq_mouse_IMC.R
   - Summary: QC, add genotype data, normalization, and integration of single nuclei data across sample pools
   - Input data: Cell ranger output from /Data/count/
   - Output data: seurat/mouse.snseq.combined.sct.RData
   - Figures: Supp Fig 1, Supp Fig 2 
   - Tables: 
-  HypoMap_snseq_mouse_IMC.R
   - Summary: Nuclei cell type annotation from HypoMap reference
   - Input data: seurat/mouse.snseq.combined.sct.RData
   - Output data: seurat/mouse.snseq.combined.sct.RData
-  Adds cell type annotation to data
   - Figures: Fig 2a, Supp Fig 3 
   - Tables: 
-  neurons_snseq_mouse_IMC.R
   - Summary: Clustering and analysis of neuronal nuclei (including WGCNA methods not used in manuscript)
   - Input data: seurat/mouse.snseq.combined.sct.RData
   - Output data: 
   - Figures: Fig 2b, Fig 2c, Fig 3c, Supp Fig 4, Supp Fig 5
   - Supplemental tables: Supp Table 1-2
-  DEG_snseq_mouse_IMC.R
   - Summary: Analysis of neuronal neuroendocrine receptors including DEG with limmatrend and RRHO with RedRibbon
   - Input data: seurat/mouse.snseq.combined.sct.RData
   - Figures: Fig 5
   - Supplemental tables: Supp Table 3, 10-20
-  DEG_celltypes_snseq_mouse_IMC.R
   - Summary: Analysis of broad cell classes including DEG with limmatrend and RRHO with RedRibbon (includes original RRHO2 not used in manuscript)
   - Input data: seurat/mouse.snseq.combined.sct.RData 
   - Figures: Fig 2d, Supp Fig 6, Supp Fig 7, Supp Fig 8
   - Supplemental tables: Supp Table 4-9
-  Neuropeptide_network_snseq_mouse_IMC.R
   - Summary: Neuropeptide network statistics and figures
   - Input data: seurat/mouse.snseq.combined.sct.RData and /gene.lists/neuropeptides.list.csv
   - Figures: Fig 3a, Fig 3d, Fig 4a, Supp Fig 9

## Figure Details
>All figures from R scripts found in directory: /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/

-  seurat_snseq_mouse_IMC.R
   - Supp Fig 1a left - QC_filtering/MaleDom/UMI.counts.filtered.male.dom.png
   - Supp Fig 1a right - QC_filtering/MaleDom/featurescatter.male.dom.filter.png
   - Supp Fig 1b left - QC_filtering/MaleSub/UMI.counts.filtered.male.sub.png
   - Supp Fig 1b right - QC_filtering/MaleSub/featurescatter.male.sub.filter.png
   - Supp Fig 1c left - QC_filtering/FemaleDom/UMI.counts.filtered.female.dom.png
   - Supp Fig 1c right - QC_filtering/FemaleDom/featurescatter.female.dom.filter.png
   - Supp Fig 1d left - QC_filtering/FemaleSub/UMI.counts.filtered.female.sub.png
   - Supp Fig 1d right - QC_filtering/FemaleSub/featurescatter.female.sub.filter.png
   - Supp Fig 2a - QC_filtering/Integrated/DomvsSub.nfeature.dimplot.comparison.all.png
   - Supp Fig 2b - QC_filtering/Integrated/DomvsSub.genotype.dimplot.comparison.all.png
   - Supp Fig 2c - QC_filtering/Integrated/DomvsSub.clusters.dimplot.comparison.all.png
-  HypoMap_snseq_mouse_IMC.R
   - Fig 2a - HypoMap/Predicted broad UMAP.png
   - Supp Fig 3 - HypoMap/Predicted UMAP projection broad unknown.png 
-  neurons_snseq_mouse_IMC.R
   - Fig 2b - neurons/UMAP/UMAP neurons recluster clusters.png
   - Fig 2c - neurons/UMAP/Clusters vs orig.ident neurons recluster genotype scaled poster.png
   - Fig 3c - neurons/neuropeptides/Bar chart overlap AVP neurons expressing OXT color.png
   - Supp Fig 4a - neurons/ncount_nfeature_by_cluster.png
   - Supp Fig 4b - neurons/Neuron cluster tree.png
   - Supp Fig 4c - neurons/composition_samples_clusters_by_percent.png
   - Supp Fig 4d - neurons/composition_orig.idents_clusters_by_number.png
   - Supp Fig 5 - neurons/neuron markers heatmap.png
-  DEG_snseq_mouse_IMC.R
   - Fig 5a - neurons/neuropeptides/Subset genes neurons percentage orig.ident genotype.png
   - Fig 5b - neurons/neuropeptides/RedRibbon/Esr1 RedRibbon_scaled.png
   - Fig 5c - neurons/neuropeptides/RedRibbon/Ar RedRibbon_scaled.png
   - Fig 5d - neurons/neuropeptides/RedRibbon/Pgr RedRibbon_scaled.png
   - Fig 5e - neurons/neuropeptides/RedRibbon/Nr3c1 RedRibbon_scaled.png
   - Fig 5f - neurons/neuropeptides/RedRibbon/Nr3c2 RedRibbon_scaled.png
-  DEG_celltypes_snseq_mouse_IMC.R
   - Fig 2d - neurons/limmatrend/RedRibbon.neurons.DvsS.sexes.png
   - Supp Fig 6a - neurons/neurons_RRHO2_number_gene_quad.png
   - Supp Fig 6b – neurons/ neurons_RRHO2_upup_cluster_markers.png
   - Supp Fig 7a left - neurons/neuron_rr_results_upup.png
   - Supp Fig 7a right- neurons/neuron_rr_results_upup_tree.png
   - Supp Fig 7b left - neurons/neuron_rr_results_downdown.png
   - Supp Fig 7b right- neurons/neuron_rr_results_ downdown _tree.png
   - Supp Fig 8a - RedRibbon/Endothelial RedRibbon.png
   - Supp Fig 8b - RedRibbon/OPCs RedRibbon.png
   - Supp Fig 8c - RedRibbon/Oligodendrocytes RedRibbon.png
   - Supp Fig 8d - RedRibbon/Immune RedRibbon.png
   - Supp Fig 8e - RedRibbon/Astrocytes RedRibbon.png
-  Neuropeptide_network_snseq_mouse_IMC.R
   - Fig 3a - neurons/neuropeptides/comparison/neuropeptides count scaled males paper.pdf
   - Fig 3a - neurons/neuropeptides/comparison/neuropeptides count scaled females paper.pdf
   - Fig 3d left - neurons/neuropeptides/comparison/neuropeptides per neuron percent genotype.png
   - Fig 3d right - neurons/neuropeptides/comparison/neuropeptides per neuron percent cumulative genotype.png
   - Fig 4a - neurons/neuropeptides/comparison/neuropeptides presence network male dom.pdf
   - Fig 4b - neurons/neuropeptides/comparison/neuropeptides presence network male sub.pdf
   - Fig 4c - neurons/neuropeptides/comparison/neuropeptides presence network female dom.pdf
   - Fig 4d - neurons/neuropeptides/comparison/neuropeptides presence network female sub.pdf
   - Supp Fig 9 – Table from ‘results.testcov.all’

## Table Details
> All tables from R scripts found in directory: /stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/

-  neurons_snseq_mouse_IMC.R
   - Supp Table 1 - neurons/neuron_seurat_cluster_markers.csv
   - Supp Table 2 - neurons/neuron_seurat_cluster_glm.csv
-  DEG_snseq_mouse_IMC.R
   - Supp Table 3 - neurons/limmatrend/neuron.limma.results.sct.df.csv
   - Supp Table 10 - neurons/neuron_seurat_genes_glm.csv
   - Supp Table 11 - neurons/neuropeptides/Ar/limmatrend/Ar.neuron.limma.results.sct.df.csv
   - Supp Table 12 - neurons/neuropeptides/Pgr/limmatrend/Pgr.neuron.limma.results.sct.df.csv
   - Supp Table 13 - neurons/neuropeptides/Esr1/limmatrend/Esr1.neuron.limma.results.sct.df.csv
   - Supp Table 14 - neurons/neuropeptides/Nr3c1/limmatrend/Nr3c1.neuron.limma.results.sct.df.csv
   - Supp Table 15 - neurons/neuropeptides/Nr3c2/limmatrend/Nr3c2.neuron.limma.results.sct.df.csv
   - Supp Table 16 - neurons/neuropeptides/RedRibbon/Ar quadrant genes.csv
   - Supp Table 17 - neurons/neuropeptides/RedRibbon/Pgr quadrant genes.csv
   - Supp Table 18 - neurons/neuropeptides/RedRibbon/Esr1 quadrant genes.csv
   - Supp Table 19 - neurons/neuropeptides/RedRibbon/Nr3c1 quadrant genes.csv
   - Supp Table 20 - neurons/neuropeptides/RedRibbon/Nr3c2 quadrant genes.csv
-  DEG_celltypes_snseq_mouse_IMC.R
   - Supp Table 4 - RedRibbon/neuron quadrant genes.csv
   - Supp Table 5 - RedRibbon/astrocytes quadrant genes.csv
   - Supp Table 6 - RedRibbon/oligodendrocytes quadrant genes.csv
   - Supp Table 7 - RedRibbon/OPCs quadrant genes.csv
   - Supp Table 8 - RedRibbon/Immune quadrant genes.csv
   - Supp Table 9 - RedRibbon/Endothelial quadrant genes.csv
