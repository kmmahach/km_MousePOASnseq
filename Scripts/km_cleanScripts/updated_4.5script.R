# R-4.3.1, Seurat v.4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# libraries 
load_packages(c("Seurat", "RRHO2", "DEsingle", "dendextend", "RedRibbon", "emmeans", "limma", "annotables",
                "pheatmap", "ggalluvial", "multcomp", "clustree", "patchwork", "tidyverse", "edgeR"),
              out_prefix = "4.5")

#### load data ####
load("data/integrated_seurat_onlyNeurons.rda")

# subset to only seurat cluster 0-9 (10-13 are too small/not represented in every subgroup)
MSCneurons.reclust = subset(MSCneurons.reclust, seurat_clusters %in% c(0:9))

# create bulk_group idents - can be done automatically in prep.for.DGE (leave uniqueID = NULL)
# but I'm picky about pretty names for readability (:
MSCneurons.reclust$bulk_group <- with(MSCneurons.reclust@meta.data,
                                      paste(paste(orig.ident, indiv_genotype, sep = "_"),
                                            paste("cluster", seurat_clusters, sep = "_"), sep = "."))

MSCneurons.reclust$bulk_group <- gsub(".data", "", MSCneurons.reclust$bulk_group)
MSCneurons.reclust$id <- sub("\\.[^.]*$", "", MSCneurons.reclust$bulk_group)
MSCneurons.reclust$Sex <- ifelse(grepl("female", MSCneurons.reclust$orig.ident), "Female", "Male")
MSCneurons.reclust$Status <- ifelse(grepl("dom", MSCneurons.reclust$orig.ident), "Dom", "Sub")

# not strictly necessary to subset here, but makes it easy to check for unique groupings
# e.g. if some clusters don't contain any glut neurons or vice versa, that will cause downstream problems w/ limma
Idents(MSCneurons.reclust) = MSCneurons.reclust@meta.data$parent_id.broad
sub.MSCneurons <- subset_by_ident(MSCneurons.reclust,
                                  c("C7-2: GABA", "C7-1: GLU"), 
                                  cluster = FALSE)

names(sub.MSCneurons) = c("GABA", "GLU")

MSCneurons.bulk <- prep.for.DGE(sub.MSCneurons, 
                                pseudo_bulk = TRUE,
                                group.by = c("Sex", "Status", 
                                             "indiv_genotype", 
                                             "seurat_clusters"),
                                uniqueID = "bulk_group")

d = apply(MSCneurons.bulk$bulk.matrix, 2, as.numeric)
dim(d)
# [1] 26007   232

d0 = DGEList(d)
dim(d0) # check

rownames(d0) <- rownames(MSCneurons.bulk$bulk.matrix)
d0 <- calcNormFactors(d0)

cutoff <- 200 # higher for snseq data due to sparsity 
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl) 
# [1]  6112 232

#### Choose best model for differential expression ####

MSCneurons.bulk$pb_metadata <- lapply(MSCneurons.bulk$pb_metadata, \(x) factor(x))
dge.dl$samples <- append(dge.dl$samples, MSCneurons.bulk$pb_metadata)
# metadata <- as.data.frame(MSCneurons.bulk$pb_metadata)

factor_list <- MSCneurons.bulk$pb_metadata[c("Sex", "Status", "seurat_clusters", "source_obj")] %>% 
  set_names(c("sex", "status", "cluster", "n_type"))

my_model <- iterative_model_gen(factor_list,
                                dge.dl$counts,
                                top_n = 10)

formula <- my_model$final_summary$model_formula[1]
topOrder_term <- substr(formula, 
                        tail(gregexpr(" ", formula)[[1]], 1) + 1,
                        nchar(formula))

design <- make_cellMeans_matrix(as.formula(paste("~0 + ", 
                                                 topOrder_term)), 
                                factor_list,
                                colnames(dge.dl$counts))

design_clean <- design[, colSums(abs(design)) > 0, drop = FALSE]
rownames(design_clean) <- rownames(design)

# check which columns were dropped
dropped <- setdiff(colnames(design), colnames(design_clean))
dropped
# if none, keep using 'design'

# correlation of expression between neuron subtypes due to individual?
MSCneurons.bulk$pb_metadata$id <- factor(paste(MSCneurons.bulk$pb_metadata$Sex, 
                                               MSCneurons.bulk$pb_metadata$Status,
                                               MSCneurons.bulk$pb_metadata$indiv_genotype,
                                               sep = "_"))

cor <- duplicateCorrelation(MSCneurons.bulk$bulk.matrix, design, 
                            block = MSCneurons.bulk$pb_metadata$id)
cor$consensus.correlation # check
# [1] 0.0109614
# not much - don't need to add random effect of id
# but probably should add weights for neuron proportions in grouped contrasts
neuron.proportions <- get.neuron.proportions(MSCneurons.reclust)
weights = neuron.proportions %>% 
  dplyr::select(mean_prop)

# limma-voom
v.dl = voom(dge.dl, design, plot = T)
vfit.dl = lmFit(v.dl, design)

categories <- colnames(design)
groups <- groups.for.contrasts(design, 
                               factor_list, 
                               list(c("sex", "status"), 
                                    c("sex", "status", "cluster"),
                                    c("sex", "status", "n_type")))

groups_filtered <- filter_groups_for_design(groups, design)
contrast_list <- make.contrast.list(groups_filtered, 
                                    weights = weights,
                                    design = design)

contrast.matrix <- makeContrasts(contrasts = contrast_list,
                                 levels = colnames(design))


vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)
efit.dl = eBayes(vfit.dl2)
p.dl.limma = efit.dl[["p.value"]]


#### Permutation analysis
R = 10000
set.seed(12345)

# to store p-values in
p.dl.rand = vector('list', length = R)

# to store "t" values (coefficients)
t.dl.rand = vector('list', length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  shuffled_factors <- shuffle.factor.design(factor_list)
  design.dl.rand <- make_cellMeans_matrix(as.formula(paste("~0 + ", 
                                                   topOrder_term)), 
                                  shuffled_factors,
                                  colnames(dge.dl$counts))
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contrast.matrix)
  efit.dl.rand = eBayes(vfit.dl.rand2)
  p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
  t.dl.rand[[g]] = efit.dl.rand[["t"]] # maybe overkill? use to compare to p
}

# perm pvals above/below observed values
# q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma) }  )) / R 

q.dl = as.data.frame(q.dl)
efit.dl[["p.value"]] <- q.dl
names(efit.dl[["p.value"]]) <- names(contrast_list)

limma_list <- vector('list', length = length(efit.dl$p.value))

for(i in seq_along(efit.dl$p.value)){
  results.table <- contrasts.fit(efit.dl, coef = i) %>% 
    topTable(, sort.by = "P", n = Inf) %>% 
    rownames_to_column('gene') 
  
  limma_list[[i]] <- results.table
  names(limma_list)[i] <- colnames(efit.dl$p.value)[i]
}


# save everything
save(efit.dl, limma_list, contrast_list, 
     file = paste0(root.dir, "/DGE_CellTypes/neurons/all_neurons/limma_perm/limma_perm_results.rda"))

# the big one, in case we want to look at distributions later... 
save(p.dl.rand, t.dl.rand, 
     file = paste0(net.dir, "/KM_limmaPerm/10kPerm_limmaNeurons_randPvalsTvals.rda"),
     compress = "xz")

# volcano plots
sex = c("Male", "Female")
status = c("Dom", "Sub")

  volcano.plot(limma_list,
               outdir = paste0(root.dir, "/DGE_CellTypes/neurons/all_neurons/limma_perm"))


#### RRHO on limma permutation results #### 
setwd(paste0(root.dir, "/DGE_CellTypes/")) 

# Reformat for RRHO
suffixes <- sub(".*sub_?", "", names(limma_list))
limma_list <- split(limma_list, suffixes)

limma_list <- lapply(limma_list, function(sublist) {
  newnames <- vapply(names(sublist), function(nm) {
    if (grepl("^F", nm)) {
      "ttFemale_Dom_vs_Female_Sub"
    } else if (grepl("^M", nm)) {
      "ttMale_Dom_vs_Male_Sub"
    } else {
      nm  # fallback if it doesn't match
    }
  }, character(1))
  
  names(sublist) <- newnames
  sublist
})


rrho_results <- get.RRHO(limma_list,
                         group.by = sex,
                         compare.across = status,
                         # new.max.log = max.log.scale,
                         outdir = "./neurons/all_neurons/limma_perm")

  save(rrho_results, file = "./neurons/all_neurons/limma_perm/rrho_results.rda")
  


  
  