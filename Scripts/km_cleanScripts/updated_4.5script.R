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
                "pheatmap", "ggalluvial", "multcomp", "clustree", "patchwork", "tidyverse", "edgeR", "tidytext"),
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
# and also: if there are large diffs in # of cells per group, use proportions as weights in model
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

cutoff <- 225
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl) 
# [1] 4961  232

#### Choose best model for differential expression ####
MSCneurons.bulk$pb_metadata <- lapply(MSCneurons.bulk$pb_metadata, \(x) factor(x))
dge.dl$samples <- append(dge.dl$samples, MSCneurons.bulk$pb_metadata)
metadata <- as.data.frame(MSCneurons.bulk$pb_metadata)

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

neuron.proportions %>% 
  dplyr::select(mean_prop) -> weights

# limma-voom
v.dl <- voom(dge.dl, design, plot = T)
  save(v.dl, file = paste0(root.dir, "/DGE_CellTypes/neurons/all_neurons/limma_perm/norm_counts.rda"))
  
vfit.dl <- lmFit(v.dl, design)

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
efit.dl <- eBayes(vfit.dl2)
p.dl.limma <- efit.dl[["p.value"]]


#### Permutation analysis
R = 10000 # 10k
set.seed(12345)

# to store p-values/t-values in
p.dl.rand = vector('list', length = R)
t.dl.rand = vector('list', length = R)

# will take a while depending on # of permutations
for(g in 1:R) {
  cat("Starting on Permutation", g, "\n")
  
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
  t.dl.rand[[g]] = efit.dl.rand[["t"]]
}

# perm pvals above/below observed values
q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma) }  )) / R 

q.dl = as.data.frame(q.dl)
efit.dl[["p.value"]] <- q.dl
names(efit.dl[["p.value"]]) <- names(contrast_list)

limma_list <- vector('list', length = length(efit.dl$p.value))

for(i in seq_along(efit.dl$p.value)) {
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
p.dl.rand = lapply(p.dl.rand, as, "dgCMatrix")
t.dl.rand = lapply(t.dl.rand, as, "dgCMatrix")


save(p.dl.rand, t.dl.rand, 
     file = paste0(net.dir, "/KM_limmaPerm/10kPerm_limmaNeurons_randPvalsTvals.rda"),
     compress = "xz")

  rm(p.dl.rand, t.dl.rand)

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
  

#### Concordant/discordant DEGs across clusters ####
setwd(paste0(root.dir, "/DGE_CellTypes/")) 
  ## if starting here:
  # load("./neurons/all_neurons/limma_perm/limma_perm_results.rda") 
  
subset_by_sex <- function(limma_list,
                          logFC_threshold = 0.2,
                          min_pval = 0.05) 
{
  sublist <- lapply(limma_list, \(lst) {
    lst <- subset(lst, abs(lst$logFC) > logFC_threshold &
                    lst$P.Value < min_pval)
    }
  )
  
  Females <- subset(sublist, grepl("^[F]", names(sublist))==TRUE)  
  Males <- subset(sublist, grepl("^[M]", names(sublist))==TRUE)  
  
  newlist <- as_named_list(Females, Males)
  
  np_DGE <- lapply(newlist, \(df) {
    
    results = bind_rows(df, .id = "contrast") %>% 
      mutate(contrast = sub(".*_", "", contrast))
    
    newdf <- list( results = results,
                   num_clust = results %>%
                     distinct(gene, contrast) %>%
                     count(gene, name = "n_clust"),      
                   num_genes = results %>% 
                     distinct(contrast, gene) %>% 
                     count(contrast, name = "n_genes") )
    return(newdf)
    }
  ) 
  return(np_DGE)
}


DGE_list <- subset_by_sex(limma_list)


# DE in > 6 clusters/subtypes, get top genes DE in the most clusters
topDEGs <- lapply(DGE_list, \(lst) {
  
  subset_genes <- subset(lst$num_clust, n_clust > 6)
  top_genes <- subset(lst$results, gene %in% subset_genes$gene) %>% 
    left_join(lst$num_clust) %>% 
    group_by(contrast) %>%
    slice_max(order_by = n_clust, n = 7)
  
  results <- subset(lst$results, gene %in% top_genes$gene)
  
  lst.df <- list( results = results,
                  num_clust = results %>%
                    distinct(gene, contrast) %>%
                    count(gene, name = "n_clust"),      
                  num_genes = results %>% 
                    distinct(contrast, gene) %>% 
                    count(contrast, name = "n_genes") )
  return(lst.df)
  }
)


# make df with contrast(cluster), sex, status, gene, logFC, num_genes, and n_clust
topDEGs <- lapply(topDEGs, \(lst){
  df = lst$results %>% 
    left_join(lst$num_clust) %>% 
    left_join(lst$num_genes)
  
  return(df)
})

topDEGs.test <- topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) 
  
topDEGs %>% 
  bind_rows(.id = "Sex") %>% 
  mutate(Status = ifelse(logFC > 0, "Dom", "Sub")) %>% 
  mutate(contrast = reorder_within(contrast, -n_genes, Sex),   
         gene = reorder_within(gene, n_clust, Sex)) %>%
  ggplot(aes(x = contrast, 
             y = gene, 
             fill = logFC)) +
  geom_tile(color = "grey90") +
  labs(x = "neuronal subgroups", 
       y = "top DEGs") +
  scale_fill_gradient2(name = bquote(~italic(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(topDEGs.test, Status, logFC),
            aes(x = 1, 
                y = 1, 
                fill = Status), 
            alpha = 0) +
  scale_fill_manual(name = paste0(sprintf("\u2191 \u2191"), " expr in"),
                    values = c("Dom" = muted("blue"), "Sub" = muted("red"))) +
  guides(fill = guide_legend(override.aes = list(values = c(muted("blue"), 
                                                            muted("red")),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  facet_wrap(~ Sex, scales = "free") +
  scale_x_reordered() +
  scale_y_reordered() +
  ggtitle("DEGs across clusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, vjust = 0.65,
                                   face = "italic", size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15, face = "italic"),
        legend.text = element_text(size = 12, face = "bold"))

ggsave('neurons/all_neurons/topDEGs_byCluster.png',
       width = 12, height = 6) 
  
#### Oxytocin and Vasopressin DE by cluster #### 
setwd(paste0(root.dir, "/DGE_CellTypes/")) 
  ## if starting here:
  # load("./neurons/all_neurons/limma_perm/norm_counts.rda") 


norm_counts <- data.frame(v.dl$E) %>% 
  rownames_to_column(var = "gene") %>% 
  filter(gene %in% c("Oxt", "Avp")) %>% 
  pivot_longer(cols = !gene, 
               values_to = "norm_expr",
               names_to = "pb_colname") %>% 
  left_join(metadata, by = "pb_colname") %>% 
  mutate(group = paste(Sex, Status, sep = " "),
         id = paste(Sex, Status, indiv_genotype, sep = "_"),
         seurat_clusters = paste("cluster", seurat_clusters, sep = " ")) %>% 
  group_by(gene, group, id, seurat_clusters) %>% 
  summarise(norm_expr = mean(norm_expr))

norm_counts %>% 
  filter(gene == "Oxt") %>% 
  ggplot(aes(group, norm_expr, 
           color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 30), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15)) + 
  ggtitle("OXT expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/oxtExpr_byCluster.png"),
              width = 10, height = 7)


norm_counts %>% 
  filter(gene == "Avp") %>% 
  ggplot(aes(group, norm_expr, 
             color = group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.4, jitter.size = 2, size = .8, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7"))+
  facet_wrap(~ seurat_clusters, nrow = 2) +
  scale_y_continuous(expand = c(0, 1)) +
  labs(y = "Normalized Expression", x = "") +
  theme_bw() +
  theme(text = element_text(size = 30), 
        axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", 
                                  size = 15)) + 
  ggtitle("AVP expr in neurons")

ggsave(paste0(root.dir, "/DGE_CellTypes/neurons/neuropeptides/avpExpr_byCluster.png"),
       width = 10, height = 7)


