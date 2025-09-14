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
  write.csv(neuron.proportions, 
            file = paste0(root.dir, "/DGE_CellTypes/neurons/cluster_stats/neuron_clusterWeights.csv"))

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

#### volcano plots ####
sex = c("Male", "Female")
status = c("Dom", "Sub")

  volcano.plot(limma_list,
               outdir = paste0(root.dir, "/DGE_CellTypes/neurons/all_neurons/limma_perm"))

femaleDEGs <- limma_list$Fdom_vs_sub_all
maleDEGs <- limma_list$Mdom_vs_sub_all

up <- "Dom"
down <- "Sub"

    dcf <- femaleDEGs %>% 
      mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
      mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
      mutate(rank = logFC * log10) %>% 
      mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                    ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))
    
    max.x <- round(max(dcf$logFC) * 2) / 2
    min.x <- round(min(dcf$logFC) * 2) / 2 
    breaks = seq(min.x - 1, max.x, 0.5)
    
    upgene_labels <- dcf %>% slice_max(order_by = rank, n = 10)
    downgene_labels <- dcf %>% slice_max(order_by = -rank, n = 10)
    
    dcf_jitters <- subset(dcf, dcf$log10 > 3 & 
                            !gene %in% c(upgene_labels$gene, downgene_labels$gene))
    
    dcf_jitters_up <- subset(dcf, dcf$gene %in% upgene_labels)
    dcf_jitters_down <- subset(dcf, dcf$gene %in% downgene_labels)
    pos_jitter <- position_jitter(width = 0.4, height = 0.1, seed = 123)

    ggplot() +
      geom_point(data = subset(dcf, !gene %in% c(dcf_jitters_up$gene,
                                                 dcf_jitters_down$gene,
                                                 dcf_jitters$gene)),
                 aes(x = logFC, y = -log10(P.Value), color = diffexpressed), 
                 alpha = 0.25, size = 3.5) +
      geom_jitter(data = dcf_jitters, 
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  height = 0.4, width = 0.1, 
                  alpha = 0.25, size = 3.5) +
      geom_jitter(data = upgene_labels,
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  position = pos_jitter,
                  alpha = 0.25, size = 3) +
      geom_jitter(data = downgene_labels,
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  position = pos_jitter,
                  alpha = 0.25, size = 3) +
      geom_text_repel(data = upgene_labels,
                      aes(x = logFC, y = -log10(P.Value), label = gene),
                      size = 5, color = "black", hjust = -1, vjust = 0.5,
                      max.overlaps = Inf) +
      geom_text_repel(data = downgene_labels,
                      aes(x = logFC, y = -log10(P.Value), label = gene),
                      size = 5, color = "black", vjust = -1.2, hjust = -0.5,
                      max.overlaps = Inf) +
      geom_hline(yintercept = 1.301, lty = 4, col = "grey", lwd = 0.8) +
      labs(x = "log2 Fold Change",
           y = bquote(~-Log[10]~italic(eFDR))) +
      annotate(geom = "text", x = 1.5, y = 0.5, 
               label = paste0("\u2191\u2191", up), 
               fontface = "bold", color = "black", size = 8) +
      annotate(geom = "text", x = -1.5, y = 0.5, 
               label = paste0("\u2191\u2191", down), 
               fontface = "bold", color = "black", size = 8) +
      scale_color_manual(values = c("UP" = "#f94449", "DOWN" = "#408D8E", "NO" = "grey")) +
      scale_x_continuous(limits = c(min.x -1, max.x), breaks = breaks) +
      ylim(0, 4) +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 1, size = 15),
            axis.text.y = element_text(hjust = 0.5, size = 20),
            axis.text = element_text(color = "#3C3C3C", size = 20),
            axis.title = element_text(size = 20),   
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            strip.background = element_blank(),
            strip.text.x = element_text(size = 20),
            text = element_text(size = 25)) +
      ggtitle("Females only")
    

    ggsave(paste0(getwd(), "/neurons/all_neurons/onlyFem_all_volcano_plot.png"),
           width = 27, height = 22, units = "cm", dpi = 600)

    
    dcm <- maleDEGs %>% 
      mutate(log10 = ifelse(P.Value==0, 4, -log10(P.Value))) %>% 
      mutate(P.Value = ifelse(P.Value==0, 0.0001, P.Value)) %>% 
      mutate(rank = logFC * log10) %>% 
      mutate(diffexpressed = ifelse(logFC > 0 & P.Value < 0.05, "UP", 
                                    ifelse(logFC < -0 & P.Value < 0.05, "DOWN", "NO")))
    
    max.x <- round(max(dcm$logFC) * 2) / 2
    min.x <- round(min(dcm$logFC) * 2) / 2 
    breaks = seq(min.x, max.x, 0.5)
    
    upgene_labels <- dcm %>% slice_max(order_by = rank, n = 10)
    downgene_labels <- dcm %>% slice_max(order_by = -rank, n = 10)
    
    dcm_jitters <- subset(dcm, dcm$log10 > 3 & 
                            !gene %in% c(upgene_labels$gene, downgene_labels$gene))
    
    dcm_jitters_up <- subset(dcm, dcm$gene %in% upgene_labels)
    dcm_jitters_down <- subset(dcm, dcm$gene %in% downgene_labels)
    pos_jitter <- position_jitter(width = 0.4, height = 0.1, seed = 123)
    
    ggplot() +
      geom_point(data = subset(dcm, !gene %in% c(dcm_jitters_up$gene,
                                                 dcm_jitters_down$gene,
                                                 dcm_jitters$gene)),
                 aes(x = logFC, y = log10, color = diffexpressed), 
                 alpha = 0.25, size = 3.5) +
      geom_jitter(data = dcm_jitters, 
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  height = 0.4, width = 0.1, 
                  alpha = 0.25, size = 3.5) +
      geom_jitter(data = upgene_labels,
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  position = pos_jitter,
                  alpha = 0.25, size = 3) +
      geom_jitter(data = downgene_labels,
                  aes(x = logFC, y = log10, color = diffexpressed), 
                  position = pos_jitter,
                  alpha = 0.25, size = 3) +
      geom_text_repel(data = upgene_labels,
                      aes(x = logFC, y = log10, label = gene),
                      size = 5, color = "black", hjust = 0.2, vjust = 0.7,
                      max.overlaps = Inf) +
      geom_text_repel(data = downgene_labels,
                      aes(x = logFC, y = log10, label = gene),
                      size = 5, color = "black", vjust = 0.5, hjust = -0.5,
                      max.overlaps = Inf) +
      geom_hline(yintercept = 1.301, lty = 4, col = "grey", lwd = 0.8) +
      labs(x = "log2 Fold Change",
           y = bquote(~-Log[10]~italic(eFDR))) +
      annotate(geom = "text", x = 1.5, y = 0.5, 
               label = paste0("\u2191\u2191", up), 
               fontface = "bold", color = "black", size = 8) +
      annotate(geom = "text", x = -1.5, y = 0.5, 
               label = paste0("\u2191\u2191", down), 
               fontface = "bold", color = "black", size = 8) +
      scale_color_manual(values = c("UP" =  "#ff7d00",  "DOWN" = "#7e38b7", "NO" = "grey")) +
      scale_x_continuous(limits = c(min.x, max.x), breaks = breaks) +
      ylim(0, 4) +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 1, size = 15),
            axis.text.y = element_text(hjust = 0.5, size = 20),
            axis.text = element_text(color = "#3C3C3C", size = 20),
            axis.title = element_text(size = 20),   
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            strip.background = element_blank(),
            strip.text.x = element_text(size = 20),
            text = element_text(size = 25)) +
      ggtitle("Males only")
    
    
    ggsave(paste0(getwd(), "/neurons/all_neurons/onlyMale_all_volcano_plot.png"),
           width = 27, height = 22, units = "cm", dpi = 600)
    
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
  # revert format
  load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")
  
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

# overlap with cluster marker genes 
neuron_clusterMarkers <- read_csv('neurons/cluster_stats/neuron_clusterMarkers_weighted.csv')

sig_maleDEGs <- maleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

sig_femaleDEGs <- femaleDEGs %>% 
  filter(abs(logFC) > 0.2, P.Value < 0.05)

conc_genes = intersect(sig_femaleDEGs$gene,
                       sig_maleDEGs$gene)

result <- limma_list %>%
  .[str_detect(names(.), "_cluster\\d+$")] %>%
  imap_dfr(~ {
    .x %>%
      filter(gene %in% conc_genes) %>%
      select(gene, logFC, P.Value) %>%
      mutate(
        sex = str_to_lower(str_sub(.y, 1, 1)),  
        cluster = str_extract(.y, "(?<=_cluster)\\d+")
      )
  }) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(logFC, P.Value),
    names_glue = "{.value}.{sex}"
  ) %>%
  filter(!is.na(logFC.f) & !is.na(logFC.m)) %>%
  group_by(cluster) %>%
  filter(all(!is.na(logFC.f)) & all(!is.na(logFC.m))) %>%
  ungroup() %>% 
  mutate(
    logFC.f   = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, logFC.f, NA),
    logFC.m   = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, logFC.m, NA),
    P.Value.f = ifelse(P.Value.f < 0.05 & abs(logFC.f) > 0.2, P.Value.f, NA),
    P.Value.m = ifelse(P.Value.m < 0.05 & abs(logFC.m) > 0.2, P.Value.m, NA),
    avg.logFC = (abs(logFC.f) + abs(logFC.m)) / 2 
  ) %>% 
  filter(!(is.na(logFC.f) | is.na(logFC.m)))

neuron_clusterMarkers_conc <- neuron_clusterMarkers %>% 
  filter(gene %in% result$gene) 

result.test = left_join(result, neuron_clusterMarkers_conc[, c("gene", "specificity")]) 

result.test %>% 
  # slice_max(avg.logFC, n = 6) %>% 
  dplyr::select(gene, cluster, logFC.f, logFC.m, avg.logFC, specificity) %>%
  complete(gene, cluster,
           fill = list(avg.logFC = NA)) %>%
  group_by(cluster) %>% 
  mutate(avg.logFC = ifelse(logFC.f * logFC.m < 0, -1 * avg.logFC, avg.logFC),
         direction = ifelse(avg.logFC > 0, "concordant", "discordant")) -> result.test

result.test %>% 
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  scale_fill_gradient2(low = "#408D8E",
                       high = "#7e38b7",
                       na.value = "gray95",
                       name = bquote(avg~abs(log[2]~FC))) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = distinct(drop_na(result.test), avg.logFC, direction),
            aes(x = 1,
                y = 1,
                fill = direction),
            alpha = 0) +
  scale_fill_manual(name = "effect of status across sex",
                    values = c("concordant" = "#7e38b7", "discordant" = "#408D8E")) +
  guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
                                                            "#408D8E"),
                                                 alpha = c(1,1),
                                                 shape = c(1,1)))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()

ggsave(file = "./neurons/cluster_marker.DEGs.png",
       width = 9.25, height = 16)


# color by enrichment quadrant from RRHO to indicate direction
result2 <- rrho_results %>%
  # keep only elements that look like clusters
  .[str_detect(names(.), "^cluster\\d+$")] %>%
  imap_dfr(~ {
    cluster_name <- .y
    df <- .x$df
    
    # iterate over quadrants
    map_dfr(names(.x$RedRibbon.quads), function(q) {
      pos <- .x$RedRibbon.quads[[q]]$positions
      
      df[pos, ] %>%
        transmute(
          gene,
          cluster = str_extract(cluster_name, "\\d+"),
          quadrant = q,
          avg.logFC = (abs(a) + abs(b)) / 2
        ) 
    })
  })
  
neuron_clusterMarkers_hiSpec <- neuron_clusterMarkers %>% 
  filter(specificity > 0.1)

result2 <- result2 %>% 
  filter(gene %in% neuron_clusterMarkers_hiSpec$gene &
           avg.logFC > 1) %>% 
  complete(gene, cluster,
           fill = list(avg.logFC = NA))


result2 %>% 
  filter(is.na(quadrant)) %>% 
  ggplot(aes(x = factor(cluster),
             y = gene,
             fill = avg.logFC)) +
  geom_tile(color = "white",
            lwd = 0.8) +
  scale_fill_viridis(na.value = "gray95") +
  geom_tile(data = subset(result2, quadrant == "upup"),
            aes(fill = avg.logFC),
    color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                       high = "#7e38b7",
                       na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'upup' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "downdown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                       high = "#408D8E",
                       na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'downdown' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "updown"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                       high = "#f94449",
                       na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'updown' quadrant") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(result2, quadrant == "downup"),
            aes(fill = avg.logFC),
            color = "white",
            lwd = 0.8) +
  scale_fill_gradient(low = "white",
                       high = "#ff7d00",
                       na.value = "gray95",
                      limits = c(0, 2)) +
  labs(fill = "avg.logFC in 'downup' quadrant") +
  labs(x = "Neuron clusters", y = "cluster marker genes") +
  # ggnewscale::new_scale_fill() +
  # geom_tile(data = distinct(drop_na(result2), avg.logFC, quadrant),
  #           aes(x = 1,
  #               y = 1,
  #               fill = quadrant),
  #           alpha = 0) +
  # scale_fill_manual(name = "effect of status across sex",
  #                   values = c("upup" = "#7e38b7",
  #                              "downdown" = "#408D8E",
  #                              "updown" = "#f94449",
  #                              "downup" = "#ff7d00")) +
  # guides(fill = guide_legend(override.aes = list(values = c("#7e38b7",
  #                                                           "#408D8E",
  #                                                           "#f94449",
  #                                                           "#ff7d00")),
  #                                                alpha = c(1,1),
  #                                                shape = c(1,1))) +
  theme_classic() +
  ggtitle("Differential expression of cluster markers") +
  theme(strip.text.y = element_blank(),
        panel.spacing = unit(1, "mm"),
        plot.title = element_text(hjust = 0, size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.65, size = 13),
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  coord_fixed()
  
ggsave(file = "./neurons/cluster_marker.quads.png",
       width = 9.25, height = 16)


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


