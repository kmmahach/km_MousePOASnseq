#### from fig S5C ####
setwd(paste0(root.dir, "/DGE_CellTypes/")) 
load("./neurons/all_neurons/limma_perm/limma_perm_results.rda")
load("./neurons/all_neurons/limma_perm/rrho_results.rda")

# pull only clusters 0-9
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
      mutate(contrast = sub(".*_", "", contrast)) %>% 
      filter(as.numeric(str_extract(contrast, "\\d+")) %in% c(0:9))
    
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

# DE in > 5 clusters/subtypes, get top genes DE in the most clusters
topDEGs <- lapply(DGE_list, \(lst) {
  
  subset_genes <- subset(lst$num_clust, n_clust > 5)
  top_genes <- subset(lst$results, gene %in% subset_genes$gene) %>% 
    left_join(lst$num_clust,
              by = "gene") 
  
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

length(unique(topDEGs$Females$results$gene))
100 * length(which(topDEGs$Females$results$t > 0))/nrow(topDEGs$Females$results) 

length(unique(topDEGs$Males$results$gene))
100 * length(which(topDEGs$Males$results$t > 0))/nrow(topDEGs$Males$results) 







