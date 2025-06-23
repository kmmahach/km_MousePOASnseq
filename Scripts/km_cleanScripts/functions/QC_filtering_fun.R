#### functions ####

# plot min cell range of genes by number of cells w/ expression
gene.min.cells.range <- function(list_of_dgCMatrices, steps_vector) { 
  lapply(list_of_dgCMatrices, (\(x) {
      temp.m <- matrix(0, nrow = length(steps), ncol = 2)
  
        for(i in 1:length(steps)) {
    
        temp <- CreateSeuratObject(counts = x,
                               project = "temp", 
                               min.cells = steps[i],
                               min.features = 700)
    
      temp.m[i,] = c(steps[i], nrow(temp@assays[["RNA"]]@meta.data))
        print(paste0("finishing step ", steps[i]))
    
        rm(temp) 
       }

      temp.df <- data.frame(temp.m)
      colnames(temp.df) <- c("min.cells", "gene.count")
      append(x, temp.df)
     }
  ))
}

gene.min.cells.plot <- function(list_with_minCells) { 
  
  lapply(list_with_minCells, (\(x) {
    
    title <- gsub('.{5}$', '', x[["group"]])
    df <- data.frame(cbind(x[["min.cells"]], x[["gene.count"]]))
    colnames(df) <- c("min.cells", "gene.count")
  
    show(ggplot(df,aes(x = min.cells,
                       y = gene.count,
                       label = paste(min.cells,
                                 gene.count,
                                 sep = ", "))) +
                geom_line() +
                geom_label() +
                ylim(0,
                max(df$gene.count)) +
                theme_classic() +
                ggtitle(paste(title)) )
        
    directory <- list.dirs(paste0(root.dir, "QC_filtering"))
  
      if (length(grep(title, directory)) > 0) {
        plotname <- paste0(root.dir, 
                           "QC_filtering/", 
                           title, 
                           "/min.cells.", 
                           title, ".png")
    } else {  
        plotname <- paste0(root.dir,
                           "QC_filtering/min.cells.", 
                           title, ".png") 
      }

    ggsave(plotname,
           units = "in", 
           width = 10, 
           height = 10, 
           bg = "white")
    
    }
  ))
}

make.seurat.obj <- function(list_of_dgCMatrices) { 
  lapply(list_of_dgCMatrices, (\(x) {
      if ({class(x)=="dgCMatrix"}) {
         CreateSeuratObject(counts = x, 
                            project = paste(deparse(substitute(x))),
                            min.cells = 3,
                            min.features = 200) -> x 
          return(x) 
    } else {
        if(require("collapse")) {
           
          CreateSeuratObject(counts = get_elem(x, \(x){class(x)=="dgCMatrix"}, 
                                         recursive = TRUE, keep.class = TRUE), 
                             project = paste(deparse(substitute(x))),
                             min.cells = 3,
                             min.features = 200) -> x 
            return(x) 
        
      } else {
         print("trying to install collapse")
         install.packages("collapse")
            if(require(collapse)){
              print("collapse installed and loaded")
          } else {
             stop("could not install collapse")
           } 
         } 
       }
     }
  ))
}

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
