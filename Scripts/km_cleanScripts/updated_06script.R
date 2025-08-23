# R-4.3.3, Seurat v4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# libraries 
load_packages(c("Seurat", "sp", "DEsingle", "pheatmap", "metan", "clustree", "patchwork", 
                "tidyverse", "ggraph", "igraph", "ecodist", "ggridges", "ggrepel"),
              out_prefix = "06")


#### load data ####
load(paste0(root.dir, "/HypoMap/data/integrated_seurat_withHypoMap_predictions.rda"))

# set idents
Idents(object = int.ldfs) <- "parent_id.broad.prob"

# subset to neurons
int.ldfs = subset(int.ldfs,
                  idents = c("C7-2: GABA", "C7-1: GLU"))

# subset with SCT data
DefaultAssay(int.ldfs) = "SCT"

#### Neuropeptide candidates ####
setwd(paste0(root.dir, "/DGE_CellTypes"))

# load neuropeptide list - where did this come from? 
neuropeptides = read.csv(paste0(net.dir, '/seurat/gene.lists/neuropeptides.list.csv'))

# genes
neuropeptides %>% 
  rename(gene = Gene.name) %>% 
  mutate(gene = str_to_title(gene)) %>% 
  filter(gene %in% rownames(int.ldfs)) -> neuropeptides.genes

# 68 NPs after filtering
unlist(as.vector(neuropeptides.genes$gene)) -> neuropeptides.genes

# subset Seurat object by neuropeptide.genes
sub.MSCneurons <- subset_by_gene(int.ldfs,
                                 neuropeptides.genes,
                                 slot = "data",
                                 min_count = 0.5)

# get presence/absence (1/0) for neuropeptide.genes
umap = data.frame(int.ldfs@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(int.ldfs@meta.data) %>%
  full_join(umap, by = "Cell_ID")

sapply(
  names(sub.MSCneurons), \(sub_name) {
    all_cells = metadata$Cell_ID
    as.integer(all_cells %in% sub.MSCneurons[[sub_name]]@meta.data$Cell_ID)
  }
) %>%
  as.data.frame() %>%
  mutate(Cell_ID = metadata$Cell_ID) -> bin_df

select_neur = full_join(metadata, bin_df, by = "Cell_ID")

sum_neur <- select_neur %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

sum_neur %>% 
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "cell_count"
  ) %>%
  ggplot(aes(x = reorder(gene, -cell_count),
             y = cell_count)) +
  geom_col() +
  theme_bw()

ggsave('./neurons/neuropeptides/neurons_countExprNeuropeptides_all.png',
       height = 5, width = 27)

# percent of neurons expr NPs by group
# only look at top 50% of expr NPs (by total sum)
sum_neur %>% 
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "cell_count"
  ) %>% 
  dplyr::select(gene, cell_count) %>% 
  filter(cell_count >= quantile(cell_count, 0.50)) -> top50


sum_neur_group <- select_neur %>%
  group_by(orig.ident) %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

sum_neur_group %>% 
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "cell_count"
  ) %>% 
  filter(gene %in% top50$gene) %>% 
  ggplot(aes(x = reorder(gene, -cell_count),
             y = cell_count,
             color = orig.ident)) +
  geom_point(size = 3,
             alpha = 0.75,
             position = position_jitter(width = 0.3, 
                                        height = 1)) +
  theme_classic() +
  ggtitle('number of neurons per group expressing top50 most abundant neuropeptide genes')

ggsave('./neurons/neuropeptides/neurons_countExpr_top50Neuropeptides_by_orig.ident.png',
       height = 5, width = 20)

# get top 25 NPs 
sum_neur %>% 
  pivot_longer(
    cols = starts_with("sum_"),
    names_to = "gene",
    names_prefix = "sum_",
    values_to = "cell_count"
  ) %>% 
  dplyr::select(gene, cell_count) %>% 
  filter(cell_count >= quantile(cell_count, 0.75)) -> top25

# percent of neurons expr select genes by individual

select_neur %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(total.cells = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(sum = ~ sum(.x),
                        pct = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop") %>% 
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "gene",
    names_prefix = "pct_",
    values_to = "percentage"
  ) %>%
  filter(gene %in% top25$gene) %>% 
  ggplot(aes(x = reorder(gene, -percentage),
             y = percentage,
             color = orig.ident)) +
  geom_boxplot(width = 0,
               position = position_dodge(0.75),
               size = 1) +
  geom_point(position = position_dodge(0.75),
             aes(shape = indiv_genotype,
                 group = orig.ident),
             size = 2) +
  theme_classic() +
  scale_color_manual(values = c("#CC5500",
                                "#FF6A00",
                                "#0077CC",
                                "#0095FF")) +
  xlab('') + ylab('% of neurons') +
  ggtitle('Percent nuclei expressing NP genes')

ggsave('neurons/neuropeptides/pctNeurons_exprTop50Neuropeptides_by_indiv_genotype.png',
       height = 5, width = 15)


