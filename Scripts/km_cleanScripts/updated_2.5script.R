# R-4.3.1, Seurat v.4.4.0

net.dir <- "/stor/work/Hofmann/All_projects/Mouse_poa_snseq"
root.dir <- "/stor/home/kmm7552/km_MousePOASnseq"
setwd(paste0(root.dir, "/Scripts/km_cleanScripts/")) 
set.seed(12345)
compression = "xz" # slower, but usually smallest compression

# functions
source("./functions/DGE_fun.R")

# packages
load_packages(c("tidyverse", "Seurat", "marginaleffects", "modglm", "emmeans", 
                "multcomp", "gridExtra", "gghalves", "patchwork"), 
              out_prefix = "2.5")

# data
load("./data/integrated_seurat_onlyNeurons.rda")

#### Plot subclass proportions ####
setwd(paste0(root.dir, "/DGE_CellTypes/neurons"))

# percent of neurons in each subclass
predicted.prob.df = MSCneurons.reclust@meta.data 

predicted.prob.df %>% 
  group_by(predicted.prob,
           orig.ident,
           indiv_genotype) %>% 
  summarise(total_cells = n()) %>% 
  ungroup() %>%
  group_by(orig.ident, 
           indiv_genotype) %>% 
  mutate(n_cells = sum(total_cells)) %>%  
  ungroup() %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident)) -> df


df %>% 
  group_by(orig.ident, 
           predicted.prob) %>% 
  summarise(group.mean = mean(total_cells/n_cells),
            total_cells = sum(total_cells)) %>% 
  ungroup() %>% 
  mutate(orig.ident = sub("\\.data", "", orig.ident)) -> df.bar


ggplot() +
  scale_y_continuous(labels = scales::percent, expand = c(0,0.001)) +
  geom_bar(aes(x = reorder(predicted.prob, -total_cells),
               y = group.mean,
               fill = orig.ident), df.bar,
           stat = "identity",
           position = "dodge") + 
  geom_point(position = position_dodge(0.9),
             aes(x = reorder(predicted.prob, -total_cells),
                 y = total_cells/n_cells,
                 shape = indiv_genotype,
                 group = orig.ident), df,
             size = 2, alpha = 0.75) +
  theme_classic() +
  scale_fill_manual(values = c("#f94449", "#408D8E", "#ff7d00", "#7e38b7")) +
  labs(x = "",
       y = "percent of neurons") +
  theme(axis.text.x = element_text(hjust = 1, 
                                   vjust = 1, 
                                   angle = 45),
        plot.margin = unit(c(1,1,1,1), "cm")) 

ggsave('all_neurons/subtypes_stats/pct_cells_bySubtype.png',
       width = 19, height = 7)

#### Prep for stats ####
# subset by cell_types
Idents(MSCneurons.reclust) <- "predicted.prob"
cell_types <- unique(MSCneurons.reclust$predicted.prob)

sub.MSCneurons <- subset_by_ident(MSCneurons.reclust,
                                  cell_types,
                                  cluster = FALSE)

# get presence/absence (1/0) for subtypess
umap = data.frame(MSCneurons.reclust@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cell_ID')

metadata = data.frame(MSCneurons.reclust@meta.data) %>%
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

# counts & percentages
sum_neur_indiv <- select_neur %>%
  group_by(orig.ident,
           indiv_genotype) %>%
  summarise(total.neuron.count = n(),
            across(all_of(names(sub.MSCneurons)),
                   list(Freq = ~ sum(.x),
                        Percent = ~ (sum(.x) / n()) * 100),
                   .names = "{fn}_{.col}"),
            .groups = "drop")

sum_neur_indiv %>%
  pivot_longer(
    cols = matches("^(Percent|Freq)_"),
    names_to = c(".value", "predicted.prob"),
    names_pattern = "(Percent|Freq)_(.*)"
  ) %>%
  dplyr::select(orig.ident,
                indiv_genotype,
                predicted.prob,
                Freq,
                total.neuron.count,
                Percent) %>% 
  mutate(Sex =
           ifelse(grepl("female", orig.ident),
                  "Female",
                  "Male")) %>%
  mutate(Status =
           ifelse(grepl("dom", orig.ident),
                  "Dom",
                  "Sub")) -> dat.for.stats

# filter to only subclasses with > 1% of neurons in each sample pool
dat.for.stats %>% 
  group_by(orig.ident, 
           predicted.prob) %>% 
  summarise(
    mean.pct = mean(Percent)
  ) %>%
  filter(mean.pct > 1) -> filtered.summary

# 21 / 50 neuron types retained
print(unique(filtered.summary$predicted.prob))

dat.for.stats %>% 
  filter(predicted.prob %in% 
           filtered.summary$predicted.prob) %>% 
  mutate(predicted.prob = gsub("-", "_", sub("^.*: ", "", predicted.prob))) -> dat.for.stats

#### GLM binomial #### 
setwd(paste0(root.dir, "/DGE_CellTypes/neurons"))

# binomial GLM on neuron proportions
neuron.sub_types.glm <- data.frame()
neuron.sub_types.pairwise <- data.frame()
neuron.sub_types.modglm <- data.frame()

for (i in unique(dat.for.stats$predicted.prob)) {
  
  # binomial GLM on sex and status w/ interaction term
  tmp = glm(cbind(Freq,
                  total.neuron.count - Freq) ~ Sex*Status,
            data = subset(dat.for.stats,
                          predicted.prob == i),
            family = binomial)
  
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  data.frame(Estimate = tmp.res$test$coefficients,
             Std.error = tmp.res$test$sigma,
             z.value = tmp.res$test$tstat,
             p.value = tmp.res$test$pvalues) %>%
    rownames_to_column('Predictor') %>%
    mutate(neuron.sub_types = i) -> tmp.res.df
  
  # FDR correction for all comparisons
  tmp.res.df %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(p.adjust.fdr =
             round(p.adjust(p.value,
                            method = 'fdr'),
                   digits = 4) ) -> neuron.sub_types.glm.fdr
  
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint) %>%
    mutate(Predictor = row.names(.)) %>%
    filter(!Predictor == "(Intercept)")
  neuron.sub_types.glm.fdr %>%
    full_join(ci_df, by = c("Predictor",
                            "Estimate")) -> neuron.sub_types.glm.fdr
  # graph results
  ci_df %>%
    ggplot(aes(x = Estimate,
               y = Predictor)) +
    geom_point() +
    geom_errorbar(aes(xmin = lwr,
                      xmax = upr),
                  width = 0.2) +
    geom_vline(xintercept = 0,
               linetype = 2) +
    labs(x = "Difference in Means",
         y = "Comparison") +
    theme_classic() +
    ggtitle(paste0(i," neurons: binomial GLM"))
  
  ggsave(paste0('all_neurons/subtypes_stats/neuron_subtypesGLM/graphs/', i, '_binomGLM_cell.count.png'),
         width = 6, height = 5)
  
  # Average marginal effects (AMEs)
  ame <- slopes(tmp) %>%
    as.data.frame() %>%
    mutate(neuron.sub_types = i,
           EffectType = "AME")
  
  # Marginal effects at representative values (MERs)
  mer <- slopes(tmp, type = "response") %>%
    as.data.frame() %>%
    mutate(neuron.sub_types = i,
           EffectType = "MER") %>%
    dplyr::select(!c(rowid,
                     predicted_lo,
                     predicted_hi)) %>%
    distinct()
  
  merplots <- plot_mer_facet(mer, group = paste(i))
  
  ggsave(paste0('all_neurons/subtypes_stats/neuron_subtypesGLM/graphs/', i, '_MER_cell.count.png'),
         plot = merplots, width = 10, height = 5)
  
  write_csv(mer,
            file = paste0('all_neurons/subtypes_stats/neuron_subtypesGLM/MER/', i, '_MERtable.csv'))
  
  dat.for.mod <- subset(dat.for.stats,
                        predicted.prob == i) %>%
    mutate(Sex = ifelse(Sex=="Female", 1, 0),
           Status = ifelse(Status=="Dom", 1, 0)) %>%
    dplyr::select(Freq, total.neuron.count, Sex, Status)
  
  tmp.for.mod <- glm(cbind(Freq,
                           total.neuron.count - Freq) ~ Sex*Status,
                     data = dat.for.mod,
                     family = binomial)
  
  # modglm interaction effects
  ints <- modglm(model = tmp.for.mod,
                 vars = c("Sex", "Status"),
                 data = dat.for.mod,
                 hyps = "means",
                 type = "dd")
  
  # Average interaction effect (AIE)
  aie_df <- data.frame(neuron.sub_types = i,
                       EffectType = "AIE",
                       term = "Sex:Status",
                       Estimate = ints$aie$aie.est,
                       Std.error = ints$aie$aie.se.delta,
                       lwr = ints$aie$aie.ll,
                       upr = ints$aie$aie.ll,
                       t.val = NA)
  
  # Interaction at hypothetical values
  inthyp_df <- data.frame(
    neuron.sub_types = i,
    EffectType = "HypotheticalInt",
    term = "Sex:Status",
    Estimate = ints$inthyp$hat,
    Std.error = ints$inthyp$se.int.est,
    lwr = ints$inthyp$hat - ints$inthyp$se.int.est,
    upr = ints$inthyp$hat + ints$inthyp$se.int.est,
    t.val = ints$inthyp$t.val)
  
  neuron.sub_types.mod <- rbind(aie_df, inthyp_df)
  
  # Tukey for all pairwise comparisons
  tmp2 = glm(cbind(Freq,
                   total.neuron.count - Freq) ~ orig.ident,
             data = subset(dat.for.stats,
                           predicted.prob == i),
             family = binomial)
  
  # use emmeans to get pairwise differences
  emm = emmeans(tmp2, ~ "orig.ident")
  emm.res = pairs(emm)
  data.frame(contrast = summary(emm.res)$contrast,
             Estimate = summary(emm.res)$estimate,
             Std.error = summary(emm.res)$SE,
             z.ratio = summary(emm.res)$z.ratio,
             adj.pval = round((summary(emm.res)$p.value),
                              digits = 4)) %>%
    mutate(neuron.sub_types = i) -> emm.res.df
  
  # combine results and save
  neuron.sub_types.glm = rbind(neuron.sub_types.glm,
                               neuron.sub_types.glm.fdr)
  neuron.sub_types.pairwise = rbind(neuron.sub_types.pairwise,
                                    emm.res.df)
  neuron.sub_types.modglm = rbind(neuron.sub_types.modglm,
                                  neuron.sub_types.mod)
  
  rm(neuron.sub_types.mod, emm.res.df, neuron.sub_types.glm.fdr)
  cat("Processed ", i, "\n")
}


write_csv(neuron.sub_types.glm,
          file = 'all_neurons/subtypes_stats/neuron_subtypesGLM/neuron_sub_typesGLM.csv')
write_csv(neuron.sub_types.pairwise,
          file = 'all_neurons/subtypes_stats/neuron_subtypesGLM/neuron_sub_types_pairwise.csv')
write_csv(neuron.sub_types.modglm,
          file = 'all_neurons/subtypes_stats/neuron_subtypesGLM/neuron_sub_types_modglm.csv')

