#### Mouse snseq genotyping analysis
### souporcell
#https://github.com/wheaton5/souporcell

### set working directory
setwd("/stor/work/Hofmann/All_projects/Mouse_poa_snseq/seurat/")


#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")
#load libraries
library(Seurat)
library(tidyverse)

###ideas
## Assignment/status vs UMI/reads 
# Do doublets have high number of reads? 
# Do unassigned have low number of UMI?
## Select filtered cells, look at assignment and status
# Do any doublets or unassigned make the cut?
# Is there a bias for a specific assignment?
## Assignment of filtered cells for cluster?
# Is there a bias for any clusters?
## compare dom to sub
# ambient RNA amount: ambient_rna.txt

#### load data ####
### load souporcell data
### load single cell data combined
load('mouse.snseq.combined.sct.RData')

## get cell barcode list
mouse.snseq.meta.data = mouse.snseq.combined.sct@meta.data %>% 
  rownames_to_column('Cell.id') 

##dom
# male
dom.male.clusters = read_tsv('../souporcell/dommale_soupercell/clusters.tsv')
#rename barcode to Cell.id
dom.male.clusters = dom.male.clusters %>%
  rename(Cell.id = barcode) %>%
  mutate(Cell.id = paste(Cell.id,
                         '1',
                         sep = '_'),
         orig.ident = 'Male.Dom')

# female
dom.female.clusters = read_tsv('../souporcell/domfemale_soupercell/clusters.tsv')
#rename barcode to Cell.id
dom.female.clusters = dom.female.clusters %>%
  rename(Cell.id = barcode) %>%
  mutate(Cell.id = paste(Cell.id,
                         '2',
                         sep = '_'),
         orig.ident = 'female.dom')

##sub
# male
sub.male.clusters = read_tsv('../souporcell/submale_soupercell/clusters.tsv')
#rename barcode to Cell.id
sub.male.clusters = sub.male.clusters %>%
  rename(Cell.id = barcode) %>%
  mutate(Cell.id = paste(Cell.id,
                         '3',
                         sep = '_'),
         orig.ident = 'Male.Sub')

# female
sub.female.clusters = read_tsv('../souporcell/subfemale_soupercell/clusters.tsv')
#rename barcode to Cell.id
sub.female.clusters = sub.female.clusters %>%
  rename(Cell.id = barcode) %>%
  mutate(Cell.id = paste(Cell.id,
                         '4',
                         sep = '_'),
         orig.ident = 'female.sub')

### combine barcodes with rownames
mouse.snseq.meta.data = dom.male.clusters %>% 
  rbind(sub.male.clusters) %>% 
  rbind(dom.female.clusters) %>% 
  rbind(sub.female.clusters) %>% 
  full_join(mouse.snseq.meta.data) %>% 
  mutate(Keep = ifelse(is.na(nCount_RNA),
                       'Remove',
                       'Keep'))

#### Graph ####
### graph souporcell clusters
##dom
#male
mouse.snseq.meta.data %>% 
  filter(orig.ident == 'Male.Dom') %>% 
  dplyr::select(status,
                assignment,
                Keep) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = reorder(assignment,
                         -Freq),
             y = Freq,
             shape = status,
             color = Keep)) +
  geom_point() +
  theme_classic() +
  ggtitle('Dom Male') +
  xlab('Genotype') 
ggsave('souporcell/Dom male genotype.png',
       height = 5,
       width = 10)

#female
mouse.snseq.meta.data %>% 
  filter(orig.ident == 'female.dom') %>% 
  dplyr::select(status,
                assignment,
                Keep) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = reorder(assignment,
                         -Freq),
             y = Freq,
             shape = status,
             color = Keep)) +
  geom_point() +
  theme_classic() +
  ggtitle('Dom Female') +
  xlab('Genotype') 
ggsave('souporcell/Dom female genotype.png',
       height = 5,
       width = 10)
  
##sub
#male
mouse.snseq.meta.data %>% 
  filter(orig.ident == 'Male.Sub') %>% 
  dplyr::select(status,
                assignment,
                Keep) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = reorder(assignment,
                         -Freq),
             y = Freq,
             shape = status,
             color = Keep)) +
  geom_point() +
  theme_classic() +
  ggtitle('Sub Male') +
  xlab('Genotype') 
ggsave('souporcell/Sub male genotype.png',
       height = 5,
       width = 10)

##sub
#female
mouse.snseq.meta.data %>% 
  filter(orig.ident == 'female.sub') %>% 
  dplyr::select(status,
                assignment,
                Keep) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = reorder(assignment,
                         -Freq),
             y = Freq,
             shape = status,
             color = Keep)) +
  geom_point() +
  theme_classic() +
  ggtitle('Sub Female') +
  xlab('Genotype') 
ggsave('souporcell/Sub female genotype.png',
       height = 5,
       width = 10)

### total count
mouse.snseq.meta.data %>% 
  filter(status == 'singlet') %>% 
  dplyr::select(Keep,
                orig.ident) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = orig.ident,
             y = Freq,
             label = Freq,
             fill = Keep)) +
  geom_bar(stat = 'identity') +
  geom_label() +
  ggtitle('Singlet count per sample') +
  theme_classic() 
ggsave('souporcell/Singlet count.png',
       height = 10,
       width = 10)





