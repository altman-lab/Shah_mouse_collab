##### Setup #####
# Data manipulation and figures
library(tidyverse)
# Multi-panel figures for ggplot
library(cowplot)
#Set seed
set.seed(4389)
#Function
source("scripts/Shah_boxplot_fxn.R")

##### Data #####
load("data_clean/Shah.clean.RData")
key <- read_tsv("data_clean/ensembl.key.txt")

#list genes of interest
genes.toPlot <- data.frame(
                    mgi_symbol = c("Xbp1", "Eif2ak3", "Ern1", "Atf6",
                    "Atf4", "Eif2a", "Ddit3", "Bid", "Hspa5", "Bax")) %>% 
  left_join(key, by = "mgi_symbol")


##### PLOTS #####
plot.all(voom.dat = dat.voom, 
         pval.dat = "results/gene_level/Shah_contrast_gene_pval.csv",
         vars = c("status","cell"),
         genes.toPlot = genes.toPlot$ensembl_gene_id,
         outdir = "figs/publication/",
         gene.key = "data_clean/ensembl.key.txt",
         fdr.toPlot = 0.3)
