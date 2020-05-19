library(tidyverse)
library(ComplexHeatmap)
library(pvclust)
#Magma color scheme
library(viridis)
dir.create("figs/heatmap", showWarnings = FALSE)

#### Data ####
counts <- read_csv("results/module_Shah_contrast_deepSplit3_minMod50/Shah_contrast_mod_voom_counts.csv") %>% 
  #Simplify module names
  mutate(module = gsub("module_Shah_contrast_","mod_",module)) 

meta <- read_csv("data_clean/Shah.metadata.csv") %>% 
  mutate(cell = factor(cell, levels=c("WT","TKO")),
      status = factor(status, levels=c("Uninfected","Infected"))) %>% 
  #Sort and format to matrix
  select(sampID, status, cell) %>% 
  arrange(status, cell) %>% 
  column_to_rownames("sampID")

#### Split mod 0 ####
#Split module 0 and calculate mean counts
load("data_clean/module_0_sorted.RData")
mod_0 <- list(mget(c(ls(pattern="mod_0_"))))
mod_0_df <- plyr::ldply(mod_0[[1]], rbind) %>% 
  pivot_longer(-'.id', values_to = "geneName") %>% 
  drop_na( ) %>% 
  select(-name) %>% 
  rename("module"='.id')

mod_0_means <- read_csv("results/gene_level/Shah_gene_voom_counts.csv") %>% 
  inner_join(mod_0_df, by="geneName") %>% 
  select(-geneName) %>% 
  group_by(module) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) 

##### Combine all count data #####
#Combine counts data
counts.all <- counts %>% 
  bind_rows(mod_0_means) %>% 
  #order columns like meta
  select(module, rownames(meta))

#### Column (sample) annotation ####
col_annot <- meta %>% 
  select(status, cell) %>% 
  HeatmapAnnotation(df=., col=list(
    "cell" = c("WT"="#a6cee3", "TKO"="#1f78b4"),
    "status" = c("Uninfected"="#b2df8a", 
                 "Infected"="#33a02c")))

row_annot_df <- data.frame(
  module = counts.all$module) %>% 
  mutate(Uninfected = 
           ifelse(module %in% c("mod_14","mod_07","mod_16",
                                "mod_09","mod_0_both_up",
                                "mod_13","mod_15","mod_0_UI_up"),
                  "up",
            ifelse(module %in% c("mod_05","mod_11","mod_17",
                                 "mod_10","mod_02","mod_0_both_down",
                                 "mod_06","mod_12","mod_0_UI_down"),
                         "down",
                         "NS")),
         Infected = 
           ifelse(module %in% c("mod_14","mod_07","mod_16",
                                "mod_09","mod_0_both_up",
                                "mod_01","mod_03", 
                                "mod_0_I_up", "mod_00"),
                  "up",
            ifelse(module %in% c("mod_05","mod_11","mod_17",
                                 "mod_10","mod_02", "mod_0_both_down",
                                 "mod_04","mod_08","mod_0_I_down"),
                         "down",
                         "NS"))) 

#### ALL mods ####
counts.sub <- counts.all %>% 
  #Remove split module 0
  filter(!grepl("mod_0_", module)) %>% 
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()

#### Tree
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")
#au = Approximately Unbiased p-value 
#bp = Bootstrap Probability
pdf(file = "figs/heatmap/heatmap_modules_tree.pdf", height=8, width=10)
plot(corr.pv)
dev.off()

#### Row (module) annotation
row_annot <- row_annot_df %>% 
  #Remove split module 0
  filter(!grepl("mod_0_", module)) %>% 
  #Format
  column_to_rownames("module") %>% 
  rowAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                              "NS"="white",
                                              "down"="#0571b0"),
                               "Uninfected" = c("up"="#ca0020",
                                                "NS"="white",
                                                "down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                        list(title = "Significant\nfold change")))

#### heatmap
pdf(file = "figs/heatmap/heatmap_modules.pdf", height=10, width=10)

draw(Heatmap(counts.sub, name = "Module log2\nexpression",
        #Expression colors
        col = magma(20),
        #Sample annot
        cluster_columns = FALSE,
        column_split = meta,
        column_gap = unit(5, "mm"),
        #Module annot
        right_annotation = row_annot,
        cluster_rows = corr.pv$hclust,
        row_split = 2, row_gap = unit(5, "mm"),
        row_dend_width = unit(3, "cm")))
dev.off()

#### SPLIT mod 0 ####
counts.sub <- counts.all %>% 
  #Remove module 0
  filter(module != "mod_00") %>% 
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()

#### Tree
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")
#au = Approximately Unbiased p-value 
#bp = Bootstrap Probability
pdf(file = "figs/heatmap/heatmap_modules.split0_tree.pdf", height=8,
    width=10)
plot(corr.pv)
dev.off()

#### Row (module) annotation
row_annot <- row_annot_df %>% 
  #Remove module 0
  filter(module != "mod_00") %>% 
  #Format
  column_to_rownames("module") %>% 
  rowAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                              "NS"="white",
                                              "down"="#0571b0"),
                               "Uninfected" = c("up"="#ca0020",
                                                "NS"="white",
                                                "down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                             list(title = "Significant\nfold change")))

#### heatmap
pdf(file = "figs/heatmap/heatmap_modules.split0.pdf", 
    height=10, width=10)

draw(Heatmap(counts.sub, name = "Module log2\nexpression",
             #Expression colors
             col = magma(20),
             #Sample annot
             cluster_columns = FALSE,
             column_split = meta,
             column_gap = unit(5, "mm"),
             #Module annot
             right_annotation = row_annot,
             cluster_rows = corr.pv$hclust,
             row_split = 2, row_gap = unit(5, "mm"),
             row_dend_width = unit(3, "cm")))
dev.off()

#### REMOVE mod 0 ####
counts.sub <- counts.all %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()

#### Tree
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")
#au = Approximately Unbiased p-value 
#bp = Bootstrap Probability
pdf(file = "figs/heatmap/heatmap_modules.remove0_tree.pdf", height=8,
    width=10)
plot(corr.pv)
dev.off()

#### Row (module) annotation
row_annot <- row_annot_df %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  #Format
  column_to_rownames("module") %>% 
  rowAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                              "NS"="white",
                                              "down"="#0571b0"),
                               "Uninfected" = c("up"="#ca0020",
                                                "NS"="white",
                                                "down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                          list(title = "Significant\nfold change")))

#### heatmap 
pdf(file = "figs/heatmap/heatmap_modules.remove0.pdf", 
    height=10, width=10)

draw(Heatmap(counts.sub, name = "Module log2\nexpression",
             #Expression colors
             col = magma(20),
             #Sample annot
             cluster_columns = FALSE,
             column_split = meta,
             column_gap = unit(5, "mm"),
             #Module annot
             right_annotation = row_annot,
             cluster_rows = corr.pv$hclust,
             row_split = 2, row_gap = unit(5, "mm"),
             row_dend_width = unit(3, "cm")))
dev.off()

##### Hallmark names #####

hallmark <- read_csv("results/GSEA/GSEA_modules_H.csv") 

#Calculate percent of genes in term
hallmark_pct <- hallmark %>% 
  #remove mod 0
  filter(group != "00") %>% 
  #Calculate proportion of genes in term
  select(group, Description, size.overlap.term, 
         size.group, p.adjust) %>% 
  mutate(pct = size.overlap.term/size.group*100) %>% 
  #Format labels
  mutate(group = paste("mod",group, sep="_"),
         Description = gsub("HALLMARK_","",Description))

#list terms with at least 1 module FDR < 0.5
hallmark_summ <- hallmark_pct %>% 
  group_by(Description) %>% 
  summarize(pct.max=max(pct), fdr.min = min(p.adjust)) %>% 
  arrange(-fdr.min)

terms.to.keep <- hallmark_summ %>% 
  filter(fdr.min<=0.3) %>% 
  select(Description) %>% unlist(use.names = FALSE)

hallmark_sub <- hallmark_pct %>% 
  filter(Description %in% terms.to.keep) %>% 
  #Wide format
  select(group, Description, pct) %>% 
  pivot_wider(names_from = Description, values_from = pct) %>% 
  arrange(group) %>% 
  #Fill NAs
  mutate_if(is.numeric, ~ifelse(is.na(.),0,.)) %>% 
  #To matrix
  column_to_rownames("group") %>% 
  as.matrix()
  
#### Row (module) annotation
row_annot <- row_annot_df %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  #Format
  column_to_rownames("module") %>% 
  rowAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                              "NS"="white",
                                              "down"="#0571b0"),
                               "Uninfected" = c("up"="#ca0020",
                                                "NS"="white",
                                                "down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                                                 list(title = "Significant\nfold change")))

#### heatmap 
pdf(file = "figs/heatmap/heatmap_hallmark.remove0.pdf", 
    height=10, width=15)

draw(Heatmap(hallmark_sub, name = "Percent genes\nin module",
             #Expression colors
             col = magma(20),
             #Sample annot
             cluster_columns = TRUE,
             #Module annot
             right_annotation = row_annot,
             cluster_rows = corr.pv$hclust,
             row_split = 2, row_gap = unit(5, "mm"),
             row_dend_width = unit(3, "cm"),
             column_names_gp = gpar(fontsize = 8)))
dev.off()

#### 4 group summary, remove mod 0 ####
counts.sub <- counts.all %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  pivot_longer(-module, names_to = "rowname") %>% 
  left_join(rownames_to_column(meta)) %>% 
  select(-rowname) %>% 
  #calculate mean w/in group
  group_by(module, status, cell) %>% 
  summarize(count.mean = mean(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  
  #wide format
  mutate(group = paste(status, cell, sep="_"),
         group = factor(group, levels=c("Uninfected_WT",
                                        "Uninfected_TKO",
                                        "Infected_WT",
                                        "Infected_TKO"))) %>% 
  arrange(group) %>% 
  select(module, group, count.mean) %>% 
  pivot_wider(names_from = group, values_from = count.mean) %>% 
  
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()

#### Tree
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")
#au = Approximately Unbiased p-value 
#bp = Bootstrap Probability
pdf(file = "figs/heatmap/heatmap_modules.4groups_tree.pdf", 
    height=8, width=10)
plot(corr.pv)
dev.off()

#### Row (module) annotation
row_annot <- row_annot_df %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  #Format
  column_to_rownames("module") %>% 
  rowAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                              "NS"="white",
                                              "down"="#0571b0"),
                               "Uninfected" = c("up"="#ca0020",
                                                "NS"="white",
                                                "down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                                                 list(title = "Significant\nfold change")))

#### heatmap 
pdf(file = "figs/heatmap/heatmap_modules.4groups.pdf", 
    height=10, width=6)

draw(Heatmap(counts.sub, name = "Module log2\nexpression",
             #Expression colors
             col = magma(20),
             #Sample annot
             cluster_columns = FALSE,
             column_split = c("Uninfected","Uninfected",
                              "Infected","Infected"),
             column_gap = unit(5, "mm"),
             #Module annot
             right_annotation = row_annot,
             cluster_rows = corr.pv$hclust,
             row_split = 2, row_gap = unit(5, "mm"),
             row_dend_width = unit(3, "cm")))
dev.off()

