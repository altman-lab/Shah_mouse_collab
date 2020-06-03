library(tidyverse)
library(ComplexHeatmap)
library(pvclust)
#Magma color scheme
library(viridis)
set.seed(546)

#### Data ####
counts <- read_csv("results/module_Shah_contrast_deepSplit3_minMod50/Shah_contrast_mod_voom_counts.csv") %>% 
  #Simplify module names
  mutate(module = gsub("module_Shah_contrast_","",module)) 

meta <- read_csv("data_clean/Shah.metadata.csv") %>% 
  mutate(cell = factor(cell, levels=c("WT","TKO")),
         status = factor(status, levels=c("Uninfected","Infected"))) %>% 
  #Sort and format to matrix
  select(sampID, status, cell) %>% 
  arrange(status, cell) %>% 
  column_to_rownames("sampID")

#### Format counts #####
counts.sub <- counts %>% 
  #Remove module 0 
  filter(module != "00") %>% 
  #Add metadata
  pivot_longer(-module, names_to = "rowname") %>% 
  left_join(rownames_to_column(meta)) %>% 
  select(-rowname) %>% 
  #calculate mean w/in group
  group_by(module, status, cell) %>% 
  summarize(count.mean = mean(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  
  #Recode treatment groups
  mutate(group = paste(status,cell, sep="_")) %>% 

  arrange(group) %>% 
  select(module, group, count.mean) %>% 
  pivot_wider(names_from = group, values_from = count.mean) %>% 
  
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()

#mutate()

#### Format meta #####
meta.sub<- meta %>% 
  #Recode treatment groups
  mutate(group = paste(status,cell, sep="_")) %>% 
  arrange(group) %>% 
  distinct() %>% 
  column_to_rownames("group") %>% 
  as.matrix()

#### MODULES ####
#### Tree ####
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")

#### column (module) annotation ####
col_annot <- data.frame(module = counts$module) %>% 
  mutate(Uninfected = 
           ifelse(module %in% c("14","07","16","09","13","15"),
                  "Up",
                  ifelse(module %in% c("05","11","17","10","02","06","12"),
                         "Down",
                         NA)),
         Infected = 
           ifelse(module %in% c("14","07","16","09","01","03","00"),
                  "Up",
                  ifelse(module %in% c("05","11","17","10","02","04","08"),
                         "Down",
                         NA))) %>% 
  #Remove module 0 
  filter(module != "00") %>% 
  #Format
  column_to_rownames("module") %>% 
  HeatmapAnnotation(df=., col=list("Infected" = c("Up"="#ca0020",
                                             
                                              "Down"="#0571b0"),
                               "Uninfected" = c("Up"="#ca0020",
                                               
                                                "Down"="#0571b0")),
                show_legend=c(TRUE,FALSE),
                annotation_legend_param = list(Uninfected = 
                                      list(title = "Significant\nfold change")),
                na_col = "white")

#### heatmap ####
rowNames <- c(expression(Infected~italic(Tollip^"-/-")),
              expression(Infected~italic(Tollip^"+/+")),
              expression(Uninfected~italic(Tollip^"-/-")),
              expression(Uninfected~italic(Tollip^"+/+")))

mod_hm <- Heatmap(t(counts.sub), name = "Mean log2\nexpression",
                  #Expression colors
                  col = magma(20),
                  #Sample annot
                  row_split = c(2,2,1,1),
                  row_gap = unit(5, "mm"),
                  row_title = " ",
                  row_labels = rowNames,
                  #Module annot
                  column_names_side = "top",
                  top_annotation = col_annot,
                  cluster_columns = corr.pv$hclust,
                  column_split = 2, column_gap = unit(5, "mm"),
                  column_dend_height = unit(2, "cm"),
                  column_names_rot = 0,
                  column_title = " ",
                  column_names_centered = TRUE,
                  #Force square
                  heatmap_height = unit(10, "cm"),
                  heatmap_width = unit(24, "cm"))

#### Save ####
pdf(file = "figs/publication/heatmap_modules.4groups.pdf", 
    height=6, width=12)

draw(mod_hm)
dev.off()

##### HALLMARK #####
#### Data ####
hallmark <- read_csv("results/GSEA/GSEA_modules_H.csv") 

FDR.cutoff <- .05

#### Format data ####
#Calculate percent of genes in term
hallmark_pct <- hallmark %>% 
  #remove mod 0
  filter(group != "00") %>% 
  #Calculate proportion of genes in term
  select(group, Description, size.overlap.term, 
         size.group, p.adjust) %>% 
  mutate(pct = size.overlap.term/size.group*100) %>% 
  #Format labels
  mutate(Description = gsub("HALLMARK_","",Description),
         Description = gsub("_", " ", Description) )

#list terms with at least 1 module FDR < 0.5
hallmark_summ <- hallmark_pct %>% 
  group_by(Description) %>% 
  summarize(pct.max=max(pct), fdr.min = min(p.adjust)) %>% 
  arrange(-fdr.min)

terms.to.keep <- hallmark_summ %>% 
  filter(fdr.min<=FDR.cutoff) %>% 
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

#### heatmap ####
hallmark_hm <- Heatmap(t(hallmark_sub), name = "Percent genes\nin module",
                       #Expression colors
                       col = magma(20),
                       #Module annot
                       column_names_side = "top",
                       top_annotation = col_annot,
                       cluster_columns = corr.pv$hclust,
                       column_split = 2, column_gap = unit(5, "mm"),
                       column_dend_height = unit(2, "cm"),
                       column_names_rot = 0,
                       column_title = " ",
                       column_names_centered = TRUE,
                       #Force square
                       heatmap_height = unit(22, "cm"),
                       heatmap_width = unit(24, "cm"))

#### save
pdf(file = paste("figs/publication/heatmap_hallmark_FDR",
                 FDR.cutoff, ".pdf", sep=""), 
    height=9, width=12)

draw(hallmark_hm)

dev.off()
