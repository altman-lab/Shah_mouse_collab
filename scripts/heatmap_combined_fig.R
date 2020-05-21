library(tidyverse)
library(ComplexHeatmap)
library(pvclust)
#Magma color scheme
library(viridis)
dir.create("figs/heatmap", showWarnings = FALSE)
set.seed(546)

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

mod_annot_df <- data.frame(
  module = counts$module) %>% 
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

#### 1. 4 group summary, remove mod 0 ####
counts.sub <- counts %>% 
  #Remove module 0 and splits
  filter(module != "mod_00") %>% 
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

meta2 <- meta %>% 
  mutate(group = paste(status, cell, sep="_")) %>% 
  arrange(group) %>% 
  distinct() %>% 
  column_to_rownames("group") %>% 
  as.matrix()

#### Tree ####
corr.pv <- pvclust(t(counts.sub), nboot=1000, 
                   method.hclust="average", method.dist="correlation")
#au = Approximately Unbiased p-value 
#bp = Bootstrap Probability
pdf(file = "figs/heatmap/heatmap_modules.4groups_tree.pdf", 
    height=8, width=10)
plot(corr.pv)
dev.off()

#### Col (module) annotation ####
col_annot <- mod_annot_df %>% 
  #Remove module 0 and splits
  filter(!grepl("mod_0_", module) & module != "mod_00") %>% 
  #Format
  column_to_rownames("module") %>% 
  HeatmapAnnotation(df=., col=list("Infected" = c("up"="#ca0020",
                                                  "NS"="white",
                                                  "down"="#0571b0"),
                                   "Uninfected" = c("up"="#ca0020",
                                                    "NS"="white",
                                                    "down"="#0571b0")),
                    show_legend=c(TRUE,FALSE),
                    annotation_legend_param = list(Uninfected = 
                       list(title = "Significant\nfold change")),
                    simple_anno_size = unit(0.25,"cm"))

#### heatmap ####
mod_hm <- Heatmap(t(counts.sub), name = "Mean module\nlog2 expression",
                  #Expression colors
                  col = magma(20),
                  #Sample annot
                  row_split = c(1,1,2,2),
                  row_gap = unit(5, "mm"),
                  #Module annot
                  column_names_side = "top",
                  top_annotation = col_annot,
                  cluster_columns = corr.pv$hclust,
                  column_split = 2, column_gap = unit(5, "mm"),
                  column_dend_height = unit(2, "cm"),
                  #Force square
                  heatmap_height = unit(10, "cm"),
                  heatmap_width = unit(24, "cm"))

pdf(file = "figs/heatmap/heatmap_modules.4groups.pdf", 
    height=4.5, width=12)

draw(mod_hm)
dev.off()

##### 2. Hallmark names #####

hallmark <- read_csv("results/GSEA/GSEA_modules_H.csv") 

FDR.cutoff <- .05

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
hallmark_hp <- Heatmap(t(hallmark_sub), name = "Percent genes\nin module",
                       #Expression colors
                       col = magma(20),
                       #Module annot
                       column_names_side = "top",
                       top_annotation = col_annot,
                       cluster_columns = corr.pv$hclust,
                       column_split = 2, column_gap = unit(5, "mm"),
                       column_dend_height = unit(2, "cm"),
                       #Force square
                       heatmap_height = unit(22, "cm"),
                       heatmap_width = unit(24, "cm"))

#### save
pdf(file = paste("figs/heatmap/heatmap_hallmark_FDR",
                 FDR.cutoff, ".pdf", sep=""), 
    height=9, width=12)

draw(hallmark_hp)

dev.off()

#### Group by function ####
hallmark_sub_group <- data.frame(name=colnames(hallmark_sub)) %>% 
  mutate(group = fct_collapse(factor(name),
                  "development" = "ADIPOGENESIS",
                  "DNA damage" = "DNA_REPAIR" ,
                  "immune" = c("INTERFERON_ALPHA_RESPONSE",
                               "INTERFERON_GAMMA_RESPONSE",
                               "INFLAMMATORY_RESPONSE"),
                  "metabolic" = c("FATTY_ACID_METABOLISM",
                                  "OXIDATIVE_PHOSPHORYLATION"),
                  "pathway" = c("PROTEIN_SECRETION",
                                "UNFOLDED_PROTEIN_RESPONSE"),
                  "proliferation" = c("MITOTIC_SPINDLE",
                                      "MYC_TARGETS_V1"),
                  "signaling" = c("MTORC1_SIGNALING",
                                  "TNFA_SIGNALING_VIA_NFKB",
                                  "KRAS_SIGNALING_UP",
                                  "PI3K_AKT_MTOR_SIGNALING",
                                  "IL2_STAT5_SIGNALING",
                                  "ESTROGEN_RESPONSE_LATE"))) %>% 
  arrange(group, name) %>% 
  #remove duplicate group names
  mutate(group2 = c("development","DNA damage",
                    "signaling", "","","","","",
                    "metabolic","",
                    "immune", "","",
                    "proliferation","",
                    "pathway",""))

#Reorder data
hallmark_sub2 <- as.data.frame(hallmark_sub) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(hallmark_sub_group$name)) %>% 
  column_to_rownames() %>% 
  as.matrix()

#### term annot ####
row_annot <-  rowAnnotation(group = anno_text(hallmark_sub_group$group2,
                                              location=0, rot=0,
                                              just="left"),
                  show_legend=FALSE)

#### heatmap ####
hallmark_hp2 <- Heatmap(t(hallmark_sub2), 
                        name = "Percent genes\nin module",
                       #Expression colors
                       col = magma(20),
                       #term annot
                       cluster_rows = FALSE,
                       row_split = as.numeric(hallmark_sub_group$group),
                       left_annotation = row_annot,
                       row_gap = unit(2.5, "mm"),
                       #Module annot
                       column_names_side = "top",
                       top_annotation = col_annot,
                       cluster_columns = corr.pv$hclust,
                       column_split = 2, column_gap = unit(5, "mm"),
                       column_dend_height = unit(2, "cm"),
                       #Force square
                       heatmap_height = unit(22, "cm"),
                       heatmap_width = unit(24, "cm"))

#### save
pdf(file = paste("figs/heatmap/heatmap_hallmark_reorder_FDR",
                 FDR.cutoff, ".pdf", sep=""), 
    height=9, width=12)

draw(hallmark_hp2)

dev.off()
