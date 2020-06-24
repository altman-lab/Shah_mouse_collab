library(tidyverse)
library(ggrepel)
#Magma color scheme
library(viridis)
set.seed(546)

##### FIGURE 6D: HALLMARK TERMS #####
#### Data ####
hallmark <- read_csv("results/GSEA/GSEA_modules_H.csv")

FDR.cutoff <- 0.05

#### Format data ####
#List terms with at least 1 significant module
terms.to.keep <- hallmark %>% 
  #Get minimum FDR per term
  group_by(Description) %>% 
  summarize(fdr.min = min(p.adjust)) %>% 
  #Apply min cutoff
  filter(fdr.min <= FDR.cutoff) %>% 
  select(Description) %>% unlist(use.names = FALSE)

#List modules with at least 1 significant term 
mods.to.keep <- hallmark %>% 
  #Get minimum FDR per term
  group_by(group) %>% 
  summarize(fdr.min = min(p.adjust)) %>% 
  #Apply min cutoff
  filter(fdr.min <= FDR.cutoff) %>% 
  select(group) %>% unlist(use.names = FALSE)

#Format data
hallmark_pct <- hallmark %>% 
  #remove mod 0
  filter(group != "00") %>% 
  #Keep terms of interest
  filter(Description %in% terms.to.keep) %>% 
  #Calculate proportion of genes in term
  select(group, Description, size.overlap.term, 
         size.group, p.adjust) %>% 
  mutate(pct = size.overlap.term/size.group*100) %>% 
  #Format labels
  mutate(Description = gsub("HALLMARK_","",Description),
         Description = gsub("_", " ", Description) ) %>%
  #Remove leading 0 in module name
  mutate(module = sub("^0+","",group))

##### Plot #####
plot <- hallmark_pct %>% 
  
  ggplot(aes(x=Description, y=-log10(p.adjust))) +
  geom_point(aes(color=size.overlap.term, size=size.group)) +
  #Label signif points
  geom_text_repel(data=filter(hallmark_pct, 
                              p.adjust <= FDR.cutoff),
                aes(label=module), show.legend = FALSE) +
  #Cutoff line
  geom_hline(yintercept = -log10(0.05)) +
  
  #Beautify
  theme_classic() +
  labs(x="", y="-log10( FDR )", size="Genes in module",
       color="Genes overlapping\nwith Hallmark term") +
  coord_flip() +
  scale_color_viridis(option="magma")

ggsave("figs/publication/Fig6D_alternate.dotplot.pdf",
       width=12, height=5)


##### Table output #####
library(kableExtra)
hallmark %>% 
  filter(group != "00") %>% 
  filter(p.adjust <= FDR.cutoff) %>% 
  #filter(Description %in% terms.to.keep) %>% 
  select(Description, group, size.group, size.overlap.term, 
         p.adjust) %>% 
  arrange(Description, group) %>% 
  
  kable() %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>% 
  collapse_rows(1, valign="top")


hallmark %>% 
  filter(group != "00") %>% 
  filter(p.adjust <= FDR.cutoff) %>% 
  #filter(Description %in% terms.to.keep) %>% 
  select(group, size.group, Description, size.overlap.term, 
         p.adjust) %>% 
  arrange(group,Description) %>% 
  
  kable() %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>% 
  collapse_rows(1:2, valign="top")
