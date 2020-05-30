library(tidyverse)
library(ggrepel)

dat <- read_csv("results/gene_level/Shah_contrast_gene_pval.csv") %>% 
  #Reorder infection status
  mutate(group = factor(group, levels=c("uninfected","infected")))

ann_text <- data.frame(group = c("infected","infected"), 
                       lab = c("FDR = 0.05", "FDR = 0.3"),
                       adj.P.Val = c(0.05,0.3),
                       logFC = c(5.8,5.8))
plot <- dat %>% 
  #create fold change color gorups
  mutate(col.group = 
           ifelse(adj.P.Val <= 0.05 & FC.group == "up", "up1",
              ifelse(adj.P.Val <= 0.05 & FC.group == "down", "down1",
                ifelse(adj.P.Val <= 0.3 & FC.group == "up", "up2",
                    ifelse(adj.P.Val <= 0.3 & FC.group == "down", "down2",
                         "none"))))) %>% 
  arrange(-adj.P.Val) %>% 
  
  #Base dot plot
  ggplot(aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=col.group)) +
  #Recolor
  scale_color_manual(values=c("#0571b0","#92c5de",
                              "grey",
                              "#ca0020","#f4a582")) +
  #Add cutoff labels
  geom_hline(yintercept = -log10(0.05)) +
  geom_hline(yintercept = -log10(0.3)) +
  geom_text(data = ann_text, aes(label = lab, vjust=-0.2)) +
  
  #Add fold change direction labels
  geom_text(aes(-4.5, 6, label="Down in Tollip-/-")) +
  geom_text(aes(4.5, 6, label="Up in Tollip-/-")) +
  
  #Label tollip
  geom_text_repel(data = filter(dat, mgi_symbol == "Tollip"),
                  aes(label=mgi_symbol),
                  vjust = 1, hjust = -1,
                  show.legend = FALSE) +
  
  #Beautify
  theme_classic() +
  labs(x="log2 fold change", y="-log10( FDR )") +
  lims(x=c(-7.5,7.5), y=c(0,6)) +
  theme(legend.position = "none") +
  #Separate by infection
  facet_wrap(~group) 

ggsave(plot, filename="figs/publication/gene_volcano_plot.png",
       width=6, height=5)
