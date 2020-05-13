"Boxplots of expression of a list of genes or modules

Saves individual plots of gene expression by variables of interest.

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu
Copyright (C) 2020 Kim Dill-Mcfarland
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
  voom.dat = Filepath to .csv or name of object in environment 
             containing voom normalized counts.
  pval.dat = Filepath to .csv or name of object in environment 
             containing limma results output by 'extract.pval.R'
  meta.dat = Filepath to .csv containing metadata. Only required if 
             voom.dat is NOT a voom object
  vars = Character vector of variables in metadata to plot
  genes.toPlot = Character vector listing genes to plot
  outdir = Filepath to directory to save results

OPTIONAL
   name = Character string to prepend to output names. Default is NULL
   gene.key = Filepath to Ensembl gene key to name genes in plots.
              Generally 'EnsemblToHGNC_GRCh38.txt'. Default is NULL
   cores = Number of parallel cores to use. Default is 1
   fdr.toPlot = Numeric cutoff for FDR values to include on plots
   
Example
  plot.all(voom.dat='P259.2_voom.counts.csv', 
           pval.dat='P259.2.gene.pval.interact.csv', 
           meta.dat='P259_all_metadata.csv', 
           genes.ToPlot=c('gene1', 'gene2'),
           vars=c('var1','var2'),
           gene.key='EnsemblToHGNC_GRCh38.txt',
           outdir='figs/gene_P259.2/', 
           name='P259.2_expression_',
           cores=3)
"

#################

plot.all <- function(voom.dat, pval.dat, meta.dat, 
                     genes.toPlot,
                     vars=NULL,
                     outdir=NULL, name=NULL, 
                     gene.key=NULL,
                     cores=1,
                     fdr.toPlot=0.05){
########## SETUP ########## 
# Data manipulation and figures
library(tidyverse)
# Multi-panel figures for ggplot
library(cowplot)
#pval annotation
library(ggpubr)
#Set seed
set.seed(333)

########## Load data ########## 
#Voom normalized counts
if(is.character(voom.dat)){
  voom.dat.loaded <- read_csv(voom.dat) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if(class(voom.dat) == "EList"){
  voom.dat.loaded <- as.data.frame(voom.dat$E) %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.toPlot)
} else if(class(voom.dat) == "data.frame"){
  voom.dat.loaded <- voom.mods %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.toPlot)
} else {
  stop("Voom data must be CSV on disk or EList/data.frame object in environment")
}

#Pvalues
if(is.character(pval.dat)){
  pval.dat.loaded <- read_csv(pval.dat) %>% 
    dplyr::select(1, adj.P.Val, group) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if(class(pval.dat) == "data.frame"){
  pval.dat.loaded <- pval.dat %>% 
    dplyr::select(1, adj.P.Val, group)%>% 
    filter(.[[1]] %in% genes.toPlot)
} else {
  stop("P-value data must be CSV on disk or data.frame in environment")
}

#Metadata
if(class(voom.dat) == "EList"){
  meta.dat.loaded <- as.data.frame(voom.dat$targets) %>% 
    dplyr::select(libID, all_of(vars))
} else if(is.character(meta.dat)){
  meta.dat.loaded <- read_csv(meta.dat) %>% 
    dplyr::select(libID, all_of(vars))
} else if(class(meta.dat) == "data.frame"){
  meta.dat.loaded <- meta.dat %>% 
    dplyr::select(libID, all_of(vars))
} else {
  stop("Metadata must be CSV on disk or part of EList/data.frame object in environment.")
}

########## Format data ########## 
#Rename 1st column to match
colnames(pval.dat.loaded)[1] <- "gene"
colnames(voom.dat.loaded)[1] <- "gene"
  
# combine logCPM, meta, and pval data
plot.dat <- voom.dat.loaded %>% 
  pivot_longer(-1, names_to = "libID", 
               values_to = "voom.count") %>% 
  left_join(meta.dat.loaded, by="libID") %>% 
  left_join(pval.dat.loaded, by="gene")

########## Plots ########## 
#List all genes/modules
to_plot <- sort(unique(plot.dat$gene))

# Setup parallel computing
library(doParallel)
library(foreach)
registerDoParallel(cores=cores)

##########  Loop through genes ########## 
foreach(i = 1:length(to_plot)) %dopar% {
  
  #Subset data to gene/module of interest
  plot.dat.sub <- plot.dat %>% 
    filter(gene == to_plot[i]) 
  
  #Interaction plot
  #Calculate mean and sd
  plot.dat.summ <- plot.dat.sub %>% 
    dplyr::select(-group, -adj.P.Val) %>% 
    distinct() %>% 
    mutate(cell = recode_factor(cell, TKO="Tollip-/-")) %>% 
    mutate(cell = factor(cell, levels = c("WT","Tollip-/-"))) %>% 
    group_by(status, cell) %>% 
    summarise(mean=mean(voom.count),
              sd=sd(voom.count)) %>% 
    ungroup()
  
  plot <- plot.dat.sub %>% 
    dplyr::select(-group, -adj.P.Val) %>% 
    distinct() %>% 
    mutate(cell = recode_factor(cell, TKO="Tollip-/-")) %>% 
    mutate(cell = factor(cell, levels = c("WT","Tollip-/-"))) %>% 
    
   ggplot(aes(x=status:cell, y=voom.count)) +
      geom_jitter(aes(shape=cell), size=2, height=0, width=0.2,
                  stroke=1, fill="white") +
      scale_shape_manual(values = c(21,22)) +
      scale_x_discrete(labels=c("Uninfected:WT" = "Uninfected", 
                                "Infected:WT" = "Infected", 
                                "Uninfected:Tollip-/-" = "", 
                                "Infected:Tollip-/-" = "")) + 
      stat_summary(fun.data=mean_sdl, 
                   fun.args = list(mult=1), 
                   geom="errorbar", color="black", width=0.1)+
      stat_summary(fun=mean, geom="errorbar", aes(ymax=..y.., ymin=..y..),
                   color="black", width=0.25) +
      theme_classic() +
      labs(y="Normalized log2 expression",
           x="", shape="") +
  
      theme(axis.text.x=element_text(hjust=-0.1, 
                                     vjust=-4, color="black"),
            axis.text.y=element_text(color="black"),
            text = element_text(size=14)) +
      coord_cartesian(xlim=c(1,4), 
                      ylim=c(min(plot.dat.sub$voom.count),
                             max(plot.dat.sub$voom.count)*1.008), 
                      clip="off")
    
#### Add signif bars for contrasts #####     
    signif.dat <- plot.dat.sub %>% 
      dplyr::select(group, adj.P.Val) %>% 
      arrange(desc(group)) %>% 
      distinct() %>% 
      mutate(FDR = round(adj.P.Val, digits = 3))  %>% 
      mutate(label = paste("FDR=", FDR, sep="")) %>% 
      mutate( group1 = c(1,3),
              group2 = c(2,4),
              y.position = max(plot.dat.sub[,"voom.count"])*1.008) %>% 
      filter(adj.P.Val <= fdr.toPlot)
      
      if(nrow(signif.dat)>0){
        plot.annot <- plot +
          stat_pvalue_manual(signif.dat, label="label",
                             tip.length = 0)
      } else{
        plot.annot <- plot
      }
  
#### Add x axis group bars#### 
  plot.ymin <- min(ggplot_build(plot.annot)$layout$coord$limits$y)
  plot.ymax <- max(ggplot_build(plot.annot)$layout$coord$limits$y)
  bar_y <- plot.ymin-((plot.ymax-plot.ymin)/11.5)
  
  plot.annot2 <- plot.annot +
    annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y,
             color="black", size=0.5) + 
    annotate("segment", x = 3, xend = 4, y = bar_y, yend = bar_y,
             color="black", size=0.5)
  
  #### Save to disk
  #### Include gene name if desired
  dir.create(path=outdir, showWarnings = FALSE)
  if (!is.null(gene.key)){
    gene.key <- read_tsv(gene.key) %>% 
      filter(ensembl_gene_id == to_plot[i])
    
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), "_",
                      unique(gene.key$mgi_symbol),
                      ".png", sep="")
    
    ggsave(filename, plot.annot2, width=5, height=5)
  } else{
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), ".png", sep="")
    ggsave(filename, plot.annot2, width=5, height=5)
  }
}

print("All plots complete.")
}
