enrich.fxn <- function(gene.list, 
                       category, 
                       subcategory=NULL,
                       genome, 
                       basename=NULL){
  require(clusterProfiler)
  require(msigdbr)
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(tidyverse)
  
  #Convert gene list to Entrez ID
  gene.entrez <- bitr(gene.list, fromType="ENSEMBL",
                           toType=c("ENTREZID","SYMBOL"),
                           OrgDb=genome)
  
  gene.entrez.list <- gene.entrez$ENTREZID
  
  #Get database of interest
  if(grepl("Hs", genome$packageName)){
    if(is.null(subcategory)){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = category))
    } else {
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                         category = category,
                                         subcategory = subcategory))}
    
  } else if(grepl("Mm", genome$packageName)){
    if(is.null(subcategory)){
    db.species <- as.data.frame(msigdbr(species = "Mus musculus",
                                        category = category)) 
    } else{
    db.species <- as.data.frame(msigdbr(species = "Mus musculus",
                                          category = category,
                                          subcategory = subcategory)) 
    }
  } else{
    stop("Function only available for human and mouse genomes.")
  }

  #run enrichment on gene list
  enrich <- enricher(gene=gene.entrez.list, 
                     TERM2GENE=dplyr::select(db.species, gs_name, entrez_gene))
  
  if (is.null(enrich)){
    print("No enriched terms.")
    
    enrich.result.clean <- data.frame(
      Description="No enriched terms",
      database=name
    )
    
    assign(name, as.data.frame(enrich.result.clean), envir=.GlobalEnv) 
  }
  else{
    #Extract results
    enrich.result <- enrich@result %>% 
      remove_rownames() %>% 
      arrange(p.adjust, Count)
    
    #Create group names for entrez+number genes ID'ed
    ## Use to separate list of entrez IDs if > 1 exist for a given term
    pivot_names <- c()
    for (i in 1:max(enrich.result$Count)){
      pivot_names[[i]] <- paste("entrez", i, sep="")
    }
    
    enrich.result.clean <- enrich.result %>% 
      #Separate entrez ID lists
      separate(geneID, into=pivot_names, sep="/") %>% 
      pivot_longer(pivot_names, names_to = "rep", values_to = "ENTREZID") %>% 
      drop_na(ENTREZID) %>% 
      #Match entrez IDs to gene IDs
      left_join(gene.entrez, by="ENTREZID") %>% 
      
      #Combine lists into single columns, sep by /
      group_by(Description, BgRatio, pvalue, p.adjust, qvalue, Count) %>% 
      mutate(ENTREZIDs = paste(ENTREZID, collapse="/"),
             SYMBOLs = paste(SYMBOL, collapse="/"),
             ENSEMBLIDs = paste(ENSEMBL, collapse="/"),) %>% 
      ungroup() %>% 
      #Keep and order vars of interest
      dplyr::select(Description, Count, GeneRatio, BgRatio, p.adjust, 
                    ENTREZIDs, SYMBOLs, ENSEMBLIDs) %>% 
      #remove repeated rows. Occur when >1 probe maps to an entrez ID
      distinct() %>% 
      #Add ID column for database name
      mutate(category=category, subcategory=subcategory)
    
    #Save results
    dir.create("results/GSEA/", showWarnings = FALSE)
    
    if(is.null(basename) & is.null(subcategory)){ 
      output.name <- category 
      filename <- paste("results/GSEA/GSEA_", output.name, ".csv", sep="")
    } else if(is.null(basename) & !is.null(subcategory)){
      output.name <- paste(category, subcategory, sep="_")
      filename <- paste("results/GSEA/GSEA_", output.name, ".csv", sep="")
    } else if(!is.null(basename) & is.null(subcategory)){
      output.name <- paste(basename, category, sep="_") 
      filename <- paste("results/GSEA/GSEA_",
                        output.name, ".csv", sep="")
    } else{ 
      output.name <- paste(basename, category, subcategory, sep="_") 
      filename <- paste("results/GSEA/GSEA_",
                        output.name, ".csv", sep="")
    }
      assign(output.name, as.data.frame(enrich.result.clean), 
             envir=.GlobalEnv) 
      write_csv(enrich.result.clean, filename)
  }}