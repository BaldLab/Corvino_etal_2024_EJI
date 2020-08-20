
# GO term enrichment in a given up or down-regulated vector of genes
GO.function <- function(markers, 
                        topn = 1000,
                        org, 
                        title.var = deparse(substitute(markers)),
                        plot = TRUE,
                        ...) {
  
  # Function is designed to take only up or downregulated genes
  # Function will find GO terms enriched in the list of up or downregulated genes
  # input requires a column "gene" and a column "FDR"
  # gene should be as SYMBOL
  
  #browser()
  pkg_name <- ifelse(org == "human", "org.Hs.eg.db", "org.Mm.eg.db")
  
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(paste("Package \"clusterProfiler\" needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  gene_list <- markers %>% 
    arrange(FDR) %>% 
    head(topn) %>% 
    pull(gene)
  
  db <- if(org == "human") org.Hs.eg.db::org.Hs.eg.db else org.Mm.eg.db::org.Mm.eg.db
  
  res <- clusterProfiler::enrichGO(gene = gene_list, 
                                   OrgDb = db, 
                                   keyType = "SYMBOL", ...)
  
  df <- as_tibble(res@result) %>% 
    arrange(p.adjust) %>% 
    head(10) %>% 
    mutate(Description = as.factor(Description)) %>% 
    mutate(Description = forcats::fct_reorder(Description, 
                                              dplyr::desc(p.adjust)))
  
  
  if(plot){
    ggplot(df, 
           mapping = aes(x = .data$Description, 
                         y = -log10(.data$p.adjust))) + 
      geom_bar(aes(fill = .data$Count), 
               stat = "identity") + 
      scale_fill_gradient2("Gene Count", 
                           low = "lightgrey", 
                           mid = "#feb24c", 
                           high = "#bd0026") + 
      coord_flip() + 
      geom_hline(yintercept = -log10(0.05), 
                 linetype = "dashed") + 
      xlab("Gene Ontology") + 
      ylab(bquote("-log"[10] ~ " adjusted p-value")) +
      theme_bw() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)) +
      ggtitle(title.var)
  }else{
    
    return(res)
    
  }
}



# Heatmap function
basic.heatmap <- function(df, 
                          Col.cluster = FALSE,
                          Row.cluster = FALSE,
                          title.var = "",
                          colsepvar = NULL,
                          row.size = 0.4,
                          dist.method = "euclidean",
                          hclust.method = "ward.D"){
  
  dissimfun <- function(x) {
    dist(x, method = paste0(dist.method)) # "euclidean", "maximum", "manhattan", "canberra", "binary"
  }
  
  clusterfun <- function(x) {
    hclust(x, method = paste0(hclust.method)) # "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
  }
  
  
  heatmap.2(as.matrix(log2(df + 1)),
            col = colorRampPalette(c("blue", "white", "red"))(100), #turquoise
            scale = c("row"),
            na.rm = TRUE, 
            trace = "none",
            Rowv = Row.cluster,
            Colv = Col.cluster, 
            distfun = dissimfun,
            hclustfun = clusterfun,
            dendrogram = "both",
            margins = c(10, 7),
            cexCol = 0.8,
            cexRow = row.size,
            main = title.var,
            key = TRUE,
            keysize = 1.3, 
            na.color = "yellow",
            # breaks = break_scale, 
            sepcolor = "black", 
            colsep = colsepvar, 
            srtCol = 90)
  
}



# Dendrogram function
dendrogram.samples <- function(data.matrix, 
                               hc_metric.var = "euclidean",
                               hc_method.var = "ward.D2", 
                               type.var = "unrooted",   #"phylogram”, “cladogram”, “fan”, “unrooted”, “radial”
                               file.path = "Output/Figures/All_samples"){
 
   # Create output file path
  if(!dir.exists(file.path)){
    dir.create(file.path)}
  
  
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop(paste("Package", "factoextra", "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop(paste("Package", "cluster", "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop(paste("Package", "ape", "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  
  
  # set plotting margins 
  par(mar=c(6,4.1,4.1,2.1))
  
  
  #Sample Dendrogram
  res.hc <- eclust(t(data.matrix), 
                   stand = TRUE,
                   FUNcluster = "hclust",
                   hc_metric = hc_metric.var,
                   hc_method = hc_method.var, 
                   #nstart = 25,
                   k = 1,
                   verbose = TRUE, 
                   graph = FALSE)
  
  
  
  x <- as.dendrogram(res.hc)
  
  pdf(paste0(file.path, "_Dendrogram.pdf"))
  print(fviz_dend(x, 
                  #k = 1, 
                  #k_colors = c("red", "blue", "green"),
                  show_labels = TRUE, 
                  color_labels_by_k = FALSE, 
                  type = "rectangle",
                  rect = TRUE,
                  rect_border = c("black"),
                  rect_lty = "solid",
                  rect_lwd = 1,
                  main = "Sample Dendrogram",
                  xlab = paste0("Dist = ", hc_metric.var, " & Clust = ",  hc_method.var), 
                  cex = 0.6))
  dev.off()
  
  
  # using package ape
  pdf(paste0(file.path, "_Dendrogram_", type.var, ".pdf"))
  
  print(plot(as.phylo(res.hc), 
             type = type.var, 
             cex = 0.6,
             no.margin = TRUE))
  
  dev.off()
  
  
  
}




# msigdb 
msigdb.function <- function(input.data,
                            species.var = "human",
                            msig.database = "H", # c("H", "C2", "C3", "C5", "C7")
                            pval.cut.val = 0.05,
                            is.gene.list = TRUE,
                            qval.cut.val = 0.2, 
                            file.path = "Output/Figures/GSEA", 
                            title.var = ""){  
  
  # If is.gene.list = TRUE, input.data = vector of genes sig up or down regulated 
  # if is.gene.list = FALSE, input.data is a named vector of genes up and down regulated genes with Foldchange and ordered
  
  #H: hallmark gene sets
  #C1: positional gene sets
  #C2: curated gene sets
  #C3: motif gene sets
  #C4: computational gene sets
  #C5: GO gene sets
  #C6: oncogenic signatures
  #C7: immunologic signatures
  
  # Create output directory
  if(!dir.exists(file.path)){
    dir.create(file.path)}
  
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(paste("Package", "msigdbr", "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(paste("Package", "clusterProfiler", "needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  
  # List of available species 
  #msigdbr_show_species() 
  
  if(species.var == "human"){
    species.var <- "Homo sapiens"
  }else{
    species.var <- "Mus musculus"}
  
  
  # Formatting database
  m_t2g <- msigdbr(species = species.var, category = msig.database) %>% 
    dplyr::select(gs_name, gene_symbol)
  
  
  
  if(is.gene.list){
    
    # enrichment of ontology terms in gene list
    e.genelist <- enricher(input.data,
                           TERM2GENE = m_t2g, 
                           pvalueCutoff = pval.cut.val, 
                           qvalueCutoff = qval.cut.val)
    
    
    # Barplot
    pdf(paste0(file.path, "/MsigDB_", msig.database, "_enrichment_barplot_", title.var, ".pdf"))
    
    print(barplot(e.genelist, 
                  showCategory = 10,
                  font.size = 7) + 
            ggtitle(paste0(title.var, " genes msigDB ", msig.database, " enrichment")))
    
    dev.off()
    
    
    # Network plot
    pdf(paste0(file.path, "/MsigDB_", msig.database, "_enrichment_cnetplot_", title.var, ".pdf"))
    
    print(cnetplot(e.genelist,
                   showCategory = 5) + 
            ggtitle(paste0(title.var, " genes msigDB ", msig.database, " enrichment")))
    
    dev.off()
    
    # Save file
    write.csv(e.genelist@result, 
              paste0(file.path, "/MsigDB_", msig.database, "_enrich_", title.var, ".csv"))
    
    
    
  }else{
    
    
    
    GSEA.msigdb <- GSEA(input.data, 
                        TERM2GENE = m_t2g,
                        eps = 0,
                        pvalueCutoff = pval.cut.val)
    
    
    if(nrow(GSEA.msigdb) != 0){
      
      # Dotplot
      pdf(paste0(file.path, "/MsigDB_", msig.database, "_GSEA_dotplot_", title.var, ".pdf"))
      
      print(dotplot(GSEA.msigdb,
                    showCategory = 10) + 
              ggtitle(paste0("GSEA msigDB ", msig.database)))
      
      dev.off()
      
      # Network plot
      pdf(paste0(file.path, "/MsigDB_", msig.database, "_GSEA_cnetplot_", title.var, ".pdf"))
      
      print(cnetplot(GSEA.msigdb,
                     showCategory = 5, 
                     foldChange = GSEA.data) + 
              ggtitle(paste0("GSEA msigDB ", msig.database)))
      
      dev.off()
      
      # Save file
      write.csv(GSEA.msigdb@result, 
                paste0(file.path, "/MsigDB_", msig.database, "_GSEA_", title.var, ".csv"))
      
      
      
    }
    
    
  }
  
}

