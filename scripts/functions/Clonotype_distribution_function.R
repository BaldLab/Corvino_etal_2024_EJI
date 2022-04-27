
Get.distribution <- function(clonotype.table,
                             query.clusters = NULL,
                             source.clusters = NULL,
                             target.clusters = NULL,
                             rename.source = NULL,
                             rename.target = NULL,
                             threshold.rowval = 3){
  
  #browser()
  # Extract clonotype by cluster info
  clonotype.overview <- table(clonotype.table$CTaa, clonotype.table$seurat_clusters)
  
  # Format table
  logic.vec <- colnames(clonotype.overview) %in% query.clusters
  clonotype.overview <- clonotype.overview[ ,logic.vec]
  
  # Create output table to store data
  output.df <- clonotype.overview
  output.df <- output.df[,1:2]
  output.df[,] <- 0
  colnames(output.df) <- c("Source", "Target")
  
  # Sum values in two source clusters together
  logic.vec <- colnames(clonotype.overview) %in% source.clusters
  
  if(length(source.clusters) > 1){
    
    output.df[,"Source"] <- rowSums(clonotype.overview[ ,logic.vec])
    
  }else{
    
    output.df[,"Source"] <- clonotype.overview[ ,logic.vec]
    
  }
 
  logic.vec <- colnames(clonotype.overview) %in% target.clusters
  
  
  if(length(target.clusters) > 1){
    
    output.df[,"Target"] <- rowSums(clonotype.overview[ ,logic.vec])
    
  }else{
    
    output.df[,"Target"] <- clonotype.overview[ ,logic.vec]
    
  }
  
  
  
  
  if(!is.null(rename.source)){
    colnames(output.df)[colnames(output.df) == "Source"] <- rename.source
  }
  
  if(!is.null(rename.target)){
    colnames(output.df)[colnames(output.df) == "Target"] <- rename.target
  }
  
  # Get rowSums
  row.values <- rowSums(output.df)
  
  # remove clonotypes with < clone.min requirement in this subsetted dataset 
  print("Summary statistics are:")
  print(summary(row.values))
  print(paste0("Before filtering clonotypes = ", nrow(output.df)))
  
  keep.clonotypes <- names(row.values)[row.values >= threshold.rowval]
  
  output.df <- output.df[rownames(output.df) %in% keep.clonotypes, ]
  
  print(paste0("After filtering by threshold.val clonotypes remaining = ", nrow(output.df)))
  
  # Get rowSums vec for new clonotypes
  row.values <- rowSums(output.df)
  
  # Convert rows to percentages 
  clonotype.freq <- output.df / row.values
  clonotype.freq <- clonotype.freq*100
  
  # Order dataset
  clonotype.freq <- clonotype.freq[order(clonotype.freq[,1], decreasing = T), ]
  
  return(clonotype.freq)}
