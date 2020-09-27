
Get.distribution <- function(clonotype.table, threshold.rowval = 3){
  
  # Extract clonotype by cluster info
  temp <- table(clonotype.table$CTaa, clonotype.table$seurat_clusters)
  
  # Format table
  logic.vec <- colnames(temp) %in% c("Stimulated_exhausted", "Exhausted_1", "Exhausted_2")
  temp <- temp[ ,logic.vec]
  
  # Sum values in two exhausted clusters together
  logic.vec <- colnames(temp) %in% c("Exhausted_1", "Exhausted_2")
  temp[,2] <- rowSums(temp[ ,logic.vec])
  colnames(temp)[2] <- "Exhausted"
  temp <- temp[,1:2]
  colnames(temp)
  head(temp)
  
  # Create variable
  clonotype.counts <- temp
  
  # Get rowSums
  row.values <- rowSums(clonotype.counts)
  
  # remove clonotypes with < clone.min requirement in this subsetted dataset 
  print("Summary statistics are:")
  print(summary(row.values))
  print(paste0("Before filtering clonotypes = ", nrow(clonotype.counts)))
  
  keep.clonotypes <- names(row.values)[row.values >= threshold.rowval]
  
  clonotype.counts <- clonotype.counts[rownames(clonotype.counts) %in% keep.clonotypes, ]
  
  print(paste0("After filtering by threshold.val clonotypes remaining = ", nrow(clonotype.counts)))
  
  # Get rowSums vec for new clonotypes
  row.values <- rowSums(clonotype.counts)
  
  # Convert rows to percentages 
  clonotype.freq <- clonotype.counts / row.values
  clonotype.freq <- clonotype.freq*100
  
  # Order dataset
  clonotype.freq <- clonotype.freq[order(clonotype.freq[,1], decreasing = T), ]
  
  return(clonotype.freq)}
