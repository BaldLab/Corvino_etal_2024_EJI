get.clonotypes <- function(input.seurat, identity.column, ident.id){
  
  # identitiy column = cluster or module 
  # ident.id = a specific module or cluster ID
  
  # This function identifies all clonotypes in a specified cluster/module subset and then pulls out all occurances of these clonotypes
  
  # set idents
  Idents(input.seurat) <- input.seurat@meta.data[ , paste0(identity.column)]
  
  # Subset by ident.id
  Ident.specific.seurat <- subset(input.seurat, idents = paste0(ident.id))
  
  # get list of clonotypes in selected cluster/module
  Ident.specific.clones <- sort(table(Ident.specific.seurat@meta.data$CTaa), decreasing = TRUE)
  
  # Find all cells that have these clonotypes
  keep.logic <- input.seurat@meta.data$CTaa %in% names(Ident.specific.clones)
  cells.with.ident.sp.clones <- rownames(input.seurat@meta.data)[keep.logic]
  
  # Subset seurat by these cell IDs
  seurat.ident.specific.clones <- subset(input.seurat, cells = cells.with.ident.sp.clones)
  
  
  return(seurat.ident.specific.clones)
  
}
