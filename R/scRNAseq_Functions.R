# Single cell RNA-sequencing functions 

# UMAP optimisation function 
UMAP.optimise <- function(input.seurat, 
                          output.dir = NULL, 
                          init.min.dist = 0.001, 
                          init.neigh = 5,
                          min.dist.step = 2,
                          neigh.step = 20){
  # Info: 
  # n.neighbor sensible values = 5 - 50 
  # min.dist sensible values = 0.001 - 0.5
  
  # Create directory if one isnt already created
  if(is.null(output.dir)){
    if(!dir.exists("output/Optimising_UMAP")){
      dir.create("output/Optimising_UMAP", 
                 recursive = TRUE)}
    
    output.dir <- "output/Optimising_UMAP/"
  }
  
  
  
  # set starting min.dist
  min.dist.val <- init.min.dist
  
  while(min.dist.val < 1){
    
    # set starting neighbor value
    n.neigh.val <- init.neigh
    
    while(n.neigh.val <= 50){
      
      print(paste0("Calculating UMAP for min.dist = ", min.dist.val, " & n.neighbor = ", n.neigh.val))
      
      # Calculate UMAP
      temp <- RunUMAP(object = input.seurat,
                      reduction = "pca",
                      dims = 1:20,
                      umap.method = "uwot",
                      n.neighbors = n.neigh.val,
                      min.dist = min.dist.val, 
                      seed.use = 42)
      
      # Plot UMAP
      print(UMAPPlot(object = temp,
                     label = TRUE, 
                     label.size = 4) + 
              ggtitle(paste0("UMAP Min dist = ", min.dist.val, " n.val = ", n.neigh.val)) + 
              NoLegend())
      
      dev.copy(pdf, paste0(output.dir, "UMAP_Min_dist_", min.dist.val, "_neighval_", n.neigh.val, ".pdf"))
      dev.off()
      
      n.neigh.val <- n.neigh.val + neigh.step
    }
    min.dist.val <- min.dist.val * min.dist.step
  }
  
  rm(temp)
  
}









