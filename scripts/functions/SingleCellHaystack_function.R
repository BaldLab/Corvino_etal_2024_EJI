scHay <- function(seurat.object, output.dir){
  
  # Create output directory
  if(!dir.exists(paste0(output.dir))){
    dir.create(paste0(output.dir), 
               recursive = T)
  }
  
  # load required libraries
  print(paste0("loading packages"))
  
  library("Seurat")
  library("ggplot2")
  library("singleCellHaystack")
  
  # Set seed for reproducibility 
  set.seed(42)
  
  ##############################
  # Set up input variables
  ##############################
  
  print(paste0("Running haystack DEG calculation"))
  
  # Run haystack DEG calculation
  res.umap <- haystack(seurat.object,
                       assay = "RNA",
                       slot = "data",
                       coord = "umap",
                       dims = NULL,
                       cutoff = 1,
                       method = "highD",
                       use.advanced.sampling = NULL)
  
  # top 10 DEGs
  print(show_result_haystack(res.haystack = res.umap,
                             n = 10))
  
  
  
  # Export data
  print(paste0("exporting results to output.dir"))
  
  write.table(res.umap$results, 
              paste0(output.dir, "DEG_singleCellHaystack.txt"), 
              sep = "\t",
              quote = FALSE)
  
  return(res.umap$results)
  
}





