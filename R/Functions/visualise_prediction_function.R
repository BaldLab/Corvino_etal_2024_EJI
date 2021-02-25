visualise.prediction <- function(input.seurat,
                                 sig.name = "IFNsig", 
                                 cluster.var = "seurat_clusters"){
  
  print(
    VlnPlot(input.seurat, 
            feature = paste0(sig.name, "1"),
            group.by = paste0(sig.name, "_Prediction"),
            pt.size = 0.1))
  
  print(
    UMAPPlot(input.seurat, 
             pt.size = 1, 
             group.by = paste0(cluster.var), 
             split.by = paste0(sig.name, "_Prediction"))
  )
  
  print(
    VlnPlot(input.seurat, 
            feature = paste0(sig.name, "1"),
            group.by = paste0(cluster.var),
            split.by = paste0(sig.name, "_Prediction"),
            pt.size = 0.1)    
  )
}



