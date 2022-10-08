
annotated.heatmap <- function(input.seurat,
                              cluster.id = "Clusters_l1",
                              assay.use = "RNA",
                              average.expression = TRUE,
                              goi, 
                              genes.to.label, 
                              col_order, 
                              col.colours, 
                              range.val = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)){
  
  #browser()
  
  
  # Load required packages
  #install_github("jokergoo/ComplexHeatmap")
  library(ComplexHeatmap)
  library(circlize)
  
  
  ####################
  # Prepare dataset
  ####################
  
  # Set Idents 
  Idents(input.seurat) <- input.seurat@meta.data[ , cluster.id]
  
  
  
  if(average.expression){
    input.seurat <- AverageExpression(input.seurat,
                                      assay = assay.use,
                                      slot = "data",
                                      verbose = TRUE,
                                      return.seurat = TRUE)
    
    input.seurat <- ScaleData(input.seurat)
    
    # Get matrix for heatmap
    heatmap.matrix <- GetAssayData(input.seurat,
                                   slot = "scale.data", 
                                   assay = assay.use)

    
    
    
      }else{
    
    # Downsample 
    input.seurat <- subset(input.seurat, downsample = 300)
    
    
    # Get matrix for heatmap
    heatmap.matrix <- GetAssayData(input.seurat,
                                   slot = "scale.data", 
                                   assay = assay.use)
    
    colnames(heatmap.matrix) <- input.seurat@meta.data[ , cluster.id]
  }
  


  
  # Keep only features of interest
  keep.logic <- rownames(heatmap.matrix) %in% goi
  heatmap.matrix <- heatmap.matrix[keep.logic, ]
  
  # Order columns
  x1 <- colnames(heatmap.matrix)
  
  f1 <- factor(x1, 
               levels = col_order)
  
  index.val <- order(f1)
  
  heatmap.matrix <- heatmap.matrix[, index.val]
  
  # Order rows
  x1 <- rownames(heatmap.matrix)
  ord <- c(goi)
  
  f1 <- factor(x1,
               levels = ord)
  
  index.val <- order(f1)
  
  heatmap.matrix <- heatmap.matrix[index.val, ]
  
  
  # Get column colour scheme 
  top.anno <- columnAnnotation(clusterID = colnames(heatmap.matrix), 
                               col = list(clusterID = col.colours))
  
  # Heatmap colour scheme
  col_fun <- colorRamp2(range.val, c(scico::scico(n= length(range.val), palette = "batlow")))
  

 # Colour scheme copied from seurat DoHeatmap which uses = PurpleAndYellow()
  
  
  # Do broad search first
  pattern.vec <- genes.to.label
  
  # only proceed with those found in dataset
  pattern.vec <- grep(paste0("^", pattern.vec, "$", collapse = "|"),
                      rownames(heatmap.matrix), 
                      value = TRUE)
  
  
  location.vec <- grep(paste0("^", pattern.vec, "$", collapse = "|"),
                       rownames(heatmap.matrix))
  
  
  row.anno.var <- rowAnnotation(foo = anno_mark(at = location.vec, labels = pattern.vec))
  
  
  col.split.var <- factor(colnames(heatmap.matrix),
                          levels = col_order)
  
  print(Heatmap(heatmap.matrix, 
                name = "Expression", 
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                show_row_names = FALSE, 
                show_column_names = FALSE,
                column_title_rot = 90,
                top_annotation = top.anno, 
                right_annotation = row.anno.var, 
                col = col_fun,
                column_split = col.split.var,
                use_raster = TRUE)) 

  
}
