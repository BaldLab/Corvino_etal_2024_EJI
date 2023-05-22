### technical setup ----------------------



rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
  
  require(matrixStats)
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  # Get the column means
  cm = rowMeans(x, na.rm = TRUE)
  # Get the column sd
  if (scale) {
    csd = rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}



calc_signature_score <- function(seurat.object,gene.signature) {
  if (is_empty(seurat.object@assays$RNA@scale.data)) stop("Please scale data")
  scaled_matrix <- Matrix::as.matrix(seurat.object@assays$RNA@scale.data)
  if (is.list(gene.signature)){
    res_enrichment =  lapply(gene.signature, function(x){
      overlap.gene.signature <- intersect(rownames(scaled_matrix),x)
      scaled_matrix_sub <- scaled_matrix[overlap.gene.signature,]
      as.vector(scale(apply(scaled_matrix_sub,2,sum)))
    })
    res_enrichment = as.data.frame(bind_cols(res_enrichment))
    rownames(res_enrichment) = colnames(seurat.object)
    colnames(res_enrichment) <- paste0("Score_",names(gene.signature))
  }
  #scaled_matrix_sub <- rowScale(scaled_matrix_sub)
  #scaled_matrix_sub <- na.omit(scaled_matrix_sub)
  #apply(scaled_matrix_sub,2,mean)
  res_enrichment
}
