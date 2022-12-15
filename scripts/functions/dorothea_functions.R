
# Functions for TF prediction out of dorothea/ SCpubr workflow
add.dorothea.data <- function(input.seurat, activities = activities){
  
  input.seurat[["dorothea"]] <- activities %>% 
    dplyr::filter(.data$statistic == "norm_wmean") %>% 
    tidyr::pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") %>% 
    tibble::column_to_rownames("source") %>% 
    Seurat::CreateAssayObject()
  
  Seurat::DefaultAssay(input.seurat) <- "dorothea"
  
  input.seurat <- Seurat::ScaleData(object = input.seurat,
                                    verbose = FALSE)
  
  return(input.seurat)
}

extract.top.TFs.dorothea <- function(input.seurat, 
                                     group.var = NULL, 
                                     return.table = FALSE,
                                     n_tfs = 10){
  
  if (is.null(group.var)) {
    input.seurat$group.var <- Seurat::Idents(input.seurat)
  }else{
    input.seurat$group.var <- input.seurat@meta.data[, group.var]
  }
  
  df <- t(as.matrix(input.seurat@assays$dorothea@scale.data)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "cell") %>% 
    dplyr::left_join(y = {
      input.seurat@meta.data[, "group.var", drop = FALSE] %>% 
        tibble::rownames_to_column(var = "cell")
    }, by = "cell") %>% 
    dplyr::select(-.data$cell) %>% 
    tidyr::pivot_longer(cols = -"group.var",
                        names_to = "source", 
                        values_to = "score") %>% 
    dplyr::group_by(.data$group.var,
                    .data$source) %>%
    dplyr::summarise(mean = mean(.data$score))
  
  if(return.table){
    return(df)
  }else{
    tfs <- df %>% 
      dplyr::ungroup() %>%
      dplyr::top_n(n_tfs, mean) %>%
      dplyr::group_by(.data$group.var) %>% 
      dplyr::arrange(group.var, mean) %>% 
      dplyr::pull(.data$source)
    
    tfs <- unique(tfs)
    
    return(tfs)
  }
}