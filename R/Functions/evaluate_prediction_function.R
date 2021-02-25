

evaluate.prediction <- function(input.seurat, 
                                sig.name = "IFNsig",
                                validation.column, 
                                validation.cluster){
  
  
  output.temp <- input.seurat@meta.data %>%
    mutate(response_var = case_when(get(validation.column) == get("validation.cluster") ~ 1,
                                    TRUE ~ 0), 
           prediction_var = case_when(get(paste0(sig.name, "_Prediction")) == "IFN_cell" ~ 1, 
                                      TRUE ~ 0))
  
  roc(response = output.temp$response_var,
      predictor = output.temp$prediction_var, 
      ci = TRUE, 
      ci.alpha = 0.9,
      stratified = FALSE,
      plot = TRUE, 
      auc.polygon = TRUE, 
      max.auc.polygon = TRUE, 
      grid = TRUE, 
      print.auc = TRUE, 
      show.thres = TRUE)
  
}