
####################################  
# Run Augur analysis on clusters
####################################  

# Metadata info 
metadata.df <- seurat.combined@meta.data

# see which cols to keep
colnames(metadata.df)

# Select just cluster and condition cols 
metadata.df <- metadata.df %>%
  dplyr::select(seurat_clusters, condition)

# Name cols for input into augur function
colnames(metadata.df) <- c("cell_type", "label")

# output expression dataset
normed.df <- seurat.combined@assays$RNA@data

# Run Augur analysis // ~45 min compute time
augur.out <- Augur::calculate_auc(normed.df,
                                  metadata.df, 
                                  n_threads = 8)

# Output is stored in AUC slot 
augur.out$AUC

# Save RDS file
saveRDS(augur.out, "Exported_RDS_files/Augur_output_clusters.rds")
rm(augur.out)



###################################
# Run Augur analysis on modules
###################################

# Metadata info 
metadata.df <- seurat.combined@meta.data

# see which cols to keep
colnames(metadata.df)

# Select just cluster and condition cols 
metadata.df <- metadata.df %>%
  dplyr::select(Module, condition)

# Name cols for input into augur function
colnames(metadata.df) <- c("cell_type", "label")


# output expression dataset
normed.df <- seurat.combined@assays$RNA@data

# Run Augur analysis // ~10 min compute time
augur.out <- Augur::calculate_auc(normed.df,
                                  metadata.df, 
                                  n_threads = 8)

# Output is stored in AUC slot 
augur.out$AUC

# Save RDS file
saveRDS(augur.out, "Exported_RDS_files/Augur_output_module.rds")
rm(augur.out)