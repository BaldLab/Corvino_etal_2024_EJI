# Complete dataset
seurat.combined <- readRDS("Exported_RDS_files/seurat_combined.rds")

# set idents and assay to use
Idents(seurat.combined) <- seurat.combined@meta.data$Clusters_l1
DefaultAssay(seurat.combined) <- "RNA"



fig.output <- c("output/figures/Presentations")

if(!dir.exists(paste0(fig.output))){
  dir.create(paste0(fig.output), 
             recursive = T)
}




cluster.groups <- unique(seurat.combined@meta.data$Clusters_l1)


cluster.list <- list(Naive = c("Naive_like_1_CM", "Naive_like_3", "Naive_like_2_SC"), 
                     Exhausted = c("Exhausted_1", "Exhausted_2", "Stimulated_exhausted"),
                     Innate = c("gd_T_g9d2", "gd_T_non_g9d2", "MAIT"),
                     Misc = c("TRM", "Proliferative"), 
                     effector = c("Cytotoxic", "Type_I_IFN", "Stimulated_1"))






for(i in seq_along(cluster.list)){
temp.cols <- rep("Grey", length(clust.cols))

logic.var <- grepl(paste0(cluster.list[[i]], collapse = "|"), names(clust.cols))

temp.cols[logic.var] <- clust.cols[logic.var]


name.val <- names(cluster.list[i])
print(UMAPPlot(seurat.combined, 
               shuffle = F,
               pt.size = 1,
               cols = temp.cols,
               group.by = "Clusters_l1")) + 
  ggtitle(paste0(name.val))


dev.copy(pdf, paste0(fig.output, "/UMAP_", name.val, ".pdf"))
dev.off()





}





for(i in seq_along(cluster.groups)){
temp.cols <- rep("Grey", length(clust.cols))
temp.cols[i] <- clust.cols[i]

print(UMAPPlot(seurat.combined, 
         shuffle = F,
         pt.size = 1,
         cols = temp.cols,
         group.by = "Clusters_l1"))




}


dev.copy(pdf, paste0(fig.output, "/UMAP_condition_coloured.pdf"))
dev.off()

small.seurat <- subset(seurat.combined, downsample = 100)

small.seurat <- ScaleData(small.seurat)


goi <- c("GZMA", "GZMB", "GZMM", "GZMK", "GNLY", "PRF1", "IFNG", "TNF", "IL2", 
         "XCL1", "XCL2", "TIGIT", "CD226", "PDCD1", "HAVCR2", "NKG7", "MKI67", "RORA", "TRAV1-2", "ITGAE", "LCK", "TGFB1", "IL32", "IL7R")

goi <- c("GZMA", "GZMB", "GZMM", "GZMK", "GNLY", "PRF1")
goi <- gene_list_plot[[1]]
DoHeatmap(small.seurat, 
          features = goi)




goi <- c("GZMA", "GZMB", "GZMM", "GZMK", "GNLY", "PRF1", "IFNG", "TNF", "IL2", 
         "XCL1", "XCL2", "TIGIT", "CD226", "PDCD1", "HAVCR2", "NKG7", "MKI67", "RORA", "TRAV1-2", "ITGAE", "LCK", "TGFB1", "IL32", "IL7R")



gene_list_plot <- list(goi.1 = c("GZMA", "GZMB", "GZMM", "GZMK", "GNLY", "PRF1"),
                       goi.2 = c("IFNG", "TNF", "XCL1", "XCL2"),
                       goi.3 = c("TCF7", "TOX", "PDCD1", "LAG3", "CD226"),
                       goi.4 = c("TRDV2", "TRAV1-2", "RORA", "ZNF683", "MKI67"))



levels(seurat.combined@meta.data$Clusters_l1)

levels(Idents(object = seurat.combined))


Idents(object = seurat.combined) <- seurat.combined@meta.data$Clusters_l1

seurat.combined <- RenameIdents(object = seurat.combined,
                                `Naive_like_1_CM` = "Naive-like-1",
                                `Naive_like_2_SC` = "Naive-like-2",
                                `Naive_like_3` = "Naive-like-3",
                                `Type_I_IFN` = "Type-I IFN",
                                `Stimulated_1` = "Stimulated-1",
                                `Stimulated_exhausted` = "Stimulated Exhausted",
                                `Exhausted_1` = "Exhausted-1",
                                `Exhausted_2` = "Exhausted-2",
                                `gd_T_g9d2` = "gd-T cell (G9D2)",
                                `gd_T_non_g9d2` = "gd-T cell (non-G9D2)")



seurat.combined@meta.data$Clusters_l1 <- Idents(object = seurat.combined)

names(clust.cols) <- c("Naive-like-1", "Naive-like-2", "Naive-like-3", 
                       "Cytotoxic", "Type-I IFN", "Stimulated-1",
                       "Stimulated Exhausted", "Exhausted-1", "Exhausted-2", 
                       "TRM", "gd-T cell (G9D2)", "gd-T cell (non-G9D2)",
                       "MAIT", "Proliferative")




for(i in seq_along(gene_list_plot)){
  
  
  print(Stacked_VlnPlot(seurat_object = seurat.combined,
                        assay = "RNA",
                        pt.size = 0,
                        features = gene_list_plot[[i]],
                        x_lab_rotate = TRUE,
                        plot_spacing = 0.15,
                        colors_use = clust.cols))
  
  
  dev.copy(pdf, paste0(fig.output, "/Stacked_VlnPlots_", i, ".pdf"))
  dev.off()
  
  
}


print(Stacked_VlnPlot(seurat_object = seurat.combined,
                      assay = "RNA",
                      pt.size = 0,
                      features = gene_list_plot[[1]],
                      x_lab_rotate = TRUE,
                      plot_spacing = 0.15,
                      colors_use = clust.cols))


UMAPPlot(seurat.combined, 
         cols = clust.cols,
         pt.size = 1) + NoLegend()



















print(Stacked_VlnPlot(seurat_object = seurat.combined,
                      assay = "RNA",
                      pt.size = 0,
                      features = gene_list_plot[[1]],
                      x_lab_rotate = TRUE,
                      plot_spacing = 0.15,
                      colors_use = clust.cols))




library("SeuratDisk")

trex.seurat <- SeuratDisk::LoadH5Seurat("Data/Trex_CD8_TCR_seurat.h5Seurat")

DimPlot(trex.seurat, 
        reduction = "tsne")


trex.seurat <- RunUMAP(trex.seurat, 
                       reduction = "Trex_UMAP",  
                       
                       reduction.key = "Trex_")


trex.seurat.V2 <- RunUMAP(object = trex.seurat,
                           reduction = "TrexPCR",
                           dims = 1:20,
                           umap.method = "uwot",
                           n.neighbors = 30, # 5 to 50
                           min.dist = 0.3, # Sensible values are in the range 0.001 to 0.5
                           seed.use = 42)



DimPlot(trex.seurat.V2, 
        reduction = "umap", 
        cols = clust.cols, label = T)




trex.seurat

for(i in seq_along(gene_list_plot)){
  
  
  print(Stacked_VlnPlot(seurat_object = seurat.combined,
                        assay = "RNA",
                        pt.size = 0,
                        features = gene_list_plot[[i]],
                        x_lab_rotate = TRUE,
                        plot_spacing = 0.15,
                        colors_use = clust.cols))
  
  
  #dev.copy(pdf, paste0(output.fig1, "/Stacked_VlnPlots_", i, ".pdf"))
  #dev.off()
  
  
}

levels(seurat.combined@meta.data$Clusters_l1)










goi <- c("GZMA", "GZMB", "GZMM", "GZMK", "GNLY", "PRF1", "IFNG", "TNF", "IL2", 
         "XCL1", "TIGIT", "CD226", "PDCD1", "HAVCR2", "NKG7", "MKI67", "RORA", "TRAV1-2", "ITGAE", "LCK", "TGFB1", "IL32", "IL7R")

DotPlot_scCustom(seurat_object = seurat.combined, 
                 features = goi,
                 flip_axes = F,
                 x_lab_rotate = TRUE)

#dev.copy(pdf, paste0(fig.output, "/Dotplot_chemokine_receptors.pdf"))
#dev.off()



goi <- c("GNLY", "PRF1", "GZMK", "GZMH", "GZMA", "TIGIT", "CD226", "PDCD1", "HAVCR2", "NKG7", "LCK")

for(i in 1:length(goi)){
  
  print(Plot_Density_Custom(seurat.combined, 
                            features = goi[i],
                            custom_palette = batlow.pal,
                            joint = FALSE, 
                            pt.size = 1, 
                            reduction = "umap"))
  
  dev.copy(pdf, paste0(fig.output, "/Nebulosa_RNA_", goi[i], ".pdf"))
  dev.off()
  
  
}
