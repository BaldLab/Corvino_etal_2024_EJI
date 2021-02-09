
# Create output directory
if(!dir.exists("output/figures/singleCellHaystack")){
  dir.create("output/figures/singleCellHaystack")
}
library("Seurat")
# Set seed for reproducibility 
set.seed(42)

##############################
# Set up input variables
##############################

# Get UMAP dims
umap.dims <- seurat.combined@reductions$umap@cell.embeddings
umap.dims <- as.data.frame(umap.dims)

# Get expression data
exp.data <- seurat.combined@assays$integrated@data
exp.data <- as.matrix(exp.data)

summary(rowSums(exp.data))

# Determine gene detection

# Above count 1
detected.genes <- exp.data > 1

# Number of detected genes after thresholding
general.detection = apply(detected.genes, 2, sum)
ggplot(umap.dims, aes(x = UMAP_1, y = UMAP_2, colour = general.detection)) + labs(x = "UMAP1", y = "UMAP2") +
  geom_point(size=2) + scale_color_gradient(low="white", high="red") + labs(color = "Det. genes")


# Run haystack DEG calculation
res.umap <- haystack(umap.dims, 
                     detection = detected.genes, 
                     method = "2D")

class(res.umap)

# top 10 DEGs
show_result_haystack(res.haystack = res.umap,
                     n = 10)


# visualize genes of interest
plot_gene_haystack(umap.dims,
                   expression = exp.data,
                   gene = "CCL4",
                   detection = detected.genes,
                   high.resolution = TRUE,
                   point.size = 1)


# set seurat params
DefaultAssay(seurat.combined) <- "integrated"
Idents(seurat.combined) <- seurat.combined@meta.data$seurat_clusters

# vis with featureplot
FeaturePlot(seurat.combined, 
            feature = "HAVCR2", 
            order = TRUE,
            pt.size = 1)

# Vis with Vln
VlnPlot(seurat.combined, 
        assay = "RNA",
        feature = "PDCD1", 
        pt.size = 0)


# get the top most significant genes, and cluster them by their distribution pattern in the 2D plot
sorted.table <- show_result_haystack(res.haystack = res.umap,
                                     p.value.threshold = 1e-10)

gene.subset <- row.names(sorted.table)
length(gene.subset) # 1,353


# Export data
write.table(sorted.table, 
            "output/tables/DEG_singleCellHaystack.txt", 
            sep = "\t",
            quote = FALSE)


# Use k-means clustering to group DEGs
km <- kmeans_haystack(umap.dims, 
                      detection = detected.genes,
                      genes = gene.subset,
                      k = 10)

km.clusters <- km$cluster


###################
# Visualise DEGs
###################

goi <- show_result_haystack(res.haystack = res.umap,
                            p.value.threshold = 1e-10, 
                            n = 20)
goi <- rownames(goi)

for(i in seq_along(goi)){
  
  # Featureplot
  print(FeaturePlot(seurat.combined, 
                    feature = goi[i], 
                    reduction = "umap",
                    order = TRUE,
                    pt.size = 1))
  
  dev.copy(pdf, paste0("output/figures/singleCellHaystack/FeaturePlot_", goi[i], ".pdf"))
  dev.off()
  
  # Vlnplot 
  print(VlnPlot(seurat.combined, 
                assay = "RNA",
                feature = goi[i], 
                pt.size = 0))
  
  dev.copy(pdf, paste0("output/figures/singleCellHaystack/VlnPlot_", goi[i], ".pdf"))
  dev.off()
  
  
}



# Vis DEG clusters

seurat.combined.small <- subset(seurat.combined, downsample = 300)


for(i in seq_len(max(km.clusters))){
  goi <- names(km.clusters[km.clusters == i])
  
  print(DoHeatmap(seurat.combined.small, 
                  features = goi) + NoLegend())
  
  dev.copy(pdf, paste0("output/figures/singleCellHaystack/DEG_heatmap_cluster_", i, ".pdf"))
  dev.off()
  
}


# Remove un-needed objects
rm(res.umap, 
   sorted.table, 
   umap.dims, 
   general.detection, 
   km.clusters, 
   detected.genes, 
   exp.data, 
   km)