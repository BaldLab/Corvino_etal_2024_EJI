
#remotes::install_github("xmc811/Scillus")

# Scillus package would not compile
# I took the functions and made a new .R file
# Some of the functions called and used in this script required modifications to work outside of the package environment


#source("~/Documents/Work/Sciebo/Scripts/Packages/Scillus_package/Scillus_functions.R")
#load("~/Documents/Work/Sciebo/Scripts/Packages/Scillus_package/pathways.hallmark.rda")


# Create output directories
if(!dir.exists("output/figures/GO_analysis/Scillus_GO")){
  dir.create("output/figures/GO_analysis/Scillus_GO", 
             recursive = T)
}


######################
# Cluster analysis
######################

# set idents to cluster
Idents(seurat.combined) <- seurat.combined@meta.data$seurat_clusters


# get a broad list of DEG markers 
Cluster.markers <- FindAllMarkers(seurat.combined,
                                  assay = "RNA",
                                  only.pos = TRUE,
                                  min.diff.pct = 0.1,
                                  logfc.threshold = 0,
                                  return.thresh = 0.01)




# Threshold marker list by P value and fold change for GO analysis 

GO.markers <- Cluster.markers[Cluster.markers$p_val_adj < 0.05, ]
GO.markers <- GO.markers[abs(GO.markers$avg_log2FC) > 0.25, ]
dim(Cluster.markers) # 8,635 genes
dim(GO.markers) # 4,890 genes

table(GO.markers$cluster) # number of genes DEG per cluster // used to justify using topn = 100

# Plot Gene Ontology analysis
# CC, BP, MF

# To plot graph for a specific cluster use below code
#plot_cluster_go(Module.markers, cluster_name = "Activated", org = "human", ont = "MF")

# CC 
pdf("output/figures/GO_analysis/Scillus_GO/GO_Clusters_CC.pdf", width = 20, height = 20)

plot_all_cluster_go(GO.markers,
                    topn = 100,
                    org = "human", 
                    ont = "CC")

dev.off()


# BP 
pdf("output/figures/GO_analysis/Scillus_GO/GO_Clusters_BP.pdf", width = 20, height = 20)

plot_all_cluster_go(GO.markers,
                    topn = 100,
                    org = "human", 
                    ont = "BP")

dev.off()

# MF 
pdf("output/figures/GO_analysis/Scillus_GO/GO_Clusters_MF.pdf", width = 20, height = 20)

plot_all_cluster_go(GO.markers,
                    topn = 100,
                    org = "human", 
                    ont = "MF")
dev.off()




# Perform GSEA 
#gsea_res <- test_GSEA(Cluster.markers, 
#                      pathway = pathways.hallmark)


#plot_GSEA(gsea_res)
#dev.copy(pdf, "output/figures/GO_analysis/Scillus_GO/GSEA_Clusters.pdf")
#dev.off()