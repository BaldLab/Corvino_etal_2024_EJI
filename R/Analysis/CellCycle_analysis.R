# Analyse if there is any enrichment for cell cycle stage

if(!dir.exists("output/figures/cell_cycle")){
  dir.create("output/figures/cell_cycle")
}

##################################
# Ensure idents is set up
##################################

Idents(seurat.combined) <- seurat.combined@meta.data$seurat_clusters

################
# Read in data
################


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
# We can segregate this list into markers of G2/M phase and markers of S phase

cc.genes$s.genes
cc.genes$g2m.genes

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


#######################
# Cell cycle analysis
#######################

seurat.combined <- CellCycleScoring(seurat.combined, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, 
                                    set.ident = FALSE)


head(seurat.combined@meta.data)


#######################
# Visualisation
#######################

UMAPPlot(seurat.combined, 
         group.by = "Phase",
         split.by = "condition",
         pt.size = 1)

dev.copy(pdf, "output/figures/cell_cycle/UMAP_CellCycle_Phase.pdf")
dev.off()


FeaturePlot(seurat.combined, 
            reduction = "umap",
            features = "S.Score",
            pt.size = 1.5, 
            label = TRUE, 
            label.size = 5,
            order = TRUE, 
            cols = c("lightgray", "orange", "red")) + NoLegend()

dev.copy(pdf, "output/figures/cell_cycle/FeaturePlot_S.pdf")
dev.off()

FeaturePlot(seurat.combined, 
            reduction = "umap",
            features = "G2M.Score",
            pt.size = 1.5, 
            label = TRUE, 
            label.size = 5,
            order = TRUE, 
            cols = c("lightgray", "orange", "red")) + NoLegend()

dev.copy(pdf, "output/figures/cell_cycle/FeaturePlot_G2M.pdf")
dev.off()

#######################
# Vln plots
#######################

VlnPlot(seurat.combined, 
        features = "S.Score", 
        pt.size = 0, 
        cols = clust.cols) + NoLegend()

dev.copy(pdf, "output/figures/cell_cycle/VlnPlot_Sscore.pdf")
dev.off()


VlnPlot(seurat.combined,
        features = "G2M.Score",
        pt.size = 0, 
        cols = clust.cols) + NoLegend()

dev.copy(pdf, "output/figures/cell_cycle/VlnPlot_G2Mscore.pdf")
dev.off()

#######################
# Ridge plots
#######################

RidgePlot(seurat.combined, 
          features = c("S.Score", "G2M.Score"), 
          log = TRUE, 
          cols = clust.cols)

dev.copy(pdf, "output/figures/cell_cycle/RidgePlot_S_and_G2Mscore.pdf")
dev.off()


table(seurat.combined@meta.data$Phase)
