setwd("/work_dir/")


# Setup -----------------------------------------------
require(tidyverse)
require(Seurat)
require(useful)
require(ggrepel)
require(clusterProfiler)
require(enrichplot)
require(org.Hs.eg.db)
require(UCell)
require(Nebulosa)
require(ggsci)
library(harmony)

source("scripts/functions.R")


# Bald dataset --------------------------------

count_dirs <- paste0(list.dirs(path = "data/Human_HNSCC_2019/Raw_input_datasets/GEX_only",
                        full.names = T,recursive = F),"/outs/filtered_feature_bc_matrix")
names(count_dirs) <- paste0("S",1:4)

counts_raw_ctrl <- Read10X(data.dir = count_dirs[1:2],
                                strip.suffix = T)
counts_raw_treated <- Read10X(data.dir = count_dirs[3:4],
                                strip.suffix = T)


seurat_ctrl <- CreateSeuratObject(counts_raw_ctrl,
                                  min.cells = 2,)
seurat_treated <- CreateSeuratObject(counts_raw_treated,
                                  min.cells = 2)

seurat_ctrl <- PercentageFeatureSet(seurat_ctrl,
                                    pattern = "^MT-",
                                    col.name = "percent.mito")

seurat_treated <- PercentageFeatureSet(seurat_treated,
                                    pattern = "^MT-",
                                    col.name = "percent.mito")

p1 <-VlnPlot(seurat_ctrl,features = c("nCount_RNA","nFeature_RNA","percent.mito"),log = T)
p2 <-VlnPlot(seurat_treated,features = c("nCount_RNA","nFeature_RNA","percent.mito"),log=T)
p1/p2


seurat_ctrl <- subset(seurat_ctrl, 
                      percent.mito < 6 &
                        nFeature_RNA >500)
seurat_treated <- subset(seurat_treated, 
                      percent.mito < 6 &
                        nFeature_RNA >500)

# nfeature > 200, <2500
# percent.mito < 10 

# Check for cells removed by higher cutoff


p1 <-VlnPlot(seurat_ctrl,features = c("nCount_RNA","nFeature_RNA","percent.mito"))
p2 <-VlnPlot(seurat_treated,features = c("nCount_RNA","nFeature_RNA","percent.mito"))
p1/p2

seurat_combined <- list(seurat_ctrl,seurat_treated)

seurat_combined <- lapply(X = seurat_combined, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = seurat_combined)
seurat.anchors <- FindIntegrationAnchors(object.list = seurat_combined,
                                         anchor.features = features)
seurat_combined <- IntegrateData(anchorset = seurat.anchors)
seurat_combined


DefaultAssay(seurat_combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_combined <- ScaleData(seurat_combined, verbose = FALSE,features = rownames(seurat_combined))
seurat_combined <- RunPCA(seurat_combined, npcs = 30, verbose = FALSE)
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:30)
seurat_combined <- FindNeighbors(seurat_combined, reduction = "pca", dims = 1:30)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)


seurat_combined$condition <- ifelse(seurat_combined$orig.ident%in%c("S1","S2"),"US","Stim")

# Visualization
p1 <- DimPlot(seurat_combined, reduction = "umap",group.by = "condition")
p2 <- DimPlot(seurat_combined, reduction = "umap", 
              label = TRUE, repel = TRUE) +
  scale_color_manual(values = pals::cols25())
p1 + p2

## Export data for velocity analysis -----------------

#load("results/saves/2021-12-08_after_signature.RData")

write.csv(Cells(seurat_object), file = "cellID_obs.csv", row.names = FALSE)



## Define signature -----------------------------------

FeaturePlot(seurat_combined,features = "IFI6")
VlnPlot(seurat_combined,features = "IFI6")

DefaultAssay(seurat_combined) <- "RNA"
inf_cluster_signature <- FindMarkers(seurat_combined,ident.1 = 3,
                                     logfc.threshold = 0,
                                     min.pct = 0) 

inf_cluster_signature <-
  inf_cluster_signature%>% 
  as_tibble(rownames = "gene_symbol") %>%
  mutate(DE=case_when(
    avg_log2FC>0 & p_val_adj <=.1 ~ "up",
    avg_log2FC<0 & p_val_adj <=.1 ~ "down",
    TRUE ~ "notDE"))

write_csv(x = inf_cluster_signature,
          file = paste0("results/exports/ifncluster_DE_result_",Sys.Date(),".csv"))


inf_cluster_signature %>%
  mutate(label=if_else(abs(avg_log2FC)>1 & p_val_adj<.1,
                       gene_symbol,"")) %>%
  ggplot(aes(avg_log2FC,y=-log10(p_val))) +
  geom_point(aes(color=avg_log2FC)) +
  geom_text_repel(aes(label=label),max.time = 2,max.overlaps = 20,segment.colour = "grey",size=3) +
  theme_minimal() +
  scale_color_distiller(palette = "RdBu",limits=c(-1,1),oob=scales::squish)


inf_cluster_signature %>%
  group_by(DE) %>% tally() 


DefaultAssay(seurat_combined) <- "RNA"
cl6_signature <- FindMarkers(seurat_combined,ident.1 = 6,
                                     logfc.threshold = 0.25,
                                     min.pct = 0) 

# All markers 

all_markers <- FindAllMarkers(seurat_combined,logfc.threshold = .4)

all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(seurat_combined, features = top10$gene,
          group.colors = pals::cols25()) + 
  NoLegend() + 
  scale_fill_viridis_c(option = "A")


write_csv(x = all_markers,
          file = paste0("results/exports/DE_result_allcluster_",Sys.Date(),".csv"))



save.image(file = paste0("results/saves/",Sys.Date(),"_after_signature.RData"))

  
## Characterize IFN cluster --------------------------

inf_cluster_signature_ea <- 
  inf_cluster_signature %>%
  mutate(metric = -log10(p_val+.0001)*avg_log2FC) %>%
  arrange(desc(metric))

gene_anno <- bitr(geneID = inf_cluster_signature_ea %>% pull(gene_symbol),
                  fromType = "SYMBOL",drop = T,
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db) %>%
  filter(!duplicated(SYMBOL))
  
inf_cluster_signature_ea <- inf_cluster_signature_ea %>%
  left_join(gene_anno %>% 
              as_tibble() %>%
              rename(SYMBOL="gene_symbol"))

  


geneList <- inf_cluster_signature_ea %>% filter(!is.na(ENTREZID)) %>% pull(metric)
names(geneList) <- inf_cluster_signature_ea %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)



# GO - GSEA
ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.1,
             verbose      = FALSE)

ego_df <- as.data.frame(ego)


ego_df %>% 
  as_tibble() %>%
  mutate(direction=if_else(NES>0,true = "up","down"),
         .before="pvalue") %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES),n=15) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(Description=fct_inorder(Description)) %>%
  ggplot(aes(x=NES,y=Description,fill=NES)) +
  geom_bar(stat="identity") +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_minimal()


# GO - EA
up_genes <- inf_cluster_signature_ea%>% filter(avg_log2FC>0 & p_val_adj < .1) %>%
  pull(ENTREZID) %>% na.omit()

ego_EA <- enrichGO(gene = up_genes,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.1,
             readable = T)

ego_EA_df <- as.data.frame(ego_EA)


ego_EA_df %>% 
  as_tibble() %>%
  slice_min(order_by = p.adjust,n=50) %>%
  mutate(ES =-log10(p.adjust) ) %>%
  arrange(ES) %>%
  mutate(Description=fct_inorder(Description)) %>%
  ggplot(aes(x=ES,y=Description,color=p.adjust,size=Count)) +
  geom_point() +
  scale_color_distiller(palette = "RdYlBu", limits=c(0,.1), oob=scales::squish,direction = 1) +
  theme_minimal()


heatplot(ego_EA, showCategory=50)


ego_EA_emap <- pairwise_termsim(ego_EA)



p1 <- emapplot(ego_EA_emap,showCategory = 100,
         cex_label_category = 1,
         cex_line = .5) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,.1),
                        oob=scales::squish,direction = 1) 

ggsave(filename = paste0("results/figures/EA_GO_emap_",Sys.Date(),".png"),
    width=20,height = 20,dpi = 600)



require(ReactomePA)

e_reactome <- enrichPathway(
                 gene = up_genes,
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.1,
                 readable = T)

e_reactome_df <- as.data.frame(e_reactome)




eractome_emap <- pairwise_termsim(e_reactome)



p1 <- emapplot(eractome_emap,showCategory = 100,
               cex_label_category = 1,
               cex_line = .5) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,.1),
                       oob=scales::squish,direction = 1) 

ggsave(plot = p1,filename = paste0("results/figures/EA_reactome_emap_",Sys.Date(),".png"),
       width=20,height = 20,dpi = 600)


inf_cluster_signature_ea %>% pull(avg_log2FC) -> fc_genes
names(fc_genes) <- inf_cluster_signature_ea %>% pull(gene_symbol)


cnetplot(e_reactome[grepl(pattern = "PTEN",e_reactome@result$Description),,asis=T],
         foldChange = fc_genes)



save.image(file = paste0("results/saves/",Sys.Date(),"_after_EA.RData"))


##




## Signature scoring - Tcell ---------------------------
load("results/saves/2021-10-29_after_EA.RData")


ifn_pos_signature <- inf_cluster_signature %>% 
  filter(avg_log2FC>=log2(2) , p_val_adj<=.1) %>%
  pull(gene_symbol)

ifn_signatures <- list(ifn_pos_signature=ifn_pos_signature)



### method 1 - ucell -------------------------------------------

DefaultAssay(seurat_combined) <- "RNA"
seurat_combined <- AddModuleScore_UCell(seurat_combined,
                                     features = ifn_signatures, 
                                     ncores = 2)


# Viz data
featnames <- paste0(names(ifn_signatures), "_UCell")
FeaturePlot(seurat_combined, features = featnames,
            pt.size = 0.1, max.cutoff = "q99", 
            ncol = 3,cols = viridis::magma(n = 50))

VlnPlot(seurat_combined, features = featnames,
        pt.size = 0, group.by = "seurat_clusters") +
  scale_fill_manual(values = pals::cols25()) +
  geom_boxplot(width=.1,fill="white")


plot_density(seurat_combined, featnames)


### method 2 - score -------------------------------------------

seurat_combined <- ScaleData(seurat_combined,features = as.vector(unlist(ifn_signatures)))
res_signatures <- calc_signature_score(seurat_combined,
                                       gene.signature = ifn_signatures)




seurat_combined <- AddMetaData(seurat_combined,metadata = res_signatures)



FeaturePlot(seurat_combined, features = names(res_signatures),
            pt.size = 0.1, max.cutoff = "q99", 
            ncol = 3,cols = viridis::magma(n = 50))



VlnPlot(seurat_combined, features = names(res_signatures),
        pt.size = 0, group.by = "seurat_clusters") +
  scale_fill_manual(values = pals::cols25())




### method 3 - AddModuleScore -------------------------------------------




DefaultAssay(seurat_combined) <- "RNA"
seurat_combined <- AddModuleScore(seurat_combined,
                                   features = ifn_signatures, 
                                   name = "srt_score_",
                                   ncores = 2)
colnames(seurat_combined@meta.data)[grepl(pattern = "srt_score_",colnames(seurat_combined@meta.data))] <-
  paste0("srt_score_",names(ifn_signatures))


FeaturePlot(seurat_combined, features = paste0("srt_score_",names(ifn_signatures)),
            pt.size = 0.1, max.cutoff = "q99", 
            ncol = 3,cols = viridis::magma(n = 50))


### method 4 - irGSEA ----------------------------

ifn_pos_signature <- inf_cluster_signature %>% 
  filter(avg_log2FC>=log2(2) , p_val_adj<=.1) %>%
  pull(gene_symbol)

ifn_signatures <- list(ifn_pos_signature=ifn_pos_signature)



require(irGSEA)
seurat_combined.final <- irGSEA.score(object = seurat_combined,
                                      assay = "RNA",slot = "data", 
                                      seeds = 123, ncores = 4,
                                      min.cells = 3, min.feature = 0,
                                      custom = T, geneset = ifn_signatures, msigdb = F, 
                                      species = "Homo sapiens", geneid = "symbol")


result.dge <- irGSEA.integrate(object = seurat_combined.final,
                               group.by = "seurat_clusters",
                               metadata = NULL, col.name = NULL)


irGSEA_methods=c("AUCell","UCell","singscore", "ssgsea")

lapply(irGSEA_methods, function(x){
  irGSEA.densityheatmap(object = seurat_combined.final,
                            method = x,
                            show.geneset = "ifn-pos-signature")
}) %>% wrap_plots(ncol = 4)


irGSEA_bp <- bind_rows(lapply(irGSEA_methods, function(x){
  as_tibble(seurat_combined.final@assays[[x]]@counts) %>%
  pivot_longer(cols = everything()) %>%
    mutate(method=x)
}))

irGSEA_bp %>%
  group_by(method) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  left_join(seurat_combined.final@meta.data %>% 
              as_tibble(rownames="name") %>% 
              dplyr::select(name,seurat_clusters)) %>%
  mutate(method=factor(method,levels = c("AUCell","UCell","singscore","ssgsea"))) %>%
  ggplot(aes(x=seurat_clusters,y=value_scaled,fill=seurat_clusters)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(facets = ~method,ncol = 4) +
  theme_bw() +
  scale_fill_manual(values = c(pals::cols25()))



seurat_combined@meta.data %>% as_tibble() %>% 
  dplyr::select(seurat_clusters,
                signature_score=Score_ifn_pos_signature,
                Seurat_AddModuleScore=srt_score_ifn_pos_signature) %>%
  pivot_longer(cols = -seurat_clusters) %>%
  group_by(name) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  ggplot(aes(x=seurat_clusters,y = value_scaled,
             fill=seurat_clusters)) +
  geom_boxplot(outlier.colour = NA) +
  facet_wrap(facets = ~name,ncol = 2)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))






### Percentages of enrichment ------------

seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>% 
  dplyr::select(cells,seurat_clusters,
                signature_score=Score_ifn_pos_signature,
                Seurat_AddModuleScore=srt_score_ifn_pos_signature) %>%
  left_join(tibble(
    cells=colnames(seurat_combined.final@assays$AUCell),
    AUCell=as.numeric(seurat_combined.final@assays$AUCell@counts),
    UCell=as.numeric(seurat_combined.final@assays$UCell@counts),
    singscore=as.numeric(seurat_combined.final@assays$singscore@counts),
    ssgsea=as.numeric(seurat_combined.final@assays$ssgsea@counts))) %>%
  pivot_longer(cols = c(-cells,-seurat_clusters)) %>%
  group_by(name) %>%
  mutate(value_scaled=scale(value)[,1]) %>%
  mutate(perc=if_else(value_scaled>1,true = "pos",false = "neg")) %>%
  ggplot(aes(x=seurat_clusters,fill=perc)) +
  geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(facets = ~name,ncol = 6)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))+
  labs(y="Percentage of cells enriched for IFN signature")





#ucell
seurat_combined$ifn_pos_signature_UCell_pos <-
  if_else(scale(seurat_combined$ifn_pos_signature_UCell)>1,"pos","neg")

p1 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  left_join(seurat_combined@reductions$umap@cell.embeddings %>% 
              as_tibble(rownames="cells")) %>%
  ggplot(aes(UMAP_1,UMAP_2,color=ifn_pos_signature_UCell_pos)) +
  geom_point(size=.1) +
  scale_color_manual(values = c("grey","black"))+
  ggtitle("UCell")

#signature_score

seurat_combined$Score_ifn_pos_signature_pos <-
  if_else(scale(seurat_combined$Score_ifn_pos_signature)>1,"pos","neg")

p2 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  left_join(seurat_combined@reductions$umap@cell.embeddings %>% 
              as_tibble(rownames="cells")) %>%
  ggplot(aes(UMAP_1,UMAP_2,color=Score_ifn_pos_signature_pos)) +
  geom_point(size=.1) +
  scale_color_manual(values = c("grey","black"))+
  ggtitle("Signature Score")

#modscore
seurat_combined$ifn_pos_signature_modscore_pos <-
  if_else(scale(seurat_combined$srt_score_ifn_pos_signature)>1,"pos","neg")

p3 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  left_join(seurat_combined@reductions$umap@cell.embeddings %>% 
              as_tibble(rownames="cells")) %>%
  ggplot(aes(UMAP_1,UMAP_2,color=ifn_pos_signature_modscore_pos)) +
  geom_point(size=.1) +
  scale_color_manual(values = c("grey","black"))+
  ggtitle("AddModuleScore")


(p1 + p2 +p3) & theme_bw() + 
  theme(legend.position = "bottom")

## bar chart 

p1 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,ifn_pos_signature_UCell_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=ifn_pos_signature_UCell_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom") +
  ggtitle("UCell")

p2 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,Score_ifn_pos_signature_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=Score_ifn_pos_signature_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom")+
  ggtitle("Signature Score")


p3 <- seurat_combined@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,ifn_pos_signature_modscore_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=ifn_pos_signature_modscore_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom")+
  ggtitle("AddModuleScore")

(p1 + p2 +p3) & theme_bw() + 
  theme(legend.position = "bottom")



## Signature scoring - Tcell wo Cluster 3 -----------------

### Score ------------------------------------

seurat_combined_minus_cl3 <- 
  subset(x = seurat_combined, 
         idents = c("3"), invert = TRUE)
seurat_combined_minus_cl3
DimPlot(seurat_combined_minus_cl3)

# ucell
DefaultAssay(seurat_combined_minus_cl3) <- "RNA"
seurat_combined_minus_cl3 <- AddModuleScore_UCell(seurat_combined_minus_cl3,
                                        features = ifn_signatures, 
                                        ncores = 2)


seurat_combined_minus_cl3$ifn_pos_signature_UCell_pos <-
  if_else(scale(seurat_combined_minus_cl3$ifn_pos_signature_UCell)>2,"pos","neg")

p1 <- seurat_combined_minus_cl3@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,ifn_pos_signature_UCell_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=ifn_pos_signature_UCell_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom") +
  ggtitle("UCell")

## Signature Score 
res_signatures_sub <- calc_signature_score(seurat_combined_minus_cl3,
                                       gene.signature = ifn_signatures)




seurat_combined_minus_cl3 <- AddMetaData(seurat_combined_minus_cl3,
                                         metadata = res_signatures_sub)

seurat_combined_minus_cl3$Score_ifn_pos_signature_pos <-
  if_else(scale(seurat_combined_minus_cl3$Score_ifn_pos_signature)>1,"pos","neg")

p2<-seurat_combined_minus_cl3@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,Score_ifn_pos_signature_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=Score_ifn_pos_signature_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom") +
  ggtitle("Signature Score")

# AddmModuleScore
seurat_combined_minus_cl3 <- AddModuleScore(seurat_combined_minus_cl3,
                                  features = ifn_signatures, 
                                  name = "srt_score_",
                                  ncores = 2)
colnames(seurat_combined_minus_cl3@meta.data)[grepl(pattern = "srt_score_",colnames(seurat_combined_minus_cl3@meta.data))] <-
  paste0("srt_score_",names(ifn_signatures))

seurat_combined_minus_cl3$ifn_pos_signature_modscore_pos <-
  if_else(scale(seurat_combined_minus_cl3$srt_score_ifn_pos_signature)>1,"pos","neg")

p3 <- seurat_combined_minus_cl3@meta.data %>% as_tibble(rownames = "cells") %>%
  group_by(seurat_clusters,ifn_pos_signature_modscore_pos) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,fill=ifn_pos_signature_modscore_pos)) +
  geom_bar(stat="identity",position="fill",width = .75)  +
  scale_fill_manual(values = c("grey","black")) +
  theme(legend.position = "bottom")+
  ggtitle("AddModuleScore")

(p1 + p2 +p3) & theme_bw() + 
  theme(legend.position = "bottom")



seurat_combined_minus_cl3.final <- irGSEA.score(object = seurat_combined_minus_cl3,
                                      assay = "RNA",slot = "data", 
                                      seeds = 123, ncores = 4,
                                      min.cells = 3, min.feature = 0,
                                      custom = T, geneset = ifn_signatures, 
                                      msigdb = F, 
                                      species = "Homo sapiens", 
                                      geneid = "symbol")




irGSEA_methods=c("AUCell","UCell","singscore", "ssgsea")

lapply(irGSEA_methods, function(x){
  irGSEA.densityheatmap(object = seurat_combined_minus_cl3.final,
                        method = x,
                        show.geneset = "ifn-pos-signature")
}) %>% wrap_plots(ncol = 4)


irGSEA_bp_woc3 <- bind_rows(lapply(irGSEA_methods, function(x){
  as_tibble(seurat_combined_minus_cl3.final@assays[[x]]@counts) %>%
    pivot_longer(cols = everything()) %>%
    mutate(method=x)
}))

irGSEA_bp_woc3 %>%
  group_by(method) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  left_join(seurat_combined_minus_cl3.final@meta.data %>% 
              as_tibble(rownames="name") %>% 
              dplyr::select(name,seurat_clusters)) %>%
  mutate(method=factor(method,levels = c("AUCell","UCell","singscore","ssgsea"))) %>%
  ggplot(aes(x=seurat_clusters,y=value_scaled,fill=seurat_clusters)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(facets = ~method,ncol = 4) +
  theme_bw() +
  scale_fill_manual(values = c(pals::cols25()))



seurat_combined_minus_cl3.final@meta.data %>% as_tibble() %>% 
  dplyr::select(seurat_clusters,
                signature_score=Score_ifn_pos_signature,
                Seurat_AddModuleScore=srt_score_ifn_pos_signature) %>%
  pivot_longer(cols = -seurat_clusters) %>%
  group_by(name) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  ggplot(aes(x=seurat_clusters,y = value_scaled,
             fill=seurat_clusters)) +
  geom_boxplot(outlier.colour = NA) +
  facet_wrap(facets = ~name,ncol = 2)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))




seurat_combined_minus_cl3.final@meta.data %>% as_tibble(rownames = "cells") %>% 
  dplyr::select(cells,seurat_clusters,
                signature_score=Score_ifn_pos_signature,
                Seurat_AddModuleScore=srt_score_ifn_pos_signature) %>%
  left_join(tibble(
    cells=colnames(seurat_combined_minus_cl3.final@assays$AUCell),
    AUCell=as.numeric(seurat_combined_minus_cl3.final@assays$AUCell@counts),
    UCell=as.numeric(seurat_combined_minus_cl3.final@assays$UCell@counts),
    singscore=as.numeric(seurat_combined_minus_cl3.final@assays$singscore@counts),
    ssgsea=as.numeric(seurat_combined_minus_cl3.final@assays$ssgsea@counts))) %>%
  pivot_longer(cols = c(-cells,-seurat_clusters)) %>%
  group_by(name) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  ggplot(aes(x=seurat_clusters,y=value_scaled,fill=seurat_clusters)) +
  geom_boxplot(outlier.colour = NA) +
  facet_wrap(facets = ~name,ncol = 6)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))






### Percentages ------------


seurat_combined_minus_cl3.final@meta.data %>% as_tibble(rownames = "cells") %>% 
  dplyr::select(cells,seurat_clusters,
                signature_score=Score_ifn_pos_signature,
                Seurat_AddModuleScore=srt_score_ifn_pos_signature) %>%
  left_join(tibble(
    cells=colnames(seurat_combined_minus_cl3.final@assays$AUCell),
    AUCell=as.numeric(seurat_combined_minus_cl3.final@assays$AUCell@counts),
    UCell=as.numeric(seurat_combined_minus_cl3.final@assays$UCell@counts),
    singscore=as.numeric(seurat_combined_minus_cl3.final@assays$singscore@counts),
    ssgsea=as.numeric(seurat_combined_minus_cl3.final@assays$ssgsea@counts))) %>%
  pivot_longer(cols = c(-cells,-seurat_clusters)) %>%
  group_by(name) %>%
  mutate(value_scaled=scale(value)[,1]) %>%
  mutate(perc=if_else(value_scaled>1,true = "pos",false = "neg")) %>%
  ggplot(aes(x=seurat_clusters,fill=perc)) +
  geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(facets = ~name,ncol = 6)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25())) +
  labs(y="Percentage of cells enriched for IFN signature")






save.image(file = paste0("results/saves/",Sys.Date(),"_after_bald_dataset_enrichment.RData"))




# External datasets -----------------------------

### cillo -----------------

# 
dirs <- str_remove(list.files(path = "data/cillo/",pattern = ".mtx.gz"),pattern = "_matrix.mtx.gz")
# run once!!
# for (i in dirs) {
#   dir.create(paste0("data/cillo/",i))
#   files <- list.files(path = "data/cillo/",pattern = i)
#   for (j in files) {
#      file.rename(from = paste0("data/cillo/",j),
#                  to = paste0("data/cillo/",i,"/",str_remove(j,paste0(i,"_"))))
#      }
# }

dirs <- list.dirs("data/cillo")[-1]
#dirs <- grep(pattern = "TIL",x = dirs,value = T)
names(dirs) <- str_remove(string = dirs,"data/cillo/")


## rename genes file 
genes_files <- list.files(path = "data/cillo/",recursive = T,pattern = "genes.tsv.gz")
for(i in genes_files){
  file.rename(from = i,
              to = str_replace(i,pattern = "genes.tsv.gz",
                               replacement = "features.tsv.gz"))
}

# After creation of the gene/barcode matrix, a cell-level filtering step was performed 
# to remove cells with either few genes per cell (<200) or many molecules per cell (>20,000). 
# Next, genes that were lowly expressed (fewer reads than 3 counts in 1% of cells, or genes 
# expressed in fewer than 1% of cells) across all samples were removed.


cillo_seurat <- Read10X(data.dir = dirs)
cillo_seurat = CreateSeuratObject(counts = cillo_seurat,
                                  min.cells = 1000)


cillo_seurat


cillo_seurat <- PercentageFeatureSet(cillo_seurat,
                                       pattern = "^MT-",
                                       col.name = "percent.mito")
cillo_seurat$grp <- "all_cells"

VlnPlot(cillo_seurat,features = c("nCount_RNA","nFeature_RNA","percent.mito"),
        group.by = "grp")



cillo_seurat <- subset(cillo_seurat, 
                      percent.mito < 15 &
                        nFeature_RNA >200 &
                      nCount_RNA<20000)
cillo_seurat



cillo_seurat <- ScaleData(cillo_seurat, verbose = FALSE)
cillo_seurat <- FindVariableFeatures(cillo_seurat)
cillo_seurat <- RunPCA(cillo_seurat, npcs = 30, verbose = FALSE)
cillo_seurat <- RunUMAP(cillo_seurat, reduction = "pca", dims = 1:30)
cillo_seurat <- FindNeighbors(cillo_seurat, reduction = "pca", dims = 1:30)
cillo_seurat <- FindClusters(cillo_seurat, resolution = 0.3)

cillo_seurat@meta.data$sample <- 
  str_split(rownames(cillo_seurat@meta.data),pattern = "_",n = 2,simplify = T)[,1]


cillo_seurat@meta.data$sample_id <- 
  paste(str_split(str_split(rownames(cillo_seurat@meta.data),
                            pattern = "_",n = 2,simplify = T)[,2],
                  pattern = "_",simplify = T)[,1],
        str_split(str_split(rownames(cillo_seurat@meta.data),
                            pattern = "_",n = 2,simplify = T)[,2],
                  pattern = "_",simplify = T)[,2],
        str_split(str_split(rownames(cillo_seurat@meta.data),
                            pattern = "_",n = 2,simplify = T)[,2],
                  pattern = "_",simplify = T)[,3],sep = "_")


cillo_seurat@meta.data$disease <- 
  str_split(cillo_seurat@meta.data$sample_id,pattern = "_",simplify = T)[,1]

cillo_seurat@meta.data$tissue <- 
  str_split(cillo_seurat@meta.data$sample_id,pattern = "_",simplify = T)[,3]

table(cillo_seurat@meta.data$disease)
table(cillo_seurat@meta.data$tissue)


p1 <- DimPlot(cillo_seurat,reduction = "umap",group.by = "tissue")+
  scale_color_jco()
p2 <- DimPlot(cillo_seurat,reduction = "umap",group.by = "disease") +
  scale_color_jama()

p1/p2


## Identify CD8 positive cells

DimPlot(cillo_seurat,reduction = "umap") +
  scale_color_manual(values = pals::cols25())

p2 <- VlnPlot(cillo_seurat,c("CD3E","CD8B"))+
  scale_fill_manual(values = pals::cols25())
p3 <- FeaturePlot(cillo_seurat,features = c("CD3E","CD8B"))
p2/p3


cillo_seurat_cd8 <- subset(cillo_seurat,idents= c(1,4,9))
cillo_seurat_cd8 <- ScaleData(cillo_seurat_cd8, verbose = FALSE)
cillo_seurat_cd8 <- FindVariableFeatures(cillo_seurat_cd8)
cillo_seurat_cd8 <- RunPCA(cillo_seurat_cd8, npcs = 30, verbose = FALSE)
cillo_seurat_cd8 <- RunUMAP(cillo_seurat_cd8, reduction = "pca", dims = 1:30)
cillo_seurat_cd8 <- FindNeighbors(cillo_seurat_cd8, reduction = "pca", dims = 1:30)
cillo_seurat_cd8 <- FindClusters(cillo_seurat_cd8, resolution = 0.5)

DimPlot(cillo_seurat_cd8,reduction = "umap") +
  scale_color_manual(values = pals::cols25())


#### score --------------------------------------------

require(irGSEA)

irGSEA_methods=c("AUCell","UCell")



cillo_seurat_cd8.final <- irGSEA.score(object = cillo_seurat_cd8,
                                       assay = "RNA",slot = "data", 
                                       seeds = 123, ncores = 4,
                                       min.cells = 3, min.feature = 0,
                                       method = irGSEA_methods,
                                       custom = T, geneset = ifn_signatures, msigdb = F, 
                                       species = "Homo sapiens", geneid = "symbol")


result.dge.cd8.final <- irGSEA.integrate(object = cillo_seurat_cd8.final,
                                         method = irGSEA_methods,
                                         group.by = "seurat_clusters",
                                         metadata = NULL, col.name = NULL)



lapply(irGSEA_methods, function(x){
  irGSEA.densityheatmap(object = cillo_seurat_cd8.final,
                        method = x,
                        show.geneset = "ifn-pos-signature")
}) %>% wrap_plots(nrow = 2)




irGSEA_bp <- bind_rows(lapply(irGSEA_methods, function(x){
  as_tibble(cillo_seurat_cd8.final@assays[[x]]@counts) %>%
    pivot_longer(cols = everything()) %>%
    mutate(method=x)
}))

irGSEA_bp %>%
  group_by(method) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  left_join(cillo_seurat_cd8.final@meta.data %>% 
              as_tibble(rownames="name") %>% 
              dplyr::select(name,seurat_clusters)) %>%
  mutate(method=factor(method,levels = c("AUCell","UCell"))) %>%
  ggplot(aes(x=seurat_clusters,y=value_scaled,fill=seurat_clusters)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(facets = ~method,nrow = 2) +
  theme_bw() +
  scale_fill_manual(values = c(pals::cols25()))





#### Percentages ------------

cillo_seurat_cd8.final@meta.data %>% as_tibble(rownames = "cells") %>% 
  dplyr::select(cells,seurat_clusters) %>%
  left_join(tibble(
    cells=colnames(cillo_seurat_cd8.final@assays$AUCell),
    AUCell=as.numeric(cillo_seurat_cd8.final@assays$AUCell@counts),
    UCell=as.numeric(cillo_seurat_cd8.final@assays$UCell@counts))) %>%
  pivot_longer(cols = c(-cells,-seurat_clusters)) %>%
  group_by(name) %>%
  mutate(value_scaled=scale(value)[,1]) %>%
  mutate(perc=if_else(value_scaled>1,true = "pos",false = "neg")) %>%
  ggplot(aes(x=seurat_clusters,fill=perc)) +
  geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(facets = ~name,ncol = 6)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))+
  labs(y="Percentage of cells enriched for IFN signature")


cillo_seurat_cd8.final$score_AUCell=as.numeric(cillo_seurat_cd8.final@assays$AUCell@counts)
cillo_seurat_cd8.final$score_UCell=as.numeric(cillo_seurat_cd8.final@assays$UCell@counts)


p1 <- FeaturePlot(object = cillo_seurat_cd8.final,
            features = "score_AUCell",pt.size = .35) +
  scale_color_viridis_c(option = "A")

p2 <- FeaturePlot(object = cillo_seurat_cd8.final,
            features = "score_UCell",pt.size = .35) +
  scale_color_viridis_c(option = "A")


p1 | p2




p1 <- cillo_seurat_cd8.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,disease,tissue) %>%
  group_by(seurat_clusters,disease,tissue) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,
             fill=disease))+
  geom_bar(stat="identity",position = "dodge") +
  theme_bw()+
  scale_fill_jco()+
  labs(y="# of cells enriched for IFN signature")


p2 <- cillo_seurat_cd8.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,disease,tissue) %>%
  group_by(seurat_clusters,disease,tissue) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,y=n,
             fill=tissue))+
  geom_bar(stat="identity",position = "fill") +
  theme_bw()+
  scale_fill_jama()+
  labs(y="# of cells enriched for IFN signature")


p1 | p2


cillo_seurat_cd8.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,disease,tissue) %>%
  group_by(seurat_clusters,disease,tissue) %>%
  tally()  %>% dplyr::filter(seurat_clusters==11)


cillo_seurat_cd8.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,disease,tissue) %>%
  group_by(seurat_clusters,disease,tissue) %>%
  tally()  %>% dplyr::filter(seurat_clusters==11)



### utility ------------------------------------------

tcell_utility <- readRDS("data/utility/data/processedData/filtered_seuratObjects_harmony.rds")
tcell_utility


p1 <- DimPlot(tcell_utility,group.by = "consensus.major") +
  scale_color_manual(values = pals::cols25())
p2 <- DimPlot(tcell_utility,group.by = "consensus.Tcell") +
  scale_color_jama(na.value="grey")
p3 <- DimPlot(tcell_utility,group.by = "functional.cluster") +
  scale_color_manual(values = c("#e5f5f9",
                                "#fde0dd","#fa9fb5","#dd3497","#ae017e","#49006a",
                                "#99d8c9","#41ae76","#006d2c"),na.value = "grey")

p1|p2|p3


cd8_types <- unique(grep(pattern = "CD8_",x = tcell_utility$functional.cluster,value = T))
tcell_utility_cd8 <- subset(tcell_utility, functional.cluster %in% cd8_types) 
table(tcell_utility_cd8$functional.cluster)

#### harmony ----------------
tcell_utility_cd8_combined <- tcell_utility_cd8 %>% 
  RunHarmony("Cohort", plot_convergence = TRUE,verbose=T)

tcell_utility_cd8_combined <- tcell_utility_cd8_combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()


DimPlot(tcell_utility_cd8_combined,
        group.by = "functional.cluster") +
  scale_color_manual(values = c("#fde0dd","#fa9fb5","#dd3497","#ae017e","#49006a")) 







p1 <- tcell_utility_cd8_combined@meta.data %>% 
  dplyr::select(seurat_clusters,Type,Tissue) %>%
  group_by(across(everything())) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,fill=Tissue,y=n))+
  geom_bar(stat="identity",position="fill") +
  scale_fill_jco() +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percentage of cells")

p2 <- tcell_utility_cd8_combined@meta.data %>% 
  dplyr::select(seurat_clusters,Type,functional.cluster) %>%
  group_by(across(everything())) %>%
  tally() %>%
  ggplot(aes(x=seurat_clusters,fill=functional.cluster,y=n))+
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values =  c("#fde0dd","#fa9fb5","#dd3497","#ae017e","#49006a")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percentage of cells")


p1/p2


#### cca -----------------------------

# tcell_utility_cd8 <- SplitObject(tcell_utility_cd8,split.by = "Cohort")
# tcell_utility_cd8 <- lapply(X = tcell_utility_cd8, function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# features <- SelectIntegrationFeatures(object.list = tcell_utility_cd8)
# seurat.anchors <- FindIntegrationAnchors(object.list = tcell_utility_cd8,
#                                          anchor.features = features)
# tcell_utility_cd8 <- IntegrateData(anchorset = seurat.anchors)
# tcell_utility_cd8

#DefaultAssay(tcell_utility_cd8_combined) <- "integrated"



#### score ---------------


DimPlot(tcell_utility_cd8_combined,group.by = "functional.cluster") +
  scale_color_manual(values = c("#e5f5f9",
                                "#fde0dd","#fa9fb5","#dd3497","#ae017e","#49006a",
                                "#99d8c9","#41ae76","#006d2c"),na.value = "grey")
DimPlot(tcell_utility_cd8_combined,group.by = "seurat_clusters") +
  scale_color_manual(values = pals::cols25(),na.value = "grey")



p1 <- DimPlot(tcell_utility_cd8_combined,group.by = "Cohort") +
  scale_color_jco()
p2 <- DimPlot(tcell_utility_cd8_combined,group.by = "Type")+
  scale_color_jco()
p3 <- DimPlot(tcell_utility_cd8_combined,group.by = "Tissue")+
  scale_color_jco()
p4 <- DimPlot(tcell_utility_cd8_combined,group.by = "Sorted")+
  scale_color_jco()



(p1+p2) / (p3+p4)


# save Robject
saveRDS(tcell_utility_cd8_combined,file = paste0("results/saves/",Sys.Date(),"_tcell_utility_cd8_combined_harmony.RDS"))



require(irGSEA)

irGSEA_methods=c("AUCell","UCell")



tcell_utility_cd8_combined.final <- irGSEA.score(object = tcell_utility_cd8_combined,
                                       assay = "RNA",slot = "data", 
                                       seeds = 123, ncores = 4,
                                       min.cells = 3, min.feature = 0,
                                       method = irGSEA_methods,
                                       custom = T, geneset = ifn_signatures, msigdb = F, 
                                       species = "Homo sapiens", geneid = "symbol")



lapply(irGSEA_methods, function(x){
  irGSEA.densityheatmap(object = tcell_utility_cd8_combined.final,
                        method = x,
                        show.geneset = "ifn-pos-signature")
}) %>% wrap_plots(nrow = 2)




irGSEA_utility <- bind_rows(lapply(irGSEA_methods, function(x){
  as_tibble(tcell_utility_cd8_combined.final@assays[[x]]@counts) %>%
    pivot_longer(cols = everything()) %>%
    mutate(method=x)
}))

irGSEA_utility %>%
  group_by(method) %>%
  mutate(value_scaled=scales::rescale(value)) %>%
  left_join(tcell_utility_cd8_combined.final@meta.data %>% 
              as_tibble(rownames="name") %>% 
              dplyr::select(name,seurat_clusters)) %>%
  mutate(method=factor(method,levels = c("AUCell","UCell"))) %>%
  ggplot(aes(x=seurat_clusters,y=value_scaled,fill=seurat_clusters)) +
  geom_boxplot(outlier.color = NA) +
  facet_wrap(facets = ~method,nrow = 2) +
  theme_bw() +
  scale_fill_manual(values = c(pals::cols25()))





#### Percentages ------------

tcell_utility_cd8_combined.final@meta.data %>% as_tibble(rownames = "cells") %>% 
  dplyr::select(cells,seurat_clusters) %>%
  left_join(tibble(
    cells=colnames(tcell_utility_cd8_combined.final@assays$AUCell),
    AUCell=as.numeric(tcell_utility_cd8_combined.final@assays$AUCell@counts),
    UCell=as.numeric(tcell_utility_cd8_combined.final@assays$UCell@counts))) %>%
  pivot_longer(cols = c(-cells,-seurat_clusters)) %>%
  group_by(name) %>%
  mutate(value_scaled=scale(value)[,1]) %>%
  mutate(perc=if_else(value_scaled>1,true = "pos",false = "neg")) %>%
  ggplot(aes(x=seurat_clusters,fill=perc)) +
  geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(facets = ~name,ncol = 6)+
  theme_bw()+
  scale_fill_manual(values = c(pals::cols25()))+
  labs(y="Percentage of cells enriched for IFN signature")


tcell_utility_cd8_combined.final$score_AUCell=as.numeric(tcell_utility_cd8_combined.final@assays$AUCell@counts)
tcell_utility_cd8_combined.final$score_UCell=as.numeric(tcell_utility_cd8_combined.final@assays$UCell@counts)


p1 <- FeaturePlot(object = tcell_utility_cd8_combined.final,
                  features = "score_AUCell",pt.size = .35) +
  scale_color_viridis_c(option = "A")

p2 <- FeaturePlot(object = tcell_utility_cd8_combined.final,
                  features = "score_UCell",pt.size = .35) +
  scale_color_viridis_c(option = "A")


p1 | p2




p1 <- tcell_utility_cd8_combined.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,Tissue) %>%
  filter(seurat_clusters==7) %>%
  group_by(Tissue) %>%
  tally() %>%
  ggplot(aes(x=Tissue,y=n,
             fill=Tissue))+
  geom_bar(stat="identity",
           position = "dodge") +
  theme_bw()+
  scale_fill_jco()+
  labs(title="Origin of Cluster 7 cells", y="# of cells") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


p2 <- tcell_utility_cd8_combined.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,functional.cluster) %>%
  filter(seurat_clusters==7) %>%
  group_by(functional.cluster) %>%
  tally() %>%
  ggplot(aes(x=functional.cluster,y=n,
             fill=functional.cluster))+
  geom_bar(stat="identity",
           position = "dodge") +
  theme_bw()+
  scale_fill_manual(values =  c("#fde0dd","#fa9fb5","#dd3497","#ae017e","#49006a")) +
  labs(title="Origin of Cluster 7 cells", y="# of cells")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


p3 <- tcell_utility_cd8_combined.final@meta.data %>% as_tibble() %>%
  dplyr::select(seurat_clusters,Type) %>%
  filter(seurat_clusters==7) %>%
  group_by(Type) %>%
  tally() %>%
  ggplot(aes(x=Type,y=n,
             fill=Type))+
  geom_bar(stat="identity",
           position = "dodge") +
  theme_bw()+
  scale_fill_jama()+
  labs(title="Origin of Cluster 7 cells", y="# of cells")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


p1 + p2 + p3


