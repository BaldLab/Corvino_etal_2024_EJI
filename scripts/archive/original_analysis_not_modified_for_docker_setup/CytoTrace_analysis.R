#####################
# Install package
#####################

#BiocManager::install("sva")

# change with path to your local download of CytoTRACE.tar.gz file
#devtools::install_local("~/Documents/Work/Sciebo/Scripts/Packages/CytoTRACE_0.3.3.tar.gz")

# run in terminal 
# pip install scanoramaCT
# pip install numpy

# Example dataset
#data.file <- read.delim("CytoTRACE_Analysis/Example_datasets/Bone_marrow_10x_matrix.txt")
#data.pheno <- read.delim("CytoTRACE_Analysis/Example_datasets/Bone_marrow_10x_metadata.txt")


##############################
# Create output directories
##############################

if(!dir.exists("output/figures/CytoTRACE")){
  dir.create("output/figures/CytoTRACE", 
             recursive = T)}


# Prepare data

# Expression matrix (Raw data)
data.file <- as.data.frame(seurat.combined@assays$RNA@counts)
head(data.file[1:10, 1:10])

# Metadata file
data.pheno <- as.data.frame(seurat.combined@meta.data)
head(data.pheno)

# Run cytoTRACE
Trace.object <- CytoTRACE(data.file,
                          enableFast = TRUE,
                          ncores = 6, 
                          subsamplesize = 1000)


# Create phenotype label vector
pheno.labels <- as.character(data.pheno$seurat_clusters)
names(pheno.labels) <- rownames(data.pheno)


# Plot data
plotCytoTRACE(Trace.object,
              emb = data.pheno[ , c("UMAP_1", "UMAP_2")],
              phenotype = pheno.labels,
              outputDir = "output/figures/CytoTRACE/")



plotCytoGenes(Trace.object, 
              numOfGenes = 50,
              outputDir = "output/figures/CytoTRACE/")




# Remove variables to save space and prevent errors
rm(data.file, data.pheno, Trace.object, pheno.labels)
