# Environment Set up
rm(list = ls()) #Clean workspace
cat("\014")     #Clean Console
gc() # Free memory

###################
# Install packages
###################

pkgs <- c("remedy", "dplyr", "rstudioapi",
          "cowplot", "ggplot2", "grid", "gridExtra",
          "styler", "stringr", "inlmisc", "RColorBrewer",
          "readxl", "devtools", "tidyverse", "hdf5r", "scales",
          "useful", "renv", "pROC")

for(i in 1:length(pkgs)){
  if(!require(pkgs[i], character.only = T)){
    install.packages(pkgs[i])
    require(pkgs[i], character.only = T)
  }else{
    require(pkgs[i], character.only = T)
  }
}

pkgs <- c("gplots", "fgsea", "biomaRt", "clusterProfiler", 
          "GSEABase", "org.Hs.eg.db", "pcaMethods",
          "SingleCellExperiment", "batchelor", 
          "DelayedArray", "DelayedMatrixStats",
          "limma", "SummarizedExperiment", "SingleR",
          #"celldex",
          "progeny", "RcisTarget",
          "doMC", "doRNG", "DT", "visNetwork", "readr", "pheatmap", "tibble")

for(i in 1:length(pkgs)){
  if(!require(pkgs[i], character.only = T)){
    BiocManager::install(pkgs[i])
    require(pkgs[i], character.only = T)
  }else{
    require(pkgs[i], character.only = T)
  }
}

#####################
# Github packages
#####################

# library("devtools")

# usethis::browse_github_pat()
# usethis::edit_r_environ()
# GITHUB_PAT = "d8207153aef7b295cdf66eb1e1b2a2ed38b0ca18"
# R_MAX_VSIZE = 30Gb


# Scillus 
#devtools::install_github("xmc811/Scillus", ref = "development")
library("Scillus")

# Nebulosa for density plotting
#devtools::install_github("powellgenomicslab/Nebulosa")
library("Nebulosa")

# Volcano
#devtools::install_github('kevinblighe/EnhancedVolcano')
library("EnhancedVolcano")

# Clustifyr
#devtools::install_github("rnabioco/clustifyr")
library("clustifyr")

# ComplexHeatmap
#install_github("jokergoo/ComplexHeatmap")
library("ComplexHeatmap")

# scPred and requirements
#devtools::install_github("immunogenomics/harmony")
#devtools::install_github("powellgenomicslab/scPred")
library("scPred")
library("magrittr")
library("doParallel")

# Cytotrace
# Download source file from - https://cytotrace.stanford.edu/
#devtools::install_local("PATH/TO/DIRECTORY/CytoTRACE_0.3.3.tar.gz")
library("CytoTRACE")

# circlize
#devtools::install_github("jokergoo/circlize")
library("circlize")


######################
# Cell-cell interaction
######################

# celltalker
#devtools::install_github("arc85/celltalker")
library("celltalker")

# iTalk
#devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library("iTALK")

# singleCellHaystack
#remotes::install_github("alexisvdb/singleCellHaystack")
library("singleCellHaystack")


######################
# Seurat
######################

#devtools::install_github('satijalab/seurat-data')
#devtools::install_github('satijalab/seurat-wrappers')
library("SeuratWrappers")


######################
# Trajectory analysis
######################

# Trajectory analysis packages are loaded when required

# Velocity 
#devtools::install_github("velocyto-team/velocyto.R")

# Dyno
#devtools::install_github("dynverse/dyno")

# Monocle
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')

#.rs.restartR()

