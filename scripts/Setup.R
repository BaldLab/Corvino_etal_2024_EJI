# pipeline variables
start.time <- Sys.time()

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

# since moving script from local to github - I want to adjust work dir to be main github dir - therefore 
setwd("..")

working.dir <- getwd()

# For reproducibility
set.seed(42)

# Explicitly load key packages
library("Seurat")
library("SeuratDisk")
library("dplyr")

# Establish colour scheme
source("scripts/variables/Colour_scheme_variable.R", local = knitr::knit_global())

# Load custom functions
source("scripts/functions/scRNAseq_function.R", local = knitr::knit_global())
source("scripts/functions/annotate_seurat_heatmap_function.R", local = knitr::knit_global())
source("scripts/functions/Enhanced_volcano_custom_defaults_function.R", local = knitr::knit_global())

# TCR specific functions
source("scripts/functions/combineMeta_function.R", local = knitr::knit_global())
source("scripts/functions/Clonotype_distribution_function.R", local = knitr::knit_global())
source("scripts/functions/get_clonotypes_function.R", local = knitr::knit_global())


###############################
# create output directories
###############################

# Common directories 
if(!dir.exists("saves")){dir.create("saves", recursive = T)}

# scRNAseq directories
if(!dir.exists("results/scRNAseq")){dir.create("results/scRNAseq", recursive = T)}
if(!dir.exists("results/scRNAseq/figures")){dir.create("results/scRNAseq/figures", recursive = T)}
if(!dir.exists("results/scRNAseq/QC")){dir.create("results/scRNAseq/QC", recursive = T)}
if(!dir.exists("results/scRNAseq/tables")){dir.create("results/scRNAseq/tables", recursive = T)}

# TCRseq directories
if(!dir.exists("results/tcr")){dir.create("results/tcr", recursive = T)}
if(!dir.exists("results/tcr/figures")){dir.create("results/tcr/figures", recursive = T)}
if(!dir.exists("results/tcr/QC")){dir.create("results/tcr/QC", recursive = T)}
if(!dir.exists("results/tcr/tables")){dir.create("results/tcr/tables", recursive = T)}

