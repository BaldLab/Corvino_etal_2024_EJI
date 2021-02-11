# pipeline variables
start.time <- Sys.time()

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

# since moving script from local to github - I want to adjust work dir to be main github dir - therefore 
setwd("..")

working.dir <- getwd()

# Establish colour scheme
source("R/Colour_scheme_variable.R", local = knitr::knit_global())

# Load custom functions
source("R/Functions/scRNAseq_function.R", local = knitr::knit_global())
source("R/Functions/annotate_seurat_heatmap_function.R", local = knitr::knit_global())
source("R/Functions/Enhanced_volcano_custom_defaults_function.R", local = knitr::knit_global())
source("R/Functions/Modified_Celltalker_function.R", local = knitr::knit_global())

# TCR specific functions
source("R/Functions/getCircles_function.R", local = knitr::knit_global())
source("R/Functions/combineMeta_function.R", local = knitr::knit_global())
source("R/Functions/Clonotype_distribution_function.R", local = knitr::knit_global())
source("R/Functions/get_clonotypes_function.R", local = knitr::knit_global())


###############################
# create output directories
###############################

# Common directories 
if(!dir.exists("Exported_RDS_files")){dir.create("Exported_RDS_files", recursive = T)}

# scRNAseq directories
if(!dir.exists("output")){dir.create("output", recursive = T)}
if(!dir.exists("output/figures")){dir.create("output/figures", recursive = T)}
if(!dir.exists("output/tables")){dir.create("output/tables", recursive = T)}
if(!dir.exists("output/QC")){dir.create("output/QC", recursive = T)}

# scTCRseq directories
if(!dir.exists("output_tcr")){dir.create("output_tcr", recursive = T)}
if(!dir.exists("output_tcr/figures")){dir.create("output_tcr/figures", recursive = T)}
if(!dir.exists("output_tcr/tables")){dir.create("output_tcr/tables", recursive = T)}
if(!dir.exists("output_tcr/QC")){dir.create("output_tcr/QC", recursive = T)}
