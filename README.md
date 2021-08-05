# HNSCC\_Descriptive

-   Human samples from HNSCC patients
-   samples 1 and 2 are unstimulated
-   samples 3 and 4 are stimulated
-   Patients in 1 match patients in 3 and sample 2 pairs with 4
-   Cells are HNSCC TILs sorted on Lymphocytes/Live/CD3+/CD4-/CD8+
-   Data acquired was transcript expression, ADT (antibody expression), and TCRA & B sequences
-   This analysis uses just Transcript and TCR data

# Data info:

-   Input data can be found in ~/Data
-   Gene expression data can be found ~/Data/10X_GEX_data
    - This folder contains the output from CellRanger (Barcodes, genes, matrix) 
    - Data is "cell by gene" matrix 
    - This serves as the raw input to seurat pipeline

-   TCR data can be found in ~/Data/TCR_data
    - two .csv files per sample (Filtered_contig_annotations and clonotypes)

-   All remaining input data are reference datasets, signature files, or other such input

-   Intermediate and processed files 
    - Can be found in ~/Exported_RDS_files
    - I have put many axillary ".rds" files into "~/Exported_RDS_files/Archive" in order to clean up the directory | this directory is not synced with github
    - seurat_combined.rds = all cells and genes following filtering and annotation
      - Contains all metadata, normalised data, imputed data, integrated US and Stim data
    - seurat_tcr_filt.rds = only cells which have a paired alpha/beta TCR associated
      - Analysis was seperated in this way to retain as many cells as possible for GEX analysis 
      - TCR analysis may also not include innate cell subsets and therefore dataset was split between GEX analysis and TCR analysis
      
      
      
# Scripts
- Knited ".html" outputs for each script can be found in "~"
- Main script for reading in CellRanger output (~/Data/10X_GEX_data) and normalising, integrating, QC etc, is found in R/01_HNSCC_Dataset_Preprocessing.Rmd
- Custom analysis functions and pipelines are located in ~/R/Functions or ~/R/Analysis
- Output Directory creation, package loading, etc is all found within setup.R and Load_packages.R scripts in ~/R
- scTCRseq analysis begins with ~/R/06_HNSCC_scTCRseq_dataset_Preprocessing.Rmd where raw CellRanger output is formatted


# Gitignore

-   all files in output/ will be ignored and not tracked in github
-   all files in output_tcr/ will be ignored and not tracked in github
-   Files in ~/Exported_RDS_files/Archive are ignored to reduce memory burden on collaborators who clone the repo

# Gitattribute

-   All .rds files in Exported\_RDS\_files are included in gitattribute file
-   Add any large files that should be tracked by github to the .gitattribute file

