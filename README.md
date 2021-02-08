# HNSCC\_Descriptive

-   Human samples from HNSCC patients
-   samples 1 and 2 are unstimulated
-   samples 3 and 4 are stimulated
-   Patients in 1 match patients in 3 and sample 2 pairs with 4
-   Cells are HNSCC TILs sorted on Lymphocytes/Live/CD3+/CD4-/CD8+
-   Data acquired was transcript expression, ADT (antibody expression), and TCRA & B sequences
-   This analysis uses just Transcript and TCR data

# Note:

-   Data/ seurat\_combined.rds = all cells and genes following filtering and annotations
-   Data/ seurat\_tcr.rds = only cells which have a paired ab TCR
-   Data/VDJ\_data/ is the raw dataframes output by 10x (Cellranger) pipeline with TCR info

# Gitignore

-   all files in output/ will be ignored and not tracked in github

# Gitattribute

-   All .rds files in Exported\_RDS\_files are included in gitattribute file
-   Add any large files that should be tracked by github to the .gitattribute file

# To Do: 

-   Overlay dataset with Cillo et al., HNSCC dataset (Any MAIT or GammaDelta cells)

-   Identify surface expressed genes in IFN cluster that can be used for sorting these cells

-   Generate IFN\_I cluster signature and validate with Cillo et al.,

-   Are IFN-I cells found in other disease settings, solid tumors, liquid tumors, autoimmune, viral, etc

# Outstanding tasks

-   Track clonotypes - for each cluster in US do clones increase/decrease and/or move between clusters

# scTCRseq tasks

-   Estimate epitope specificities with vdjdb tools
-   GLIPH analysis - note there is now a Gliph V2 algorithm
-   GLIPH visualisation
-   Logo TCR motif generation

# scRNAseq tasks

-   Transcription factor network analysis
-   Establish a signature for TypeI-IFN cluster
-   Apply TypeI-IFN signature to other TIL datasets and look for presence/absence of these cells across HNSCC and other cancer entities
-   Velocity analysis
