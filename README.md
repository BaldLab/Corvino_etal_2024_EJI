# HNSCC_Descriptive

- Human samples from HNSCC patients
- samples 1 and 2 are unstimulated 
- samples 3 and 4 are stimulated 
- Patients in 1 match patients in 3 and sample 2 pairs with 4
- Cells are HNSCC TILs sorted on Lymphocytes/Live/CD3+/CD4-/CD8+
- Data acquired was transcript expression, ADT (antibody expression), and TCRA & B sequences
- This analysis uses just Transcript and TCR data

# Note: 
- Dir/file paths are not written in a github compatible structure - these file paths will need to be amended 
- Data/ seurat_combined.rds = all cells and genes following filtering and annotations 
- Data/ seurat_tcr.rds = only cells which have a paired ab TCR 
- Data/VDJ_data/ is the raw dataframes output by 10x (Cellranger) pipeline with TCR info
- R/ HNSC_scRNAseq_Descriptive_paper_V*.Rmd = full script with entire analysis so far - including many work in progress sections of code 
- R/ HNSCC_scTCRseq_analysis.Rmd = just TCRseq relevant code chunks taken from the full script


# Outstanding tasks
- Track clonotypes - for each cluster in US do clones increase/decrease and/or move between clusters 

# scTCRseq tasks
- Estimate epitope specificities with vdjdb tools
- GLIPH analysis - note there is now a Gliph V2 algorithm 
- GLIPH visualisation 
- Logo TCR motif generation 

# scRNAseq tasks
- Transcription factor network analysis
- Establish a signature for TypeI-IFN cluster
- Apply TypeI-IFN signature to other TIL datasets and look for presence/absence of these cells across HNSCC and other cancer entities
- Velocity analysis

