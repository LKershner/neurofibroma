This repository contains scripts and output files from "Multiple Nf1 Schwann cell populations reprogram the plexiform neurofibroma tumor microenvironment" (add link to manuscript here once published).

These scripts are intended only to show the specific parameters used in this manuscript's analyses, and as such they may not be the most up to date or relevant packages and parameters if used for your specific study. Certain variable names are edited from their original form for clarity. (For example, the orig.ident for samples here versus the original RDS file may differ slightly.) Also, while I have made an effort to double check and remove unnecessary dependencies, some dependencies may not be necessary to run the example script. In these instances, the extra dependencies were likely used for additional downstream analyses not included in the shortened script.

FASTQ and h5 files can be found at GSE accession GSE181985

The following is a short description of the files in this repository:

SCRIPTS

cellRanger_count: script to run cellRanger; full documentation available at https://github.com/10XGenomics/cellranger

scrublet: script to run scrublet; full documentation available at https://github.com/AllonKleinLab/scrublet 

seurat-3.0.R: script to run basic Seurat v3 workflow; full documentation available at https://satijalab.org/seurat/index.html

Seurat_integration: script to run standard Seurat v3 integration (Stuart T, Butler A, et al., 2019) on 4 mouse conditions; full documentation available at https://satijalab.org/seurat/archive/v3.1/integration.html

cellphoneDB: script to run cellphoneDB; full documentation available at https://github.com/Teichlab/cellphonedb

SUPPLEMENTAL FILES

SupplementalInfoMouseSamples.xlsx: Supplemental table showing mouse sample information and cellRanger quality metrics output

Top50GenesCombined.txt: Top 50 marker genes for each cluster of the integrated mouse object as determined by Seurat's FindAllMarkers function

barcode_to_groups.txt: group assignments for each barcode in the mouse integrated Seurat object

barcode_to_clusters.txt: cluster assignments for each barcode in the mouse integrated Seurat object

networks_1.2-fold: PDFs of differentially expressed gene networks for each cluster in tumor vs 7 month age-matched control (CellHarmony), fold change > 1.2

