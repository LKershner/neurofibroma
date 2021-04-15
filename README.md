This repository contains scripts and output files from "Transcriptional reprogramming of Schwann cell populations and the tumor microenvironment during plexiform neurofibroma formation" (add link to manuscript here once published).

These scripts are intended only to show the specific parameters used in this manuscript's analyses, and as such they may not be the most up to date or relevant packages and parameters if used for your specific study. Certain variable names are edited from their original form for clarity. (For example, the orig.ident for samples here versus the original RDS file may differ slightly.) Also, while I have made an effort to double check and remove unnecessary dependencies, some dependencies may not be necessary to run the example script. In these instances, the extra dependencies were likely used for additional downstream analyses not included in the shortened script.

The following is a short description of the files in this repository:

SCRIPTS

cellRanger_count: script to run cellRanger; full documentation available at https://github.com/10XGenomics/cellranger
scrublet: script to run scrublet; full documentation available at https://github.com/AllonKleinLab/scrublet 
cellphoneDB: script to run cellphoneDB; full documentation available at https://github.com/Teichlab/cellphonedb

SUPPLEMENTAL FILES

SupplementalInfoMouseSamples.xlsx: Supplemental table (update) showing mouse sample information and cellRanger quality metrics output
Top50GenesCombined.txt: Top 50 marker genes for each cluster of the integrated mouse object as determined by Seurat's FindAllMarkers function


