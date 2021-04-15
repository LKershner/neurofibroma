#for complete documentation see https://satijalab.org/seurat/archive/v3.1/integration.html

#note that this analysis was performed using the standard workflow. 

#note that this script is different from that used to integrate human and mouse, using mouse as a reference, which is uploaded separately

#Lastly, I have commented out when the script has been modified for clarity from its exact original form. 
#Some of the original directory and variable names corresponded to sample designations given by the gene expression and sequencing cores. 
#These names often include extra information (e.g., sample submission date, sequencing date) and thus can be lengthy and confusing. 
#Here, variable and directory names are shortened and simplified to the note only the sample group with no extraneous information.

#before running this, will need to take following steps in terminal
#module load python/2.7.13
#module load rstudio/1.1.463
#rstudio

#in rstudio (interactively)

library(Seurat)
library(cowplot)
library(reticulate) #necessary if using python tools
library(devtools)
library(dplyr)
library(Matrix)
library(base)
library(methods)
library(utils)
library(stats)
library(gdata)
library(graphics)
library(grDevices)
library(ggplot2)
library(ggridges) #not necessary unless making ridge plots

import("numpy") #python numpy
import("umap")  #python umap

#load in Seurat objects for integration

setwd("/data/tumorDirectory") #modified for clarity
tumor <- readRDS("tumor.rds") #modified for clarity
setwd("/data/pretumorDirectory") #modified for clarity
pretumor <- readRDS("pretumor.rds") #modified for clarity
 setwd("/data/twomonthcontrol") #modified for clarity
twomoctrl <- readRDS("twomonthcontrol.rds") #modified for clarity
setwd("/data/sevenmonthcontrol") #modified for clarity
sevenmoctrl <-readRDS("sevenmonthcontrol.rds") #modified for clarity

#add samples to list

samples.list = list(pretumor, tumor, twomoctrl, sevenmoctrl)

#use standard integration normalization process with 2000 features/anchors

samples.list <- lapply(X = samples.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#find anchors and integrate data

samples.anchors <- FindIntegrationAnchors(object.list = samples.list, dims = 1:30)
samples.combined <- IntegrateData(anchorset = samples.anchors, dims = 1:30)

#set default assay to integrated for downstream steps

DefaultAssay(samples.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
samples.combined <- ScaleData(samples.combined, verbose = TRUE)
samples.combined <- RunPCA(samples.combined, npcs = 30, verbose = TRUE)

# t-SNE and Clustering
samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:30)
samples.combined <- FindNeighbors(samples.combined, reduction = "pca", dims = 1:30)
samples.combined <- FindClusters(samples.combined, resolution = 0.5)

# Visualization (note this is the standard workflow, actual script for final manuscript plots may vary in their paramters)
p1 <- DimPlot(samples.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(samples.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(samples.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

## Finding differentially expressed features (cluster biomarkers)

samples.markers <- FindAllMarkers(object = samples, only.pos = TRUE, min.pct=0.25, logfc.threshold = 0.25)

top50 <- samples.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

write.table(top50,file="Top50GenesCombined.txt",sep="\t",col.names= NA)
write.table(x = Idents(object = samples.combined),"Cluster.txt",sep="\t",col.names= NA) #this table is available as an output file in this GitHub
dev.off()
