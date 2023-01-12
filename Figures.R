### R4.0.2
### Script for manuscript figures
### Isabelle Seufert
### Department of Chromatin Networks, DKFZ
### 12.01.2023

########## Load libraries #########
library(ArchR)
library(gUtils)
library(Repitools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)

devtools::install_github("RippeLab/RWire-IFN/RWire10x")
library(RWire10x)

### For installation:
# if (!requireNamespace(“devtools”, quietly = TRUE)) install.packages(“devtools”)
# if (!requireNamespace(“BiocManager”, quietly = TRUE)) install.packages(“BiocManager”)
# devtools::install_github(“GreenleafLab/ArchR”, ref=“master”, repos = BiocManager::repositories())
# library(ArchR)
# ArchR::installExtraPackages()

# devtools::install_github(‘mskilab/gUtils’)
# BiocManager::install(“Repitools”)
# BiocManager::install(“BSgenome.Mmusculus.UCSC.mm10")
# BiocManager::install("GenomicRanges")


########## Set up analysis ##########
addArchRThreads(threads = 19)
addArchRGenome("mm10")

### Data available from GEO (GSE160764)
proj <- loadArchRProject("../01_ArchR_basics/ESC_IFNb")

### Default ArchR workflow (see Methods section)
# Creation of Arrow files and ArchR projects, doublet detection
# Cell filtering criteria ESCs: log10(nFrags) >= 3.5; log10(nFrags) < 5; TSS ratio > 4; Blacklist ratio <= 0.0225
# Cell filtering criteria MEFs: log10(nFrags) >= 3.5; log10(nFrags) < 5; TSS ratio > 4; Blacklist ratio <= 0.0165; Doublet filtering ratio = 1.35

### Iterative LSI
proj_ESC <- addIterativeLSI(ArchRProj = proj_ESC,
                            useMatrix = "TileMatrix", 
                            name = "IterativeLSI", 
                            iterations = 2, 
                            clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10), 
                            varFeatures = 25000, 
                            dimsToUse = 1:30)
                            
proj_MEF <- addIterativeLSI(ArchRProj = proj_MEF,
                            useMatrix = "TileMatrix", 
                            name = "IterativeLSI", 
                            iterations = 2, 
                            clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10), 
                            varFeatures = 25000, 
                            dimsToUse = 2:12) # Relevant components determined from SVD deviation. LSI component 1 is highly correlated with sequencing depth.

### UMAP embedding
proj_ESC <- addUMAP(ArchRProj = proj_ESC, 
                    reducedDims = "IterativeLSI", 
                    name = "UMAP", 
                    nNeighbors = 30, 
                    minDist = 0.5, 
                    metric = "cosine")
                    
proj_MEF <- addUMAP(ArchRProj = proj_MEF, 
                    reducedDims = "IterativeLSI", 
                    name = "UMAP", 
                    nNeighbors = 30, 
                    minDist = 0.5, 
                    metric = "cosine")

### Clustering
proj_MEF <- addClusters(input = proj_MEF, 
                        reducedDims = "IterativeLSI", 
                        method = "Seurat", 
                        name = "Cluster", 
                        resolution = 0.2)
                        
### Integration with scRNA-seq data
proj_MEF <- addGeneIntegrationMatrix(ArchRProj = proj_MEF, 
                                     useMatrix = "GeneScoreMatrix",
                                     matrixName = "GeneIntegrationMatrix",
                                     reducedDims = "IterativeLSI",
                                     seRNA = proj_RNA,
                                     addToArrow = FALSE,
                                     groupRNA = "cell_type",
                                     nameCell = "predictedCell",
                                     nameGroup = "predictedGroup",
                                     nameScore = "predictedScore")

### Figures 4 A-C
plotUmaps(proj_ESC, embeddings = "UMAP", coloringLayer = "Sample")
plotUmaps(proj_MEF, embeddings = "UMAP", coloringLayer = "Sample")
plotUmaps(proj_MEF, embeddings = "UMAP", coloringLayer = "Cluster")
plotUmaps(proj_MEF, embeddings = "UMAP", coloringLayer = "predictedGroup_Un_cellType")


########## Print time and success status ##########
print("DONE!")
print("End Time:")
Sys.time()
