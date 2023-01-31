### R4.0.2
### Exemplary script for co-accessibility analysis of IFNb-treated mouse cells
### Isabelle Seufert
### Department of Chromatin Networks, DKFZ
### 27.05.2020


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
proj <- loadArchRProject("ArchRProject")

### Default ArchR workflow (see Methods section)
# Creation of Arrow files and ArchR projects, doublet detection
# Cell filtering criteria ESCs: log10(nFrags) >= 3.5; log10(nFrags) < 5; TSS ratio > 4; Blacklist ratio <= 0.0225
# Cell filtering criteria MEFs: log10(nFrags) >= 3.5; log10(nFrags) < 5; TSS ratio > 4; Blacklist ratio <= 0.0165; Doublet filtering ratio = 1.35


########## Create peak set ##########
### Call peaks
proj <- addGroupCoverages(proj, groupBy = "Sample")
proj <- addReproduciblePeakSet(proj, groupBy = "Sample", pathToMacs2 = "/path/to/macs2", 
                               extendSummits = 1000, reproducibility = "1")
peaks <- getPeakSet(proj)

### Add 2 kb regions around STAT1/2 bound sites (Table_S4) and TSSs of ISGs (Table_S5) to peak set
peaks <- c(peaks, gUtils::gr.mid(STAT_sites)+1000, gUtils::gr.mid(ISG_TSSs)+1000)

### Annotate peaks overlapping with STAT1/2 bound sites
peaks$STAT <- FALSE
peaks$STAT[IRanges::findOverlaps(peaks, STAT_sites) %>% IRanges::queryHits(.)] <- TRUE

### Add peak set
proj <- addPeakSet(proj, peakSet = peaks, force = TRUE)
proj <- addPeakMatrix(proj)


########## Define analysis conditions and compensate biases in nCells and nFrags ##########
### Conduct analysis for individual cell types and conditions, seperately
samples <- c("ESC_IFNb0h", "ESC_IFNb6h", 
                "MEF_IFNb0h_mesenchymal", "MEF_IFNb1h_mesenchymal", "MEF_IFNb6h_mesenchymal",
                "MEF_IFNb0h_epithelial", "MEF_IFNb1h_epithelial", "MEF_IFNb6h_epithelial")

### Define reference sample and nCells
n <- 2700
ref <- "MEF_IFNb0h_epithelial"

### Define cells to keep
keepCells <- sample(proj$cellNames[proj$Sample == ref], n)
ref_nFrags <- proj$nFrags[proj$cellNames == keepCells[[ref]]]
keepCells <- lapply(samples[-ref], function(sample){query_nFrags <- proj$nFrags[proj$Sample == sample];
                                                    names(query_nFrags) <- proj$cellNames[proj$Sample == sample];
                                                    keep <- names(query_nFrags)[ Closest(x=query_nFrags, ref_nFrags[1], which = TRUE)[1] ];
                                                    for (a in ref_nFrags[-1]){
                                                            lookup <- query_nFrags[-keep]; names(lookup) <- names(query_nFrags)[-keep];
                                                            keep <- c(keep, names(lookup)[ Closest(x=lookup, a, which = TRUE)[1] ]);
                                                    }
                                                    return(keep) }) %>% unlist(.) %>% c(keepCells, .)

### Filter cells
proj <- subsetCells(ArchRProj = proj, cellNames = keepCells)


########## Compute co-accessibility ##########
### Split up ArchRProjects for samples
proj_list <- sapply(samples, function(sample){proj_x <- saveArchRProject(ArchRProj = proj, outputDirectory = paste0("ArchhRProj_", sample), load = TRUE); 
                                               subsetCells(ArchRProj = proj_x, cellNames = proj_x$cellNames[which(proj_x$Sample == sample)]) })

### Get co-accessibility for unaggregated cells
proj_list <- lapply(proj_list, function(x){RWire10x::addCoAccessibility(x, reducedDims = "IterativeLSI", maxDist = 2e+06, 
                                                                        cellAggregation = "duplicated", knnIteration = n, k=1) })
coacc_list <- lapply(proj_list, function(x){RWire10x::getCoAccessibility(x, corCutOff = -1, resolution = 1, returnLoops = TRUE)[[1]] })

### Get background co-accessibility
background_list <- lapply(proj_list, function(x){RWire10x::getBackgroundCoAccessibility(x, reducedDims = "IterativeLSI", maxDist=2e+06, 
                                                                                        cellAggregation = "duplicated", knnIteration = n, k=1) })

### Determine 99th percentile maximum of cell- and feature-shuffled background co-accessibility
background_cutoff <- lapply(background_list, function(x){max(quantile(x$featShuffle$correlation, seq(0.01,1,by=0.01), na.rm=T)[99],
                                                             quantile(x$cellShuffle$correlation, seq(0.01,1,by=0.01), na.rm=T)[99]) }) %>% max(.) 

### Filter links with correlation coefficient above background cutoff, p-value below 0.01, and overlap with STAT1/2 bound site
coacc_list_filtered <- lapply(samples, function(sample){coacc_filtered <- coacc_list[[sample]][ coacc_list[[sample]]$correlation > background_cutoff ];
                                                        coacc_filtered <- coacc_filtered[ coacc_filtered$Pval < 0.01 ];
                                                        coacc_filtered <- coacc_filtered[ coacc_filtered$STAT ];
                                                        return(coacc_filtered) })

### Visualize co-accessibility scores around exemplary STAT1/2 bound site
bs <- GRanges(seqnames="chr1", IRanges(start=35274688, end=35275064))
plotBrowserTrack(proj, groupBy = "Sample", region = gr.mid(bs)+500000, loops = coacc_list_filtered, normMethod = "nFrags")


########## Print time and success status ##########
print("DONE!")
print("End Time:")
Sys.time()
