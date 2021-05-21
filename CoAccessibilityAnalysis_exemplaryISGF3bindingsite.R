### R4.0.2
### Exemplary script for co-accessibility analysis of IFNb-treated mouse cells
### Isabelle Lander
### Department of Chromatin Networks, DKFZ
### 21.05.2020


########## Load libraries #########
library(magrittr)
library(ArchR)
library(qs)
library(gUtils)
library(Repitools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(RWire10x)


########## Set up analysis ##########
addArchRThreads(threads = 19)
addArchRGenome("mm10")

### Data available from GEO (GSE160764)
proj <- loadArchRProject("../01_ArchR_basics/ESC_IFNb")

### Exemplary ISGF3 binding site
bs <- GRanges(seqnames="chr1", IRanges(35274688, end=35275064))

### Generate peak base set over all chromosomes (otherwise bug in ArchR addCoAccessibility function for chr16 only analyses)
peaks_base <- GRanges(
        seqnames=Rle(paste0("chr", c(1:19, "X")), 1),
        ranges=IRanges(rep(10000, 20), end=rep(11000, 20)))

### Get chromosome information for mm10
chrInfo <- getChromInfoFromUCSC("mm10")

### Conditions
conditions <- c("ESC_IFNb0h", "ESC_IFNb6h")


########## Co-accessibility analysis for one exemplary binding site ##########
### Generate tile set (tile size 2 kb) for a 1 Mb window around the binding site
if (start(gr.mid(bs)) > chrInfo$size[chrInfo$chrom==seqnames(bs)@values]-500000){
    tiles <- tile(GRanges(seqnames=seqnames(bs), IRanges(start(gr.mid(bs))-499000, end=start(gr.mid(bs))+1000+floor((chrInfo$size[chrInfo$chrom==seqnames(bs)@values]-start(gr.mid(bs))-1000) / 2000)*2000)), n=floor((chrInfo$size[chrInfo$chrom==seqnames(bs)@values]-start(gr.mid(bs))-1000+500000) / 2000))[[1]]
} else if (start(gr.mid(bs)) < 500000) {
    stop("Code border region at chr start.")
} else {
    tiles <- tile(gr.mid(bs)+499000, n=499)[[1]]
}

### Add GC content for tiles
tiles$GC <- Repitools::gcContentCalc(tiles, organism=Mmusculus)

### Add peak set
proj <- addPeakSet(proj, peakSet = c(peaks_base, tiles), force = TRUE)

### Generate separate projects for conditions
proj_list <- sapply(conditions, function(x){proj_x <- saveArchRProject(ArchRProj = proj, outputDirectory = paste0("01_ArchRProjs/", x), load = TRUE); subsetCells(ArchRProj = proj_x, cellNames = proj$cellNames[which(proj$Sample_cluster==x)])})

### Add peak/tile set matrix
proj_list <- lapply(proj_list, addPeakMatrix)

### Add co-accessibility
proj_list <- lapply(proj_list, function(x){addCoAccessibility(x, reducedDims = "IterativeLSI", maxDist = 1e+06)})

### Get co-accessibility scores and accessible cell ratio (ACR) for all tiles
coacc_list <- lapply(proj_list, function(x){RWire10x::getCoAccessibility(x, corCutOff = 0, resolution = 1, returnLoops = TRUE)[[1]]})

### Get background co-accessibility by randomly shuffling accessibility in cells and tiles
background_list <- lapply(proj_list, function(x){RWire10x::getBackgroundCoAccessibility(x, reducedDims = "IterativeLSI_filtered_25000varFeatures", maxDist=1e+06)})

### Determine 99th percentile maximum of cell- and feature-shuffled background co-accessibility as threshold
background_cutoff <- max(quantile(background_list$featShuffle$correlation, seq(0.01,1,by=0.01), na.rm=T)[99],
                         quantile(background_list$cellShuffle$correlation, seq(0.01,1,by=0.01), na.rm=T)[99])

### Filter out links with background co-accessibility above 0.13 and accessible cell ratios of genomic tiles below 0.01
for (condition in conditions){
    coacc_list_filtered[[condition]] <- coacc_list[[condition]][coacc_list[[condition]]$correlation>background_cutoff]
    coacc_list_filtered[[condition]] <- coacc_list_filtered[[condition]][(coacc_list_filtered[[condition]]$accessCellRatio1>0.01 & coacc_list_filtered[[condition]]$accessCellRatio2>0.01)]
}

### Visualize co-accessibility scores around ISGF3 binding site
plotBrowserTrack(proj, groupBy = "Sample", region = gr.mid(bs)+500000, loops = coacc_list[1], normMethod = "nFrags")


########## Print time and success status ##########
print("DONE!")
print("End Time:")
Sys.time()