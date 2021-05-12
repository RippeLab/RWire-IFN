#' plotCoAccessibilityMap
#' 
#' This function returns co-accessibility map of specified genomic region.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param region GRanges object of length 1 with region of interest for which co-accessibility should be calculated.
#' @param resolution Tile size for which co-accessibility should be calculated. Defaults to 1 kb.
#' @param ArchRGenome Used genome and gene annotation.
#' @return Co-accessibility map.
#' @keywords co-accessibility map coAccessibilityMap archr rwire archrwire
#' @examples Coming soon.
#' @export


plotCoAccessibilityMap <- function(
    ArchRProj = NULL,
    reducedDims = names(ArchRProj@reducedDims)[1],
    region = NULL,
    resolution = 1000,
    ArchRGenome = getArchRGenome()
    ){
    
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    .validInput(input = region, name = "region", valid = c("GRanges"))
    .validInput(input = resolution, name = "resolution", valid = c("numeric"))
    .validInput(input = ArchRGenome, name = "ArchRGenome", valid = c("character"))
    
    ### Create temporary ArchRProject
    dir.create('tmp', showWarnings = FALSE)
    ArchRProj <- saveArchRProject(ArchRProj = ArchRProj, outputDirectory = 'tmp/ArchRProj_plotCoAccessibilityMap', load = TRUE)
    
    ### Generate tile set
    if (ArchRGenome == 'mm10'){ # base tile set on all chr (otherwise bug for chr16) 
        tiles_base <- GRanges(seqnames=Rle(paste0("chr", c(1:19, "X")), 1), ranges=IRanges(rep(10000,20), end=rep(11000,20)))
    } else if (ArchRGenome == 'Hg38'){
        tiles_base <- GRanges(seqnames=Rle(paste0("chr", c(1:22, "X")), 1), ranges=IRanges(rep(10000,23), end=rep(11000,23)))
    } else {
        stop('Genome built not included yet.')
    }

    numTiles <- (end(region)-start(region))/resolution
    tiles_region <- tile(region, n=numTiles)[[1]]
    tiles <- c(tiles_base, tiles_region)

    if (ArchRGenome == 'mm10'){ # GC content for addCoAccessibility function (does not work without GC content for normalization)
        library(BSgenome.Mmusculus.UCSC.mm10)
        tiles$GC <- Repitools::gcContentCalc(tiles, organism=Mmusculus)
    } else if (ArchRGenome == 'Hg38'){
        library(BSgenome.Hsapiens.UCSC.hg38)
        tiles$GC <- Repitools::gcContentCalc(tiles, organism=Hsapiens)
    }
    
    ### Calculate co-accessibility scores in region
    ArchRProj <- addPeakSet(ArchRProj, peakSet = tiles, force=TRUE)
    ArchRProj <- addPeakMatrix(ArchRProj)
    
    ArchRProj <- addCoAccessibility(ArchRProj, reducedDims = reducedDims, maxDist = 2*(end(region)-start(region)))
    coacc_tiles <- ArchRWire::getCoAccessibility(ArchRProj, corCutOff = 0, resolution = 1, returnLoops = TRUE)[[1]]
      
    ### Co-accessibility matrix visualization with RWire    
    df_background <- expand.grid(tiles_region %>% gUtils::gr.mid(.) %>% start(.), tiles_region %>% gUtils::gr.mid(.) %>% start(.))
    colnames(df_background) <- c('start', 'end')
    df_background$correlation <- rep(0, nrow(df_background))
    
    df_coacc_tiles <- data.frame(start=start(coacc_tiles), end=end(coacc_tiles), correlation=coacc_tiles$value)
    df_coacc_tiles <- rbind(df_coacc_tiles, df_background)
    df_coacc_tiles <- df_coacc_tiles[!duplicated(df_coacc_tiles[, c('start', 'end')]),]
    
    mtx_coacc_tiles <- xtabs(correlation ~ start + end, data = df_coacc_tiles, sparse = TRUE)
    
    loop_max <- max(df_coacc_tiles$correlation)
    palette <- rgb(1, 1 - (0:255)/255, 1 - (0:255)/255)
    
    fields::image.plot(as.matrix(mtx_coacc_tiles), zlim=c(0, loop_max), axes = FALSE, main=as.character(region), col=palette)
    abline(a=0, b=1, col='grey')
    axis(side=1, at=0, labels=start(region), tick=FALSE)
    axis(side=3, at=1, labels=end(tiles[length(tiles)]), tick=FALSE)
    
    p <- grDevices::recordPlot()
    
    return(list(CoAccessibilityMap = p, CoAccessibility = coacc_tiles))
}