#' Get the peak co-accessibility from an ArchRProject
#' 
#' This function obtains co-accessibility data from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the co-accessibility signal as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export


getCoAccessibility <- function(
  ArchRProj = NULL, 
  corCutOff = 0, 
  resolution = 1, 
  returnLoops = TRUE
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){

    return(NULL)
  }

  
  if(is.null(metadata(ArchRProj@peakSet)$CoAccessibility)){
  
    return(NULL)
  
  }else{
   
    coA <- metadata(ArchRProj@peakSet)$CoAccessibility
    coA <- coA[coA$correlation >= corCutOff,,drop=FALSE]
      
    ### ADDED IL ###
    # Get peak count matrix and calculate ACR per peak
    PeakMatrix <- getMatrixFromProject(ArchRProj, "PeakMatrix")
    PeakMatrix@elementMetadata$numFrags <- Matrix::rowSums(PeakMatrix@assays@data$PeakMatrix)
    PeakMatrix@elementMetadata$numCells <- Matrix::rowSums(PeakMatrix@assays@data$PeakMatrix > 0)
    PeakMatrix@elementMetadata$ACR <- PeakMatrix@elementMetadata$numCells/ncol(PeakMatrix@assays@data$PeakMatrix)
    
    # Find corresponding peaks in co-accessibility matrix and peak count matrix
    PeakSet <- metadata(coA)$peakSet
    PeakMatrix@elementMetadata$chr <- as.character(seqnames(PeakSet)) # all(which(PeakMatrix@elementMetadata$idx == PeakSet$idx)) TRUE
    PeakMatrix@elementMetadata$peakSummit <- resize(PeakSet, 1, "center") %>% start(.)
    # queryHits and subjectHits in coA correspond to row number in peak set and peak count matrix
    minReg_rowID <- ifelse(PeakMatrix@elementMetadata$peakSummit[coA$queryHits] == pmin(PeakMatrix@elementMetadata$peakSummit[coA$queryHits], 
               PeakMatrix@elementMetadata$peakSummit[coA$subjectHits]), 
                           coA$queryHits, coA$subjectHits)
    maxReg_rowID <- ifelse(PeakMatrix@elementMetadata$peakSummit[coA$queryHits] == pmax(PeakMatrix@elementMetadata$peakSummit[coA$queryHits], 
               PeakMatrix@elementMetadata$peakSummit[coA$subjectHits]), 
                           coA$queryHits, coA$subjectHits)
    
    # Add information on numFrags and ACR to coA
    coA$reg1_numFrags <- PeakMatrix@elementMetadata$numFrags[minReg_rowID]
    coA$reg2_numFrags <- PeakMatrix@elementMetadata$numFrags[maxReg_rowID]
    coA$reg1_ACR <- PeakMatrix@elementMetadata$ACR[minReg_rowID]
    coA$reg2_ACR <- PeakMatrix@elementMetadata$ACR[maxReg_rowID]
    coA$ACR <- apply(coA[, c("reg1_ACR", "reg2_ACR")], 1, mean)
    ###
      
    if(returnLoops){
      
      peakSummits <- resize(metadata(coA)$peakSet, 1, "center")

      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
      }
    
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[coA[,1]]),
        start = summitTiles[coA[,1]],
        end = summitTiles[coA[,2]]
      )
      metadata(coA) <- list()
      mcols(loops) <- coA[,-c(1:3)]
      mcols(loops)$value <- coA$correlation
        
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))

      loops <- SimpleList(CoAccessibility = loops)

      return(loops)

    }else{

      return(coA)

    }

  }

}


.constructGR <- function(
  seqnames = NULL, 
  start = NULL, 
  end = NULL, 
  ignoreStrand = TRUE
  ){
  .validInput(input = seqnames, name = "seqnames", valid = c("character", "rleCharacter"))
  .validInput(input = start, name = "start", valid = c("integer"))
  .validInput(input = end, name = "end", valid = c("integer"))
  .validInput(input = ignoreStrand, name = "ignoreStrand", valid = c("boolean"))
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignoreStrand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-" 
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)
}