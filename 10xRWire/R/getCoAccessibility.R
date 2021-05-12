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
      
    ### Add information on percent of accessible cells per linked peak
    PeakMatrix <- getMatrixFromProject(ArchRProj, 'PeakMatrix')  
    PeakMatrix@elementMetadata$numFrags <- Matrix::rowSums(PeakMatrix@assays@data$PeakMatrix)
    PeakMatrix@elementMetadata$numCells <- Matrix::rowSums(PeakMatrix@assays@data$PeakMatrix>0)
    PeakMatrix@elementMetadata$accessCellRatio <- PeakMatrix@elementMetadata$numCells/ncol(PeakMatrix@assays@data$PeakMatrix)
    
    minReg_idx <- getPeakSet(ArchRProj)$idx[pmin(coA$queryHits, coA$subjectHits)]
    maxReg_idx <- getPeakSet(ArchRProj)$idx[pmax(coA$queryHits, coA$subjectHits)]
    coA$numFrags1 <- PeakMatrix@elementMetadata$numFrags[match(minReg_idx, PeakMatrix@elementMetadata$idx)]
    coA$numFrags2 <- PeakMatrix@elementMetadata$numFrags[match(maxReg_idx, PeakMatrix@elementMetadata$idx)]
    coA$accessCellRatio1 <- PeakMatrix@elementMetadata$accessCellRatio[match(minReg_idx, PeakMatrix@elementMetadata$idx)]
    coA$accessCellRatio2 <- PeakMatrix@elementMetadata$accessCellRatio[match(maxReg_idx, PeakMatrix@elementMetadata$idx)]

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