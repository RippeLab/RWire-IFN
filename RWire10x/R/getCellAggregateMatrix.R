#' getBackgroundCoAccessibility
#' 
#' This function returns background co-accessibility scores for peaks of a given ArchRProject. The background is simulated by shuffling of features and cells.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses. Number of cells per aggregate.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`. Number of aggregates.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @keywords aggregates matrix archr rwire archrwire
#' @examples Coming soon.
#' @export


getCellAggregateMatrix <- function(
    ArchRProj = NULL, 
    reducedDims = "IterativeLSI",
    k = 100, 
    knnIteration = 500, 
    corCutOff = 0.75, 
    dimsToUse = 1:30, 
    overlapCutoff = 0.8, 
    scaleTo = 10^4
    ){
  
    set.seed(1)
    
    #Get Reduced Dims
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
    if(!is.null(NULL)){
        rD <- rD[NULL, ,drop=FALSE]
    }

    #Subsample
    idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

    knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

    #Determin Overlap
    keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

    #Keep Above Cutoff
    knnObj <- knnObj[keepKnn==0,]

    #Convert To Names List
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
        rownames(rD)[knnObj[x, ]]
    }) %>% SimpleList

    #Check Chromosomes
    chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
    chrj <- gtools::mixedsort(unique(paste0(seqnames(tiles))))
    stopifnot(identical(chri,chrj))  

    #Features
    featureDF <- mcols(tiles)
    featureDF$seqnames <- chri[1]

    #Group Matrix
    groupMat <- .getGroupMatrix(
        ArrowFiles = getArrowFiles(ArchRProj), 
        featureDF = featureDF, 
        groupList = knnObj, 
        useMatrix = "PeakMatrix",
        verbose = FALSE
    )

    #Scale
    #Peak Matrix ColSums
    cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
    gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))
    groupMat <- t(t(groupMat) / gS) * scaleTo

    return(groupMat)
}