#' Add Peak Co-Accessibility to an ArchRProject
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param cellAggregation Choose cell aggregation method from "duplicated" or "unique". I no cell aggregation should be performed, set k to 1.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export

addCoAccessibility <- function (ArchRProj = NULL, reducedDims = "IterativeLSI", dimsToUse = 1:30, 
    scaleDims = NULL, corCutOff = 0.75, cellsToUse = NULL, cellAggregation = "duplicated", k = 100, 
    knnIteration = 500, overlapCutoff = 0.8, maxDist = 1e+05, 
    scaleTo = 10^4, log2Norm = TRUE, seed = 1, threads = getArchRThreads(), 
    verbose = TRUE, logFile = createLogFile("addCoAccessibility")) 
{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", 
        "null"))
    .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", 
        "null"))
    .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", 
        "null"))
    .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", 
        "null"))
    .validInput(input = k, name = "k", valid = c("integer"))
    .validInput(input = knnIteration, name = "knnIteration", 
        valid = c("integer"))
    .validInput(input = overlapCutoff, name = "overlapCutoff", 
        valid = c("numeric"))
    .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
    .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
    .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    tstart <- Sys.time()
    .startLogging(logFile = logFile)
    .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "addCoAccessibility Input-Parameters", logFile = logFile)
    set.seed(seed)
    peakSet <- getPeakSet(ArchRProj)
    
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, 
        corCutOff = corCutOff, dimsToUse = dimsToUse)
    if (!is.null(cellsToUse)) {
        rD <- rD[cellsToUse, , drop = FALSE]
    }
    
    ### Select most distant (randomly distributed) cells in low dimensional embedding
    if(cellAggregation == "unique"){
        # From: https://stackoverflow.com/questions/22152482/choose-n-most-distant-points-in-r
        bestavg <- 0
        bestSet <- NA
        for (i in 1:nrow(rD)){
            subset <- rD[sample(1:nrow(rD),knnIteration),]
            avg <- mean(dist(subset))
            if (avg > bestavg) {
                bestavg <- avg
                bestSet <- subset
            }
        }
        idx <- match(rownames(bestSet), rownames(rD))
    } else if (cellAggregation == "duplicated"){
        idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= 
        knnIteration)
    }
    
    .logDiffTime(main = "Computing KNN", t1 = tstart, verbose = verbose, 
        logFile = logFile)
    
    ### Select every cell only once
    if(cellAggregation == "unique"){
        query1 <- rD[idx[1], ] %>% as.matrix(.) %>% t(.)
        rownames(query1) <- rownames(rD)[1]
        y <- .computeKNN(data = rD[-idx,], query = query1, k = k-1)
        names <- rownames(rD[-idx,])[as.integer(y)]
        knnObj <- c(idx[1], match(names, rownames(rD))) %>% as.matrix(.) %>% t(.)
        
        for (i in 2:knnIteration){
            query1 <- rD[idx[i], ] %>% as.matrix(.) %>% t(.)
            rownames(query1) <- rownames(rD)[i]
            y <- .computeKNN(data = rD[-unique(c(idx, knnObj %>% as.integer())),], query = query1, k = k-1)
            names <- rownames(rD[-unique(c(idx, knnObj %>% as.integer())),])[as.integer(y)]
            knnObj <- rbind(knnObj, c(idx[i], match(names, rownames(rD))))
        }
    } else if (cellAggregation == "duplicated"){
        knnObj <- .computeKNN(data = rD, query = rD[idx, ], k = k)
    }
    
    .logDiffTime(main = "Identifying Non-Overlapping KNN pairs", 
        t1 = tstart, verbose = verbose, logFile = logFile)
    keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * 
        k))
    knnObj <- knnObj[keepKnn == 0, ]
    .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), 
        t1 = tstart, verbose = verbose, logFile = logFile)
    #### convert integer to matrix:array
    if (k == 1){
        knnObj <- matrix(knnObj, nrow = knnIteration)
    }
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x) {
        rownames(rD)[knnObj[x, ]]
    }) %>% SimpleList
    chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), 
        subGroup = "PeakMatrix"))
    chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
    stopifnot(identical(chri, chrj))
    peakSummits <- resize(peakSet, 1, "center")
    peakWindows <- resize(peakSummits, maxDist, "center")
    o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
    o <- o[o[, 1] != o[, 2], ]
    o$seqnames <- seqnames(peakSet)[o[, 1]]
    o$idx1 <- peakSet$idx[o[, 1]]
    o$idx2 <- peakSet$idx[o[, 2]]
    o$correlation <- -999.999
    o$Variability1 <- 0
    o$Variability2 <- 0
    cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, 
        useMatrix = "PeakMatrix")
    gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], 
        na.rm = TRUE)))
    ### Add pseudo-count to gS for non-accessible cell aggregates in all regions of interest
    gS <- gS + 1
    ###
    for (x in seq_along(chri)) {
        .logDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", 
            chri[x], x, length(chri)), t1 = tstart, verbose = verbose, 
            logFile = logFile)
        featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == 
            chri[x]), ]
        featureDF$seqnames <- chri[x]
        groupMat <- .getGroupMatrix(ArrowFiles = getArrowFiles(ArchRProj), 
            featureDF = featureDF, groupList = knnObj, useMatrix = "PeakMatrix", 
            threads = threads, verbose = FALSE)
        groupMat <- t(t(groupMat)/gS) * scaleTo
        if (log2Norm) {
            groupMat <- log2(groupMat + 1)
        }
        idx <- BiocGenerics::which(o$seqnames == chri[x])
        corVals <- rowCorCpp(idxX = o[idx, ]$idx1, idxY = o[idx, 
            ]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
        .logThis(head(corVals), paste0("SubsetCorVals-", x), 
            logFile = logFile)
        rowVars <- as.numeric(matrixStats::rowVars(groupMat))
        o[idx, ]$correlation <- as.numeric(corVals)
        o[idx, ]$Variability1 <- rowVars[o[idx, ]$idx1]
        o[idx, ]$Variability2 <- rowVars[o[idx, ]$idx2]
        .logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
        .logThis(o[idx, ], paste0("SubsetCoA-", x), logFile = logFile)
    }
    o$idx1 <- NULL
    o$idx2 <- NULL
    o <- o[!is.na(o$correlation), ]
    o$TStat <- (o$correlation/sqrt((pmax(1 - o$correlation^2, 
        1e-17, na.rm = TRUE))/(length(knnObj) - 2)))
    o$Pval <- 2 * pt(-abs(o$TStat), length(knnObj) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")
    o$VarQuantile1 <- .getQuantiles(o$Variability1)
    o$VarQuantile2 <- .getQuantiles(o$Variability2)
    mcols(peakSet) <- NULL
    o@metadata$peakSet <- peakSet
    metadata(ArchRProj@peakSet)$CoAccessibility <- o
    .endLogging(logFile = logFile)
    ArchRProj
}