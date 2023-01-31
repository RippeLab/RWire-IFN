####################################################################
# Hidden Helper Utils for Arrow Files
####################################################################

.validArrow <- function(ArrowFile = NULL){
  o <- h5closeAll()
  if(h5read(ArrowFile,"Class")!="Arrow"){
    warning(
      "This file is not a valid ArrowFile, this most likely is a bug with previous function where the class was not added.\n",
      "To fix your ArrowFiles :\n",
      "\tlapply(getArrowFiles(ArchRProj), function(x) h5write(obj = 'Arrow', file = x, name = 'Class'))",
      "\nThis will be an error in future versions."
    )
  }
  o <- h5closeAll()
  return(ArrowFile)
}

.isProtectedArray <- function(matrixName = NULL, exclude = NULL){
  protectedArrays <- tolower(c("peakmatrix", "tilematrix", "genescorematrix"))
  if(!is.null(exclude)){
    protectedArrays <- protectedArrays[protectedArrays %ni% tolower(exclude)]
  }
  if(tolower(matrixName) %in% protectedArrays){
    stop(sprintf("Error %s cannot be used as this conflicts with another predefined matrix function!", matrixName))
  }
  matrixName
}

.availableArrays <- function(ArrowFiles = NULL, threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  availableArrays <- .safelapply(seq_along(ArrowFiles), function(x){
    groups <- h5ls(ArrowFiles[x]) %>% {.[.$group=="/" & .$otype=="H5I_GROUP","name"]}
    groups <- groups[!grepl("Fragments|Metadata", groups)]
    groups
  }, threads = threads) %>% Reduce("intersect", .)
  o <- h5closeAll()
  return(availableArrays)
}

.availableSeqnames <- function(ArrowFiles = NULL, subGroup = "Fragments", threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  seqList <- .safelapply(seq_along(ArrowFiles), function(x){
    seqnames <- h5ls(ArrowFiles[x]) %>% {.[.$group==paste0("/",subGroup),]$name}
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  }, threads = threads)
  if(!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],seqList[[1]]))))){
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}

.availableChr <- function(ArrowFiles = NULL, subGroup = "Fragments"){
  seqnames <- .availableSeqnames(ArrowFiles, subGroup)
  # if(getArchRChrPrefix()){
  #   seqnames <- seqnames[grep("chr", seqnames, ignore.case = TRUE)]
  # }
  if(length(seqnames) == 0){
    stop("No Chr Found in ArrowFiles!")
  }
  return(seqnames)
}

.availableCells <- function(ArrowFile = NULL, subGroup = NULL, passQC = TRUE){
  if(is.null(subGroup)){
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, "Metadata/CellNames")
    if(passQC){
      passQC <- tryCatch({
        h5read(ArrowFile, "Metadata/PassQC")
      }, error = function(x){
        rep(1, length(cellNames))
      })
      cellNames <- cellNames[which(passQC==1)]
    }
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }else{
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, paste0(subGroup, "/Info/CellNames"))
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  return(paste0(sampleName,"#",cellNames))
}

.sampleName <- function(ArrowFile = NULL){
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}

.summarizeArrowContent <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  h5DF <- h5ls(ArrowFile)

  #Re-Organize Content Info
  h5DF <- h5DF[-which(h5DF$group == "/"),]
  groups <- stringr::str_split(h5DF$group, pattern = "/", simplify=TRUE)[,2]
  groupList <- split(h5DF, groups)

  #Split Nested Lists
  groupList2 <- lapply(seq_along(groupList), function(x){
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]),"",groupDFx$group)
    if(all(groupx=="")){
      groupDFx
    }else{
      subDF <- groupDFx[-which(groupx == ""),]
      split(subDF, stringr::str_split(subDF$group, pattern = "/", simplify=TRUE)[,3])
    }
  })
  names(groupList2) <- names(groupList)


  o <- h5closeAll()

  return(groupList2)

}

.getMetadata <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- .summarizeArrowContent(ArrowFile)$Metadata

  #Which are same dimensions as cell names
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name=="CellNames"]),]

  #Load these into a S4 DataFrame
  md <- lapply(seq_len(nrow(arrowMD)), function(x){
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x],"/",arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)

  #Correct CellNames
  md$CellNames <- paste0(sampleName,"#",md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md)=="CellNames")]
  md <- md[,order(colnames(md))]

  o <- h5closeAll()

  return(md)
}

.getFeatureDF <- function(ArrowFiles = NULL, subGroup = "TileMatrix", threads = getArchRThreads()){

  threads <- min(threads, length(ArrowFiles))
  
  .helpFeatureDF <- function(ArrowFile = NULL, subGroup = NULL){
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup,"/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }

  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)

  if(length(ArrowFiles) > 1){
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- .safelapply(seq_along(ArrowFiles), function(x){
        fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
        identical(fdfx, fdf)
    }, threads = threads) %>% unlist %>% all
    if(!checkIdentical){
      stop("Error not all FeatureDF for asssay is the same!")
    }
  }

  #Re-Order for Split Check!
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% {lapply(seq_along(.), function(x) .[[x]])} %>% Reduce("c", .)
  fdf[newOrder,]
  
}

#####################################################################
# Dropping Group From Hdf5 File
#####################################################################
.createArrowGroup <- function(
  ArrowFile = NULL, 
  group = "GeneScoreMatrix", 
  force = FALSE, 
  verbose = FALSE,
  logFile = NULL
  ){
  
  ArrowInfo <- .summarizeArrowContent(ArrowFile)
  if(group == "Fragments"){ #This shouldnt happen but just in case
    .logMessage(".createArrowGroup : Cannot create Group over Fragments in Arrow!", logFile = logFile)
    stop("Cannot create Group over Fragments in Arrow!")
  }

  if(group %in% names(ArrowInfo)){
    #We Should Check How Big it is if it exists
    ArrowGroup <- ArrowInfo[[group]]
    ArrowGroup <- ArrowGroup[names(ArrowGroup) %ni% c("Info")]
    if(length(ArrowGroup) > 0){
      if(!force){
        .logMessage(".createArrowGroup : Arrow Group already exists! Set force = TRUE to continue!", logFile = logFile)
        stop("Arrow Group already exists! Set force = TRUE to continue!")
      }else{
        .logMessage(".createArrowGroup : Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!", logFile = logFile)
        if(verbose) message("Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!")
        o <- .dropGroupsFromArrow(ArrowFile = ArrowFile, dropGroups = group, verbose = verbose, logFile = logFile)
        tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
        invisible(return(0))
      }
    }
  }else{
    tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
    invisible(return(0))
  }

}

.dropGroupsFromArrow <- function(
  ArrowFile = NULL, 
  dropGroups = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL
  ){

  tstart <- Sys.time()

  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(ArrowFile)

  .logMessage(".dropGroupsFromArrow : Initializing Temp ArrowFile", logFile = logFile)

  #We need to transfer first
  outArrow <- .tempfile(fileext = ".arrow")
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")

  #1. Metadata First
  .logMessage(".dropGroupsFromArrow : Adding Metadata to Temp ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)

  mData <- ArrowInfo[[groupName]]
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name)
  }

  #2. Other Groups
  .logMessage(".dropGroupsFromArrow : Adding SubGroups to Temp ArrowFile", logFile = logFile)

  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% "Metadata"]
  if(!is.null(dropGroups)){
    groupsToTransfer <- groupsToTransfer[tolower(groupsToTransfer) %ni% tolower(dropGroups)]
  }

  for(k in seq_along(groupsToTransfer)){

    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)

    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }
    
    for(j in seq_along(seqOrder)){

      if(verbose) message(j, " ", appendLF = FALSE)

      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])
      o <- h5createGroup(outArrow, groupJ)

      #Sub mData
      mDataj <- mData[[seqOrder[j]]]

      #Transfer Components
      for(i in seq_len(nrow(mDataj))){
        h5name <- paste0(groupJ, "/", mDataj$name[i])
        .suppressAll(h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name, level = level))
      }

    }

    gc()
    
    if(verbose) message("")

  }

  .logMessage(".dropGroupsFromArrow : Move Temp ArrowFile to ArrowFile", logFile = logFile)

  rmf <- file.remove(ArrowFile)
  out <- .fileRename(from = outArrow, to = ArrowFile)
  
  .logDiffTime("Completed Dropping of Group(s)", tstart, logFile = logFile, verbose = verbose)

  ArrowFile

}

.copyArrows <- function(
  inArrows = NULL,
  outArrows = NULL,
  cellsKeep = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL,
  threads = 1
  ){

  stopifnot(length(inArrows) == length(outArrows))

  unlist(.safelapply(seq_along(inArrows), function(x){
    .copyArrowSingle(
      inArrow = inArrows[x],
      outArrow = outArrows[x],
      cellsKeep = cellsKeep,
      level = level,
      verbose = verbose,
      logFile = logFile
    )
  }, threads = threads))

}

.copyArrowSingle <- function(
  inArrow = NULL,
  outArrow = NULL,
  cellsKeep = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL
  ){

  tstart <- Sys.time()

  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(inArrow)
  sampleName <- .sampleName(inArrow)

  .logMessage(".copyArrow : Initializing Out ArrowFile", logFile = logFile)

  #We need to transfer first
  o <- .suppressAll(file.remove(outArrow))
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")

  #1. Metadata First
  .logMessage(".copyArrow : Adding Metadata to Out ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)

  mData <- ArrowInfo[[groupName]]
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    h5write(.h5read(inArrow, h5name), file = outArrow, name = h5name)
  }

  #2. scATAC-Fragments
  .logDiffTime(paste0("Transferring Fragments"), tstart, verbose = verbose, logFile = logFile)

  #Create Group
  groupName <- "Fragments"
  o <- h5createGroup(outArrow, groupName)
  
  #Sub Data
  mData <- ArrowInfo[[groupName]]
  
  #Get Order Of Sub Groups (Mostly Related to Seqnames)
  seqOrder <- sort(names(mData))
  if(any(grepl("chr", seqOrder))){
    seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
  }
  
  for(j in seq_along(seqOrder)){

    if(verbose) message(j, " ", appendLF = FALSE)

    #Create Group
    groupJ <- paste0(groupName, "/", seqOrder[j])
    o <- h5createGroup(outArrow, groupJ)

    #Sub mData
    mDataj <- mData[[seqOrder[j]]]

    #Read In Fragments
    RGLengths <- .h5read(inArrow, paste0(groupJ, "/RGLengths"))
    RGValues <- .h5read(inArrow, paste0(groupJ, "/RGValues"))
    RGRle <- Rle(paste0(sampleName, "#", RGValues), RGLengths)
    
    #Determine Which to Keep
    idx <- BiocGenerics::which(RGRle %bcin% cellsKeep)
    RGRle <- RGRle[idx]
    RGLengths <- RGRle@lengths

    #print(head(RGRle@values))
    RGValues <- stringr::str_split(RGRle@values, pattern = "#", simplify = TRUE)[,2]

    #Create Data Sets
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/Ranges"), storage.mode = "integer", dims = c(length(RGRle), 2), level = level))
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/RGLengths"), storage.mode = "integer", dims = c(length(RGRle), 1), level = level))
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/RGValues"), storage.mode = "character", dims = c(length(RGRle), 1), level = level, 
    #         size = max(nchar(RGValues) + 1)))

    #Write Barcodes
    o <- .suppressAll(h5write(RGLengths, file = outArrow, name = paste0(groupJ, "/RGLengths"), level = level))
    o <- .suppressAll(h5write(RGValues, file = outArrow, name = paste0(groupJ, "/RGValues"), level = level))

    #Write Ranges
    o <- .suppressAll(
      h5write(
        obj = .h5read(inArrow, paste0(groupJ, "/Ranges"))[idx, ], 
        file = outArrow, 
        name = paste0(groupJ, "/Ranges"), 
        level = level
      )
    )

  }
  
  if(verbose) message("")

  #3. Other Matrices
  .logMessage(".copyArrow : Adding SubMatrices to Out ArrowFile", logFile = logFile)
  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% c("Metadata", "Fragments")]

  for(k in seq_along(groupsToTransfer)){

    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)

    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }

    cellNames <- paste0(sampleName, "#", .h5read(inArrow, paste0(groupName, "/Info/CellNames")))
    featureDF <- .getFeatureDF(ArrowFile = inArrow, subGroup = groupName)

    for(j in seq_along(seqOrder)){

      if(verbose) message(j, " ", appendLF = FALSE)

      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])

      if(seqOrder[j] == "Info"){

        o <- h5createGroup(outArrow, groupJ)

        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        idxCL <- which(mDataj$dim == mDataj$dim[mDataj$name=="CellNames"])
        idxCL <- idxCL[mDataj$name[idxCL] %ni% "FeatureDF"]
        idxKeep <- which(cellNames %in% cellsKeep)

        #Transfer Components
        for(i in seq_len(nrow(mDataj))){

          h5name <- paste0(groupJ, "/", mDataj$name[i])
          
          if(i %in% idxCL){
            .suppressAll(h5write(.h5read(inArrow, h5name)[idxKeep], file = outArrow, name = h5name))
          }else{
            .suppressAll(h5write(.h5read(inArrow, h5name), file = outArrow, name = h5name))
          }

        }

      }else{

        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        addAnalysis <- mDataj[mDataj$name %ni% c("i", "jLengths", "jValues", "x"), "name"]

        mat <- .getMatFromArrow(
          ArrowFile = inArrow,
          featureDF = featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqOrder[j]),],
          useMatrix = groupName,
          cellNames = cellNames[cellNames %in% cellsKeep]
        )

        o <- .addMatToArrow(
          mat = mat, 
          ArrowFile = outArrow, 
          Group = paste0(groupName, "/", seqOrder[j]), 
          binarize = all(mat@x == 1),
          addColSums = "colSums" %in% addAnalysis,
          addRowSums = "rowSums" %in% addAnalysis,
          addRowMeans = "rowMeans" %in% addAnalysis,
          addRowVars = "rowVars" %in% addAnalysis,
          addRowVarsLog2 = "rowVarsLog2" %in% addAnalysis
        )

        rm(mat)

      }

    }

    gc()
    
    if(verbose) message("")

  }

  .logDiffTime("Completed Copying ArrowFile", tstart, logFile = logFile, verbose = verbose)

  outArrow

}
                                                                  
                                                                  
                                                                  ####################################################################
# Reading fragments from Arrow Files
####################################################################

#' Get the fragments from an ArrowFile 
#' 
#' This function retrieves the fragments from a given ArrowFile as a GRanges object.
#'
#' @param ArrowFile The path to the ArrowFile from which fragments should be obtained.
#' @param chr A name of a chromosome to be used to subset the fragments `GRanges` object to a specific chromsome if desired.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  cellNames = NULL, 
  verbose = TRUE,
  logFile = createLogFile("getFragmentsFromArrow")
  ){

  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = chr, name = "chr", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getFragmentsFromArrow Input-Parameters", logFile = logFile)

  ArrowFile <- .validArrow(ArrowFile)

  if(is.null(chr)){
    chr <- .availableSeqnames(ArrowFile, subGroup = "Fragments")
  }

  if(any(chr %ni% .availableSeqnames(ArrowFile, subGroup = "Fragments"))){
    stop("Error Chromosome not in ArrowFile!")
  }
  
  out <- lapply(seq_along(chr), function(x){
    .logDiffTime(sprintf("Reading Chr %s of %s", x, length(chr)), t1 = tstart, verbose = verbose, logFile = logFile)
    .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = chr[x], 
      out = "GRanges", 
      cellNames = cellNames, 
      method = "fast"
    )
  })


  .logDiffTime("Merging", tstart, t1 = tstart, verbose = verbose, logFile = logFile)

  out <- tryCatch({

    o <- .suppressAll(unlist(GRangesList(out, compress = FALSE)))

    if(.isGRList(o)){
      stop("Still a GRangesList")
    }

    o

  }, error = function(x){
    
    o <- c()
    
    for(i in seq_along(out)){
      if(!is.null(out[[i]])){
        if(i == 1){
          o <- out[[i]]
        }else{
          o <- c(o, out[[i]])
        }
      }
    }
    
    o

  })

  out

}

.getFragsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  out = "GRanges", 
  cellNames = NULL, 
  method = "fast"
  ){

  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }

  o <- h5closeAll()
  ArrowFile <- .validArrow(ArrowFile)
  
  avSeq <- .availableSeqnames(ArrowFile)
  if(chr %ni% avSeq){
    stop(paste0("Chromosome ", chr ," not in ArrowFile! Available Chromosomes are : ", paste0(avSeq, collapse=",")))
  }

  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)

  o <- h5closeAll()
  nFrags <- sum(.h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method))

  if(nFrags==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), RG = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }

  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }

  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){

      stop("None of input cellNames are in ArrowFile availableCells!")

    }else{

      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
          {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }

  }
  
  o <- h5closeAll()

  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }

  return(output)
}

####################################################################
# Reading Matrices/Arrays from Arrow Files
####################################################################

#' Get a data matrix stored in an ArchRProject
#' 
#' This function gets a given data matrix from an `ArchRProject`.
#'
#' @param ArchRProj An `ArchRProject` object to get data matrix from.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromProject <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromProject Input-Parameters", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)

  cellNames <- ArchRProj$cellNames


  seL <- .safelapply(seq_along(ArrowFiles), function(x){

    .logDiffTime(paste0("Reading ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
      t1 = tstart, verbose = FALSE, logFile = logFile)

    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]

    if(length(allCells) != 0){

      o <- getMatrixFromArrow(
        ArrowFile = ArrowFiles[x],
        useMatrix = useMatrix,
        useSeqnames = useSeqnames,
        cellNames = allCells, 
        ArchRProj = ArchRProj,
        verbose = FALSE,
        binarize = binarize,
        logFile = logFile
      )

      .logDiffTime(paste0("Completed ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
        t1 = tstart, verbose = FALSE, logFile = logFile)

      o

    }else{

      NULL
      
    }

  }, threads = threads) 

  #ColData
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, logFile = logFile)
  cD <- lapply(seq_along(seL), function(x){
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  #RowData
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, logFile = logFile)
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x){
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if(!rD){
    stop("Error with rowData being equal for every sample!")
  }

  #Assays
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i){
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", i, length(nAssays)), t1 = tstart, verbose = verbose, logFile = logFile)
    m <- lapply(seq_along(seL), function(j){
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, verbose = verbose, logFile = logFile)
  se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)  
  rm(seL)
  gc()

  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, logFile = logFile)

  se

}

#' Get a data matrix stored in an ArrowFile
#' 
#' This function gets a given data matrix from an individual ArrowFile.
#'
#' @param ArrowFile The path to an ArrowFile from which the selected data matrix should be obtained.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells from
#' the provided ArrowFile using `getCellNames()`.
#' @param ArchRProj An `ArchRProject` object to be used for getting additional information for cells in `cellColData`.
#' In some cases, data exists within the `ArchRProject` object that does not exist within the ArrowFiles. To access this data, you can
#' provide the `ArchRProject` object here.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromArrow <- function(
  ArrowFile = NULL, 
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  cellNames = NULL, 
  ArchRProj = NULL,
  verbose = TRUE,
  binarize = FALSE,
  logFile = createLogFile("getMatrixFromArrow")
  ){

  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromArrow Input-Parameters", logFile = logFile)

  ArrowFile <- .validArrow(ArrowFile)
  sampleName <- .sampleName(ArrowFile)

  seqnames <- .availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- .getFeatureDF(ArrowFile, subGroup = useMatrix)
  .logThis(featureDF, paste0("featureDF ", sampleName), logFile = logFile)

  if(!is.null(useSeqnames)){
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }

  if(length(seqnames) == 0){
    stop("No seqnames available!")
  }

  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]

  .logDiffTime(paste0("Getting ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
    t1 = tstart, verbose = verbose, logFile = logFile)

  if(!is.null(cellNames)){
    allCells <- .availableCells(ArrowFile = ArrowFile, subGroup = useMatrix)
    if(!all(cellNames %in% allCells)){
      stop("cellNames must all be within the ArrowFile!!!!")
    }
  }

  mat <- .getMatFromArrow(
    ArrowFile = ArrowFile, 
    featureDF = featureDF, 
    cellNames = cellNames, 
    useMatrix = useMatrix,
    binarize = binarize,
    useIndex = FALSE
  )
  .logThis(mat, paste0("mat ", sampleName), logFile = logFile)

  .logDiffTime(paste0("Organizing SE ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
    t1 = tstart, verbose = verbose, logFile = logFile)
  matrixClass <- h5read(ArrowFile, paste0(useMatrix, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    rownames(mat) <- paste0(featureDF$name)
    splitIdx <- split(seq_len(nrow(mat)), featureDF$seqnames)
    mat <- lapply(seq_along(splitIdx), function(x){
      mat[splitIdx[[x]], , drop = FALSE]
    }) %>% SimpleList
    names(mat) <- names(splitIdx)
    featureDF <- featureDF[!duplicated(paste0(featureDF$name)), ,drop = FALSE]
    featureDF <- featureDF[,which(colnames(featureDF) %ni% "seqnames"), drop=FALSE]
    rownames(featureDF) <- paste0(featureDF$name)
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix    
  }

  colData <- .getMetadata(ArrowFile)
  colData <- colData[colnames(mat[[1]]),,drop=FALSE]

  if(!is.null(ArchRProj)){
    projColData <- getCellColData(ArchRProj)[rownames(colData), ]
    colData <- cbind(colData, projColData[ ,colnames(projColData) %ni% colnames(colData)])
  }

  rowData <- tryCatch({
    makeGRangesFromDataFrame(featureDF, keep.extra.columns = TRUE)
  }, error = function(x){
    featureDF
  })

  se <- SummarizedExperiment(
    assays = mat,
    rowData = rowData,
    colData = colData
  )
  .logThis(se, paste0("se ", sampleName), logFile = logFile)

  se

}

.getMatFromArrow <- function(
  ArrowFile = NULL, 
  featureDF = NULL, 
  binarize = NULL, 
  cellNames = NULL,
  useMatrix = "TileMatrix", 
  useIndex = FALSE,
  threads = 1
  ){

  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }

  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }

  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))

  o <- h5closeAll()

  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix", "Sparse.Assays.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }

  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }

  seqnames <- unique(featureDF$seqnames)

  mat <- .safelapply(seq_along(seqnames), function(x){

    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx

    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
      )

    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]

    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }

    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)

    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)

    return(mat)

  }, threads = threads) %>% Reduce("rbind", .)

  o <- h5closeAll()

  colnames(mat) <- matColNames[idxCols]

  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL

  if(!is.null(cellNames)){
    mat <- mat[,cellNames,drop=FALSE]
  }

  return(mat)

}

####################################################################
# Helper read functioning
####################################################################
.getGroupMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  groupList = NULL,
  threads = 1, 
  useIndex = FALSE, 
  verbose = TRUE, 
  useMatrix = "TileMatrix",
  asSparse = FALSE,
  tstart = NULL
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################
  seqnames <- unique(featureDF$seqnames)
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  cellNames <- unlist(groupList, use.names = FALSE) ### UNIQUE here? doublet check QQQ

  allCellsList <- lapply(seq_along(ArrowFiles), function(x){
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if(length(allCells) != 0){
      allCells
    }else{
      NULL
    }
  })

  mat <- .safelapply(seq_along(seqnames), function(x){

    .logDiffTime(sprintf("Constructing Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)

    #Construct Matrix
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex), ]

    matChr <- matrix(0, nrow = nrow(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDFx)

    for(y in seq_along(ArrowFiles)){

      allCells <- allCellsList[[y]]
      
      if(!is.null(allCells)){

        maty <- .getMatFromArrow(
          ArrowFile = ArrowFiles[y], 
          useMatrix = useMatrix,
          featureDF = featureDFx, 
          cellNames = allCells, 
          useIndex = useIndex
        )

        for(z in seq_along(groupList)){

          #Check Cells In Group
          cellsGroupz <- groupList[[z]]
          idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)

          #If In Group RowSums
          if(length(idx) > 0){
            matChr[,z] <- matChr[,z] + Matrix::rowSums(maty[,idx,drop=FALSE])
          }

        }

        rm(maty)

      }
     

      if(y %% 20 == 0 | y %% length(ArrowFiles) == 0){
        gc()
      } 

    }

    if(asSparse){
      matChr <- as(matChr, "dgCMatrix")
    }

    .logDiffTime(sprintf("Finished Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    matChr

  }, threads = threads) %>% Reduce("rbind", .)

  mat <- mat[rownames(featureDF), , drop = FALSE]
  
  .logDiffTime("Successfully Created Group Matrix", tstart, verbose = verbose)

  gc()

  return(mat)
  
}

.getPartialMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  cellNames = NULL, 
  progress = TRUE, 
  threads = 1, 
  useMatrix = "TileMatrix",
  doSampleCells = FALSE, 
  sampledCellNames = NULL, 
  tmpPath = .tempfile(pattern = paste0("tmp-partial-mat")), 
  useIndex = FALSE,
  tstart = NULL,
  verbose = TRUE
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################

  mat <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(sprintf("Getting Partial Matrix %s of %s", x, length(ArrowFiles)), tstart, verbose = verbose)

    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]

    if(length(allCells) == 0){
      if(doSampleCells){
        return(list(mat = NULL, out = NULL))
      }else{
        return(NULL)
      }
    }

    o <- h5closeAll()
    matx <- .getMatFromArrow(
      ArrowFile = ArrowFiles[x], 
      featureDF = featureDF, 
      cellNames = allCells,
      useMatrix = useMatrix, 
      useIndex = useIndex
    )
   
    if(doSampleCells){

      #Save Temporary Matrix
      outx <- paste0(tmpPath, "-", .sampleName(ArrowFiles[x]), ".rds")
      saveRDS(matx, outx, compress = FALSE)     

      #Sample Matrix 
      matx <- matx[, which(colnames(matx) %in% sampledCellNames),drop = FALSE]
      
      return(list(mat = matx, out = outx))

    }else{
      
      return(matx)

    }

  }, threads = threads)

  gc()
  

  if(doSampleCells){

    matFiles <- lapply(mat, function(x) x[[2]]) %>% Reduce("c", .)
    mat <- lapply(mat, function(x) x[[1]]) %>% Reduce("cbind", .)
    mat <- mat[,sampledCellNames]

    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(list(mat = mat, matFiles = matFiles))

  }else{

    mat <- Reduce("cbind", mat)
    mat <- mat[,cellNames]
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(mat)

  }


}

########################################################################
# Compute Summary Statistics!
########################################################################

.getRowSums <- function(
  ArrowFiles = NULL,
  useMatrix = NULL,
  seqnames = NULL,
  verbose = TRUE,
  tstart = NULL,
  filter0 = FALSE,
  threads = 1,
  addInfo = FALSE
  ){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
    
  if(is.null(seqnames)){
    seqnames <- .availableSeqnames(ArrowFiles, useMatrix)
  }

  #Compute RowSums
  summaryDF <- .safelapply(seq_along(seqnames), function(x){
    o <- h5closeAll()
    for(y in seq_along(ArrowFiles)){
      if(y == 1){
        sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
      }else{
        sumy1 <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
        if(length(sumy1) != length(sumy)){
          stop("rowSums lengths do not match in ArrowFiles for a seqname!")
        }else{
          sumy <- sumy + sumy1
        }
      }
    }
    #Return Setup In Feature DF Format (seqnames, idx columns)
    DataFrame(seqnames = Rle(seqnames[x], lengths = length(sumy)), idx = seq_along(sumy), rowSums = as.vector(sumy))
  }, threads = threads) %>% Reduce("rbind", .)
  
  if(addInfo){
    featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, "_", featureDF$idx)
    rownames(summaryDF) <- paste0(summaryDF$seqnames, "_", summaryDF$idx)
    featureDF <- featureDF[rownames(summaryDF), , drop = FALSE]
    featureDF$rowSums <- summaryDF[rownames(featureDF), "rowSums"]
    summaryDF <- featureDF
    rownames(summaryDF) <- NULL
    remove(featureDF)
  }

  if(filter0){
    summaryDF <- summaryDF[which(summaryDF$rowSums > 0), ,drop = FALSE]
  }

  return(summaryDF)

}

.getRowVars <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  threads = 1
  ){
  
  .combineVariances <- function(dfMeans = NULL, dfVars = NULL, ns = NULL){

  #https://rdrr.io/cran/fishmethods/src/R/combinevar.R

  if(ncol(dfMeans) != ncol(dfVars) | ncol(dfMeans) != length(ns)){
    stop("Means Variances and Ns lengths not identical")
  }

  combinedMeans <- rowSums(t(t(dfMeans) * ns)) / sum(ns)
  summedVars <- rowSums(t(t(dfVars) * (ns - 1)) + t(t(dfMeans^2) * ns))
  combinedVars <- (summedVars - sum(ns)*combinedMeans^2)/(sum(ns)-1)

  data.frame(combinedVars = combinedVars, combinedMeans = combinedMeans)

  }

  featureDF <- .getFeatureDF(ArrowFiles, useMatrix)

  if(!is.null(seqnames)){
    featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames),]
  }

  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  fnames <- rownames(featureDF)

  featureDF <- S4Vectors::split(featureDF, as.character(featureDF$seqnames))

  ns <- lapply(seq_along(ArrowFiles), function(y){
    length(.availableCells(ArrowFiles[y], useMatrix))
  }) %>% unlist

  #Compute RowVars
  summaryDF <- .safelapply(seq_along(featureDF), function(x){
    
    o <- h5closeAll()
    seqx <- names(featureDF)[x]
    meanx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    varx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))

    for(y in seq_along(ArrowFiles)){
      meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeans"))
      varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVars"))
    }

  cbind(featureDF[[x]], DataFrame(.combineVariances(meanx, varx, ns)))

  }, threads = threads) %>% Reduce("rbind", .)

  summaryDF <- summaryDF[fnames, , drop = FALSE]
  
  return(summaryDF)

}

.getColSums <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  verbose = TRUE,
  tstart = NULL,
  threads = 1
  ){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Compute ColSums
  cS <- .safelapply(seq_along(seqnames), function(x){
    
    lapply(seq_along(ArrowFiles), function(y){
      
      o <- h5closeAll()
      cSy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/colSums"))
      rownames(cSy) <- .availableCells(ArrowFiles[y], useMatrix)
      cSy
      
    }) %>% Reduce("rbind", .)
    
  }, threads = threads) %>% Reduce("cbind", .) %>% rowSums

  .logDiffTime("Successfully Computed colSums", tstart, verbose = verbose)

  return(cS)

}

# h5read implementation for optimal reading
.h5read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
  ){

  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}

.getMatrixClass <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
  ){

  threads <- min(length(ArrowFiles), threads)

  matrixClass <- .safelapply(seq_along(ArrowFiles), function(i){
    h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Class"))
  }, threads = threads) %>% unlist %>% unique

  if(length(matrixClass) != 1){
    stop("Not all matrix classes are the same!")
  }

  matrixClass

}

.getMatrixUnits <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
  ){

  threads <- min(length(ArrowFiles), threads)

  matrixUnits <- .safelapply(seq_along(ArrowFiles), function(i){
    tryCatch({ #This handles backwards compatibility!
      h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Units"))
    }, error = function(x){
      "None"
    })
  }, threads = threads) %>% unlist %>% unique

  if(length(matrixUnits) != 1){
    stop("Not all matrix units are the same!")
  }

  matrixUnits

}