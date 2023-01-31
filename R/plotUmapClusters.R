#' plotUmapClusters
#' 
#' This function returns a list of UMAPs with different single cell clusterings.
#' @param ArchRProj ArchRProject with priorly calculated reduced dimensions and single-cell embedding.
#' @param embedding Optional name of calculated single-cell embedding (e. g. calculated with addUMAP function). Defaults to first embedding in ArchRProj@embeddings.
#' @param reducedDims Optional name of calculated reduced dimensions (e. g. calculated with addIterativeLSI function). Defaults to first reduced dimensions in ArchRProj@reducedDims.
#' @return List of three UMAPs with varying resolution of single cell clustering (0.2, 0.8, and 1.4).
#' @keywords umap cluster archr rwire archrwire
#' @examples Coming soon.
#' @export


plotUmapClusters <- function (ArchRProj, embedding = names(ArchRProj@embeddings)[1], reducedDims = names(ArchRProj@reducedDims)[1], resolutions = c(0.2, 0.8, 1.4)){
    #### Clusters are not returned to ArchRProject!!
    
    umaps <- lapply(resolutions, function(resolution){
        ArchRProj <- addClusters(input = ArchRProj, reducedDims = reducedDims, method = "Seurat", name = paste0("Clusters_resolution", resolution, "_", reducedDims), resolution = resolution);
        plotEmbedding(ArchRProj, colorBy = "cellColData", 
        name = paste0("Clusters_resolution", resolution, "_", reducedDims), embedding = embedding, 
        size = 0.3) + ggtitle(paste0("Clustering (resolution ", resolution, ")")) + 
        theme(plot.title = element_text(size = 30), axis.title = element_text(size = 28), 
            axis.text = element_text(size = 28), legend.title = element_text(size = 0), 
            legend.text = element_text(size = 28))
    })
    return(umaps)
}