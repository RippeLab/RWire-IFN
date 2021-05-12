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



plotUmapClusters <- function(ArchRProj, embedding=names(ArchRProj@embeddings)[1], reducedDims=names(ArchRProj@reducedDims)[1]){
    #### Clusters are not returned to ArchRProject!!
    
    ArchRProj <- addClusters(
    input = ArchRProj,
    reducedDims = reducedDims,
    method = "Seurat",
    name = paste0("Clusters0.2_", reducedDims),
    resolution = 0.2
    )
    ArchRProj <- addClusters(
    input = ArchRProj,
    reducedDims = reducedDims,
    method = "Seurat",
    name = paste0("Clusters0.8_", reducedDims),
    resolution = 0.8
    )
    ArchRProj <- addClusters(
    input = ArchRProj,
    reducedDims = reducedDims,
    method = "Seurat",
    name = paste0("Clusters1.4_", reducedDims),
    resolution = 1.4
    )
    
    umap_cluster_02 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = paste0("Clusters0.2_", reducedDims), embedding = embedding, size = 0.3) +
    ggtitle("Clustering (resolution 0.2)") +
    theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=28))
    
    umap_cluster_08 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = paste0("Clusters0.8_", reducedDims), embedding = embedding, size = 0.3) +
    ggtitle("Clustering (resolution 0.8)") +
    theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=28))
    
    umap_cluster_14 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = paste0("Clusters1.4_", reducedDims), embedding = embedding, size = 0.3) +
    ggtitle("Clustering (resolution 1.4)") +
    theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=28))
    

    return(list(umap_cluster_02, umap_cluster_08, umap_cluster_14))
}