#' plotUmaps
#' 
#' This function returns a list of UMAPs of different single-cell embeddings.
#' @param ArchRProj ArchRProject with priorly calculated reduced dimensions and single-cell embedding.
#' @param embeddings Optional character vector of calculated single-cell embeddings (e. g. calculated with addUMAP function). Defaults to all ArchRProj@embeddings.
#' @param coloringLayer Optional name of cellColData column to use as coloring layer of single cells. Defaults to 'Sample' column.
#' @param coloringPalette Optional coloring palette for split levels of single cells.
#' @param size Point size of cells in UMAP embedding.
#' @return List of UMAPs with different embeddings of single cells.
#' @keywords umap embeddings archr rwire archrwire
#' @examples Coming soon.
#' @export



plotUmaps <- function(ArchRProj, embeddings=names(ArchRProj@embeddings), coloringLayer='Sample', coloringPalette=NULL, size=0.3){
    k <- unique(ArchRProj$Sample) %>% length(.)
    umap_list <- list()
    
    for (i in 1:length(embeddings)) {
        umap_list[[i]] <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = coloringLayer, embedding = embeddings[i], size = size, pal = coloringPalette) +
                            ggtitle(embeddings[i]) + guides(color=guide_legend(nrow=ceiling(k/3), byrow=TRUE, override.aes = list(size=3))) +
                            theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.spacing.y = unit(0.5, "cm"))
    }
    
    return(umap_list)
}