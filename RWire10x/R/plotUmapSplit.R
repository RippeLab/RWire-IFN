#' plotUmapSplit
#' 
#' This function returns a list of UMAPs with cells of single split levels.
#' @param ArchRProj ArchRProject with priorly calculated reduced dimensions and single-cell embedding.
#' @param embedding Optional name of calculated single-cell embedding (e. g. calculated with addUMAP function). Defaults to first embedding in ArchRProj@embeddings.
#' @param splitLayer Optional name of cellColData column to use as split layer of single cells. Defaults to 'Sample' column.
#' @param selectLevels Optional selection of certain levels in split layer to visualize. Defaults to all levels in cellColData 'Sample' column.
#' @param coloringPalette Optional coloring palette for split levels of single cells.
#' @return List of UMAPs with cells of single split levels.
#' @keywords umap split archr rwire archrwire
#' @examples Coming soon.
#' @export



plotUmapSplit <- function(ArchRProj, embedding=names(ArchRProj@embeddings)[1], splitLayer='Sample', selectLevels=unique(ArchRProj$Sample), coloringPalette=NULL){
    split_size <- selectLevels %>% length(.)
    
    if(is.null(coloringPalette)){
        coloringPalette <- paletteDiscrete(selectLevels)
    }
    
    umap_split <- list()
    
    for (i in 1:split_size){
        pal <- coloringPalette
        pal[-i] <- '#FFFFFF'
        
        umap_split[[i]] <- plotEmbedding(ArchRProj, colorBy = 'cellColData', name = splitLayer, embedding = embedding, size = 0.3, pal = pal) +
                            ggtitle(names(coloringPalette)[i]) + guides(color=guide_legend(nrow=ceiling(split_size/3), byrow=TRUE, override.aes = list(size=3))) +
                            theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.spacing.y = unit(0.5, 'cm'))
    }
    
    return(umap_split)
}