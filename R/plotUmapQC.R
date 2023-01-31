#' plotUmapQC
#' 
#' This function returns a list of UMAPs with cells colored by various quality control meassures (total number of fragments, blacklist ratio, TSS enrichment, promoter ratio, nucleosome ratio).
#' @param ArchRProj ArchRProject with priorly calculated reduced dimensions and single-cell embedding.
#' @param embedding Optional name of calculated single-cell embedding (e. g. calculated with addUMAP function). Defaults to first embedding in ArchRProj@embeddings.
#' @param coloringLayer Optional name of cellColData column to use as comparative coloring layer of single cells. Defaults to 'Sample' column.
#' @param title Optional title of quality control UMAPs.
#' @param coloringPalette Optional coloring palette for reference coloring of single cells.
#' @param doubletScore Optional addition of doublet score as QC metric.
#' @return List of UMAPs. First with comparative coloring layer for quality assessment of single cells. Other UMAPs colored by total number of fragments, blacklist ratio, TSS enrichment, promoter ratio, nucleosome ratio, and if set doublet score.
#' @keywords umap quality qc archr rwire archrwire
#' @examples Coming soon.
#' @export



plotUmapQC <- function(ArchRProj, embedding=names(ArchRProj@embeddings)[1], coloringLayer='Sample', title='', coloringPalette=NULL, doubletScore=FALSE){
    k <- unique(ArchRProj@cellColData[,coloringLayer]) %>% length(.)
    
    umap_samples <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = coloringLayer, embedding = embedding, size = 0.3, pal = coloringPalette) +
    ggtitle(title) + guides(color=guide_legend(nrow=ceiling(k/3), byrow=TRUE, override.aes = list(size=3))) +
    theme(plot.title = element_text(size=33), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.spacing.y = unit(0.5, "cm"))
    
    umap_nfrags <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "nFrags", embedding = embedding, size = 0.3) +
    ggtitle("Number of fragments") + 
    theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3,"cm"))

    umap_blacklist <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "BlacklistRatio", embedding = embedding, size = 0.3) + 
    ggtitle("Blacklist ratio") +
    theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3,"cm"))

    umap_tss <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "TSSEnrichment", embedding = embedding, size = 0.3) + 
    ggtitle("TSS enrichment") +
    theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3,"cm"))

    umap_promoter <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "PromoterRatio", embedding = embedding, size = 0.3) + 
    ggtitle("Promoter ratio") +
    theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3, "cm"))

    umap_nuc <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "NucleosomeRatio", embedding = embedding, size = 0.3) + 
    ggtitle("Nucleosome ratio") +
    theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3, "cm"))
    
    if(doubletScore){
        umap_doub <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "DoubletScore", embedding = embedding, size = 0.3) + 
        ggtitle("Doublet score") +
        theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3, "cm"))
        
        umap_doub2 <- plotEmbedding(ArchRProj, colorBy = "cellColData", name = "DoubletEnrichment", embedding = embedding, size = 0.3) + 
        ggtitle("Doublet enrichment") +
        theme(plot.title = element_text(size=28), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=0), legend.text = element_text(size=15), legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3, "cm"))
        
        return(list(umap_samples, umap_nfrags, umap_blacklist, umap_tss, umap_promoter, umap_nuc, umap_doub, umap_doub2))
    } else {
        return(list(umap_samples, umap_nfrags, umap_blacklist, umap_tss, umap_promoter, umap_nuc))
    }
}