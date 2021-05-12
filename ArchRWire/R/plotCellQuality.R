#' plotCellQuality
#' 
#' This function returns a ggplot object with a scatter plot of TSS Enrichment over log10 Unique Fragments per cell. Default filtering criteria of TSS Enrichment of 4 and minimal number of unique fragments of 1000 are visualized.
#' @param ArchRProj ArchRProject.
#' @param sample Sample name.
#' @param cutoff_frags Optional additional cutoff for number of log10 Unique Fragments to visualize.
#' @param cutoff_frags Optional additional cutoff for TSS Enrichment to visualize.
#' @return GGplot object of scatter plot (TSS Enrichment over log10 Unique Fragemnts).
#' @keywords cell quality qc archr rwire archrwire
#' @examples Coming soon.
#' @export



plotCellQuality <- function(ArchRProj, sample, cutoff_frags=NULL, cutoff_tss=NULL){
    p <- ggPoint(x = getCellColData(ArchRProj[ArchRProj$cellNames[which(ArchRProj$Sample==sample)], ], select = "log10(nFrags)")[,1], 
                 y = getCellColData(ArchRProj[ArchRProj$cellNames[which(ArchRProj$Sample==sample)], ], select = "TSSEnrichment")[,1], 
                 colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment",
                 xlim = c(log10(500), quantile(getCellColData(ArchRProj[ArchRProj$cellNames[which(ArchRProj$Sample==sample)], ], select = "log10(nFrags)")[,1], probs = 1)),
                 ylim = c(0, quantile(getCellColData(ArchRProj[ArchRProj$cellNames[which(ArchRProj$Sample==sample)], ], select = "TSSEnrichment")[,1], probs = 1))) + 
            geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") + 
            ggtitle(paste0(sample, " (", sum(ArchRProj$Sample==sample), " cells)")) + 
            theme(plot.title = element_text(size=30), axis.title = element_text(size=28), axis.text = element_text(size=28), legend.title = element_text(size=28), legend.text = element_text(size=28), 
                  legend.key.height = unit(0.75, "cm"), legend.key.width = unit(3,"cm"))
    
    if(!is.null(cutoff_frags)){
        p <- p + geom_vline(xintercept = cutoff_frags, lty = "dashed", color="red")
    }
    if(!is.null(cutoff_tss)){
        p <- p + geom_hline(yintercept = cutoff_tss, lty = "dashed", color="red")
    }
    
    return(p)
}