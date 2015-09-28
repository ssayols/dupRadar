#' Barplot showing the cumulative read counts fraction
#'
#' \code{cumulativeDuprateBarplot} Barplot showing the cumulative read counts
#'   fraction
#'
#' This function makes a barplot showing the cumulative read counts fraction
#' from the duplication matrix calculated by \code{analyzeDuprates}.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param stepSize The size of the windows used for plotting
#' @param ... Other params sent to barplot
#' @return nothing 
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # call the plot
#' cumulativeDuprateBarplot(DupMat=dm)
cumulativeDuprateBarplot <- function(DupMat, stepSize=0.05, ...) {
    ## Attention 1:
    ## Usually we have a certain amount of genes in our tables that have
    ## no reads at all.
    ## In these cases the duplication rate is NaN, so no points will be
    ## plotted for these. Thus the resulting plot seems to end at low
    ## expression rates.
    ## We can see from the lower percentile for which there is a point,
    ## which fraction of genes is non-expressed in the current dataset
    ##
    ## Attention 2:
    ## The global duplication rate can increase with decreasing RPKM percentiles.
    ## 
    ##
    
    vec <- seq(0, 1, stepSize)
    
    barplot(sapply(vec,
                   getRpkCumulativeReadCountFraction,
                   DupMat),
            names=as.character(vec*100),
            las=2,
            ylim=c(0,1),
            ylab="duplication rate",
            xlab="percentile RPK",
            main="global dup rate <= RPK percentile",
            ...)
    
    
}

