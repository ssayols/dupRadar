#' Duplication rate ~ total reads per kilobase (RPK) boxplot
#'
#' \code{duprateExpBoxplot} Duplication rate ~ total reads per kilobase (RPK) boxplot
#'
#' This function makes a boxplot showing the distribution of per gene duplication rate versus
#' the reads per kilobase (RPK) inside gene expression bins.
#
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param stepSize Expression bin seze for the boxplot
#' @param ... Other params sent to boxplot
#' @return nothing
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # duprate boxplot
#' duprateExpBoxplot(DupMat=dm)
duprateExpBoxplot <- function(DupMat, stepSize=0.05, ...) {
  
    ## define expression bins for boxplot
    binVec <- seq(stepSize, 1, stepSize)
  
    ## get mean duplication rate per bin
    ## for axis annotation

    binMeanExpression <- sapply(binVec, getBinRpkMean, stepSize, DupMat)


    ## get list with duplication rate vectors per bin
    ## for boxplot
    binDupDistribution <- sapply(binVec, getBinDuplication, stepSize, DupMat)
  
    names(binDupDistribution) <- paste0(as.character(100 * (binVec - stepSize)),
                                        " - ",as.character(100 * binVec)," %",
                                        " / ", binMeanExpression)
 
    opar <- par(mar=c(6,4,4,2))
    boxplot(binDupDistribution,
            ylim=c(0,1),
            las=3,
            ylab="duplication (%)",
            axes=FALSE,
            ...)
    axis(2)
    axis(1,,las=2,cex.axis=.6,
         labels=names(binDupDistribution),
         at=1:length(names(binDupDistribution)))
    title(xlab="mean expression (reads/kbp)", line=5)

    par(opar)
}
