#' Barplot of percentage of reads falling into expression bins
#'
#' \code{readcountExpBoxplot} Barplot of percentage of reads falling into 
#' expression bins
#'
#' This function makes a barplot of percentage of reads falling into expression
#' bins
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param stepSize The number of bars to be shown
#' @param ... Other parameters sent to barplot()
#' @return nothing Other parameters sent to barplot()
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#'
#' # barplot of percentage of reads falling into expression bins
#' readcountExpBoxplot(DupMat=dm)
readcountExpBoxplot <- function(DupMat, stepSize=0.05, ...) {

    ## define bin upper bounds
    binVec <- seq(stepSize, 1, stepSize)
  
    ## read counts in each RPK bin
    RpkBinReadCounts <- sapply(binVec,
                               getRpkBinReadCountFraction,
                               DupMat=DupMat,
                               stepSize=stepSize)
    
    names(RpkBinReadCounts) <- paste(as.character(binVec - stepSize),
                                     as.character(binVec), sep="-")
  
    barplot(RpkBinReadCounts,
            ylim=c(0,1),
            las=3,
            ylab="% of total reads in expression bin",
            main="Percentage of Reads in Expression Bin.",
            ...)
    
    title(xlab="expression bin", line=8)
  
}

