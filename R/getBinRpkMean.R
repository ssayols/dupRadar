#' Helper function used in \code{duprateExpBoxplot}
#'
#' \code{getBinRpkMean} get mean duplication rate per bin
#'
#' @param p The vector of bins
#' @param stepSize The window size
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The averaged RPK per bin
getBinRpkMean <- function(p, stepSize, DupMat) {
    ## get subset of DupMat 
    dupMatBin <- getDupMatBin(p=p,
                              stepSize=stepSize,
                              value="RPK",
                              DupMat=DupMat)
    
    binRpkMean <- round(mean(dupMatBin[,"RPK"]), 1)
    
    return(binRpkMean)
}
