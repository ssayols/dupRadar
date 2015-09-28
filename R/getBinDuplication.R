#' Helper function used in \code{duprateExpBoxplot}
#'
#' \code{getBinDuplication} get duplication rate for a subset of the
#' duplication matrix
#'
#' @param p The vector of bins
#' @param stepSize The window size
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The duplication rate per bin
getBinDuplication <- function(p, stepSize, DupMat) {
    ## get subset of DupMat
    dupMatBin <- getDupMatBin(p=p,
                              stepSize=stepSize,
                              value="RPK",
                              DupMat=DupMat)
    
    binDuprates <- dupMatBin[,"dupRate"]
  
    return(binDuprates)
}
