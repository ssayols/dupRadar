#' Helper function used in \code{getBinDuplication} and
#'   \code{getBinRpkMean}
#'
#' \code{getDupMatBin} get a subset of the matrix for values in a specific bin
#' defined by the upper bound p and stepSize
#'
#' @param p The vector of bins
#' @param stepSize The window size
#' @param value The column to be subset
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The subseted matrix
getDupMatBin <- function(p, stepSize=0.05, value="allCounts", DupMat) {
    ## subset of the matrix for values in a specific bin
    ## defined by the upper bound p and stepSize
    dupMatBin <- DupMat[DupMat[,value] <
                        quantile(DupMat[,value], p) &
                        DupMat[,value] >=
                        quantile(DupMat[,value], p - stepSize),]
    
    return(dupMatBin)
}
