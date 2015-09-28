#' Helper function used in \code{readcountExpressionBoxplot}
#'
#' \code{readcountExpressionBoxplot} Calculates the fraction of total reads in
#'   a vector of bins
#'
#' @param p The vector of bins
#' @param stepSize The window size
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The fraction of total reads in a vector of bins
getRpkBinReadCountFraction <- function(p, stepSize=stepSize, DupMat=DupMat) {

    ## generate a subset of the matrix for RPK values given bin
    ## defined by the upper bound p and the lower bound p-step
    tm<-DupMat[DupMat[,"RPK"] <  quantile(DupMat[,"RPK"], p) &
               DupMat[,"RPK"] >= quantile(DupMat[,"RPK"], p-stepSize),]
    
    ## calculate the fraction of total reads in given bin
    BinReadFraction <- sum(tm[,"allCounts"])/sum(DupMat[,"allCounts"])
    return(BinReadFraction)
}
