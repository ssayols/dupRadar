#' Helper function used in \code{readcountExpressionBoxplot}
#'
#' \code{getRpkCumulativeReadCountFraction} get the cumulative read count
#'   fraction
#'
#' @param p The vector of bins
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The cumulative read count fraction
getRpkCumulativeReadCountFraction <- function(p, DupMat=DupMat) {

    ## generate a subset of the matrix for RPK values < than RPK percentile p
    tm <- DupMat[DupMat[,"RPK"] < quantile(DupMat[,"RPK"], p),]
    ## calculate the global duplication rate for the remaining values
    ReadFractionBelowP <- (sum(tm[,"dupsPerId"])/sum(tm[,"allCounts"]))
    
    return(ReadFractionBelowP)
}
