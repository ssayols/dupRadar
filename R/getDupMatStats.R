#' Report duplication stats on regions
#'
#' \code{getDupMatStats} Report duplication stats based on the data calculated
#'   in the duplication matrix
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return A data.frame containing the stats about the number of genes covered
#'   (1+ tags) and the number of genes containing duplicates (1+)
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # call the plot and identify genes
#' getDupMatStats(DupMat=dm)
getDupMatStats <- function(DupMat) {

    nRegions <- nrow(DupMat)
    nRegionsCovered <- nrow(DupMat[DupMat[,"allCounts"]>0,])
    fRegionsCovered <- nRegionsCovered/nRegions
    nRegionsDuplication <- nrow(DupMat[DupMat[,"dupRate"]>0,])
    fRegionsDuplication <- nRegionsDuplication/nRegions
    fCoveredRegionsDuplication <- nRegionsDuplication/nRegionsCovered

    return(c(nRegions=nRegions,
             nRegionsCovered=nRegionsCovered,
             fRegionsCovered=fRegionsCovered,
             nRegionsDuplication=nRegionsDuplication,
             fRegionsDuplication=fRegionsDuplication,
             fCoveredRegionsDuplication=fCoveredRegionsDuplication))
}
