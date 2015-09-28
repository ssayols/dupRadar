#' Draw histogram with the expression values
#'
#' \code{expressionHist} Draw histogram with the expression values
#'
#' This function draws a histogram of the expression values from the
#' duplication matrix calculated by \code{analyzeDuprates}.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param value The column from the duplication matrix containing the
#'   expression values
#' @param ... Other parameters sent to hist()
#' @return nothing
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # histogram of expression values for annotation
#' expressionHist(DupMat=dm)
expressionHist <- function(DupMat, value="RPK", ...) {
    hist(log10(DupMat[,value]),
         axes=FALSE,
         breaks=100,
         xlab="reads per kilobase (RPK)",
         main="",
         ...)
    r <- round(range(log10(DupMat[,"RPK"]), finite=TRUE))
    s <- seq(r[1], r[2], 1)
    axis(2)
    axis(1, labels=10^s, at=s)
    abline(v=3, col="red")
}
