#' Identify genes plotted by \code{duprateExpPlot}
#'
#' \code{duprateExpIdentify} Identify genes plotted by \code{duprateExpPlot}
#'
#' This function makes a barplot showing the cumulative read counts fraction
#' from the duplication matrix calculated by \code{analyzeDuprates}.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param idCol The column from the duplication matrix containing the labels
#' @return The identified points. \code{x} and \code{y} values match the ones 
#'  from \code{duprateExpPlot}
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # call the plot and identify genes
#' duprateExpPlot(DupMat=dm)
#' duprateExpIdentify(DupMat=dm)
duprateExpIdentify <- function(DupMat, idCol="ID") {
    marked <- identify(log10(DupMat[,"RPK"]),
                       100*DupMat[,"dupRate"],
                       labels=DupMat[,idCol])
    return(marked)
}

