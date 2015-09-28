#' Dynamic range
#'
#' \code{getDynamicRange} Calculate the dynamic range of the RNAseq experiment
#'
#' This function calculates the dynamic range of the RNAseq eperiment
#'
#' @param dm The duplication matrix calculated by \code{analyzeDuprates}
#' @return A list with 2 elements, containing the dynamic range counting all
#' reads and the dynamic range after removing duplicates.
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # calculate the dynamic range
#' getDynamicRange(dm)
getDynamicRange <- function(dm) {

  step <- 0.01
  qvec <- seq(step, 1, step)
  
  medianrpkm<- function(p, dm, col) {
    ## generate a subset of the matrix for read counts in a given bin
    tm <- dm[dm[,col] <  quantile(dm[,col], p) &
             dm[,col] >= quantile(dm[,col], p-step),]
    return(median(tm[,col]))
  }

  ## calculate dynamic range:
  ## (median highest 1% expressed / median lowest 1% expressed)
  
  ## median bin counts (all)
  mbc.a <- sapply(qvec, medianrpkm, dm, "RPKM")
  dr.a <- mbc.a[length(mbc.a)]/mbc.a[mbc.a>0 & !is.na(mbc.a)][1]

  ## median bin counts (duprm)
  dm$RPKMU <- (1 - dm$dupRate) * dm$RPKM
  dm$RPKMU <- ifelse(is.na(dm$RPKMU),dm$RPKM,dm$RPKMU)
  mbc.u <- sapply(qvec, medianrpkm, dm, "RPKMU");
  dr.u  <- mbc.u[length(mbc.u)]/mbc.u[mbc.u>0 & !is.na(mbc.u)][1]

  return(list("dynrange.all"=dr.a,
              "dynrange.duprm"=dr.u))
  
}

