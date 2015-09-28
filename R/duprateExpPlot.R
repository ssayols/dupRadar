#' Duplication rate ~ total read count plot
#'
#' \code{duprateExpPlot} Duplication rate ~ total read count plot
#'
#' This function makes a smooth scatter plot showing the per gene duplication 
#' rate versus the total read count.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param tNoAlternative Display threshold of 1000 reads per kilobase
#' @param tRPKM Display threshold at a given RPKM level
#' @param tRPKMval The given RPKM level
#' @param addLegend Whether to add a legend to the plot
#' @param ... Other parameters sent to smoothScatter()
#' @return nothing
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # duprate plot
#' duprateExpPlot(DupMat=dm)
duprateExpPlot <- function(DupMat,
                           tNoAlternative=TRUE,
                           tRPKM=TRUE,
                           tRPKMval=0.5,
                           addLegend=TRUE,
                           ...) {
    smoothScatter(log10(DupMat[,"RPK"]),
                  100*DupMat[,"dupRate"],
                  xlab="expression level (reads/kbp)",
                  ylab="duplication level (% duplicate reads)",
                  axes=FALSE,
                  ylim=c(0,100),
                  ...)
    axis(2)

    ## FIX: axis annotation is too short if we take the range of RPKM
    ## instead of the range of RPK
    ## r <- round(range(log10(DupMat[,"RPK"]/DupMat[,"geneLength"]), finite=TRUE))
    r <- round(range(log10(DupMat[,"RPK"]), finite=TRUE))
    s <- seq(r[1], r[2], 1)
    axis(1, labels=10^s, at=s)
  
    ## threshold of 1000 reads per kilobase
    ## ( log10(1000)=3 )
    ## -> no way to put more than 1000 reads on 1000bp
    ## without duplication
    if (tNoAlternative) {
        ## 
        abline(v=log10(1000), col="red", lty=2)
    }
  
    ## threshold at given RPKM level
    ## default: 0.5
    if (tRPKM) {
        ## get the RPK for the closest greater value of RPKM in question 
        RPKgt <- sort(DupMat[DupMat[,"RPKM"]>=tRPKMval,"RPK"])[1]
        ## get the RPK for the closest greater value of RPKM in question
        RPKlt <- sort(DupMat[DupMat[,"RPKM"]<=tRPKMval,"RPK"], decreasing=TRUE)[1]
        RPKMthreshold <- mean(RPKgt, RPKlt)
        abline(v=log10(RPKMthreshold), col="green", lty=2)
    }
  
    if (addLegend) {
        text <- c()
        col  <- c()
        
        if (tNoAlternative) {
            text <- c(text, "1 read/bp")
            col  <- c(col, "red")
        }
    
        if (tRPKM) {
            text <- c(text, paste(tRPKMval, "RPKM"))
            col  <- c(col, "green")
        }
        legend("bottomright",
               legend=text,
               col=col,
               lty=2)
    }
}

