#' Duplication rate ~ total read count plot
#'
#' \code{duprateExpDensPlot} Duplication rate ~ total read count plot
#'
#' This function makes a scatter plot showing the per gene duplication rate
#' versus the total read count.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @param pal The color palette to use to display the density
#' @param tNoAlternative Display threshold of 1000 reads per kilobase
#' @param tRPKM Display threshold at a given RPKM level
#' @param tRPKMval The given RPKM level
#' @param tFit Whether to fit the model
#' @param addLegend Whether to add a legend to the plot
#' @param ... Other parameters sent to plot()
#' @return nothing
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # duprate plot
#' duprateExpDensPlot(DupMat=dm)
duprateExpDensPlot <- function(DupMat,
                               pal=c("cyan","blue","green","yellow","red"), 
                               tNoAlternative=TRUE, tRPKM=TRUE, tRPKMval=0.5, 
                               tFit=TRUE, addLegend=TRUE,
                               ...) {
  
  x <- log10(DupMat$RPK)
  y <- 100*DupMat$dupRate
  
  ## density-coloured scatterplot
  cols  <- colorRampPalette(pal)
  dcols <- densCols(x=x,y=y,colramp=cols,nbin=500)
  xaxs  <- seq(min(x[!is.infinite(x)]),max(x),length.out=5)
  yaxs  <- seq(0,100,length.out=5)
  plot(x=x,y=y,
       col=dcols,
       pch=20,cex=.25,
       axes=FALSE,xlab="expression (reads/kbp)",
       ylab="% duplicate reads",
       ylim=c(0,100),
       ...)
  Axis(side=2,at=yaxs,labels=as.character(round(yaxs,2)), las=2)

  ## r <- round(range(log10(DupMat[,"RPK"]/DupMat[,"geneLength"]), finite=TRUE))

  r <- round(range(log10(DupMat[,"RPK"]), finite=TRUE))
  s <- seq(r[1], r[2], 1)
  Axis(side=1, labels=10^s, at=s)
  
  ## fit logistic regression model
  if(tFit) {
      fit <- duprateExpFit(DupMat)
      if(class(fit)[1] != "try-error") { 
          lines(sort(x),
                100*predict(fit$glm,data.frame(x=sort(x)),type="response"),
                col="black",lwd=2,lty=3)
          legend("topleft", legend=c(paste("Int:", round(fit$intercept,2)),
                                paste("Sl:" , round(fit$slope,2))))
          
      }
  }
  
  ## threshold of 1000 reads per kilobase ( log10(1000)=3 )
  ## -> no way to put more than 1000 reads on 1000bp without duplication
  if (tNoAlternative) {
    abline(v=log10(1000), col="red", lty=2)
  }
  
  ## threshold at given RPKM level default: 0.5
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
