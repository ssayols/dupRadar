#' Duplication rate ~ total read count fit model
#'
#' \code{duprateExpDensPlot} Duplication rate ~ total read count fit model
#'
#' Fit a Generalized Linear Model using a logit function between thegene 
#' duplication rate and the total read count.
#'
#' @param DupMat The duplication matrix calculated by \code{analyzeDuprates}
#' @return The GLM and the coefficients of the fitted logit function
#' @examples
#' # dm is a duplication matrix calculated by analyzeDuprates:
#' # R> dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' attach(dupRadar_examples)
#'
#' # duprate plot
#' duprateExpFit(DupMat=dm)
duprateExpFit <- function(DupMat) {
  
	x <- log10(DupMat$RPK)
	y <- DupMat$dupRate
  
	fit <- try(suppressWarnings(glm(y ~ x,family=binomial(link="logit"))))
	
	## alternative:
	## fit without the 0% duprate genes (mostly loci with only 1 read)
	## Note: result deviates more strongly from our empirical distribution 
	## than taking all
	## yclean <- y[y!=0]
	## xclean <- x[y!=0]
	## fit <- try(suppressWarnings(glm(yclean ~ xclean,
	##                                 family=binomial(link="logit"))))
	## newx <- sort(xclean)
	## newy <- predict(fit,data.frame(xclean=sort(xclean)),type="response")
	if(class(fit)[1] != "try-error") {
		x <- exp(coef(fit))
		fit <- list(glm=fit,
                            intercept=x[1],
                            slope=x[2])
            }
        
	return(fit)
    }
