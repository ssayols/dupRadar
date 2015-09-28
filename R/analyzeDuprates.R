#' Read in a BAM file and count the tags falling on the features described
#' in the GTF file
#'
#' \code{analyzeDuprates} returns a data.frame with tag counts
#'
#' This function makes use of the Rsubread package to count tags on the GTF
#' features in different scenarios. The scenarios are the 4 possible
#' combinations of allowing multimappers (yes/no) and duplicate reads (yes/no).
#'
#' @param bam The bam file containing the duplicate-marked reads
#' @param gtf The gtf file describing the features
#' @param stranded Whether the reads are strand specific
#' @param paired Paired end experiment?
#' @param threads The number of threads to be used for counting
#' @param verbose Whether to output Rsubread messages into the console
#' @param ... Other params sent to featureCounts
#' @return A data.frame with counts on features, with and without taking into
#'   account multimappers/duplicated reads
#' @examples
#' bam <- system.file("extdata",
#'                    "wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2_duprm.bam",
#'                    package="dupRadar")
#' gtf <- system.file("extdata","genes.gtf",package="dupRadar")
#' stranded <- 2    # '0' (unstranded), '1' (stranded) and '2' (reverse)
#' paired   <- FALSE
#' threads  <- 4
#'
#' # call the duplicate marker and analyze the reads
#' dm <- analyzeDuprates(bam,gtf,stranded,paired,threads)
analyzeDuprates <- function(bam,gtf,stranded=0,paired=FALSE,threads=1,
                            verbose=FALSE,...) {

    # check input parameters
    if(!file.exists(bam))     stop("file",bam,"not found!")
    if(!file.exists(gtf))     stop("file",gtf,"not found!")
    if(!is.logical(paired))   stop("paired has to be either TRUE/FALSE")
    if(!is.numeric(stranded)) stop("stranded has to be a number [0-2]")
    if(stranded < 0 || stranded > 2) stop("stranded has to be a number [0-2]")
    if(!is.numeric(threads)) stop("threads has to be a number")

    # featureCounts simplified call
    count <- function(mh,dup) {
        Rsubread::featureCounts(files=bam,annot.ext=gtf,
			isGTFAnnotationFile=TRUE,GTF.featureType="exon",
			GTF.attrType="gene_id",nthreads=threads, isPairedEnd=paired,
			strandSpecific=stranded,ignoreDup=dup,countMultiMappingReads=mh,...)
    }

    ## HK: status info for each call of featureCounts()? In silent mode it's
    ## not obvious any more what is happening. This would require assigning the
    ##results to individual variables first and generating the list afterwards.
    if (verbose) {
        counts <- list(mhdup    =count(mh=TRUE ,dup=FALSE), 
                       mhnodup  =count(mh=TRUE ,dup=TRUE ), 
                       nomhdup  =count(mh=FALSE,dup=FALSE), 
                       nomhnodup=count(mh=FALSE,dup=TRUE ))
    }
    else { 
        ## assign output on STDOUT from featureCounts to variable silencer,
        ## which we ignore afterwards. 
        silencer <- capture.output(
            counts <- list(mhdup    =count(mh=TRUE ,dup=FALSE), 
                           mhnodup  =count(mh=TRUE ,dup=TRUE ), 
                           nomhdup  =count(mh=FALSE,dup=FALSE), 
                           nomhnodup=count(mh=FALSE,dup=TRUE ))
            )
    }
        
    # calculate RPK & RPKM per gene
    x <- lapply(counts,function(x) { 
        # total number of mapped reads
        N <- sum(x$stat[,2]) - x$stat[x$stat$Status == "Unassigned_Unmapped",2]

        x <- data.frame(gene=rownames(x$counts),
                        width=x$annotation$Length[match(rownames(x$counts),
                                                  x$annotation$GeneID)],
                        counts=x$counts[,1],
                        RPK=0,
                        RPKM=0)
        x$RPK  <- x$counts * (10^3 / x$width)
        x$RPKM <- x$RPK * ( 10^6 / N)

        return(x)
    })

    # calculate duprates per gene, count duplicates per gene
    x <- data.frame(ID=x[[1]]$gene,
                    geneLength=x[[1]]$width,
                    allCountsMulti=x[[1]]$counts,
                    filteredCountsMulti=x[[2]]$counts,
                    dupRateMulti=(x[[1]]$counts-x[[2]]$counts) / x[[1]]$counts,
                    dupsPerIdMulti=x[[1]]$counts - x[[2]]$counts,
                    RPKMulti=x[[1]]$RPK,
                    RPKMMulti=x[[1]]$RPKM,
                    allCounts=x[[3]]$counts,
                    filteredCounts=x[[4]]$counts,
                    dupRate=(x[[3]]$counts - x[[4]]$counts) / x[[3]]$counts,
                    dupsPerId=x[[3]]$counts - x[[4]]$counts,
                    RPK=x[[3]]$RPK,
                    RPKM=x[[3]]$RPKM)
}
