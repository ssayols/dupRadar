#' Program dispatchers to mark duplicated reads from a BAM file
#'
#' \code{markDuplicates} Mark duplicated reads from a BAM file by calling
#'   widely used tools.
#'
#' This function works as a wrapper for several tools widely adopted tr mark
#' duplicated reads in a BAM file. Currently, it supports PICARD and BamUtils.
#'
#' @param dupremover The tool to be called. Currently, "picard" and "bamutils"
#'   are supported
#' @param bam The bam file to mark duplicates from
#' @param out Regular expression describing the transformation on the original
#'   filename to get the output filename. By default, a "_duprm" suffix is 
#'   added before the bam extension
#' @param rminput Whether to keep the original, non duplicate-marked, bam file
#' @param path Path to the duplicate marker binaries
#' @param verbose Redirect all the program output to the R console
#' @param ... Other parameters sent to the caller function
#' @return The output filename
#' @examples
#' \dontrun{
#' bam <- system.file("extdata","sample1Aligned.out.bam",package="dupRadar")
#' gtf <- "genes.gtf"
#' stranded <- 2    # '0' (unstranded), '1' (stranded) and '2' (reverse)
#' paired   <- FALSE
#' threads  <- 4
#' 
#' # call the duplicate marker and analyze the reads
#' bamDuprm <- markDuplicates(dupremover="bamutil",bam,
#'                            path="/opt/bamUtil-master/bin",rminput=FALSE)
#' dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)
#' }
markDuplicates <- function(dupremover="bamutil",bam=NULL,
                           out=gsub("\\.bam$","_duprm.bam",bam),rminput=TRUE,
                           path=".",verbose=TRUE,...) {

    # check input parms
    if(!file.exists(bam))  stop("file",bam,"not found!")
    if(!file.exists(path)) stop("dir",path,"not found!")
    if(!file.exists(dirname(out))) {
        if(verbose) cat("creating dir:",dirname(out),fill=TRUE)
        dir.create(dirname(out),recursive=TRUE)
    }

    # dispatch the function
    if(dupremover == "picard") {
        picardMarkDuplicates (bam=bam,out=out,path=path,verbose=verbose,...)
    }
    else if(dupremover == "bamutil") {
        bamutilMarkDuplicates(bam=bam,out=out,path=path,verbose=verbose,...)
    }
    else {
        stop("dupremover must be either picard or bamutil!")
    }

    # remove non duplicate-marked bam file
    if(rminput) file.remove(bam)

    invisible(out)
}

#' Mark duplicates using Picard tools
#'
#' \code{picardMarkDuplicates} Mark duplicated reads from a BAM file by calling
#'   picard tools
#'
#' This function is supposed to be called through the \code{markDuplicates}
#'   wrapper
#'
#' @param bam The bam file to mark duplicates from
#' @param out Regular expression describing the transformation on the original
#'   filename to get the output filename. By default, a "_duprm" suffix is
#'   added before the bam extension
#' @param path Path to the duplicate marker binaries
#' @param verbose Redirect all the program output to the R console
#' @param threads Number of threads to use
#' @param maxmem Max memory assigned to the jvm
#' @return The return code of the system call
picardMarkDuplicates <- function(bam,out,path,verbose,threads=1,maxmem="4g") {

    # call picard MarkDuplicates
    cmd <- paste0(" -XX:ParallelGCThreads=",threads,
                  " -Xmx",maxmem,
                  " -jar ",file.path(path,"MarkDuplicates.jar"),
                  " INPUT=",bam,
                  " OUTPUT=",out,
                  " METRICS_FILE=",gsub("\\.bam$","_picard_metrics.txt",bam),
                  " REMOVE_DUPLICATES=false",
                  " ASSUME_SORTED=true",
                  " PROGRAM_RECORD_ID=\'null\'")

    if(verbose) cat(paste("Running",cmd),fill=TRUE)
    rc <- system2("java",cmd,stdout=ifelse(verbose,"","FALSE"),
                  stderr=ifelse(verbose,"","FALSE"))

    # check RC
    if(rc == 127) stop("command java",cmd,"could not be run!")
    if(rc != 0)   stop("picard returned error code",rc,"in command: java",cmd)
    invisible(rc)
}


#' Mark duplicates using bamutil
#'
#' \code{bamutilMarkDuplicates} Mark duplicated reads from a BAM file by
#'   calling bamutil
#'
#' This function is supposed to be called through the \code{markDuplicates}
#'   wrapper
#'
#' @param bam The bam file to mark duplicates from
#' @param out Regular expression describing the transformation on the original
#'   filename to get the output filename. By default, a "_duprm" suffix is
#'   added before the bam extension
#' @param path Path to the duplicate marker binaries
#' @param verbose Redirect all the program output to the R console
#' @return The return code of the system call
bamutilMarkDuplicates <- function(bam,out,path,verbose) {

    # call bamutil dedup
    cmd <- paste("dedup --in",bam,"--out",out)

    if(verbose) cat(paste("Running",cmd),fill=TRUE)
    rc <- system2(file.path(path,"bam"),cmd,stdout=ifelse(verbose,"","FALSE"),
                  stderr=ifelse(verbose,"","FALSE"))

    # check RC
    if(rc == 127) stop("command bam",cmd,"could not be run!")
    if(rc != 0)   stop("bamutil returned error code",rc,"in command: bam",cmd)
    invisible(rc)
}
