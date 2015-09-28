#' dupRadar.
#'
#' @name dupRadar
#' @docType package
#' @import Rsubread 
NULL

#' Example data containing precomputed matrices for two RNASeq experiments
#'
#' Precomputed duplication matrices for two RNASeq experiments used as examples
#' of a good and a failed (in terms of high redundancy of reads) experiments. 
#' The experiments come from the ENCODE project, as a source of a broad variety
#' of protocols, library types and sequencing facilities.
#'
#' @docType data
#' @keywords datasets
#' @name dupRadar_examples
#' @usage data(dupRadar_examples)
#' @format A list with two example duplication matrices
NULL

#' Duplication matrix of a good RNASeq experiment
#'
#' A dataset containing the duplication matrix of a good RNASeq experiment, in
#' terms of duplicates. Comes from the GM12878 SR1x75 replicate 2 from Caltech
#' (UCSC's table Browser name: wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2)
#'
#' @docType data
#' @keywords datasets
#' @name dm
#' @usage data(dupRadar_examples)
#' @format A data frame with 23228 rows and 14 variables
NULL

#' Duplication matrix of a failed RNASeq experiment
#'
#' A dataset containing the duplication matrix of a failed RNASeq experiment,
#' containing unusual duplication rate. Comes from the HCT116 PE2x75 replicate 1
#' from Caltech (UCSC's table Browser name: 
#' wgEncodeCaltechRnaSeqHct116R2x75Il200AlignsRep1V2)
#'
#' @docType data
#' @keywords datasets
#' @name dm.bad
#' @usage data(dupRadar_examples)
#' @format A data frame with 23228 rows and 14 variables
NULL
