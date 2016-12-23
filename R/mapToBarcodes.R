#' Map case and control fastq files to unique oligo barcodes and generate counts
#' table
#' 
#' This function allows you to map your sequencing files back to the oligo pool
#' and compile a table describing counts for each oligo. For a read to map to an
#' oligo, 2 conditions must be satisfied: (1) The first 10 nucleotides of the
#' read should perfectly map one barcode in the pool and (2) The remaining
#' nucleotides should have no more than 2 mismatches to the variable sequence
#' linked with the barcode.
#' 
#' @param fastqCases character vector with the name of 'case' fastq files.
#' @param fastqControls character vector with the name of 'control' fastq files.
#' @param conditionLabels character vector of length two which contains the 
#'   condition labels for the two conditions that are being compared. Default 
#'   value is c("case", "control")
#' @param oligoMap the name of the tab seperated file which provides a link 
#'   between oligo name, barcode and variable sequence. The file should have 
#'   columns for the name, barcode sequence and variable sequence of the oligo.
#' @param oligoOut the name of the directory where the output file will be 
#'   written. If the output directory does not exist, one will be created. 
#'   Defaults to oligoOut
#' @import rPython
#' @keywords fastq mapping
#' @export

# @examples
# fastqCases <- c(system.file("extdata", "case1.fastq.gz", package = "oligoGames"), system.file("extdata", "case2.fastq.gz", package = "oligoGames"))
# fastqControls <- c(system.file("extdata", "control1.fastq.gz", package = "oligoGames"), system.file("extdata", "control2.fastq.gz", package = "oligoGames"))
# oligoMap <- system.file("extdata", "oligoMap.fa", package = "oligoGames")
# oligoOut <- system.file("extdata", "oligoOut", package = "oligoGames")
# conditionLabels <- c("nuclear", "total")
# mapToBarcodes(fastqCases, fastqControls, conditionLabels, oligoMap, oligoOut)

mapToBarcodes <- function(fastqCases, fastqControls, conditionLabels=c("case", "control"),
                          oligoMap="oligoMap.fa", oligoOut="oligoOut") {
  python.load(system.file("exec", "mapToBarcodes.py", package = "oligoGames"))
  python.call("barcodeCounts", fastqCases, fastqControls, conditionLabels, oligoMap, oligoOut)
}
