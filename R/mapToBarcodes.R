#' A function to map to barcodes
#'
#' This function allows you to map fastq files to barcodes and get counts table
#' @param fastqCases string of comma seperated fastq files for cases
#' @param fastqControls string of comma seperated fastq files for controls
#' @param labels string of comma seperated labels for case and control. Defaults
#'   to nuclear,total
#' @param oligoMap string with the path to a file linking each barcode to a
#'   variable sequence and transcript ID. Defaults to oligoMap.tsv
#' @param oligoOut string with the path to the directory where all output files
#'   will be written. If the output directory does not exist, one will be
#'   created. Defaults to oligoOut
#' @import rPython
#' @keywords fastq mapping
#' @export
# @examples
# fastqCases <- system.file("extdata", "case.fastq.gz", package = "oligoGames")
# fastqCases <- system.file("extdata", "control.fastq.gz", package = "oligoGames")
# oligoMap <- system.file("extdata", "oligoMap.fa", package = "oligoGames")
# oligoOut <- system.file("extdata", "oligoOut", package = "oligoGames")
# labels <- "nuclear,total"
#

mapToBarcodes <- function(fastqCases, fastqControls, labels="nuclear,total",
                          oligoMap="oligoMap.fa", oligoOut="oligoOut") {
  python.load(system.file("exec", "mapToBarcodes.py", package = "oligoGames"))
  python.call("barcodeCounts", fastqCases, fastqControls, labels, oligoMap, oligoOut)
}
