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
#' @keywords fastq mapping
#' @export

mapToBarcodes <- function(fastqCases,fastqControls,labels,oligoMap,oligoOut) {
  path <- paste(system.file(package="oligoGames"), 'exec', "mapToBarcodes.py", sep="/")
  options <- paste("--labels", labels, "--oligoMap", oligoMap, "--out", oligoOut)
  command <- paste("python", path, options, fastqCases, fastqControls)
  response <- system(command, intern=T)
  print(response)
}