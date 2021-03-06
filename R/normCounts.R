#' A function to normalize counts
#' 
#' This function performs library normalization for the experiment counts
#' @param rawCounts the name of the file with counts for each oligo. Contains 
#'   one row per oligo counts with each sample its own column. Please see 
#'   \code{system.file("inst", "extdata", "allTranscriptsCounts_Raw.tsv",
#'   package = "oligoGames")} for an example. Make sure that the name of the
#'   transcript is in the same format as described in this file. Briefly, it is
#'   "GeneName_oligoNum_Barcode". It is critical that oligoNum is second to last
#'   and separated by "_". Defaults to allTranscriptsCounts_Raw.tsv.
#' @param normType The method of library normalization the user prefers. Two 
#'   options: median or quantile. Defaults to median
#' @param quantile The quantile the user wishes to use. Should be a numeric 
#'   between 0 and 1.
#' @keywords normalization
#' @return normCounts matrix of normalized counts
#' @export
#' @examples 
#' rawCounts = read.table(system.file("extdata", "allTranscriptsCounts_Raw.tsv", package = "oligoGames"), header=T, sep='\t')
#' normalizedCounts <- normCounts(rawCounts, normType='median')

normCounts <- function(rawCounts='allTranscriptsCounts_Raw.tsv',normType='median',quantile=0.5) {
  #counts <- utils::read.table(rawCounts, header=TRUE, stringsAsFactors=FALSE)
  counts <- rawCounts
  counts2 <- counts[,c(2:ncol(counts))]
  if (normType=='median') {
    sizes <- EBSeq::MedianNorm(counts2)
  } else {
    if (normType=='quantile'){
      sizes <- EBSeq::QuantileNorm(counts,quantile)
    }
  }
  counts2 <- EBSeq::GetNormalizedMat(counts2, sizes)
  counts[,c(2:ncol(counts))] <- counts2
  rm(counts2)
  normalizedCounts <- counts
  return(normalizedCounts)
}
