#' A function to normalize counts
#' 
#' This function performs library normalization for the experiment counts
#' @param countsFile the name of the file with counts for each oligo. Contains 
#'   one row per oligo counts with each sample its own column. Please see 
#'   \code{system.file("data", "allTranscriptsCounts_Raw.tab", package =
#'   "oligoGames")} for an example. Defaults to allTranscriptsCounts_Raw.tsv
#' @param normType The method of library normalization the user prefers. Two 
#'   options: median or quantile. Defaults to median
#' @param quantile The quantile the user wishes to use. Should be a numeric 
#'   between 0 and 1.
#' @keywords normalization
#' @return normCounts matrix of normalized counts
#' @export

normCounts <- function(countsFile='allTranscriptsCounts_Raw.tsv',normType='median',quantile=0.5) {
  counts <- utils::read.table(countsFile, header=TRUE, stringsAsFactors=FALSE)
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
