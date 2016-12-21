#' A function to normalize counts
#'
#' This function performs library normalization for the experiment counts
#' @param normalizedCounts matrix of normalized counts. Defaults to
#'   allTranscripts_Norm.tab
#' @param metaData A file with the meta data of the experiment. Should include
#'   details such as the name of genes tiled, window size, length of oligo.
#'   Please see oligoMeta.tab in extdata for reference.
#' @param modelMethod The modeling method used to get nucleotide counts from
#'   oligo counts. Must be either "median", "sum", or "pgm".  Defaults to 
#'   "median".
#' @param oligoLen numeric value for the length of the oligos if all oligos have
#'   the same value and this information is not included in the metadata file.   
#' @import reshape2
#' @import dplyr
#' @import tidyr
#' @importFrom stats median
#' @importFrom utils read.table
#' @keywords modeling
#' @return modeledNucs matrix of modeled nucleotide counts
#' @export

##################################################
# Things to Do:
#   1). Implement a PGM to model nucleotide counts
#   2). Extend the function to recognize the # of
#       columns for "summarise"
##################################################

modelNucCounts <- function(normalizedCounts, metaData, 
                           modelMethod=c("median", "sum", "pgm"),
                           oligoLen=NULL){
  #edit the normalized counts to include oligo number & oligoID
  oligoID <- sapply(normalizedCounts$Transcript, function(x) unlist(strsplit(x, "_")))
  normalizedCounts$oligoNum <- as.numeric(sapply(oligoID, function(x) x[(length(x)-1)])) + 1
  oligoID <- sapply(oligoID, function(x) x[-c((length(x)-1):length(x))])
  oligoID <- sapply(oligoID, function(x) paste(x, collapse="_"))
  normalizedCounts$oligoID <- oligoID

  meta <- read.table(metaData, header=TRUE, stringsAsFactors=FALSE)
  # change index to 1-based
  meta$startBarcodeIndex <- meta$startBarcodeIndex + 1

  # add oligo index - unlike barcode index, no overlaps
  meta$startOligoIndex <- NA
  meta$endOligoIndex <- cumsum(meta$numOfOligos)
  meta$startOligoIndex  <- meta$endOligoIndex - meta$numOfOligos + 1

  # create a mapping of which row of meta dataframe (which lncRNA) corresponds to
  # each row of the normalizedCounts dataframe
  x <- match(normalizedCounts$oligoID, meta$name2)

  # add bp start and end positions of each oligo assuming
  # that the gap (bp) between this oligo and the next is
  # usually 10, except for one really long lncRNA
  # and the last oligo of each lncRNA can be shorter than 10 so
  # that the last position of the last oligo corresponds to the
  # last position of the lncRNA
  normalizedCounts$Gap <- meta$window[x]
  if ("oligoLen" %in% colnames(meta)){
     normalizedCounts$oligoLen <- meta$oligoLen[x]
  }else if (!is.null(oligoLen)){
     normalizedCounts$oligoLen <- oligoLen
     meta$oligoLen <- oligoLen
  }else{
    stop("Error: Need to specify the length of the oligos either in the
          metadata file or using the oligoLen parameter if they are all equal")
  }
  # add bp index (within lncRNA) to normalizedCounts
  normalizedCounts$bpStart <- (normalizedCounts$oligoNum-1)*normalizedCounts$Gap + 1
  normalizedCounts$bpEnd <- normalizedCounts$bpStart + normalizedCounts$oligoLen - 1

  last <- first <- second <- penultimate <- rep(NA, length(meta$name2))
  for (ol in 1:length(meta$name2)){
    last[ol] <- which(normalizedCounts$oligoNum == meta$numOfOligos[ol] &
                        normalizedCounts$oligoID == meta$name2[ol])
    first[ol] <- which(normalizedCounts$oligoNum == 1 & normalizedCounts$oligoID == meta$name2[ol])

    penultimate[ol] <- which(normalizedCounts$oligoNum == (meta$numOfOligos[ol]-1) &
                               normalizedCounts$oligoID == meta$name2[ol])
    second[ol] <- which(normalizedCounts$oligoNum == 2 & normalizedCounts$oligoID == meta$name2[ol])
  }
  normalizedCounts$bpStart[last] <- (meta$seqLen - normalizedCounts$oligoLen[last] + 1)
  normalizedCounts$bpEnd[last] <- meta$seqLen


  # check that there are no missing oligos
  numOfOligos1 <- 0
  numOfOligos2 <- 0
  for (j in unique(normalizedCounts$oligoID)){
    nums <- normalizedCounts$oligoNum[normalizedCounts$oligoID==j]
    numOfOligos1 <- numOfOligos1 + length(unique(nums))
    numOfOligos2 <- numOfOligos2 + length(min(nums):max(nums))
  }
  numOfOligos1 == sum(meta$numOfOligos)  # (should be true)
  numOfOligos2 == sum(meta$numOfOligos)
  
  if (!numOfOligos1 | !numOfOligos2){
    stop("Error: missing oligos in metadata file")
  }

  # convert one row per oligo to one row per basepair
  normalizedCounts$bps <- sapply(1:nrow(normalizedCounts), function(x)
    paste(normalizedCounts$bpStart[x]:normalizedCounts$bpEnd[x], collapse=","))
  bps <- (normalizedCounts %>%
            transform(bps=strsplit(bps,",")) %>%
            unnest(bps))
  bps$bps <- as.numeric(bps$bps)

  # for bps covered by more than one Oligo, add up the reads for each one
  # can change median() calls in the next line to sum() or some other function

  # Generalize the number of columns here

  if (modelMethod=="median") {
    bps_model <- bps %>%
      group_by(oligoID, bps) %>%
      summarise(N1=median(Nuclei_BioRep1), N2=median(Nuclei_BioRep2),
                N3=median(Nuclei_BioRep3), N4=median(Nuclei_BioRep4),
                T1=median(Total_BioRep1), T2=median(Total_BioRep2),
                T3=median(Total_BioRep3), T4=median(Total_BioRep4))
  }else if (modelMethod=="sum") {
      bps_model <- bps %>%
        group_by(oligoID, bps) %>%
        summarise(N1=sum(Nuclei_BioRep1), N2=sum(Nuclei_BioRep2),
                  N3=sum(Nuclei_BioRep3), N4=sum(Nuclei_BioRep4),
                  T1=sum(Total_BioRep1), T2=sum(Total_BioRep2),
                  T3=sum(Total_BioRep3), T4=sum(Total_BioRep4))
  } else {
    message("pgm to be implemented here...")
  }

  # correct for the oligo that has 1/4 the overlapping oligos
  # because they are spaced out at lower resolution
  # this needs to be robustified
  bps_model[bps_model$oligoID=="NR_002728",c(3:10)] <-
    bps_model[bps_model$oligoID=="NR_002728",c(3:10)]*4

  # make one object with all lncRNA data,
  # discard basepairs with same signal by design
  IDs <- unique(meta$name1)
  modeledNucs <- NULL
  chr <- pos <- NULL
  ct <- 1
  for (id in IDs){

    geneID <- meta$name2[meta$name1==id]
    gap <- meta$window[meta$name1==id]
    fish <- meta$fishClass[meta$name1==id]
    NUCS <- data.frame(bps_model[bps_model$oligoID==geneID,-1])
    rownames(NUCS) <- NUCS[,1]
    NUCS <- as.matrix(NUCS[,-1])

    diffMat <- apply(NUCS, 2, diff)
    chgpts <- unique(c(seq(1,nrow(NUCS), by=gap) - 1, nrow(NUCS) - meta[meta$name1==id,]$oligoLen))
    if (meta[meta$name1==id,]$oligoLen %% gap != 0){
      starting <- seq(1,nrow(NUCS), by=gap) - 1
      chgpts <- unique(sort(c(chgpts, starting+meta[meta$name1==id,]$oligoLen)))
      rmv <- which(chgpts > nrow(NUCS))
      if (length(rmv)>0){
        chgpts <- chgpts[-rmv]
      }
    }
    diffMatO <- diffMat[-chgpts,]
    diffMatC <- diffMat[chgpts,]

    if(sum(diffMatO != 0) != 0 | sum(diffMatC != 0) == 0){
      stop("Error pulling out basepairs with same signal by design")
    }
    chgpts <- chgpts + 1

    rownames(NUCS) <- rep(id, nrow(NUCS))
    modeledNucs <- rbind(modeledNucs, NUCS[chgpts,])

    chr <- c(chr, rep(id, length(chgpts)))
    pos <- c(pos, chgpts)
    ct <- ct + 1
  }

  modeledNucs <- cbind(pos, modeledNucs)
  return(modeledNucs)

  ## modeledNucs is a matrix with one row per 10bp gap (or 40bp gap) (since all
  ## basepairs in the gap between oligos will have the identical signal
  ## rownames indicate the oligo name
  ## first column (pos) indicates the bp position in that lncRNA
  ## column names of col 2 through 13 indicate condition
  ## (C = cytosol, N = nuclear, T = total), reps 1-4
}
