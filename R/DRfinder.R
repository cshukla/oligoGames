#' Perform inference to detect differential regions in oligo experiments
#'
#' This is the main inference function that aims to find regions with
#' differential signal between two conditions.  The two main steps of the 
#' procedure are (1) detect candidate regions from the nucleotide level 
#' signal, optionally smoothing to combat loss of power due to low coverage, 
#' and (2) evaluate the test statistic of condition difference across each
#' candidate region by comparing to a global null distribution generated 
#' by permuting sample labels.
#'
#' @param OligoSignal a data frame in the format returned by the function
#'   \code{modelNucCounts}. Contains one row per nucleotide count.
#'    The first column contains the basepair positions of the nucleotides and
#'    the remaining columns hold the counts themselves (one column per sample).
#' @param conditionLabels character vector of length two which contains the 
#'    condition labels for the two conditions that are being compared.
#' @param minInSpan positive integer that represents the minimum number of
#'    nucleotides in a smoothing span window if \code{smooth} is TRUE.  
#'    Default value is 10.
#' @param minNumRegion positive integer that represents the minimum number of
#'    nucleotides to consider for a candidate region. Default value is 5.
#' @param cutoff scalar value that represents the absolute value (or a vector 
#'    of two numbers representing a lower and upper bound) for the cutoff of 
#'    the single nucleotide condition coefficient that is used to discover 
#'    candidate regions. 
#' @param smooth logical value that indicates whether or not to smooth the 
#'    nucleotide level signal when discovering candidate regions.
#'    Defaults to FALSE.
#' @param bpSpan a positive integer that represents the length in basepairs
#'    of the smoothing span window if \code{smooth} is TRUE.  Default value is 
#'    100
#' @param verbose logical value that indicates whether progress messages
#'    should be printed to stout. Defaults value is TRUE.
#' @param workers positive integer that represents the number of cores to 
#'    use if parallelization is desired of the smoothing step. 
#' @param logT logical value that indicates whether to model the log2 
#'    transformed signal (plus a pseudocount of 1).  Default is TRUE.  Only
#'    set to false if transformation has been done prior to running this 
#'    function, or if distribution of raw values looks relatively symmetric.
#' @param sampleSize positive integer that represents the number of samples
#'    in each condition.  Defaults to \code{(ncol(OligoSignal)-1)/2}.
#' @param maxPerms a positive integer that represents the maximum number 
#'    of permutations that will be used to generate the global null 
#'    distribution of test statistics.
#' @return a data.frame that contains the results of the inference. The
#'    data.frame contains one row for each candidate region, and 
#'    9 columns, in the following order: 1. chr = 
#'    region level labels such as chromosome, gene, or lncRNA, 2. start = 
#'    start basepair position of the region, 3. end = end basepair position
#'    of the region,
#'    4. indexStart = the index of the region's starting nucleotide, 
#'    5. indexEnd = the index of the region's ending nucleotide,
#'    6. length = the number of nucleotides contained in the region,
#'    7. stat = the test statistic for the condition difference,
#'    8. pval = the permutation p-value for the significance of the test
#'    statistic, and 9. qval = the q-value for the test statistic (adjustment
#'    for multiple comparisons to control false discovery rate).
#' @keywords inference
#' @inheritParams bumphunt
#' @importFrom stats model.matrix
#' @importFrom stats p.adjust
#' @importFrom utils combn
#' @export
#' @examples 
#' \dontrun{
#' normalizedCounts <- normalize(rawCounts = system.file("extdata", 
#' "allTranscriptsCounts_Raw.tsv", package = "oligoGames"))
#' metaData <- system.file("extdata", "oligoMeta.tsv", package = "oligoGames")
#' oligoLen <- 110
#' conditionLabels <- c("Nuclei", "Total")
#' modeledNucs <- modelNucCounts(normalizedCounts, metaData, 
#' conditionLabels, modelMethod = "median", oligoLen = 110)
#' DRregions <- DRfinder(modeledNucs, conditionLabels, minInSpan = 5, 
#' bpSpan = 50, minNumRegion = 3, cutoff = 0.05, smooth = TRUE, verbose = TRUE, 
#' workers = 1, sampleSize = 4, maxPerms = 50)
#' }

DRfinder <- function(OligoSignal, 
                     conditionLabels=c("condition1", "condition2"),
                     minInSpan=10,
                     bpSpan=100, 
                     minNumRegion=5, 
                     cutoff = NULL, 
                     smooth = FALSE,
                     verbose = TRUE,
                     workers = NULL, 
                     sampleSize=(ncol(OligoSignal)-1)/2,
                     maxPerms=50, logT=TRUE,coef=2){
  
  cond <- c(rep(conditionLabels[1], sampleSize), 
            rep(conditionLabels[2], sampleSize))
  design = model.matrix(~cond)
  
  nlocs = nrow(OligoSignal) 
  oligo.mat = OligoSignal[,-1]
  
  chr = rownames(OligoSignal)
  pos = OligoSignal[,1]
  rm(OligoSignal)
  
  maxGap <- maxGapSmooth <-max(diff(pos))+1
  
  # get observed stats
  
  res = bumphunt(oligo.mat, 
                 minInSpan=minInSpan,
                 design = design, chr = chr, pos = pos, 
                 coef = coef, minNum=minInSpan, maxGapSmooth=maxGapSmooth, 
                 minNumRegion=minNumRegion,
                 cutoff = cutoff, maxGap = maxGap, 
                 smooth = smooth, 
                 verbose = verbose,
                 workers = workers, logT=logT) 
  
  
  if (nrow(design)%%2==0){
    perms <- combn(seq(1, nrow(design)), sampleSize)
    perms <- perms[, 2:(ncol(perms)/2)]
    res.flip <- NULL
    
    if (maxPerms < ncol(perms)){
      # subset on 'balanced perms'
      if (sampleSize > 3 & sampleSize < 6){ 
        sg <- apply(perms, 2, function(x) sum(x > sampleSize))
        perms <- perms[, sg < (sampleSize-1) & sg >= 2 ]
        maxPerms <- min(maxPerms, ncol(perms))
      }else if(sampleSize >= 6){
        sg <- apply(perms, 2, function(x) sum(x > sampleSize))
        perms <- perms[, sg >= floor(sampleSize/2) & sg <= ceiling(sampleSize/2) ]
      }
      perms <- perms[,sort(sample(1:ncol(perms), maxPerms, replace=FALSE))]
    }
    
    # Now rerun on flipped designs and concatenate results
    for (j in 1:ncol(perms)){
      reorder <- perms[,j]
      designr <- design
      designr[,2] <- 0
      designr[reorder,2] <- 1
      
      res.flip.p = bumphunt(oligo.mat, 
                            minInSpan=minInSpan,
                            design = designr, chr = chr, pos = pos, 
                            coef = coef, minNum=minInSpan, 
                            maxGapSmooth=maxGapSmooth, 
                            minNumRegion=minNumRegion,
                            cutoff = cutoff,  maxGap = maxGap, 
                            smooth = smooth,
                            verbose = verbose,
                            workers = workers, logT=logT) 
      res.flip.p$permNum <- j 
      
      res.flip <- rbind(res.flip, res.flip.p)
      print(paste0(j, " out of ", ncol(perms), " permutations completed"))
    }	
  }else{
    stop("Error: Currently only balanced designs supported")
  }
  rm(res.flip.p)
  
  perm.ordered <- c(sort(abs(res$stat), method="quick"), Inf)
  pval <- rep(NA, nrow(res))
  pval[!is.na(res$stat)] <- (1 + 
                               sapply(abs(res$stat[!is.na(res$stat)]), 
                                      function(x) length(perm.ordered) - 
                                        min(which(x <= perm.ordered)))) /
    (1 + sum(!is.na(res$stat)))							
  res$pval <- pval
  res$qval <- p.adjust(pval, method="BH")
  return(res)
}


#' Detect and score candidate regions
#'
#' This is an internal workhorse function called by \code{DRfinder} that
#' calculates the nucleotide-level signal, and calls the \code{regionFinder} 
#' function to determine candidate regions and score them.
#'
#' @param oligo.mat a matrix that contains the nucleotide level counts that 
#'    has one row per nucleotide and
#'    one column per sample.
#' @param design a model matrix with one row per sample and one column per 
#'    independent covariate.  
#' @param chr a character vector of labels for region-level characteristics, 
#'   with length equal to the number of rows in \code{oligo.mat} (and in the
#'   same order).  This can indicate the chromosome, gene, lncRNA, etc.
#' @param pos a numeric vector of basepair positions for each nucleotide in
#'   \code{oligo.mat} (and in the same order).
#' @param coef positive integer that indicates which column of the design
#'   matrix in \code{design} contains the condition covariate of interest
#' @param maxGap positive integer that indicates the maximum number of basepairs
#'   that can separate two nucleotides before they will be divided into two 
#'   separate candidate regions.  Defaults to 50.
#' @param maxGapSmooth positive integer that indicates the maximum number of 
#'   basepairs
#'   that can separate two nucleotides before they will be divided into two 
#'   separate smoothing regions.  Defaults to 50.
#' @param minNum positive integer that represents the minimum number of 
#'    nucleotides overall in a region to be smoothed (if \code{smooth} is TRUE).
#'    Default value is 10
#' @return a data.frame that contains the results of region detection. The
#'    data.frame contains one row for each candidate region, and 
#'    7 columns, in the following order: 1. chr = 
#'    region level labels such as chromosome, gene, or lncRNA, 2. start = 
#'    start basepair position of the region, 3. end = end basepair position
#'    of the region,
#'    4. indexStart = the index of the region's starting nucleotide, 
#'    5. indexEnd = the index of the region's ending nucleotide,
#'    6. length = the number of nucleotides contained in the region,
#'    and 7. stat = the test statistic for the condition difference.
#' @inheritParams DRfinder
#' @import bumphunter
#' @import foreach
#' @import doParallel
#' @keywords inference
  
##The WGBS version of the 'BumphunterEngine' function with 'LM' and 'GLM'
bumphunt = function(oligo.mat, design, 
                    chr = NULL, pos, coef = 2, minInSpan=10, minNum=10, 
                    minNumRegion=5,
                    cutoff = NULL, maxGap = 50, maxGapSmooth=50,
                    smooth = FALSE, bpSpan=100,
                    verbose = TRUE, workers=NULL, logT=TRUE, ...)
{
  if (!is.matrix(oligo.mat))
    stop("'oligo.mat' must be a matrices.")
  if (ncol(oligo.mat) != nrow(design))
    stop("Total number of columns in 'oligo.mat' must  
         match number of rows of 'design'")
  if (!(is.null(cutoff) || length(cutoff) %in% 1:2))
    stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
  if (length(cutoff) == 2)
    cutoff <- sort(cutoff)
  
  if (!getDoParRegistered())
    registerDoSEQ()
  registerDoParallel()
  if (is.null(workers)) { workers <- 1 }
  backend <- getDoParName()
  version <- getDoParVersion()
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose) {
    if (workers == 1) {
      mes <- "Using a single core (backend: %s, version: %s)."
      message(sprintf(mes, backend, version))
    }else {
      mes <- "Parallelizing using %s workers/cores (backend: %s, version: %s)."
      message(sprintf(mes, workers, backend, version))
    }
  }
  if (is.null(chr))
    chr <- rep("Unspecified", length(pos))
  if (is.factor(chr))
    chr <- as.character(chr)
  
  cluster <- bumphunter::clusterMaker(chr, pos, maxGap = maxGap)
  
  if (verbose)
    message("Computing coefficients.")
  
  if(logT){
    est <- bumphunter::.getEstimate(log2(oligo.mat+1), 
                                     design, coef, full=TRUE) 
  }else{
    est <- bumphunter::.getEstimate(oligo.mat, 
                                     design, coef, full=TRUE) 
  }
  rawBeta <- est$coef
  
  if (smooth) {
    sd.raw = est$sigma*est$stdev.unscaled
    weights <- sd.raw 
    weights[sd.raw < 10^-5] = mean(sd.raw[sd.raw > 10^-5]) 
    rm(sd.raw)
  }
  rm(est)
  gc()
  
  if (smooth) {
    if (verbose)
      message("Smoothing coefficients.")
    
    beta <- vector("list", 3)
    beta[[1]] <- beta[[2]] <- rep(NA, length(pos))
    for (chromosome in unique(chr)){
      beta.tmp <- smoother(y = rawBeta[chr==chromosome],
                           x = pos[chr==chromosome],
                           workers=workers, chr=chr[chr==chromosome],
                           maxGapSmooth=maxGapSmooth, weights = weights,
                           minNum = minNum, minInSpan = minInSpan, bpSpan = bpSpan,
                           verbose = verbose)
      beta[[1]][chr==chromosome] <- beta.tmp[[1]]
      beta[[2]][chr==chromosome] <- beta.tmp[[2]]
    }
    beta[[3]] <- beta.tmp[[3]]
    names(beta) <- names(beta.tmp)
    rm(beta.tmp)
    Index <- which(beta$smoothed)
    beta <- beta$fitted
    
  }else {
    beta <- rawBeta
    Index <- seq(along = beta)
  }
  rm(rawBeta)
  gc()
  if (verbose)
    message("Finding regions.")
  tab <- regionFinder(x = beta, chr = chr, pos = pos, cluster = cluster,
                      cutoff = cutoff, ind = Index, minNumRegion = minNumRegion,
                      oligo.mat = oligo.mat,
                      verbose = TRUE,
                      design = design, workers=workers, logT=logT)
  rm(beta);
  rm(chr);
  rm(pos);
  gc()
  if (nrow(tab) == 0) {
    if (verbose)
      message("No candidate regions found!")
    return(NA)
  }else {
    if (verbose)
      message(sprintf("Found %s candidate regions.", nrow(tab)))
  }
  
  return(table = tab)
}



#' Break up nucleotide level signal into candidate regions and score them 
#'
#' This is an internal workhorse function for \code{bumphunt} that takes the
#' nucleotide-level signal and parses it into contigous regions that pass
#' the threshold and form
#' the candidates, and then scores each one based on a test statistic
#' of the difference.
#' @param x a vector of condition coefficients (for the covariate of interest)
#'  for each nucletide
#' @param cluster a vector of cluster membership values for each nucleotide
#'  determined by the \code{clusterMaker} function in the \code{bumphunter}
#'  package
#' @param ind a vector if indices of \code{x} which are non-NULL. Defaults to 
#'  all indices of \code{x}.
#' @param order logical that indicates whether or not to order the candidate
#'  regions by the test statistic magnitude (largest to smallest).
#'  Defaults to TRUE.
#' @param assumeSorted logical that indicates whether the nucleotides are 
#'  sorted in ascending order.  Defaults to FALSE.
#' @inheritParams DRfinder
#' @inheritParams bumphunt
#' @return a data.frame that contains the results of region detection. The
#'    data.frame contains one row for each candidate region, and 
#'    7 columns, in the following order: 1. chr = 
#'    region level labels such as chromosome, gene, or lncRNA, 2. start = 
#'    start basepair position of the region, 3. end = end basepair position
#'    of the region,
#'    4. indexStart = the index of the region's starting nucleotide, 
#'    5. indexEnd = the index of the region's ending nucleotide,
#'    6. length = the number of nucleotides contained in the region,
#'    and 7. stat = the test statistic for the condition difference.
#' @import bumphunter
#' @import reshape2
#' @importFrom stats lm
#' @importFrom splines ns
#' @importFrom stats quantile
#' @importFrom parallel mclapply
#' @keywords inference

regionFinder <- function(x, chr, pos, cluster=NULL,
                         ind=seq(along=x),order=TRUE, minNumRegion=5,
                         maxGap, cutoff=quantile(abs(x), 0.99),
                         assumeSorted = FALSE, oligo.mat=oligo.mat,
                         verbose = TRUE,
                         design=design, workers=workers, logT=TRUE){
  if(any(is.na(x[ind]))){
    warning("NAs found and removed. ind changed.")
    ind <- intersect(which(!is.na(x)),ind)
  }
  if(is.null(cluster))
    cluster <- bumphunter::clusterMaker(chr, pos, maxGap=maxGap, 
                                         assumeSorted = assumeSorted)
  Indexes <- bumphunter::getSegments(x = x[ind], f = cluster[ind], 
                                      cutoff = cutoff, 
                                      assumeSorted = assumeSorted, 
                                      verbose = verbose)
  clusterN <- table(cluster)[as.character(cluster)]

  # only keep up and down indices
  Indexes <- Indexes[1:2]

  # only keep clusters that have at least minNumRegion elements
  for(i in 1:2){
    lns <- sapply(Indexes[[i]], length)
    Indexes[[i]] <- Indexes[[i]][lns >= minNumRegion]
  }


  stat.OligoSignal.spline <- function(ix){
    X = as.vector(sapply(design[,2], function(x) rep(x, nrow(oligo.mat[ix,]))))
    L = as.vector(rep(pos[ix], nrow(design)))
    if(logT){
      Y = log2(oligo.mat[ix,] + 1)
    }else{
      Y = oligo.mat[ix,]
    }
    Y <- melt(Y)

    Y$X <- X
    Y$L <- L
    int.knots <- ceiling((max(pos[ix]) - min(pos[ix])) / 100) + 1
    return(summary(lm(value~ X + ns(L, df=int.knots), data=Y))$coef[2,3])
  }

  cov.upper <- quantile(oligo.mat, 0.9)
  t1 <- proc.time()
  res <- vector("list",2)
  for(i in 1:2){
    res[[i]]<-
      data.frame(chr=sapply(Indexes[[i]], function(Index) chr[ind[Index[1]]]),
              start=sapply(Indexes[[i]], function(Index) min(pos[ind[Index]])),
              end=sapply(Indexes[[i]], function(Index) max(pos[ind[Index]])),
              indexStart=sapply(Indexes[[i]], function(Index) min(ind[Index])),
              indexEnd = sapply(Indexes[[i]], function(Index) max(ind[Index])),
              length = sapply(Indexes[[i]], length), stringsAsFactors=FALSE)

    res[[i]]$stat = unlist(mclapply(Indexes[[i]], function(Index)
        stat.OligoSignal.spline(ind[Index]),
        mc.cores=workers))
  }
  t2 <- proc.time()
  t2 - t1
  names(res) <- c("up","dn")
  res <- rbind(res$up,res$dn)
  if(order & nrow(res)>0) res <- res[order(-res$stat),]

  return(res)
}


#' A smoother function for coefficient values
#'
#' This is an internal workhorse function that (optionally) smooths the
#' coefficient values for the covariate of interest.
#' @param weights a numeric vector that contains the smoothing weights.  
#' @inheritParams DRfinder
#' @inheritParams bumphunt
#' @inheritParams regionFinder
#' @return a list object of length four that contains the following objects 
#'   in order: 1. fitted = a numeric matrix containing the coefficient 
#'   estimates after smoothing those that satisify the smoothing conditions 
#'   for each nucleotide,
#'   2. smoothed = a logical vector that indicates whether each nucleotide 
#'   was smoothed, 3. smoother = a character that represents the type 
#'   of smoothing that was performed, and 4. idx a numeric vector that 
#'   indicates which nucleotides are non-NULL.
#' @import bumphunter
#' @import locfit
#' @import doRNG
#' @importFrom iterators iter
#' @importFrom stats preplot
#' @keywords inference
smoother <- function(y, x=NULL, weights=NULL, chr=chr, 
                     maxGapSmooth=maxGapSmooth,
                     verbose=TRUE, workers=1, minNum = minNum,
                     minInSpan = minInSpan, bpSpan = bpSpan) {


  locfitByCluster2 <- function(ix) {

    ## if y is vector change to matrix
    yi = matrix(y[ix], ncol=1)
    xi = x[ix]
    weightsi = matrix(weights[ix], ncol=1)
    clusteri = clusterC[ix]


    if(is.null((yi)))
      stop("y (rawBeta) is missing")
    if(is.null(xi))
      stop("x (pos) is missing")
    if(is.null(clusteri))
      stop("cluster is missing")
    if(is.null(weightsi))
      weightsi <- matrix(1, nrow = nrow(yi), ncol = ncol(yi))

    Indexes <- split(seq(along = clusteri), clusteri)
    clusterL <- sapply(Indexes, length)
    smoothed <- rep(TRUE, nrow(yi))

    for(i in seq(along=Indexes)) {
      #if(verbose) if(i %% 1e4 == 0) cat(".") # print the cluster index
      Index <- Indexes[[i]]
      if(clusterL[i] >= minNum & 
         sum(rowSums(is.na(yi[Index,,drop=FALSE]))==0) >= minNum)  {
        nn <- minInSpan / length(Index)
        for(j in 1:ncol(yi)) {
          sdata <- data.frame(posi = xi[Index],
                              yi = yi[Index, j],
                              weightsi = weightsi[Index,j])
          fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan), data = sdata,
                        weights = weightsi, family = "gaussian", maxk = 10000)
          pp <- preplot(fit, where = "data", band = "local",
                        newdata = data.frame(posi = xi[Index]))
          yi[Index,j] <- pp$trans(pp$fit)
        }
      } else {
        yi[Index,] <- NA
        smoothed[Index] <- FALSE
      }
    }
    return(list(fitted=yi, smoothed=smoothed, smoother="locfit"))
  }

  if(is.null(dim(y)))
    y <- matrix(y, ncol=1) 
  if(!is.null(weights) && is.null(dim(weights)))
    weights <- matrix(weights, ncol=1)
  if (!getDoParRegistered())
    registerDoSEQ()

  ret.all <- NULL
  cores <- workers
  # loop through each chromosome to mitigate memory spikes
  for (chromosome in unique(chr)){

    clusterC <- bumphunter::clusterMaker(chr[chr==chromosome], 
                                          x[chr==chromosome], 
                                          maxGap = maxGapSmooth)


    Indexes <- split(seq(along=clusterC), clusterC)
    IndexesChunks <- vector("list", length = cores)
    baseSize <- length(Indexes) %/% cores
    remain <- length(Indexes) %% cores
    done <- 0L
    for(ii in 1:cores) {
      if(remain > 0) {
        IndexesChunks[[ii]] <- done + 1:(baseSize + 1)
        remain <- remain - 1L
        done <- done + baseSize + 1L
      } else {
        IndexesChunks[[ii]] <- done + 1:baseSize
        done <- done + baseSize
      }
    }


    IndexesChunks <- lapply(IndexesChunks, function(idxes) {
      do.call(c, unname(Indexes[idxes]))
    })

    idx <- NULL ## for R CMD check
    ret <- foreach(idx = iter(IndexesChunks)) %dorng% {
      sm <- locfitByCluster2(idx)
      c(sm, list(idx = seq(1, length(chr))[chr==chromosome][idx]))
    }

    attributes(ret)[["rng"]] <- NULL
    ## Paste together results from different workers
    ret <- bumphunter::reduceIt(ret)

    if(is.null(ret.all)){
      ret.all <- ret
    }else{
      ret.all$smoother <- c(ret.all$smoother, ret$smoother)
      ret.all$fitted <- c(ret.all$fitted, ret$fitted)
      ret.all$smoothed <- c(ret.all$smoothed, ret$smoothed)
      ret.all$idx <- c(ret.all$idx, ret$idx)
    }
    rm(ret)
  }

  ## Now fixing order issue
  revOrder <- ret.all$idx
  names(revOrder) <- seq_along(ret.all$idx)
  revOrder <- sort(revOrder)
  revOrder <- as.integer(names(revOrder))

  ret.all$smoother <- ret.all$smoother[1]
  ret.all$fitted <- ret.all$fitted[revOrder,,drop=FALSE]
  ret.all$smoothed <- ret.all$smoothed[revOrder]
  ret.all$idx <- NULL

  return(ret.all)
}
